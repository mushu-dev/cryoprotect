# Task 2.3: Implement RLS Implementation Fix Module

## Objective
Implement the RLS (Row Level Security) Implementation fix module in the maintenance utility to enable automated fixing of security policy issues.

## Context
Row Level Security is a critical component for data protection in the Supabase database. The `fix_rls_implementation.py` script exists as a standalone solution, but integrating it into the maintenance utility will provide a more consistent approach to database maintenance. RLS verification reports indicate that while progress has been made, the implementation may still be incomplete.

## Acceptance Criteria
- The `fix_rls_implementation` function in the maintenance utility is fully implemented
- All tables have appropriate RLS policies enabled
- Service role access is properly configured
- Verification functionality confirms RLS is working correctly
- The fix can be run via the command line and interactive menu

## Implementation Steps

1. Examine the existing RLS fix script:
   ```bash
   cat fix_rls_implementation.py
   ```

2. Implement the `fix_rls_implementation` function in `maintenance_utils.py`:
   ```python
   def fix_rls_implementation(args):
       """
       Fix RLS implementation issues in the database.
       
       This function:
       1. Enables RLS on all tables
       2. Creates appropriate RLS policies
       3. Configures service role access
       
       Args:
           args: Command line arguments (may include dry_run, verify, etc.)
       
       Returns:
           bool: True if successful, False otherwise
       """
       logger.info("Running RLS Implementation Fix...")
       
       try:
           # Get list of tables
           supabase = get_supabase_client()
           response = supabase.rpc("get_all_tables").execute()
           
           if response.error:
               logger.error(f"Error getting tables: {response.error.message}")
               return False
           
           tables = response.data
           
           # Enable RLS on all tables
           for table in tables:
               table_name = table['table_name']
               
               # Skip certain system tables
               if table_name.startswith('pg_') or table_name in ['schema_migrations', 'spatial_ref_sys']:
                   continue
               
               logger.info(f"Enabling RLS on table {table_name}")
               
               # Enable RLS
               sql = f"ALTER TABLE {table_name} ENABLE ROW LEVEL SECURITY;"
               response = supabase.rpc("execute_sql", {"sql": sql}).execute()
               
               if response.error:
                   logger.error(f"Error enabling RLS on {table_name}: {response.error.message}")
                   continue
               
               # Check if table has created_by column for user-based policies
               has_created_by = False
               response = supabase.rpc("column_exists", {
                   "p_table_name": table_name,
                   "p_column_name": "created_by"
               }).execute()
               
               if response.error:
                   logger.error(f"Error checking for created_by column: {response.error.message}")
               else:
                   has_created_by = response.data
               
               # Create appropriate policies
               if has_created_by:
                   # User access policy - users can view their own data
                   sql = f"""
                   CREATE POLICY "{table_name}_user_select" ON {table_name}
                   FOR SELECT
                   USING (auth.uid() = created_by);
                   """
                   response = supabase.rpc("execute_sql", {"sql": sql}).execute()
                   
                   if response.error and "already exists" not in response.error.message:
                       logger.error(f"Error creating select policy on {table_name}: {response.error.message}")
                   
                   # User modify policy - users can update/delete their own data
                   sql = f"""
                   CREATE POLICY "{table_name}_user_modify" ON {table_name}
                   FOR UPDATE
                   USING (auth.uid() = created_by);
                   """
                   response = supabase.rpc("execute_sql", {"sql": sql}).execute()
                   
                   if response.error and "already exists" not in response.error.message:
                       logger.error(f"Error creating update policy on {table_name}: {response.error.message}")
                   
                   # User delete policy
                   sql = f"""
                   CREATE POLICY "{table_name}_user_delete" ON {table_name}
                   FOR DELETE
                   USING (auth.uid() = created_by);
                   """
                   response = supabase.rpc("execute_sql", {"sql": sql}).execute()
                   
                   if response.error and "already exists" not in response.error.message:
                       logger.error(f"Error creating delete policy on {table_name}: {response.error.message}")
                   
                   # User insert policy - users can insert data with their own user ID
                   sql = f"""
                   CREATE POLICY "{table_name}_user_insert" ON {table_name}
                   FOR INSERT
                   WITH CHECK (auth.uid() = created_by);
                   """
                   response = supabase.rpc("execute_sql", {"sql": sql}).execute()
                   
                   if response.error and "already exists" not in response.error.message:
                       logger.error(f"Error creating insert policy on {table_name}: {response.error.message}")
               
               # Create service role policy for all tables
               sql = f"""
               CREATE POLICY "{table_name}_service_role" ON {table_name}
               FOR ALL
               USING (auth.role() = 'service_role');
               """
               response = supabase.rpc("execute_sql", {"sql": sql}).execute()
               
               if response.error and "already exists" not in response.error.message:
                   logger.error(f"Error creating service role policy on {table_name}: {response.error.message}")
           
           # Create RLS verification report
           timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
           report_path = f"reports/auth/rls_verification_report_{timestamp}.json"
           os.makedirs(os.path.dirname(report_path), exist_ok=True)
           
           report = {
               "timestamp": timestamp,
               "tables_processed": len(tables),
               "rls_enabled": True,
               "policies_created": True
           }
           
           with open(report_path, 'w') as f:
               json.dump(report, f, indent=2)
           
           # Copy as latest report
           latest_path = "reports/auth/rls/latest.json"
           os.makedirs(os.path.dirname(latest_path), exist_ok=True)
           shutil.copy(report_path, latest_path)
           
           logger.info("RLS Implementation fix completed successfully!")
           return True
       except Exception as e:
           logger.error(f"RLS Implementation fix failed: {str(e)}", exc_info=True)
           return False
   ```

3. Update the function mapping in the maintenance utility:
   ```python
   # Ensure the function is listed in the FIX_FUNCTIONS dictionary
   FIX_FUNCTIONS = {
       # Existing mappings...
       "rls_implementation": fix_rls_implementation,
       # Other mappings...
   }
   ```

4. Add verification function:
   ```python
   def verify_rls_implementation():
       """Verify that RLS is properly implemented on all tables."""
       try:
           # Get list of tables
           supabase = get_supabase_client()
           response = supabase.rpc("get_all_tables").execute()
           
           if response.error:
               logger.error(f"Error getting tables: {response.error.message}")
               return False
           
           tables = response.data
           all_good = True
           
           # Check RLS is enabled on each table
           for table in tables:
               table_name = table['table_name']
               
               # Skip certain system tables
               if table_name.startswith('pg_') or table_name in ['schema_migrations', 'spatial_ref_sys']:
                   continue
               
               # Check if RLS is enabled
               sql = f"""
               SELECT relrowsecurity FROM pg_class
               WHERE relname = '{table_name}'
               AND relkind = 'r';
               """
               response = supabase.rpc("execute_sql", {"sql": sql}).execute()
               
               if response.error:
                   logger.error(f"Error checking RLS on {table_name}: {response.error.message}")
                   all_good = False
               elif not response.data or not response.data[0]['relrowsecurity']:
                   logger.error(f"RLS is not enabled on table {table_name}")
                   all_good = False
               
               # Check if policies exist
               sql = f"""
               SELECT COUNT(*) FROM pg_policy
               WHERE tablename = '{table_name}';
               """
               response = supabase.rpc("execute_sql", {"sql": sql}).execute()
               
               if response.error:
                   logger.error(f"Error checking policies on {table_name}: {response.error.message}")
                   all_good = False
               elif not response.data or response.data[0]['count'] == 0:
                   logger.error(f"No policies found for table {table_name}")
                   all_good = False
           
           if all_good:
               logger.info("RLS implementation verification passed!")
           else:
               logger.error("RLS implementation verification failed!")
           
           return all_good
       except Exception as e:
           logger.error(f"RLS implementation verification failed: {str(e)}", exc_info=True)
           return False
   ```

5. Add verification to the main verification function:
   ```python
   def verify_all_fixes():
       """Run all verification functions."""
       # Existing verifications...
       success = success and verify_rls_implementation()
       # More verifications...
       return success
   ```

## Files to Modify
- `maintenance_utils.py` - Add the RLS implementation fix

## Verification
1. Run the utility with the RLS implementation fix:
   ```bash
   python maintenance_utils.py --fix rls_implementation
   ```

2. Check that the fix succeeds without errors

3. Verify the RLS implementation:
   ```bash
   python maintenance_utils.py --verify
   ```

4. Test access with different user roles:
   ```bash
   # Regular user query
   curl -H "Authorization: Bearer $USER_TOKEN" http://localhost:5000/api/v1/molecules
   
   # Service role query
   curl -H "Authorization: Bearer $SERVICE_ROLE_TOKEN" http://localhost:5000/api/v1/molecules
   
   # Anonymous query (should be restricted)
   curl http://localhost:5000/api/v1/molecules
   ```

## Notes for Roo Code Agent
- This implementation creates appropriate RLS policies for all tables
- Special consideration is given to tables with a created_by column for user-specific policies
- The service role policy allows full access for maintenance operations
- Make sure to test with different user roles to ensure policies are working correctly
- Consider additional policies for team-based access if needed
- Document the RLS implementation thoroughly for security audits