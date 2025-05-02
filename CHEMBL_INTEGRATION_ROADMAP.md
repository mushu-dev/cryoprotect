# ChEMBL Integration Implementation Roadmap

## Overview
This document provides a comprehensive roadmap for resolving the ChEMBL data acquisition issues and successfully integrating ChEMBL data into the CryoProtect v2 database. It includes specific file paths, code references, and implementation guidance for ROO agents.

## Critical Resources

### MCP Configuration
- **MCP Location**: `.roo/mcp.json` - Primary MCP configuration file for all ROO agents
- **Access Pattern**: Use the MCP SQL execution capabilities for all database operations
- **Authentication**: MCP already has the necessary permissions for successful database operations

### Key Implementation Files

#### ChEMBL Integration
1. **Primary Integration Script**: 
   - **Path**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ChEMBL_CryoProtectants_Supabase.py`
   - **Key Sections**:
     - Lines 65-125: Main data acquisition routine
     - Lines 180-240: Data transformation logic
     - Lines 320-380: Error handling (needs enhancement)

2. **ChEMBL Client**: 
   - **Path**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/chembl/client.py`
   - **Key Sections**:
     - Lines 45-70: API request handling
     - Lines 110-140: Rate limiting implementation

3. **Cache Management**:
   - **Path**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/chembl/cache.py`
   - **Cache Storage**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/cache/chembl/`

#### Database Schema
1. **Core Table Structure**:
   - **Path**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/001_initial_schema.sql`
   - **Key Tables**:
     - Lines 25-50: `molecules` table definition
     - Lines 55-70: `molecular_properties` table definition

2. **Schema Verification**:
   - **Path**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/verify_database_integrity.py`
   - **Key Functions**: 
     - `check_table_structure()`: Validates table schema
     - `check_foreign_key_integrity()`: Verifies relationships

#### RLS Policy Management
1. **RLS Implementation**:
   - **Path**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/006_rls_policies.sql`
   - **Path**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/007_service_role_rls.sql`

2. **Service Role Authentication**:
   - **Path**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/service_role_helper.py`

## Issue Resolution Strategy

### Phase 1: Diagnostic Investigation

1. **MCP SQL Execution Test**
   Execute a simple query using the correct MCP path (`.roo/mcp.json`):
   ```sql
   SELECT current_setting('role'), current_user;
   ```
   This will verify the role being used for database access.

2. **RLS Policy Inspection**
   Execute this SQL through MCP to list all policies affecting the molecules table:
   ```sql
   SELECT tablename, policyname, permissive, roles, cmd, qual, with_check 
   FROM pg_policies 
   WHERE tablename = 'molecules';
   ```

3. **Schema Validation**
   Execute this SQL through MCP to verify table structure:
   ```sql
   SELECT column_name, data_type, is_nullable
   FROM information_schema.columns
   WHERE table_name = 'molecules'
   ORDER BY ordinal_position;
   ```

### Phase 2: Critical Fixes

1. **RLS Bypass Implementation**
   Execute this SQL through MCP to temporarily disable RLS for data import:
   ```sql
   -- Option 1: Temporarily disable RLS completely
   ALTER TABLE molecules DISABLE ROW LEVEL SECURITY;
   
   -- Option 2: Create a permissive policy for service role
   CREATE POLICY "Service role has full access to molecules"
   ON molecules
   FOR ALL
   USING (auth.jwt() ->> 'role' = 'service_role')
   WITH CHECK (auth.jwt() ->> 'role' = 'service_role');
   ```

2. **Permission Adjustment**
   Execute this SQL through MCP to ensure proper permissions:
   ```sql
   GRANT ALL PRIVILEGES ON TABLE molecules TO authenticated;
   GRANT ALL PRIVILEGES ON TABLE molecular_properties TO authenticated;
   GRANT USAGE, SELECT ON SEQUENCE molecules_id_seq TO authenticated;
   ```

3. **Service Role Authentication Enhancement**
   Modify the `ChEMBL_CryoProtectants_Supabase.py` script:
   - Replace direct Supabase client operations with MCP SQL execution
   - Implement proper error handling with detailed error messages
   - Add transaction management for batch operations

### Phase 3: Modified Implementation Approach

1. **ChEMBL Data Transformation**
   Implement a robust transformation function:
   ```python
   def transform_chembl_to_molecule(chembl_compound):
       """Transform ChEMBL compound data to match our molecules table schema."""
       return {
           "name": chembl_compound.get("pref_name") or chembl_compound.get("molecule_chembl_id"),
           "smiles": chembl_compound.get("canonical_smiles"),
           "inchi": chembl_compound.get("std_inchi"),
           "inchikey": chembl_compound.get("std_inchi_key"),
           "formula": chembl_compound.get("molecular_formula"),
           "molecular_weight": chembl_compound.get("full_mwt"),
           "created_by": user_profile_id,
           "data_source": "ChEMBL Integration",
           "version": 1
       }
   ```

2. **Direct SQL Insert via MCP**
   Replace Supabase client operations with direct SQL:
   ```python
   # Instead of client.table('molecules').insert(...)
   molecule_json = json.dumps(molecule_data)
   query = f"""
   INSERT INTO molecules (name, smiles, inchi, inchikey, formula, molecular_weight, created_by, data_source, version)
   VALUES (
       '{molecule_data["name"]}',
       '{molecule_data["smiles"]}',
       '{molecule_data["inchi"]}',
       '{molecule_data["inchikey"]}',
       '{molecule_data["formula"]}',
       {molecule_data["molecular_weight"]},
       '{molecule_data["created_by"]}',
       '{molecule_data["data_source"]}',
       {molecule_data["version"]}
   )
   ON CONFLICT (inchikey) DO NOTHING
   RETURNING id;
   """
   
   # Execute using MCP
   result = execute_sql_through_mcp(query)
   ```

3. **Batch Processing with Transactions**
   Implement batch processing with proper transaction handling:
   ```python
   def batch_insert_molecules(molecules_batch):
       transaction_queries = ["BEGIN;"]
       
       for molecule in molecules_batch:
           # Create escaped query for each molecule
           query = f"""
           INSERT INTO molecules (name, smiles, inchi, inchikey, formula, molecular_weight, created_by, data_source, version)
           VALUES (
               '{escape_sql_string(molecule["name"])}',
               '{escape_sql_string(molecule["smiles"])}',
               '{escape_sql_string(molecule["inchi"])}',
               '{escape_sql_string(molecule["inchikey"])}',
               '{escape_sql_string(molecule["formula"])}',
               {molecule["molecular_weight"] or 'NULL'},
               '{molecule["created_by"]}',
               '{escape_sql_string(molecule["data_source"])}',
               {molecule["version"]}
           )
           ON CONFLICT (inchikey) DO NOTHING
           RETURNING id;
           """
           transaction_queries.append(query)
       
       transaction_queries.append("COMMIT;")
       combined_query = "\n".join(transaction_queries)
       
       # Execute using MCP
       return execute_sql_through_mcp(combined_query)
   ```

### Phase 4: Implementation Verification

1. **Progressive Testing**
   Implement a step-by-step validation process:
   ```python
   def verify_chembl_integration():
       """Verify ChEMBL integration by checking database state."""
       # Check molecules table
       molecules_count = execute_sql_through_mcp("SELECT COUNT(*) FROM molecules WHERE data_source LIKE '%ChEMBL%';")
       print(f"ChEMBL molecules count: {molecules_count[0]['count']}")
       
       # Check properties table
       properties_count = execute_sql_through_mcp("""
           SELECT COUNT(*) FROM molecular_properties mp
           JOIN molecules m ON mp.molecule_id = m.id
           WHERE m.data_source LIKE '%ChEMBL%';
       """)
       print(f"ChEMBL properties count: {properties_count[0]['count']}")
       
       # Verify a sample molecule
       sample = execute_sql_through_mcp("""
           SELECT * FROM molecules
           WHERE data_source LIKE '%ChEMBL%'
           LIMIT 1;
       """)
       print(f"Sample ChEMBL molecule: {sample[0] if sample else 'None'}")
   ```

2. **Re-Enable RLS After Import**
   Execute after successful data import:
   ```sql
   -- If RLS was completely disabled
   ALTER TABLE molecules ENABLE ROW LEVEL SECURITY;
   
   -- If a temporary policy was created
   DROP POLICY IF EXISTS "Service role has full access to molecules" ON molecules;
   ```

3. **Document Import Statistics**
   Generate a comprehensive import report:
   ```python
   def generate_import_report():
       """Generate a report of the ChEMBL integration results."""
       stats = {}
       
       # Get total counts
       stats["total_molecules"] = execute_sql_through_mcp("SELECT COUNT(*) FROM molecules;")[0]["count"]
       stats["chembl_molecules"] = execute_sql_through_mcp("SELECT COUNT(*) FROM molecules WHERE data_source LIKE '%ChEMBL%';")[0]["count"]
       stats["total_properties"] = execute_sql_through_mcp("SELECT COUNT(*) FROM molecular_properties;")[0]["count"]
       
       # Get property distribution
       property_counts = execute_sql_through_mcp("""
           SELECT pt.name, COUNT(mp.id) as count
           FROM molecular_properties mp
           JOIN property_types pt ON mp.property_type_id = pt.id
           JOIN molecules m ON mp.molecule_id = m.id
           WHERE m.data_source LIKE '%ChEMBL%'
           GROUP BY pt.name
           ORDER BY count DESC;
       """)
       stats["property_distribution"] = property_counts
       
       return stats
   ```

## Performance Optimization Guidelines

1. **ChEMBL API Rate Limiting**
   - Configure proper rate limiting in `chembl/rate_limiter.py`
   - Set `MAX_REQUESTS_PER_SECOND = 5` to avoid API throttling
   - Implement exponential backoff for failed requests

2. **Batch Processing Configuration**
   - Optimal batch size for ChEMBL integration: 50-100 compounds
   - Use transaction batching to minimize database round trips
   - Implement checkpointing to resume interrupted operations

3. **Error Handling & Resilience**
   - Catch and log specific error types separately
   - Implement retry logic for transient failures
   - Use cache to avoid redundant API calls

## Next Steps After Successful Integration

1. **Data Quality Verification**
   - Run validation scripts to verify scientific accuracy
   - Check for missing values in critical fields
   - Verify property value distributions 

2. **Derived Data Generation**
   - Generate mixtures based on chemical properties
   - Create experimental protocols using integrated compounds
   - Generate predictions using computational models

3. **Performance Testing**
   - Test API response times with populated database
   - Verify query performance for common operations
   - Ensure RLS policies don't introduce performance issues

## Success Criteria

1. **Data Volume Targets**
   - At least 1,000 molecules from ChEMBL
   - At least 10,000 molecular properties
   - All critical properties populated for each molecule

2. **Integration Quality**
   - No structural errors in imported data
   - Proper normalization of chemical identifiers
   - Consistent property value units

3. **System Performance**
   - Database queries complete in < 500ms
   - API endpoints handle the data volume efficiently
   - No RLS policy conflicts with normal operation

## Troubleshooting Guide

### Common Issues and Solutions

1. **"Error checking if molecule exists"**
   - **Cause**: Usually permissions or RLS policy issues
   - **Solution**: Execute SQL to verify and fix permissions:
     ```sql
     -- Check permissions
     SELECT grantee, privilege_type 
     FROM information_schema.role_table_grants 
     WHERE table_name = 'molecules';
     
     -- Grant if missing
     GRANT SELECT, INSERT ON molecules TO authenticated;
     ```

2. **"Error inserting molecule"**
   - **Cause**: Schema mismatch or constraint violation
   - **Solution**: Verify the transformation function aligns with table constraints

3. **"SMTP connection timeout"**
   - **Cause**: ChEMBL API rate limiting or network issues
   - **Solution**: Implement API call retry with exponential backoff

4. **"Transaction aborted due to statement timeout"**
   - **Cause**: Large transaction taking too long
   - **Solution**: Reduce batch size and increase statement timeout:
     ```sql
     SET statement_timeout = '300s';
     ```