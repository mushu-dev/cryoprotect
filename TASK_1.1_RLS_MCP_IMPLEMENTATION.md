# Task 1.1: RLS Policy Implementation via MCP

This document provides detailed instructions for implementing missing RLS policies via Supabase MCP for the CryoProtect v2 project.

## Objective

Implement missing Row Level Security (RLS) policies for the following tables and views:
- experiment_with_results
- migrations
- mixture_with_components
- molecule_with_properties

Additionally, create a scientific data audit system with performance optimizations for RLS policies.

## Implementation Files

1. **SQL Script**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/supabase_rls_policies.sql`
   - Contains all RLS policies as a single atomic transaction
   - Creates performance indexes, audit system, and verification functions

2. **MCP Execution Script**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/execute_rls_sql_via_mcp.py`
   - Executes the SQL script via Supabase MCP
   - Captures and logs execution results
   - Generates verification report

3. **MCP Tools**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/supabase_mcp_tools.py`
   - Contains utility functions for interacting with Supabase via MCP
   - Handles SQL execution, table listing, and error handling

## Implementation Steps

### 1. Environment Setup

Ensure the following environment variables are set:
- `SUPABASE_PROJECT_ID` - Your Supabase project ID (default: "tsdlmynydfuypiugmkev")
- Access token is configured in the MCP tools file

### 2. Execute RLS Implementation

Run the script to implement RLS policies:

```bash
python execute_rls_sql_via_mcp.py
```

The script will:
1. Read the SQL file
2. Execute it via Supabase MCP
3. Capture verification results
4. Generate a detailed execution report
5. Display implementation status

### 3. Verify Implementation

The implementation includes self-verification that checks:
- Tables with RLS enabled
- Views with SECURITY INVOKER set
- Audit triggers created
- Performance indexes created

A detailed verification report will be saved to:
`reports/security/rls_execution_report_YYYYMMDD_HHMMSS.json`

### 4. Implementation Details

The RLS implementation includes:

1. **Performance Indexes**:
   - Project membership indexes (e.g., `idx_user_profile_project_id`)
   - Foreign key indexes for relationship traversal (e.g., `idx_molecular_property_molecule_id`)

2. **View Security**:
   - SECURITY INVOKER for all views to respect underlying table RLS

3. **Scientific Data Audit System**:
   - `scientific_data_audit` table for tracking all data operations
   - Audit triggers on molecule, mixture, and experiment tables
   - Bypass option for bulk data loading

4. **Service Role Policies**:
   - Access policies for the service role to manage scientific data
   - Required for database population

## Troubleshooting

If the implementation fails:

1. **Check Logs**:
   - See `logs/rls_implementation_YYYYMMDD_HHMMSS.log` for error details
   - Error report in `reports/security/rls_execution_error_YYYYMMDD_HHMMSS.json`

2. **Common Issues**:
   - MCP tools connection issues
   - SQL syntax errors
   - Permissions problems
   - Table/view doesn't exist

3. **Recovery Options**:
   - Correct any issues in the SQL file
   - Re-run the implementation script (it's designed to be idempotent)

## Next Steps After Implementation

Once RLS policies are successfully implemented:

1. Proceed to Task 1.2: Database Population Enhancement
2. Test RLS effectiveness with different user roles
3. Configure audit logging for scientific operations
4. Use `set_audit_bypass(true)` before bulk operations

## Expected Outcome

After successful implementation, you should have:
1. RLS enabled on all required tables and views
2. Performance indexes for RLS policies
3. Scientific data audit system
4. Service role policies for data management

These security foundations ensure proper data isolation and access control, which are essential for the scientific database in CryoProtect v2.