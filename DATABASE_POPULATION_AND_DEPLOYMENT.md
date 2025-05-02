# DATABASE POPULATION AND DEPLOYMENT FOR CRYOPROTECT V2

## Overview

The CryoProtect v2 project has completed CI/CD implementation but lacks data in the Supabase database. This prompt focuses on two critical tasks:

1. Implementing RLS policies via Supabase MCP
2. Populating the database with scientific data

## Critical Context

- All required code exists but needs proper execution
- Supabase requires both security (RLS) and data before deployment
- The CI/CD pipeline is ready but needs configured environment variables

## Task 1: Implement RLS Policies

The first step is to implement Row Level Security (RLS) policies for missing tables and views:

1. Execute `supabase_rls_policies.sql` using Supabase MCP
2. The SQL file includes:
   - RLS for experiment_with_results, migrations, mixture_with_components, molecule_with_properties
   - Performance indexes for RLS columns
   - Scientific data audit system
   - Verification functions

### Implementation Steps

1. Examine `/mnt/c/Users/1edwa/Documents/CryoProtect v2/supabase_rls_policies.sql`
2. Create a focused script to execute this SQL using `supabase_mcp_tools.py`
3. Run the script and verify RLS implementation through Supabase
4. Ensure all tables have appropriate RLS policies

## Task 2: Populate Database with Scientific Data

After securing the database with RLS, populate it with scientific data:

1. Review and execute `/mnt/c/Users/1edwa/Documents/CryoProtect v2/populate_database_supabase.py`
2. This script contains:
   - Complete dataset of scientifically accurate cryoprotectants
   - Comprehensive molecular properties
   - Scientifically valid mixtures and experiments
   - Functionality to bypass audit for bulk loading

### Implementation Steps

1. Examine the existing `populate_database_supabase.py` script
2. Enhance for transaction support if needed
3. Execute and verify data population
4. Generate a verification report

## Technical Details

### RLS Implementation

- The Supabase MCP tool in `supabase_mcp_tools.py` provides SQL execution functions
- Use `execute_sql_on_supabase()` to run the SQL script
- Set audit bypass for bulk loading: `set_audit_bypass(project_id, True)`

### Database Population

- User authentication via service role in `service_role_helper.py`
- Transaction support through `execute_transaction()` function
- The script contains data for:
  - 11+ cryoprotectant molecules with detailed properties
  - 8+ scientifically accurate cryoprotectant mixtures
  - 5+ detailed experiments with protocols and results
  - 20+ property types and calculation methods

## Environment Setup

Ensure these environment variables are configured:

```
SUPABASE_URL=your_supabase_url
SUPABASE_KEY=your_service_role_key
SUPABASE_USER=admin_email
SUPABASE_PASSWORD=admin_password
SUPABASE_PROJECT_ID=your_project_id
```

## Verification

After completing both tasks:

1. Query tables in Supabase to confirm data exists
2. Verify RLS policies are active
3. Test application connectivity to Supabase
4. Ensure CI/CD pipeline is correctly configured with environment variables

## Next Steps After Completion

Once the database is populated and secured with RLS:

1. Run end-to-end tests with the populated database
2. Configure environment variables in GitHub Actions
3. Deploy to staging environment
4. Verify application functionality
5. Prepare for production deployment

## Success Criteria

1. RLS policies implemented for all tables and views
2. Database populated with comprehensive scientific data
3. No errors in database verification report
4. Application connects properly to populated database