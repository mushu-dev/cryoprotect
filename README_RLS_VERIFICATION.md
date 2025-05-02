# CryoProtect v2 Database RLS Verification

This directory contains scripts and templates for verifying the Row Level Security (RLS) policies, access controls, query performance, and data relationships in the CryoProtect v2 Supabase database.

## Overview

The verification process tests:

1. **RLS Effectiveness**: Verifies that RLS policies correctly restrict access to data based on user roles and ownership.
2. **Access Control Validation**: Tests access patterns with different user roles (anonymous, authenticated, service role).
3. **Performance Impact Assessment**: Measures the performance impact of RLS on query execution time.
4. **Scientific Data Relationship Verification**: Validates the integrity of scientific data relationships.

## Files

- `verify_rls_policies.py`: Main verification script that tests RLS policies via the Supabase REST API.
- `execute_rls_verification_via_mcp.py`: Contains functions for testing RLS policies via direct SQL execution.
- `execute_rls_sql_via_mcp.py`: Helper script for executing SQL queries via the MCP server.
- `run_rls_verification.py`: Entry point for running the verification tests using the previous scripts.
- `run_mcp_verification.py`: Alternative entry point that uses the MCP server directly.
- `RLS_Verification_Report_Template.md`: Template for the comprehensive verification report.

## Prerequisites

- Python 3.7+
- Access to the CryoProtect v2 Supabase project
- Supabase MCP server connection

## Configuration

The scripts are pre-configured with the following settings:

- **Project ID**: tsdlmynydfuypiugmkev
- **Project URL**: https://tsdlmynydfuypiugmkev.supabase.co
- **Anon Key**: eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InRzZGxteW55ZGZ1eXBpdWdta2V2Iiwicm9sZSI6ImFub24iLCJpYXQiOjE3NDQ3ODAxODMsImV4cCI6MjA2MDM1NjE4M30.T6udmD-3lhlTy4bY2p0y2lX5-11yvyn425PWrlnIPLU

## Usage

### Option 1: Using the MCP Server

This is the recommended approach as it allows direct SQL execution with different user roles.

1. Run the verification script:

```bash
python run_mcp_verification.py
```

This will:
- Execute SQL queries via the MCP server
- Test RLS policies with different user roles
- Measure query performance
- Validate data relationships
- Generate a comprehensive verification report

### Option 2: Using the REST API

This approach uses the Supabase REST API to test RLS policies.

1. Run the verification script:

```bash
python run_rls_verification.py
```

## Output

The verification process generates the following output files:

- `RLS_Verification_Results.json`: Raw results of the verification tests in JSON format.
- `CryoProtect_RLS_Verification_Report.md`: Comprehensive verification report in Markdown format.

## Verification Tests

### 1. RLS Effectiveness Tests

Tests the effectiveness of RLS policies by:
- Attempting to access private data as an anonymous user
- Attempting to modify data without proper permissions
- Verifying that RLS policies are correctly applied

### 2. Access Control Tests

Tests access patterns with different user roles:
- Anonymous user: Should only be able to access public data
- Authenticated user: Should be able to access own data and public data
- Service role: Should be able to access all data

### 3. Performance Impact Tests

Measures the performance impact of RLS by:
- Executing queries with RLS enabled
- Executing the same queries with RLS disabled
- Calculating the performance difference

### 4. Data Relationship Tests

Validates scientific data relationships by:
- Verifying referential integrity between tables
- Checking that foreign keys reference valid records
- Ensuring that RLS doesn't interfere with proper data access patterns

## Customization

You can customize the verification process by:

1. Modifying the SQL queries in the scripts
2. Adding additional tests for specific tables or access patterns
3. Adjusting the performance measurement parameters
4. Customizing the report template

## Troubleshooting

If you encounter issues:

1. Verify that the Supabase project is accessible
2. Check that the MCP server is properly connected
3. Ensure that you have the necessary permissions to execute SQL queries
4. Check the console output for error messages

## Security Considerations

- The scripts include the project's anonymous key, which is safe to share as it only provides public access.
- For service role testing, you would need to provide a service role key, which should be kept secure.
- The scripts do not modify any data by default, but be cautious when enabling modification tests.

## Contributing

To contribute to the verification process:

1. Add additional tests for specific tables or access patterns
2. Improve the report template with more detailed analysis
3. Enhance the performance measurement methodology
4. Add support for additional user roles or access patterns