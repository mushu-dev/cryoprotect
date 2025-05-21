#!/bin/bash

# Run RLS Policy Verification Tests
# This script runs the RLS policy verification tests to check the optimized RLS implementation

# Set up environment
echo "Setting up environment variables..."
if [ -f .env ]; then
    echo "Loading .env file"
    export $(grep -v '^#' .env | xargs)
else
    echo "No .env file found, please make sure database connection parameters are set in the environment"
fi

# Check if required variables exist
if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_NAME" ] || [ -z "$SUPABASE_DB_USER" ] || [ -z "$SUPABASE_DB_PASSWORD" ]; then
    echo "Error: Required database connection variables are not set"
    echo "Please set SUPABASE_DB_HOST, SUPABASE_DB_NAME, SUPABASE_DB_USER, and SUPABASE_DB_PASSWORD"
    exit 1
fi

# Create Python virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python -m venv venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Install required packages
echo "Installing required packages..."
pip install -r requirements.txt
pip install python-dotenv psycopg2-binary

# Run the tests
echo "Running RLS policy verification tests..."
python -m tests.test_rls_policies_verification "$@"

# Generate report
echo "Generating test report..."
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
REPORT_FILE="reports/rls_verification_${TIMESTAMP}.md"

mkdir -p reports

# Create a markdown report
cat > "$REPORT_FILE" << EOF
# RLS Policy Verification Report

Date: $(date)
Database: $SUPABASE_DB_NAME
Host: $SUPABASE_DB_HOST

## Test Results

The following RLS policy verification tests were performed:

1. Security Definer Functions (TC-RLS-001 to TC-RLS-007)
2. Table Access Policies (TC-RLS-010 to TC-RLS-017)
3. Relationship Policies (TC-RLS-020 to TC-RLS-025)
4. Service Role Access (TC-RLS-030 to TC-RLS-033)

See the log file for detailed results: rls_verification_tests.log

## Summary

The RLS Policy verification tests have been completed. Any failures indicate issues with the RLS implementation that need to be addressed.

### Recommended Next Steps

1. Review all test failures and fix the corresponding RLS policies
2. Test RLS policies with actual API endpoints to ensure proper integration
3. Verify performance with the optimized RLS implementation using real-world queries

## Test Log

\`\`\`
$(cat rls_verification_tests.log 2>/dev/null || echo "Log file not found")
\`\`\`
EOF

echo "Test report generated: $REPORT_FILE"
echo "RLS policy verification tests completed!"

# Deactivate virtual environment
deactivate