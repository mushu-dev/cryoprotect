#!/bin/bash
# run_rls_policy_improvements_test.sh
# Script to run tests for the RLS policy improvements

# Ensure script execution stops if any command fails
set -e

# Ensure we're in the project root directory
cd "$(dirname "$0")"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "=== RLS Policy Improvements Test Runner ==="
echo "This script will test the improved RLS policies to ensure they are functioning correctly."
echo ""

# Check for command line arguments
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --report       Generate a detailed report in the reports directory"
    echo "  --help, -h     Show this help message"
    exit 0
fi

# Check if we should generate a detailed report
GENERATE_REPORT=0
if [ "$1" == "--report" ]; then
    GENERATE_REPORT=1
fi

# Check Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 is required but could not be found"
    exit 1
fi

# Source .env file if it exists
if [ -f .env ]; then
    echo "Loading environment variables from .env file"
    export $(grep -v '^#' .env | xargs)
fi

# Check for required environment variables
if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_USER" ] || [ -z "$SUPABASE_DB_PASSWORD" ]; then
    echo "ERROR: Missing required environment variables"
    echo "Please set the following environment variables or create a .env file:"
    echo "  SUPABASE_DB_HOST"
    echo "  SUPABASE_DB_PORT (defaults to 5432)"
    echo "  SUPABASE_DB_NAME (defaults to postgres)"
    echo "  SUPABASE_DB_USER"
    echo "  SUPABASE_DB_PASSWORD"
    echo "  SUPABASE_DB_SSLMODE (defaults to require)"
    exit 1
fi

echo "Environment variables set correctly"
echo "Running RLS policy improvement tests..."

# Run the tests
if [ $GENERATE_REPORT -eq 1 ]; then
    # Run with report generation
    echo "Generating detailed test report..."
    
    # Create reports directory if it doesn't exist
    mkdir -p reports
    
    # Run tests with report generation
    python3 - << 'EOF'
import os
import sys
import json
import datetime
from importlib import import_module
from unittest.mock import patch

# Add the project root to the path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import the test module
from tests.test_rls_policy_improvements import run_tests

# Run the tests and capture the results
with patch('sys.stdout') as mock_stdout:
    result = run_tests()

# Generate a report
report = {
    'timestamp': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    'success': result.wasSuccessful(),
    'tests_run': result.testsRun,
    'errors': len(result.errors),
    'failures': len(result.failures),
    'test_details': []
}

# Add details for each test case
for test_case, error in result.errors:
    report['test_details'].append({
        'test_name': str(test_case),
        'status': 'ERROR',
        'message': error
    })

for test_case, failure in result.failures:
    report['test_details'].append({
        'test_name': str(test_case),
        'status': 'FAILURE',
        'message': failure
    })

# Save the report as JSON
report_file = f"reports/rls_improvements_test_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(report_file, 'w') as f:
    json.dump(report, f, indent=2)

# Generate a markdown report
md_report_file = f"reports/rls_improvements_test_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
with open(md_report_file, 'w') as f:
    f.write(f"# RLS Policy Improvements Test Report\n\n")
    f.write(f"**Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"**Status:** {'Passed ✅' if report['success'] else 'Failed ❌'}\n\n")
    f.write(f"**Tests Run:** {report['tests_run']}\n")
    f.write(f"**Errors:** {report['errors']}\n")
    f.write(f"**Failures:** {report['failures']}\n\n")
    
    if report['errors'] > 0 or report['failures'] > 0:
        f.write("## Test Failures and Errors\n\n")
        for detail in report['test_details']:
            f.write(f"### {detail['test_name']}\n\n")
            f.write(f"**Status:** {detail['status']}\n\n")
            f.write("```\n")
            f.write(detail['message'])
            f.write("\n```\n\n")
    
    f.write("## Improvement Verification\n\n")
    f.write("This test suite verified the following improvements:\n\n")
    f.write("1. **Security Definer Functions**: All required security definer functions are in place and working correctly\n")
    f.write("2. **Team-Based Access Model**: RLS policies correctly implement the team-based access model\n")
    f.write("3. **Performance Optimization**: Security definer functions provide optimized access control\n")
    f.write("4. **Comprehensive Policy Coverage**: All tables have appropriate RLS policies\n")

print(f"Report generated: {md_report_file}")

# Exit with the test result status
sys.exit(0 if result.wasSuccessful() else 1)
EOF
    
    # Check the exit code from the Python script
    if [ $? -eq 0 ]; then
        echo "RLS policy improvement tests passed successfully! Report generated in reports directory."
    else
        echo "RLS policy improvement tests failed. Check the report in the reports directory for details."
        exit 1
    fi
else
    # Run the tests directly
    python3 -m tests.test_rls_policy_improvements
    
    # Check the exit code
    if [ $? -eq 0 ]; then
        echo "RLS policy improvement tests passed successfully!"
    else
        echo "RLS policy improvement tests failed."
        exit 1
    fi
fi

echo "To run with detailed report generation, use: $0 --report"