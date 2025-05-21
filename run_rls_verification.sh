#!/bin/bash
# Script to run RLS verification tests

echo "Running RLS verification tests..."

# Set up Python environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cryoprotect || {
    echo "Error: Failed to activate conda environment 'cryoprotect'. Make sure it exists."
    exit 1
}

# Run the tests
echo "Running tests..."
python -m unittest tests/test_rls_policies_verification.py

# Check the exit status
if [ $? -eq 0 ]; then
    echo "Tests completed successfully."
else
    echo "Some tests failed. Check the log for details."
fi

# Generate report
echo "Generating verification report..."
timestamp=$(date "+%Y%m%d_%H%M%S")
report_file="reports/rls_verification_report_${timestamp}.md"

echo "# RLS Verification Report" > "$report_file"
echo "## Generated on $(date)" >> "$report_file"
echo "## Test Results" >> "$report_file"

# Append test log if it exists
if [ -f "rls_verification_tests.log" ]; then
    echo '```' >> "$report_file"
    cat rls_verification_tests.log >> "$report_file"
    echo '```' >> "$report_file"
    echo "Log file appended to report."
else
    echo "No log file found."
fi

echo "Report saved to: $report_file"
echo "RLS verification completed."