#!/bin/bash

echo "==================================================="
echo "CryoProtect v2 - Database Performance Validation"
echo "==================================================="
echo
echo "This script will run performance tests to validate:"
echo "1. Query performance on tables with RLS enabled"
echo "2. Query performance on tables with foreign key relationships"
echo "3. Query performance on junction tables"
echo
echo "Press Ctrl+C to cancel or Enter to continue..."
read

echo
echo "Setting up environment..."
if [ -f "./setup_environment.sh" ]; then
    source ./setup_environment.sh
    if [ $? -ne 0 ]; then
        echo "Failed to set up environment. Please run setup_environment.sh manually first."
        exit 1
    fi
else
    echo "setup_environment.sh not found. Please ensure environment is set up correctly."
fi

echo
echo "Running performance validation tests..."
python test_database_performance_remediation.py
if [ $? -ne 0 ]; then
    echo "Performance validation failed."
    exit 1
fi

echo
echo "Performance validation completed successfully."
echo "Results are available in:"
echo "- database_performance_validation_report.txt"
echo "- database_performance_validation_report.json"
echo
echo "Press Enter to exit..."
read