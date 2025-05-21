#!/bin/bash

echo "==============================================================="
echo "CryoProtect v2 - Test Database Remediation"
echo "==============================================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python is not installed or not in the PATH."
    echo "Please install Python and try again."
    exit 1
fi

echo "This script will test the database remediation process in a safe environment."
echo "It will:"
echo "1. Create a test schema with sample data that mimics the issues in production"
echo "2. Run the remediation process on this test schema"
echo "3. Verify that the remediation was successful"
echo "4. Clean up the test schema"
echo
echo "This is a safe way to test the remediation process without affecting your production data."
echo

read -p "Do you want to continue? (y/n): " CONFIRM
if [[ $CONFIRM != "y" && $CONFIRM != "Y" ]]; then
    echo "Operation cancelled."
    exit 0
fi

echo
echo "Select an option:"
echo "1. Run test and clean up afterward"
echo "2. Run test and keep the test schema for inspection"
echo

read -p "Enter option (1-2): " OPTION

if [[ $OPTION == "1" ]]; then
    echo
    echo "Running test and cleaning up afterward..."
    python3 test_database_remediation.py
elif [[ $OPTION == "2" ]]; then
    echo
    echo "Running test and keeping the test schema..."
    python3 test_database_remediation.py --keep-schema
else
    echo "Invalid option."
    exit 1
fi

echo
echo "Test process completed."
echo "Check the log file for details."
echo

read -p "Press Enter to exit..."