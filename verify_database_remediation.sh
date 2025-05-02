#!/bin/bash

echo "==============================================================="
echo "CryoProtect v2 - Verify Database Remediation"
echo "==============================================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python is not installed or not in the PATH."
    echo "Please install Python and try again."
    exit 1
fi

echo "This script will verify that the database remediation was successful."
echo "It will check:"
echo
echo "1. RLS is enabled on all tables"
echo "2. Anonymous access is properly restricted"
echo "3. All tables have appropriate RLS policies"
echo "4. All foreign keys have corresponding indexes"
echo "5. Application roles are created correctly"
echo "6. Table names are standardized to plural form"
echo "7. Junction tables are created to fix fan traps"
echo

read -p "Do you want to continue? (y/n): " CONFIRM
if [[ $CONFIRM != "y" && $CONFIRM != "Y" ]]; then
    echo "Operation cancelled."
    exit 0
fi

echo
echo "Running verification..."
python3 verify_database_remediation.py

echo
echo "Verification process completed."
echo "Check the log file for details."
echo

read -p "Press Enter to exit..."