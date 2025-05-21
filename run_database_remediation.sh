#!/bin/bash

echo "==============================================================="
echo "CryoProtect v2 - Database Remediation"
echo "==============================================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python is not installed or not in the PATH."
    echo "Please install Python and try again."
    exit 1
fi

echo "This script will run the complete database remediation process."
echo "It will address the following critical issues:"
echo
echo "1. SECURITY: Enable Row Level Security (RLS)"
echo "2. STRUCTURE: Standardize Schema & Fix Relationships"
echo "3. PERFORMANCE: Add Missing Indexes"
echo "4. ROLES: Create Application-Specific Roles"
echo "5. DATA: Consolidate Duplicate Tables"
echo
echo "WARNING: This will make significant changes to your database."
echo "Make sure you have a backup before proceeding."
echo

read -p "Do you want to continue? (y/n): " CONFIRM
if [[ $CONFIRM != "y" && $CONFIRM != "Y" ]]; then
    echo "Operation cancelled."
    exit 0
fi

echo
echo "Select an option:"
echo "1. Run all phases"
echo "2. Run in dry-run mode (no changes will be made)"
echo "3. Run a specific phase"
echo

read -p "Enter option (1-3): " OPTION

if [[ $OPTION == "1" ]]; then
    echo
    echo "Running all phases..."
    python3 complete_database_remediation.py
elif [[ $OPTION == "2" ]]; then
    echo
    echo "Running in dry-run mode..."
    python3 complete_database_remediation.py --dry-run
elif [[ $OPTION == "3" ]]; then
    echo
    echo "Select a phase to run:"
    echo "1. SECURITY: Enable Row Level Security (RLS)"
    echo "2. STRUCTURE: Standardize Schema & Fix Relationships"
    echo "3. PERFORMANCE: Add Missing Indexes"
    echo "4. ROLES: Create Application-Specific Roles"
    echo "5. DATA: Consolidate Duplicate Tables"
    echo
    
    read -p "Enter phase (1-5): " PHASE
    
    if [[ $PHASE -ge 1 && $PHASE -le 5 ]]; then
        echo
        echo "Running phase $PHASE..."
        python3 complete_database_remediation.py --phase $PHASE
    else
        echo "Invalid phase number."
        exit 1
    fi
else
    echo "Invalid option."
    exit 1
fi

echo
echo "Database remediation process completed."
echo "Check the log file for details."
echo

read -p "Press Enter to exit..."