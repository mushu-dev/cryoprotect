#!/bin/bash

echo
echo "============================================================"
echo "CryoProtect Analyzer - Apply Enhanced RLS Policies"
echo "============================================================"
echo "This script will apply improved RLS policies to missing tables/views:"
echo "1. experiment_with_results"
echo "2. migrations"
echo "3. mixture_with_components"
echo "4. molecule_with_properties"
echo
echo "Enhancements include:"
echo "- Transaction support for atomic operations"
echo "- Enhanced verification and effectiveness testing"
echo "- Performance benchmarking for RLS impact"
echo "- Comprehensive audit trail for scientific data"
echo "- Optimized RLS policy conditions for better performance"
echo
echo "Prerequisites:"
echo "- Python installed with psycopg2 package"
echo "- Supabase connection configured in .env file"
echo

# Check for Python
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Python is not installed or not in PATH."
    echo "Please install Python and try again."
    exit 1
fi

# Ensure backup before proceeding
echo "Creating database backup before applying RLS policies..."
python3 create_database_backup.py --description "Pre-RLS-enhancement backup"

if [ $? -ne 0 ]; then
    echo
    echo "WARNING: Could not create database backup. This is recommended before proceeding."
    echo "         Do you want to continue anyway? (Y/N)"
    read -p "> " CONTINUE
    if [[ ! "$CONTINUE" =~ ^[Yy]$ ]]; then
        echo "Operation cancelled by user."
        exit 0
    fi
fi

# Check for required packages
if ! python3 -c "import psycopg2" 2>/dev/null; then
    echo "Installing required packages..."
    pip3 install psycopg2-binary python-dotenv
fi

# Make script executable
chmod +x apply_missing_rls_policies_improved.py

echo "Running enhanced RLS policy application..."
python3 apply_missing_rls_policies_improved.py

if [ $? -ne 0 ]; then
    echo
    echo "ERROR: RLS policy application failed. Check the logs for details."
    exit 1
else
    echo
    echo "SUCCESS: Enhanced RLS policies have been applied successfully."
    echo "         A detailed report has been generated in the reports/security directory."
    echo
fi

# Optional: Open the report
echo "Do you want to view the latest report? (Y/N)"
read -p "> " OPEN_REPORT
if [[ "$OPEN_REPORT" =~ ^[Yy]$ ]]; then
    LATEST_REPORT=$(ls -t reports/security/rls_implementation_summary_*.md | head -n1)
    if command -v xdg-open &> /dev/null; then
        xdg-open "$LATEST_REPORT"
    elif command -v open &> /dev/null; then
        open "$LATEST_REPORT"
    else
        echo "Report is available at: $LATEST_REPORT"
        echo "Use your preferred text editor to open it."
    fi
fi

exit 0