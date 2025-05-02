#!/bin/bash

echo
echo "============================================================"
echo "CryoProtect Analyzer - Apply Missing RLS Policies"
echo "============================================================"
echo "This script will apply RLS policies to missing tables/views:"
echo "1. experiment_with_results"
echo "2. migrations"
echo "3. mixture_with_components"
echo "4. molecule_with_properties"
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

# Check for required packages
if ! python3 -c "import psycopg2" 2>/dev/null; then
    echo "Installing required packages..."
    pip3 install psycopg2-binary python-dotenv
fi

# Make script executable
chmod +x apply_missing_rls_policies.py

echo "Running RLS policy application..."
python3 apply_missing_rls_policies.py

if [ $? -ne 0 ]; then
    echo
    echo "ERROR: RLS policy application failed. Check the logs for details."
    exit 1
else
    echo
    echo "SUCCESS: RLS policies have been applied successfully."
    echo
fi

exit 0