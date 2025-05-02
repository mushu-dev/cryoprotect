#!/bin/bash

echo "==================================================="
echo "CryoProtect v2 - Service Role RLS Fix Application"
echo "==================================================="
echo ""
echo "This script will apply the service role RLS fix to your Supabase project."
echo ""

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python is not installed or not in PATH."
    echo "Please install Python and try again."
    exit 1
fi

# Check if the migration file exists
if [ ! -f "migrations/007_service_role_rls.sql" ]; then
    echo "Error: Migration file not found."
    echo "Please make sure migrations/007_service_role_rls.sql exists."
    exit 1
fi

echo "Step 1: Applying RLS policies..."
echo ""
echo "Please apply the RLS policies manually through the Supabase dashboard:"
echo "1. Log in to the Supabase Dashboard (https://app.supabase.com/)"
echo "2. Select your project"
echo "3. Go to the SQL Editor"
echo "4. Copy the contents of apply_service_role_rls_manual.sql"
echo "5. Paste into the SQL Editor and run the query"
echo ""
echo "Press Enter when you have completed this step..."
read -p ""

echo ""
echo "Step 2: Testing the populate_molecules.py script..."
echo ""
python3 populate_molecules.py
if [ $? -ne 0 ]; then
    echo ""
    echo "Error: The populate_molecules.py script failed."
    echo "Please check the error message above and try again."
    exit 1
fi

echo ""
echo "==================================================="
echo "Service Role RLS Fix Applied Successfully!"
echo "==================================================="
echo ""
echo "The service role RLS fix has been applied to your Supabase project."
echo "You can now run scripts like populate_molecules.py without RLS policy violations."
echo ""
echo "For more information, see README_SERVICE_ROLE_RLS_FIX.md."
echo ""
read -p "Press Enter to exit..."