#!/bin/bash

echo "==================================================="
echo "CryoProtect v2 Database Remediation Verification"
echo "==================================================="
echo

# Activate virtual environment if it exists
if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
elif [ -f "venv/bin/activate" ]; then
    source venv/bin/activate
else
    echo "WARNING: Virtual environment not found. Using system Python."
fi

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Python not found. Please install Python 3.8 or higher."
    exit 1
fi

# Check if required packages are installed
echo "Checking required packages..."
python3 -c "import supabase" &> /dev/null
if [ $? -ne 0 ]; then
    echo "Installing supabase package..."
    pip install supabase
fi

python3 -c "import requests" &> /dev/null
if [ $? -ne 0 ]; then
    echo "Installing requests package..."
    pip install requests
fi

python3 -c "import dotenv" &> /dev/null
if [ $? -ne 0 ]; then
    echo "Installing python-dotenv package..."
    pip install python-dotenv
fi

# Create reports directory if it doesn't exist
mkdir -p reports

echo
echo "Running database verification tests..."
echo

# Run the verification script
python3 verify_database_remediation.py --verbose

echo
echo "Verification complete. Check the reports directory for detailed results."
echo