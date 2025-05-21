#!/bin/bash

echo "==================================================="
echo "CryoProtect Database Schema Standardization"
echo "==================================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Python 3 is not installed or not in PATH."
    echo "Please install Python 3.6+ and try again."
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2)
echo "Detected Python version: $PYTHON_VERSION"
echo

# Create a virtual environment if it doesn't exist
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv .venv
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to create virtual environment."
        exit 1
    fi
fi

# Activate the virtual environment
echo "Activating virtual environment..."
source .venv/bin/activate
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate virtual environment."
    exit 1
fi

# Install required packages
echo "Installing required packages..."
pip install -q requests
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to install required packages."
    exit 1
fi

# Run the test script first
echo
echo "==================================================="
echo "Testing Supabase connection..."
echo "==================================================="
python3 test_supabase_connection.py
if [ $? -ne 0 ]; then
    echo
    echo "Connection test failed. Please fix the issues before continuing."
    echo
    read -p "Do you want to continue anyway? (y/n): " CONTINUE
    if [[ ! "$CONTINUE" =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo
echo "==================================================="
echo "IMPORTANT: This script will standardize the database schema"
echo "by converting singular table names to plural and adding"
echo "proper constraints. This operation is potentially destructive."
echo "==================================================="
echo
echo "Please make sure you have:"
echo " 1. Backed up your database"
echo " 2. Reviewed the standardize_schema.py script"
echo " 3. Tested your connection with test_supabase_connection.py"
echo
read -p "Are you sure you want to continue? (y/n): " CONFIRM

if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
    echo "Operation cancelled by user."
    exit 1
fi

echo
echo "==================================================="
echo "Running schema standardization..."
echo "==================================================="
python3 standardize_schema.py

if [ $? -eq 0 ]; then
    echo
    echo "==================================================="
    echo "Schema standardization completed successfully!"
    echo "==================================================="
else
    echo
    echo "==================================================="
    echo "Schema standardization completed with errors."
    echo "Please check the log file for details."
    echo "==================================================="
fi

# Deactivate virtual environment
deactivate

echo
echo "Press Enter to exit..."
read