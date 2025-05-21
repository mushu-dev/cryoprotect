#!/bin/bash
# CryoProtect v2 - Test Runner Script
# This script runs all the tests and generates a comprehensive test report.

echo "CryoProtect v2 - Running Tests"
echo "============================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH"
    exit 1
fi

# Check if the virtual environment exists
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv .venv
    
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create virtual environment"
        exit 1
    fi
fi

# Activate the virtual environment
if [ -f ".venv/bin/activate" ]; then
    echo "Activating virtual environment..."
    source .venv/bin/activate
elif [ -f ".venv/Scripts/activate" ]; then
    echo "Activating virtual environment..."
    source .venv/Scripts/activate
else
    echo "Error: Could not find activation script for virtual environment"
    exit 1
fi

# Install required packages
echo "Installing required packages..."
pip install -r requirements.txt

if [ $? -ne 0 ]; then
    echo "Error: Failed to install required packages"
    exit 1
fi

# Run the test runner
echo "Running tests..."
python run_all_tests.py

# Check the exit code
if [ $? -eq 0 ]; then
    echo "All tests passed!"
    echo "See TEST_RESULTS_REPORT_*.md for detailed results"
else
    echo "Some tests failed. See TEST_RESULTS_REPORT_*.md for detailed results"
fi

# Deactivate the virtual environment
deactivate

echo
echo "Testing complete"