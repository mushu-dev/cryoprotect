#!/bin/bash

echo "Running CryoProtect v2 API Integration Tests..."
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Python 3 is not installed or not in the PATH."
    echo "Please install Python 3 and try again."
    exit 1
fi

# Check if the virtual environment exists
if [ ! -d ".venv" ]; then
    echo "Virtual environment not found."
    echo "Creating virtual environment..."
    python3 -m venv .venv
    if [ $? -ne 0 ]; then
        echo "Failed to create virtual environment."
        exit 1
    fi
fi

# Activate the virtual environment
source .venv/bin/activate

# Install required packages
echo "Installing required packages..."
pip install requests
if [ $? -ne 0 ]; then
    echo "Failed to install required packages."
    exit 1
fi

# Run the test script
echo
echo "Starting API integration tests..."
echo
python3 test_api_integration.py
TEST_RESULT=$?

# Deactivate the virtual environment
deactivate

# Return the test result
exit $TEST_RESULT