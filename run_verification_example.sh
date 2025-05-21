#!/bin/bash

# Run Database Verification Example
# This script runs the database verification example

# Set up environment
echo "Setting up environment variables..."
if [ -f .env ]; then
    echo "Loading .env file"
    export $(grep -v '^#' .env | xargs)
else
    echo "No .env file found, please make sure database connection parameters are set in the environment"
    exit 1
fi

# Create Python virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python -m venv venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Install required packages
echo "Installing required packages..."
pip install -r requirements.txt
pip install python-dotenv psycopg2-binary

# Run the example
echo "Running database verification example..."
python examples/run_database_verification_example.py "$@"

# Deactivate virtual environment
deactivate