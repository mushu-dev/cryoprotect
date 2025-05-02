#!/bin/bash
# Script to implement enhanced RLS policies for CryoProtect v2

# Ensure we're in the project root directory
cd "$(dirname "$0")"

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Python 3 is not installed. Please install Python 3 and try again."
    exit 1
fi

# Check if virtual environment exists
if [ -d ".venv" ]; then
    echo "Activating virtual environment..."
    source .venv/bin/activate
fi

# Check if .env file exists
if [ ! -f ".env" ]; then
    echo "Warning: .env file not found. Using default configuration."
    if [ -f ".env.template" ]; then
        echo "Consider copying .env.template to .env and updating the values."
    fi
fi

# Create logs directory if it doesn't exist
mkdir -p logs
mkdir -p reports/security

# Run the Python script
echo "Implementing enhanced RLS policies..."
python3 implement_enhanced_rls_policies.py "$@"

# Check if the script executed successfully
if [ $? -eq 0 ]; then
    echo "✅ Enhanced RLS policies implemented successfully!"
    echo "Check the reports/security directory for verification reports."
else
    echo "❌ Error implementing enhanced RLS policies. Check logs for details."
    exit 1
fi