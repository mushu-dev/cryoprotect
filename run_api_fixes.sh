#!/bin/bash

echo "Running API integration fixes for CryoProtect v2..."
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH."
    echo "Please install Python 3 and try again."
    exit 1
fi

# Make the script executable
chmod +x fix_api_integration.py

# Run the fix script
python3 fix_api_integration.py

if [ $? -ne 0 ]; then
    echo
    echo "Error: Failed to apply API fixes."
    echo "Please check the logs for details."
    exit 1
else
    echo
    echo "API integration fixes applied successfully!"
    echo "The API now works with the standardized database schema."
    echo
    echo "Changes made:"
    echo "1. Updated API resource classes to use plural table names"
    echo "2. Fixed endpoint duplication issues"
    echo "3. Implemented consistent error handling"
    echo "4. Added retry logic for API calls"
    echo
    echo "Backups of the original files were created with .bak extension."
fi

read -p "Press Enter to continue..."