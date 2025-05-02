#!/bin/bash

echo "==============================================================="
echo "CryoProtect v2 - Apply Critical Performance Improvements"
echo "==============================================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in the PATH."
    echo "Please install Python 3 and try again."
    exit 1
fi

# Check if the script exists
if [ ! -f "apply_performance_improvements.py" ]; then
    echo "Error: apply_performance_improvements.py not found."
    echo "Please make sure the script is in the current directory."
    exit 1
fi

echo "Applying critical performance improvements..."
echo

# Make the script executable
chmod +x apply_performance_improvements.py

# Run the script with backup option
python3 apply_performance_improvements.py --backup

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo
    echo "Performance improvements applied successfully!"
    exit 0
elif [ $EXIT_CODE -eq 1 ]; then
    echo
    echo "Performance improvements applied with warnings."
    echo "Please check the log file for details."
    exit 1
else
    echo
    echo "Error: Failed to apply performance improvements."
    echo "Please check the log file for details."
    exit 2
fi