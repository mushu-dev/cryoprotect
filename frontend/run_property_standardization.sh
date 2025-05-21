#!/bin/bash
# Script to run property format standardization

set -e

echo "Starting property format standardization process..."

# Check if we're requesting analyze-only mode
if [ "$1" == "--analyze-only" ]; then
    echo "Running in analyze-only mode..."
    python3 standardize_property_formats.py --analyze-only
    exit $?
fi

# Check if we're requesting dry-run mode
if [ "$1" == "--dry-run" ]; then
    echo "Running in dry-run mode..."
    python3 standardize_property_formats.py --dry-run
    exit $?
fi

# Run the full standardization
echo "Running full standardization process..."
python3 standardize_property_formats.py

# Check the exit status
if [ $? -eq 0 ]; then
    echo "Property standardization completed successfully!"
    exit 0
else
    echo "Property standardization failed!"
    exit 1
fi