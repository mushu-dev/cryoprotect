#!/bin/bash
# CryoProtect Analyzer API - Standardization Runner (Unix/Linux/macOS)
#
# This script runs the API audit and applies standardization to all API endpoints.
# It provides a convenient way to standardize the API in one step.

echo "Running API standardization process..."

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Python 3 is not installed or not in PATH. Please install Python 3 and try again."
    exit 1
fi

# Parse command line arguments
DRY_RUN=""
VERBOSE=""
SKIP_AUDIT=""
SKIP_APPLY=""

for arg in "$@"; do
    case $arg in
        --dry-run)
            DRY_RUN="--dry-run"
            ;;
        --verbose)
            VERBOSE="--verbose"
            ;;
        --skip-audit)
            SKIP_AUDIT="--skip-audit"
            ;;
        --skip-apply)
            SKIP_APPLY="--skip-apply"
            ;;
    esac
done

# Make script executable if it's not already
chmod +x api_audit.py
chmod +x apply_api_standardization.py
chmod +x run_api_standardization.py

# Run the standardization process
python3 run_api_standardization.py $DRY_RUN $VERBOSE $SKIP_AUDIT $SKIP_APPLY

if [ $? -ne 0 ]; then
    echo "API standardization process failed."
    exit 1
fi

echo "API standardization process completed successfully."
exit 0