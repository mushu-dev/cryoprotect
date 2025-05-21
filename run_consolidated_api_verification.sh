#!/bin/bash
# Run the consolidated API verification script

# Make script executable
chmod +x ./verify_consolidated_api.py

# Run the verification script
python verify_consolidated_api.py

# Check exit code
if [ $? -eq 0 ]; then
    echo "Consolidated API verification completed successfully!"
else
    echo "Consolidated API verification failed!"
    exit 1
fi