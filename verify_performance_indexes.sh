#!/bin/bash

echo "==============================================================="
echo "CryoProtect v2 - Verify Performance Indexes"
echo "==============================================================="
echo

python3 verify_performance_indexes.py

if [ $? -ne 0 ]; then
    echo
    echo "Error: Performance indexes verification failed."
    echo "Please check the log file for details."
    exit 1
fi

echo
echo "Performance indexes verification completed successfully!"
echo