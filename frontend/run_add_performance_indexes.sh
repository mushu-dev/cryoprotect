#!/bin/bash
# Script to add performance-optimizing indexes to the database

set -e

echo "Starting performance index creation process..."

# Check if we're requesting dry-run mode
if [ "$1" == "--dry-run" ]; then
    echo "Running in dry-run mode (only checking existing indexes)..."
    python3 apply_performance_indexes.py --dry-run
    exit $?
fi

# Check if we're requesting simulation mode
if [ "$1" == "--simulate" ]; then
    echo "Running with simulation of query improvements..."
    python3 apply_performance_indexes.py --simulate
    exit $?
fi

# Run the full index creation
echo "Running full index creation process..."
python3 apply_performance_indexes.py

# Check the exit status
if [ $? -eq 0 ]; then
    echo "Performance index creation completed successfully!"
    exit 0
else
    echo "Performance index creation failed!"
    exit 1
fi