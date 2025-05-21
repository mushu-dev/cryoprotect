#!/bin/bash

echo "==============================================================="
echo "CryoProtect v2 - Apply Performance Indexes Migration"
echo "==============================================================="
echo

# Check if Node.js is installed
if ! command -v node &> /dev/null; then
    echo "Error: Node.js is not installed or not in the PATH."
    echo "Please install Node.js and try again."
    exit 1
fi

echo "Applying performance indexes migration..."
echo

# Run the migration script
node migrations/apply_performance_indexes_migration.js

if [ $? -ne 0 ]; then
    echo
    echo "Error: Migration failed."
    echo "Please check the error messages above."
    exit 1
fi

echo
echo "Migration completed successfully!"
echo