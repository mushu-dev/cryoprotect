#!/bin/bash

echo "==============================================================="
echo "CryoProtect v2 - Create Production Database Backup"
echo "==============================================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in the PATH."
    echo "Please install Python 3 and try again."
    exit 1
fi

# Check if the script exists
if [ ! -f "create_production_backup.py" ]; then
    echo "Error: create_production_backup.py not found."
    echo "Please make sure the script is in the current directory."
    exit 1
fi

echo "Creating production database backup..."
echo "This may take several minutes depending on the database size."
echo

# Make the script executable
chmod +x create_production_backup.py

# Run the backup script
python3 create_production_backup.py --format both

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo
    echo "Production database backup created successfully!"
    echo "Please check the production_backups directory for the backup files."
    exit 0
elif [ $EXIT_CODE -eq 1 ]; then
    echo
    echo "Production database backup created with warnings."
    echo "Please check the log file for details."
    exit 1
else
    echo
    echo "Error: Failed to create production database backup."
    echo "Please check the log file for details."
    exit 2
fi