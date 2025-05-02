#!/bin/bash
# CryoProtect v2 - Database Migration Runner Script
# This script provides a convenient way to run the migration script with different options

echo "CryoProtect v2 - Database Migration to Plural Tables"
echo "===================================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH"
    exit 1
fi

# Check if the migration script exists
if [ ! -f "migrate_to_plural_tables.py" ]; then
    echo "Error: migrate_to_plural_tables.py not found"
    exit 1
fi

# Parse command line arguments
DRY_RUN=""
ROLLBACK=""
REPORT=""
TEST=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN="--dry-run"
            shift
            ;;
        --rollback)
            ROLLBACK="--rollback"
            shift
            ;;
        --report)
            REPORT="--report $2"
            shift 2
            ;;
        --test)
            TEST=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: ./run_migration.sh [--dry-run] [--rollback --report REPORT_FILE] [--test]"
            exit 1
            ;;
    esac
done

# Run the test script if requested
if [ -n "$TEST" ]; then
    echo "Running migration test script..."
    echo
    python3 test_migration_script.py
    if [ $? -ne 0 ]; then
        echo
        echo "Migration test failed."
        exit 1
    fi
    echo
    echo "Migration test completed successfully."
    exit 0
fi

# Run the migration script with the specified options
if [ -n "$ROLLBACK" ]; then
    if [ -z "$REPORT" ]; then
        echo "Error: --report is required with --rollback"
        echo "Usage: ./run_migration.sh --rollback --report migration_report_YYYYMMDD_HHMMSS.json"
        exit 1
    fi
    echo "Running migration rollback..."
    echo
    python3 migrate_to_plural_tables.py --rollback $REPORT
elif [ -n "$DRY_RUN" ]; then
    echo "Running migration in dry run mode..."
    echo
    python3 migrate_to_plural_tables.py --dry-run
else
    echo
    echo "WARNING: You are about to run the actual migration."
    echo "This will modify the database structure and should be done during a maintenance window."
    echo
    read -p "Are you sure you want to continue? (y/n): " CONFIRM
    if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
        echo "Migration cancelled."
        exit 0
    fi
    
    echo
    echo "Running migration..."
    echo
    python3 migrate_to_plural_tables.py
fi

if [ $? -ne 0 ]; then
    echo
    echo "Migration failed. Check the log file for details."
    exit 1
fi

echo
echo "Migration completed successfully."
exit 0