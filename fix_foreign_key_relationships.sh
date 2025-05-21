#!/bin/bash

echo "CryoProtect v2 - Foreign Key Relationship Fixer"
echo "==============================================="
echo

# Activate the virtual environment if it exists
if [ -f "./.venv/bin/activate" ]; then
    source ./.venv/bin/activate
else
    echo "Warning: Virtual environment not found. Make sure dependencies are installed."
fi

echo
echo "Available options:"
echo "1. Run with dry-run (show what would be done without making changes)"
echo "2. Verify only (check for issues without making changes)"
echo "3. Apply all fixes"
echo "4. Rollback to previous state"
echo

read -p "Enter option (1-4): " option

case $option in
    1)
        echo
        echo "Running in dry-run mode..."
        python fix_foreign_key_relationships.py --dry-run
        ;;
    2)
        echo
        echo "Verifying foreign key constraints..."
        python fix_foreign_key_relationships.py --verify-only
        ;;
    3)
        echo
        echo "Applying all foreign key fixes..."
        python fix_foreign_key_relationships.py
        ;;
    4)
        echo
        echo "Rolling back to previous state..."
        python fix_foreign_key_relationships.py --rollback
        ;;
    *)
        echo
        echo "Invalid option. Please run the script again and select a valid option."
        exit 1
        ;;
esac

echo
echo "Script execution completed."
read -p "Press Enter to continue..."