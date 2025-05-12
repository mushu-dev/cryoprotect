#!/bin/bash
# Script to apply consolidated molecule triggers and test them

set -e

echo "Applying consolidated molecule triggers migration..."
PGPASSWORD=$DB_PASSWORD psql -h $DB_HOST -U $DB_USER -d $DB_NAME -f migrations/023_create_consolidated_molecule_triggers.sql

echo ""
echo "Running trigger tests..."
python test_consolidated_molecule_triggers.py

echo ""
echo "Done!"