#!/bin/bash
# Script to apply consolidated molecule constraints and indexes

set -e

echo "Applying consolidated molecule constraints and indexes..."
PGPASSWORD=$DB_PASSWORD psql -h $DB_HOST -U $DB_USER -d $DB_NAME -f migrations/026_consolidated_molecule_constraints_indexes.sql

echo ""
echo "Done!"