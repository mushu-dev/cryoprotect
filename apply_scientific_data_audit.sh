#!/bin/bash
# Script to create the scientific_data_audit table

set -e

echo "Creating scientific_data_audit table..."
PGPASSWORD=$DB_PASSWORD psql -h $DB_HOST -U $DB_USER -d $DB_NAME -f migrations/027_create_scientific_data_audit.sql

echo ""
echo "Done!"