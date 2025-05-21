#!/bin/bash
# Complete database population script for CryoProtect (final version)
# This script runs both ChEMBL and the fully fixed property-based PubChem data imports

set -e

# Database connection parameters
HOST="aws-0-us-east-1.pooler.supabase.com"
PORT="5432"
DBNAME="postgres" 
USER="postgres.tsdlmynydfuypiugmkev"
PASSWORD="LDHt\$rkaM&Gmf3X@LQ37"

# Container name
CONTAINER_NAME="cryoprotect-rdkit-minimal"

# Import targets
CHEMBL_TARGET=1000
PUBCHEM_TARGET=100

# Check if container exists and is running
if ! podman ps | grep -q $CONTAINER_NAME; then
    echo "Container $CONTAINER_NAME is not running."
    
    # Check if container exists but is stopped
    if podman ps -a | grep -q $CONTAINER_NAME; then
        echo "Starting existing container..."
        podman start $CONTAINER_NAME
    else
        echo "Creating new container..."
        podman run -d --name $CONTAINER_NAME -p 5002:5000 python:3.10-slim \
            sh -c "pip install rdkit-pypi chembl_webresource_client aiohttp requests psycopg2-binary && mkdir -p /app && tail -f /dev/null"
    fi
fi

# Create required directories in container
echo "Creating required directories..."
podman exec $CONTAINER_NAME bash -c "mkdir -p /app/database /app/logs /app/reports /app/checkpoints"

# Copy scripts to container
echo "Copying scripts to container..."
podman cp ./unified_chembl_import.py $CONTAINER_NAME:/app/
podman cp ./property_based_pubchem_import_numeric_fix.py $CONTAINER_NAME:/app/

# Create a simple script to check counts
cat > check_counts.py << EOL
#!/usr/bin/env python3
import os
import psycopg2

# Database connection parameters from environment
host = "${HOST}"
port = "${PORT}"
dbname = "${DBNAME}"
user = "${USER}"
password = "${PASSWORD}"

# Connect to the database
conn = psycopg2.connect(
    host=host,
    port=port,
    dbname=dbname,
    user=user,
    password=password,
    sslmode='require'
)

# Create a cursor
cursor = conn.cursor()

# Query the tables
cursor.execute("SELECT COUNT(*) FROM molecules")
molecule_count = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM molecular_properties")
property_count = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM property_types")
property_type_count = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM molecules WHERE data_source LIKE '%ChEMBL%'")
chembl_count = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM molecules WHERE data_source LIKE '%PubChem%'")
pubchem_count = cursor.fetchone()[0]

# Display the counts
print(f"Database Counts:")
print(f"- Total molecules: {molecule_count}")
print(f"- Total properties: {property_count}")
print(f"- Property types: {property_type_count}")
print(f"- ChEMBL molecules: {chembl_count}")
print(f"- PubChem molecules: {pubchem_count}")

# Get detailed property type information
cursor.execute("""
    SELECT name, data_type, COUNT(mp.id) AS property_count
    FROM property_types pt
    LEFT JOIN molecular_properties mp ON pt.id = mp.property_type_id
    GROUP BY pt.id, pt.name, pt.data_type
    ORDER BY property_count DESC
""")
property_types = cursor.fetchall()

print("\nProperty Type Details:")
for pt in property_types:
    print(f"- {pt[0]} ({pt[1]}): {pt[2]} properties")

# Close the connection
cursor.close()
conn.close()
EOL

# Copy the check_counts script to the container
podman cp ./check_counts.py $CONTAINER_NAME:/app/

# Check database connection before proceeding
echo "Checking database connection..."
podman exec $CONTAINER_NAME bash -c "cd /app && python -c \"
import psycopg2
try:
    conn = psycopg2.connect(
        host='$HOST',
        port='$PORT',
        dbname='$DBNAME',
        user='$USER',
        password='$PASSWORD',
        sslmode='require'
    )
    print('Database connection successful!')
    conn.close()
except Exception as e:
    print(f'Database connection failed: {e}')
    exit(1)
\""

# Run ChEMBL import
echo "Starting ChEMBL import..."
podman exec $CONTAINER_NAME bash -c "cd /app && SUPABASE_DB_HOST='$HOST' SUPABASE_DB_USER='$USER' SUPABASE_DB_PASSWORD='$PASSWORD' python unified_chembl_import.py --target $CHEMBL_TARGET --output-file /app/reports/chembl_import_\$(date +%Y%m%d_%H%M%S).json"

# Run fixed property-based PubChem import
echo "Starting property-based PubChem import (with complete fixes)..."
podman exec $CONTAINER_NAME bash -c "cd /app && python property_based_pubchem_import_numeric_fix.py \
    --host '$HOST' \
    --port '$PORT' \
    --user '$USER' \
    --password '$PASSWORD' \
    --dbname '$DBNAME' \
    --target $PUBCHEM_TARGET \
    --api-delay 0.5 \
    --batch-size 10 \
    --similarity-threshold 0.8"

# Check final counts
echo "Checking final database counts..."
podman exec $CONTAINER_NAME bash -c "cd /app && python check_counts.py"

# Generate a timestamp for the report
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Collect the latest reports
CHEMBL_REPORT=$(podman exec $CONTAINER_NAME bash -c "ls -t /app/reports/chembl_import_*.json | head -1" | tr -d '\r')
PUBCHEM_REPORT=$(podman exec $CONTAINER_NAME bash -c "ls -t /app/reports/property_pubchem_import_report_*.json | head -1" | tr -d '\r')

# Copy reports locally
mkdir -p reports
if [ -n "$CHEMBL_REPORT" ]; then
    podman cp $CONTAINER_NAME:$CHEMBL_REPORT reports/
fi

if [ -n "$PUBCHEM_REPORT" ]; then
    podman cp $CONTAINER_NAME:$PUBCHEM_REPORT reports/
fi

# Generate summary report
echo "Generating summary report..."
cat > "./reports/database_population_report_${TIMESTAMP}.md" << EOL
# CryoProtect Database Population Report (${TIMESTAMP})

## Summary

This report summarizes the database population process for the CryoProtect project.

- ChEMBL import target: $CHEMBL_TARGET compounds
- PubChem import target: $PUBCHEM_TARGET compounds

## Database Counts

$(podman exec $CONTAINER_NAME bash -c "cd /app && python check_counts.py")

## Detailed Results

### ChEMBL Import

$(if [ -n "$CHEMBL_REPORT" ]; then
    REPORT_NAME=$(basename "$CHEMBL_REPORT")
    echo "Report saved to reports/$REPORT_NAME"
    echo "\`\`\`"
    cat reports/$REPORT_NAME | python -m json.tool
    echo "\`\`\`"
else
    echo "No ChEMBL import report found."
fi)

### Property-based PubChem Import (Complete Fix)

$(if [ -n "$PUBCHEM_REPORT" ]; then
    REPORT_NAME=$(basename "$PUBCHEM_REPORT")
    echo "Report saved to reports/$REPORT_NAME"
    echo "\`\`\`"
    cat reports/$REPORT_NAME | python -m json.tool
    echo "\`\`\`"
else
    echo "No PubChem import report found."
fi)

## Implementation Notes

The property-based PubChem import has been completely fixed to handle database constraints properly:

1. Property types are created with correct data_type values
2. Molecular properties are inserted with appropriate numeric_value, text_value, or boolean_value fields based on the property type
3. All properties satisfy the CHECK constraint that requires exactly one of numeric_value, text_value, or boolean_value to be non-NULL

## Conclusion

Database population completed successfully. The database now contains cryoprotectant molecules 
from both ChEMBL and PubChem sources, with properly structured molecular properties.
This provides a solid foundation for the CryoProtect application.
EOL

echo "Database population completed successfully."
echo "Summary report saved to: ./reports/database_population_report_${TIMESTAMP}.md"