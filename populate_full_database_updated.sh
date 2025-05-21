#!/bin/bash
# Complete database population script for CryoProtect
# This script runs both ChEMBL and property-based PubChem data imports

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
podman cp ./property_based_pubchem_import.py $CONTAINER_NAME:/app/
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

# Run property-based PubChem import
echo "Starting property-based PubChem import..."
podman exec $CONTAINER_NAME bash -c "cd /app && python property_based_pubchem_import.py \
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
if [ -n "$CHEMBL_REPORT" ]; then
    mkdir -p reports
    podman cp $CONTAINER_NAME:$CHEMBL_REPORT reports/
fi

if [ -n "$PUBCHEM_REPORT" ]; then
    mkdir -p reports
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

### Property-based PubChem Import

$(if [ -n "$PUBCHEM_REPORT" ]; then
    REPORT_NAME=$(basename "$PUBCHEM_REPORT")
    echo "Report saved to reports/$REPORT_NAME"
    echo "\`\`\`"
    cat reports/$REPORT_NAME | python -m json.tool
    echo "\`\`\`"
else
    echo "No PubChem import report found."
fi)

## Conclusion

Database population completed successfully. The database now contains cryoprotectant molecules 
from both ChEMBL and PubChem sources, providing a solid foundation for the CryoProtect application.
EOL

echo "Database population completed successfully."
echo "Summary report saved to: ./reports/database_population_report_${TIMESTAMP}.md"