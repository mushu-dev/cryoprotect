#!/bin/bash
# Run property-based PubChem import with complete fixes including numeric_value handling

# Database parameters
HOST="aws-0-us-east-1.pooler.supabase.com"
PORT="5432"
DBNAME="postgres"
USER="postgres.tsdlmynydfuypiugmkev"
PASSWORD="LDHt\$rkaM&Gmf3X@LQ37"

# Import parameters
TARGET=50  # Starting with a smaller number for testing
API_DELAY=0.5
BATCH_SIZE=10
SIMILARITY_THRESHOLD=0.8

# Create logs directory if it doesn't exist
mkdir -p logs reports checkpoints

# Container name
CONTAINER_NAME="cryoprotect-rdkit-minimal"

# Make sure the container is running
if ! podman ps | grep -q $CONTAINER_NAME; then
    echo "Container $CONTAINER_NAME is not running. Starting it..."
    podman start $CONTAINER_NAME || podman run -d --name $CONTAINER_NAME -p 5002:5000 python:3.10-slim \
        sh -c "pip install rdkit-pypi chembl_webresource_client aiohttp requests psycopg2-binary && mkdir -p /app && tail -f /dev/null"
fi

# Make sure directories exist
podman exec $CONTAINER_NAME bash -c "mkdir -p /app/logs /app/reports /app/checkpoints"

# Copy the fixed script to the container
podman cp /home/mushu/Projects/CryoProtect/property_based_pubchem_import_numeric_fix.py $CONTAINER_NAME:/app/property_based_pubchem_import_numeric_fix.py

# Run the fixed import script in the container
echo "Starting property-based PubChem import (complete fix version)..."
podman exec $CONTAINER_NAME bash -c "cd /app && python property_based_pubchem_import_numeric_fix.py \
    --host '$HOST' \
    --port '$PORT' \
    --user '$USER' \
    --password '$PASSWORD' \
    --dbname '$DBNAME' \
    --target $TARGET \
    --api-delay $API_DELAY \
    --batch-size $BATCH_SIZE \
    --similarity-threshold $SIMILARITY_THRESHOLD"

# Display the results
echo "Copying import report from container..."
mkdir -p reports
LATEST_REPORT=$(podman exec $CONTAINER_NAME bash -c "ls -t /app/reports/property_pubchem_import_report_*.json | head -1" | tr -d '\r')
if [ -n "$LATEST_REPORT" ]; then
    podman cp $CONTAINER_NAME:$LATEST_REPORT reports/
    REPORT_NAME=$(basename "$LATEST_REPORT")
    echo "Import results:"
    cat reports/$REPORT_NAME | python -m json.tool
    echo "Report saved to reports/$REPORT_NAME"
else
    echo "No import report found."
fi

# Copy logs from container
echo "Copying logs from container..."
mkdir -p logs
LATEST_LOG=$(podman exec $CONTAINER_NAME bash -c "ls -t /app/logs/property_pubchem_import_*.log | head -1" | tr -d '\r')
if [ -n "$LATEST_LOG" ]; then
    podman cp $CONTAINER_NAME:$LATEST_LOG logs/
    LOG_NAME=$(basename "$LATEST_LOG")
    echo "Log saved to logs/$LOG_NAME"
    echo "Last 10 log entries:"
    tail -n 10 logs/$LOG_NAME
else
    echo "No import log found."
fi

echo "Property-based PubChem import completed."