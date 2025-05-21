#!/bin/bash
# Run property-based PubChem import tests in the container

# Container name
CONTAINER_NAME="cryoprotect-rdkit-minimal"

# Database parameters
HOST="aws-0-us-east-1.pooler.supabase.com"
PORT="5432"
DBNAME="postgres"
USER="postgres.tsdlmynydfuypiugmkev"
PASSWORD="LDHt\$rkaM&Gmf3X@LQ37"

# Make sure the container is running
if ! podman ps | grep -q $CONTAINER_NAME; then
    echo "Container $CONTAINER_NAME is not running. Starting it..."
    podman start $CONTAINER_NAME || podman run -d --name $CONTAINER_NAME -p 5002:5000 python:3.10-slim \
        sh -c "pip install rdkit-pypi chembl_webresource_client aiohttp requests psycopg2-binary && mkdir -p /app && tail -f /dev/null"
fi

# Make sure directories exist
podman exec $CONTAINER_NAME bash -c "mkdir -p /app/logs /app/reports /app/checkpoints"

# Copy the test script and implementation to the container
podman cp /home/mushu/Projects/CryoProtect/property_based_pubchem_import_fixed.py $CONTAINER_NAME:/app/
podman cp /home/mushu/Projects/CryoProtect/test_property_based_import_fixed.py $CONTAINER_NAME:/app/

# Run the test in the container with environment variables
echo "Running property-based PubChem import tests (fixed)..."
podman exec -e TEST_DB_HOST="$HOST" -e TEST_DB_PORT="$PORT" -e TEST_DB_NAME="$DBNAME" -e TEST_DB_USER="$USER" -e TEST_DB_PASSWORD="$PASSWORD" \
    $CONTAINER_NAME bash -c "cd /app && python -m unittest test_property_based_import_fixed.py"

echo "Tests completed."