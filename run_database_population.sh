#!/bin/bash
# Database population script for CryoProtect
# This script runs both ChEMBL and PubChem data imports

set -e

CONTAINER_NAME="cryoprotect-rdkit-minimal"

# Make sure the container is running
if ! podman ps | grep -q $CONTAINER_NAME; then
    echo "Container $CONTAINER_NAME is not running. Starting it..."
    podman start $CONTAINER_NAME || podman run -d --name $CONTAINER_NAME -p 5002:5000 python:3.10-slim \
        sh -c "pip install rdkit-pypi chembl_webresource_client aiohttp requests psycopg2-binary && tail -f /dev/null"
fi

# Database connection parameters
DB_HOST="${SUPABASE_DB_HOST:-localhost}"
DB_PORT="${SUPABASE_DB_PORT:-5432}"
DB_NAME="${SUPABASE_DB_NAME:-postgres}"
DB_USER="${SUPABASE_DB_USER:-postgres}"
DB_PASSWORD="${SUPABASE_DB_PASSWORD:-postgres}"

# Check if we need environment variables
if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_USER" ] || [ -z "$SUPABASE_DB_PASSWORD" ]; then
    echo "Please set the SUPABASE_DB_HOST, SUPABASE_DB_USER, and SUPABASE_DB_PASSWORD environment variables."
    echo "Example:"
    echo "export SUPABASE_DB_HOST='your_host'"
    echo "export SUPABASE_DB_USER='postgres'"
    echo "export SUPABASE_DB_PASSWORD='your_password'"
    exit 1
fi

# Run ChEMBL import in background
echo "Starting ChEMBL import in the background..."
podman exec -d $CONTAINER_NAME bash -c "cd /app && python unified_chembl_import.py --target 1000 --output-file /app/reports/chembl_import_\$(date +%Y%m%d_%H%M%S).json" &
CHEMBL_PID=$!

# Run PubChem import in background
echo "Starting PubChem import in the background..."
podman exec -d $CONTAINER_NAME bash -c "cd /app && SUPABASE_DB_HOST='$DB_HOST' SUPABASE_DB_USER='$DB_USER' SUPABASE_DB_PASSWORD='$DB_PASSWORD' python simple_pubchem_import.py --target 100 --api-delay 0.5 --batch-size 10" &
PUBCHEM_PID=$!

echo "Both import processes are running in the background."
echo "Monitor ChEMBL import: podman exec $CONTAINER_NAME bash -c 'tail -f /app/logs/unified_chembl_import.log'"
echo "Monitor PubChem import: podman exec $CONTAINER_NAME bash -c 'tail -f /app/logs/pubchem_import_simple_*.log'"

# Wait for both processes to complete
wait $CHEMBL_PID
wait $PUBCHEM_PID

echo "Database population complete. Check the logs and reports directories for details."