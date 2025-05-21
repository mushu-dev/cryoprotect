#!/bin/bash
# Run fixed PubChem import with correct database parameters

# Database parameters from unified_chembl_import.py
HOST="aws-0-us-east-1.pooler.supabase.com"
PORT="5432"
DBNAME="postgres"
USER="postgres.tsdlmynydfuypiugmkev"
PASSWORD="LDHt\$rkaM&Gmf3X@LQ37"

# Number of compounds to import
TARGET=30

# Run in container
podman exec cryoprotect-rdkit-minimal bash -c "cd /app && python fixed_pubchem_import.py --host '$HOST' --port '$PORT' --user '$USER' --password '$PASSWORD' --dbname '$DBNAME' --target $TARGET --api-delay 0.5 --batch-size 10"