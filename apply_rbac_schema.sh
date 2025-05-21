#!/bin/bash
echo "Applying RBAC schema migration..."
python3 apply_rbac_schema.py
if [ $? -ne 0 ]; then
    echo "Error applying RBAC schema migration"
    exit 1
fi
echo "RBAC schema migration applied successfully"