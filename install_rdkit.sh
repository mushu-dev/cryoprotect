#!/bin/bash
# Script to install RDKit and other dependencies for unified ChEMBL import

echo "Installing RDKit and dependencies for unified ChEMBL import..."

# Install pip if not already installed
if ! command -v pip &> /dev/null; then
    echo "Installing pip..."
    python -m ensurepip --upgrade
fi

# Install ChEMBL client
echo "Installing ChEMBL WebResource Client..."
pip install chembl_webresource_client

# Install RDKit - this is the lightweight version without 3D capabilities
echo "Installing RDKit..."
pip install rdkit-pypi

# Install psycopg2 for PostgreSQL connectivity
echo "Installing PostgreSQL adapter..."
pip install psycopg2-binary

echo "Installation complete!"
echo "You can now run the unified ChEMBL import with RDKit enabled."