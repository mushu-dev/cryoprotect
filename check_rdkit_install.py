#!/usr/bin/env python3
try:
    import rdkit
    from rdkit import Chem
    print(f'RDKit is installed. Version: {rdkit.__version__}')
except Exception as e:
    print(f'RDKit is NOT installed: {e}')