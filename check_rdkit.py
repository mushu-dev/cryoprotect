#!/usr/bin/env python

# Simple script to check if RDKit is installed
try:
    import rdkit
    print(f"RDKit is installed (version {rdkit.__version__})")
    exit(0)
except ImportError:
    print("RDKit is NOT installed")
    exit(1)
except Exception as e:
    print(f"Error checking RDKit: {e}")
    exit(2)