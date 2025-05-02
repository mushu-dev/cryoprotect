# File Structure Guide for CryoProtect v2

This guide clarifies which files exist, which need to be created, and the proper directory structure to ensure agents can correctly reference and create files.

## Existing Directories

```
chembl/               ✅ EXISTS
  __init__.py         ✅ EXISTS
  cache.py            ✅ EXISTS
  client.py           ✅ EXISTS
  error_handler.py    ✅ EXISTS (but needs implementation)
  logging.py          ✅ EXISTS
  rate_limiter.py     ✅ EXISTS
  utils.py            ✅ EXISTS

pubchem/              ✅ EXISTS
  __init__.py         ✅ EXISTS
  cache.py            ✅ EXISTS
  client.py           ✅ EXISTS
  rate_limiter.py     ✅ EXISTS
  utils.py            ✅ EXISTS
  worker_pool.py      ✅ EXISTS

logs/                 ✅ EXISTS

checkpoints/          ✅ EXISTS
```

## Files to Create

```
chembl/
  checkpoint.py       🔴 CREATE
  property_mapper.py  🔴 CREATE
  structure_normalizer.py 🔴 CREATE
  metadata_handler.py 🔴 CREATE
  reference_compounds.py 🔴 CREATE
  worker.py           🔴 CREATE
  progress_tracker.py 🔴 CREATE

run_chembl_import.py  🔴 CREATE
verify_chembl_data.py 🔴 CREATE
```

## Reference Files

These files provide useful reference implementations:

```
PubChem_CryoProtectants_Supabase.py
PubChem_CryoProtectants_Supabase_Enhanced.py
PubChem_CryoProtectants_Supabase_Enhanced_MCP.py
ChEMBL_Integrated_Import.py
```

## File Access Guidelines

1. **Always Check If File Exists First**: Before reading or modifying files, check if they exist.

2. **Create Directory Structure First**: Ensure required directories exist before creating files.

3. **Use Absolute Paths**: Always use absolute paths when referencing files.

4. **For New Files**:
   - First check if the file already exists
   - If it doesn't exist, create it with initial structure
   - Then implement the required functionality

5. **For Existing Files**:
   - Read the file first to understand its current state
   - Implement changes while preserving existing structure

## Creating New Files

When creating new files:

```python
# Example for creating a new Python module
import os

# Ensure directory exists
os.makedirs("chembl", exist_ok=True)

# Check if file exists
file_path = "chembl/checkpoint.py"
if not os.path.exists(file_path):
    # Create initial structure
    with open(file_path, "w") as f:
        f.write('''"""
Checkpoint management for ChEMBL data import.

This module provides functionality for saving and loading progress 
during ChEMBL data import operations.
"""

import os
import json
from datetime import datetime
from typing import Dict, Any, List, Optional
import logging

logger = logging.getLogger(__name__)

# Implementation goes here
''')
    print(f"Created {file_path}")
else:
    print(f"{file_path} already exists")
```

## Reference Implementation Pattern Example

When implementing a file based on a reference:

```python
# Step 1: Check if the file exists
import os

file_path = "chembl/worker.py"
if not os.path.exists(file_path):
    print(f"File {file_path} does not exist. Creating...")
    
    # Step 2: Read reference implementation
    reference_path = "pubchem/worker_pool.py"
    if os.path.exists(reference_path):
        with open(reference_path, "r") as f:
            reference_content = f.read()
    else:
        print(f"Reference file {reference_path} does not exist")
        reference_content = ""
    
    # Step 3: Create initial structure based on reference
    with open(file_path, "w") as f:
        f.write("""
# Implementation based on reference: pubchem/worker_pool.py
# Adapted for ChEMBL-specific functionality

# Your implementation here
""")
else:
    print(f"File {file_path} already exists")
```

By following these guidelines, you'll avoid wasting tokens on file not found errors and ensure proper implementation.