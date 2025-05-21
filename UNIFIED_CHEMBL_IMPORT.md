# Unified ChEMBL Import

This document provides an overview of the unified ChEMBL import solution. This streamlined approach combines what were previously multiple scripts into a single, efficient implementation.

## Overview

The unified ChEMBL import solution:

1. Imports cryoprotectant molecules from ChEMBL
2. Calculates all required molecular properties using RDKit
3. Resolves PubChem cross-references using multiple methods
4. Stores everything with proper data integrity
5. Verifies import quality against defined thresholds

## Key Features

- **All-in-one Solution:** No more juggling multiple scripts and dependencies
- **Enhanced Property Calculation:** Complete RDKit-based property calculation for all molecules
- **Cross-Reference Resolution:** Find PubChem CIDs using ChEMBL ID, InChIKey, and SMILES
- **Consistent Data Storage:** Properties stored in both DB tables and JSONB field
- **Verification Reporting:** Built-in quality checks with detailed reports
- **Resumable Processing:** Checkpoint-based operation for resilience
- **Batch Processing:** Efficient handling of large datasets
- **Detailed Logging:** Comprehensive logging for troubleshooting

## Usage

You can use the unified import script directly or through the wrapper shell script.

### Basic Usage

```bash
# Run the script directly
python unified_chembl_import.py --limit 1000 --batch-size 50

# Or use the wrapper script
./run_unified_chembl_import.sh
```

### Common Options

- `--limit NUMBER`: Maximum number of compounds to import (default: 5000)
- `--batch-size NUMBER`: Batch size for import operations (default: 50)
- `--dry-run`: Don't actually insert data, just simulate
- `--verify-only`: Only verify existing data, don't import
- `--verify-limit NUMBER`: Maximum number of molecules to verify (default: all)

### Verification Only

To check the quality of existing data without importing anything new:

```bash
python unified_chembl_import.py --verify-only
```

This will analyze your database and produce a report showing property completeness and cross-reference coverage.

## Quality Thresholds

The verification process uses the following quality thresholds:

- **Property Completeness:** At least 95% of molecules should have all properties
- **JSONB Property Completeness:** At least 95% of molecules should have all properties in the JSONB field
- **PubChem Cross-Reference Coverage:** At least 90% of molecules should have PubChem CIDs

If any of these thresholds are not met, the script will display a warning.

## Requirements

- Python 3.8+
- RDKit (optional but recommended for property calculation)
- chembl_webresource_client (required for ChEMBL API access)
- psycopg2 (required for database access)

## Code Organization

The unified script integrates the following components:

1. **ChEMBL Client:** Fetches compounds from the ChEMBL API
2. **Progress Tracker:** Manages checkpoints and progress reporting
3. **PubChem Resolver:** Finds PubChem CIDs using multiple methods
4. **Property Calculator:** Calculates and stores molecular properties
5. **Database Import:** Handles batch database operations
6. **Verification:** Checks import quality against thresholds

## Next Steps

- Add more test coverage for the unified script
- Enhance error handling for better resilience
- Optimize performance for large datasets
- Add support for more property types