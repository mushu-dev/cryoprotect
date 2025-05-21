# ChEMBL Import Guide

This guide provides instructions for importing chemical compound data from the ChEMBL database into CryoProtect.

## Overview

The ChEMBL import process fetches compounds with cryoprotectant properties from the ChEMBL database, transforms them to match our database schema, calculates additional properties, resolves cross-references with PubChem, and validates the imported data.

## Prerequisites

Before running the ChEMBL import, ensure you have:

1. Database connection credentials set in the environment variables
2. RDKit installed (optional but recommended for property calculation)
3. Sufficient database permissions to insert/update records
4. Property types table populated with standardized property definitions

## Import Process

The ChEMBL import process consists of several steps:

1. **Data Acquisition**: Fetch compounds from ChEMBL with cryoprotectant properties
2. **Data Transformation**: Convert ChEMBL data format to CryoProtect schema
3. **Property Calculation**: Calculate missing properties using RDKit when available
4. **Cross-reference Resolution**: Link ChEMBL IDs to PubChem CIDs
5. **Data Validation**: Verify structure and property integrity
6. **Database Integration**: Insert validated compounds and properties into the database

## Running the Import

### Test Import

To run a small test import:

```bash
# Run test import (10 compounds)
./run_chembl_test_import.sh
```

This will:
- Import 10 compounds from ChEMBL
- Validate the imported data
- Generate a detailed validation report

### Full Import

For a full import of all cryoprotectant compounds:

```bash
# Run full import
python unified_chembl_import.py --limit 1000 --batch-size 50 --checkpoint-interval 100
```

Command line arguments:

- `--limit`: Maximum number of compounds to import (default: 1000)
- `--batch-size`: Number of compounds to process in each database batch (default: 50)
- `--checkpoint-interval`: Number of compounds after which to save checkpoint (default: 100)
- `--dry-run`: Simulate import without writing to database
- `--output-dir`: Directory to save import logs and reports
- `--resume`: Resume from previous checkpoint if available
- `--test-mode`: Enable test mode with additional validation

## Validation

After importing data, validate the import using:

```bash
python test_chembl_import_validation.py --output validation_report.json
```

The validation tool:
- Verifies structural integrity of imported molecules
- Validates property data against established ranges
- Checks for missing required properties
- Ensures cross-references are properly resolved
- Generates a detailed validation report

## Troubleshooting

### Common Issues

1. **API Rate Limiting**: The ChEMBL API imposes rate limits. The import script includes retry logic, but imports may take longer if rate limiting occurs.

2. **Missing Structures**: Some compounds in ChEMBL may have missing structural information. These will be flagged during validation.

3. **Property Standardization**: ChEMBL property names and units may differ from our standards. The import script maps these properties to our standardized format.

4. **Database Performance**: Large imports can affect database performance. Consider running imports during off-peak hours.

### Solutions

1. **Resumable Imports**: Use the `--resume` flag to continue from the last checkpoint if an import is interrupted.

2. **Offline Validation**: If database access is limited, use the `--offline` flag with validation to simulate validation.

3. **Missing Properties**: Use the `populate_rdkit_properties.py` script to calculate missing molecular properties.

4. **Cross-references**: Run `resolve_missing_pubchem_references.py` to fill in missing PubChem CIDs.

## Import Verification

After running a full import, verify that:

1. Data integrity is maintained (run validation script)
2. Required properties are present for all compounds
3. PubChem cross-references are resolved where possible
4. Performance indexes are in place for efficient querying

## Post-Import Tasks

After completing the import:

1. Update the molecule consolidation view if it exists
2. Rebuild any performance indexes that may need optimization
3. Generate an import summary report for documentation
4. Review validation results and address any critical issues

## Related Documentation

- [Property Standardization Guide](PROPERTY_STANDARDIZATION_GUIDE.md)
- [Database Performance Guide](DATABASE_PERFORMANCE_GUIDE.md)
- [RDKit Integration Guide](../RDKIT_INTEGRATION_GUIDE.md)
- [ChEMBL API Documentation](https://chembl.gitbook.io/chembl-interface-documentation/web-services)