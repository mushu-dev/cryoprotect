# CryoProtect Data Validation

This document describes the data validation script for the CryoProtect Analyzer system, which verifies the integrity of cryoprotectant data imported into the Supabase database.

## Overview

The `validate_cryoprotectants_data.py` script performs comprehensive validation checks on the molecule and molecular_property tables to ensure data quality and integrity. It identifies issues such as missing required fields, invalid numeric values, duplicate molecules, and inconsistent relationships.

## Prerequisites

- Python 3.6+
- Supabase project with the CryoProtect schema applied
- Required Python packages:
  - supabase-py (`pip install supabase`)
  - python-dotenv (`pip install python-dotenv`)
- `.env` file with Supabase credentials:
  ```
  SUPABASE_URL=your_supabase_url
  SUPABASE_KEY=your_supabase_key
  SUPABASE_USER=your_supabase_user (optional)
  SUPABASE_PASSWORD=your_supabase_password (optional)
  ```

## Usage

Run the validation script with:

```bash
python validate_cryoprotectants_data.py [--project-id PROJECT_ID] [--report REPORT_PATH]
```

### Arguments

- `--project-id`: (Optional) Validate data for a specific project ID only
- `--report`: (Optional) Path to save the validation report (default: `validation_report.json`)

## Validation Checks

The script performs the following validation checks:

### Molecule Table Validation

1. **Required Fields**: Verifies that all required fields are present and non-null:
   - name
   - smiles
   - inchi
   - inchikey
   - formula
   - molecular_weight

2. **Numeric Value Ranges**: Checks that numeric values are within reasonable ranges:
   - molecular_weight: 0-2000

3. **Duplicate Detection**: Identifies duplicate molecules based on InChIKey

### Molecular Property Validation

1. **Required Fields**: Verifies that all required fields are present and non-null:
   - molecule_id
   - property_type
   - value

2. **Numeric Value Ranges**: Checks that numeric values are within reasonable ranges based on property type:
   - LogP: -10 to 10
   - TPSA: 0 to 500
   - H-Bond Donors: 0 to 20
   - H-Bond Acceptors: 0 to 20
   - Total Score: 0 to 200

3. **Relationship Integrity**: Ensures that properties reference valid molecules

### Property Consistency Validation

1. **Expected Properties**: Checks that molecules have all expected properties:
   - LogP
   - TPSA
   - H-Bond Donors
   - H-Bond Acceptors
   - Total Score
   - PubChem CID

## Validation Report

The script generates a JSON validation report with the following structure:

```json
{
  "validation_time": "ISO-formatted timestamp",
  "summary": {
    "molecules_checked": 15,
    "properties_checked": 87,
    "missing_required_fields": 2,
    "invalid_numeric_values": 3,
    "duplicate_molecules": 1,
    "orphaned_properties": 0,
    "invalid_property_references": 0,
    "total_issues": 6
  },
  "issues": [
    {
      "type": "issue_type",
      "description": "Detailed description of the issue",
      "entity_id": "UUID of the affected entity",
      "severity": "ERROR or WARNING"
    },
    ...
  ]
}
```

### Issue Types

- `missing_required_fields`: Required field is missing or null
- `invalid_numeric_values`: Numeric value is outside the expected range
- `duplicate_molecules`: Multiple molecules with the same InChIKey
- `orphaned_properties`: Property references a non-existent molecule
- `invalid_property_references`: Property references an invalid entity

### Severity Levels

- `ERROR`: Critical issue that should be fixed
- `WARNING`: Potential issue that should be reviewed

## Console Output

The script also prints a summary of the validation results to the console:

```
==================================================
VALIDATION SUMMARY
==================================================
Validation completed at: 2025-04-17 00:25:30.123456
Molecules checked: 15
Properties checked: 87
Total issues found: 6

Issue breakdown:
  - missing_required_fields: 2
  - invalid_numeric_values: 3
  - duplicate_molecules: 1
  - orphaned_properties: 0
  - invalid_property_references: 0
==================================================
```

## Logging

The script logs detailed information to `data_validation.log`, which can be useful for debugging and tracking the validation process.

## Exit Codes

- `0`: Validation completed successfully with no issues
- `1`: Validation completed with issues or errors

## Example

```bash
# Validate all data
python validate_cryoprotectants_data.py

# Validate data for a specific project
python validate_cryoprotectants_data.py --project-id 15e7f0ae-b909-42ff-b980-73a10ebbfcca

# Save report to a custom location
python validate_cryoprotectants_data.py --report custom_report.json
```

## Integration with Import Process

While the validation script is designed to run independently of the import process, it can be integrated into the import workflow:

1. Run the import script: `python PubChem_CryoProtectants_Supabase.py`
2. Run the validation script: `python validate_cryoprotectants_data.py`
3. Review the validation report and fix any issues

## Troubleshooting

### Authentication Issues

If you encounter authentication errors, ensure:
- Your `.env` file contains valid Supabase credentials
- Your user has appropriate permissions to access the tables

### Row Level Security (RLS) Issues

If the script cannot access certain data due to RLS policies:
- Ensure you're providing valid `SUPABASE_USER` and `SUPABASE_PASSWORD` in the `.env` file
- Verify that the authenticated user has appropriate permissions

### Empty Results

If the script reports zero molecules or properties:
- Check that the import process completed successfully
- Verify that you're using the correct project ID (if specified)
- Ensure your database connection is working properly