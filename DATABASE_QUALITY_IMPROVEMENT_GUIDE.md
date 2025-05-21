# Database Quality Improvement Guide

This guide describes the tools and processes implemented to improve data quality in the CryoProtect database.

## Overview

The CryoProtect database integrates data from multiple sources (PubChem, ChEMBL, etc.), which can lead to various data quality issues:

1. **Duplicate molecules** with the same InChIKey
2. **Inconsistent SMILES notation** formats
3. **Users without team assignments**

This document outlines the tools developed to address these issues and provides guidance on using them.

## Tools Overview

| Tool | Purpose | Description |
|------|---------|-------------|
| `consolidate_duplicate_molecules.py` | Deduplicate molecules | Identifies molecules with identical InChIKeys and consolidates them |
| `standardize_smiles_notation.py` | Standardize chemical formats | Uses RDKit to convert all SMILES strings to canonical format |
| `assign_users_to_teams.py` | Fix user permissions | Assigns users without a team to a default team |

## 1. Molecule Consolidation

Duplicate molecules can cause confusion in the UI, inconsistent property assignment, and inefficient storage.

### Usage

```bash
# Dry run (preview changes without executing)
python run_duplicate_molecule_consolidation.py --dry-run

# Execute consolidation
python run_duplicate_molecule_consolidation.py
```

### Process

1. Identifies molecules with identical InChIKeys
2. Selects the primary molecule from each group (prioritizing those with PubChem CIDs)
3. Updates all dependent tables to reference only the primary molecule
4. Creates audit records for all changes

### Verification

Check the `duplicate_consolidation_summary.json` file for a detailed report of the consolidation process.

## 2. SMILES Standardization

Inconsistent SMILES notation formats can lead to issues with similarity calculations, substructure searches, and data exchange.

### Usage

```bash
# Dry run (preview changes without executing)
python standardize_smiles_notation.py --dry-run

# Execute standardization (process all molecules)
python standardize_smiles_notation.py

# Process with batch size of 100
python standardize_smiles_notation.py --batch-size 100

# Process only 1000 molecules
python standardize_smiles_notation.py --max-molecules 1000
```

### Process

1. Retrieves molecules with non-null SMILES strings
2. Uses RDKit to convert each SMILES to its canonical form
3. Updates the database with standardized notation
4. Creates audit records for all changes

### Requirements

This tool requires RDKit to be installed:

```bash
# Install with conda
conda install -c conda-forge rdkit

# Or with pip
pip install rdkit
```

## 3. Team Assignment

Users without team assignments may have limited access to functionality or incorrectly applied permissions.

### Usage

```bash
# Dry run (preview changes without executing)
python assign_users_to_teams.py --dry-run

# Assign users as members (default)
python assign_users_to_teams.py

# Assign users as admins
python assign_users_to_teams.py --role admin
```

### Process

1. Creates a "Default Team" if it doesn't exist
2. Identifies users without any team assignment
3. Assigns those users to the default team with the specified role
4. Logs all assignments for verification

## Running the Complete Workflow

For a comprehensive data quality improvement process, run the tools in this order:

1. **Consolidate duplicate molecules** to establish clean primary molecular entities
2. **Standardize SMILES notation** to ensure consistent chemical representations
3. **Assign users to teams** to fix permission issues

This sequence ensures that standardization is applied to the consolidated molecules and that all users have proper access to view the improved data.

## Audit and Reporting

All operations performed by these tools create entries in the `scientific_data_audit` table, which can be queried to understand the changes made:

```sql
-- View molecule consolidation audit records
SELECT * FROM scientific_data_audit
WHERE operation = 'consolidated'
ORDER BY timestamp DESC;

-- View SMILES standardization audit records
SELECT * FROM scientific_data_audit
WHERE operation = 'standardize_smiles'
ORDER BY timestamp DESC;
```

## Monitoring and Maintenance

These tools should be run periodically to maintain data quality:

- Run consolidation after major data imports
- Run SMILES standardization when adding new molecule sources
- Run team assignment after user onboarding processes

The tools are designed to be idempotent, meaning they can be run multiple times safely.