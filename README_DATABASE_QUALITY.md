# CryoProtect Database Quality Tools

This repository contains tools to improve and maintain the quality of data in the CryoProtect database.

## Quick Start

Run all quality improvement tools in dry-run mode:
```bash
python run_database_quality_improvements.py --dry-run
```

Execute all quality improvements:
```bash
python run_database_quality_improvements.py
```

## Available Tools

### 1. Molecule Consolidation

Consolidates duplicate molecules with the same InChIKey, selecting the best representative as the primary molecule and updating all related tables.

```bash
python run_duplicate_molecule_consolidation.py [--dry-run]
```

### 2. SMILES Standardization

Standardizes SMILES notation using RDKit's canonical SMILES format.

```bash
python standardize_smiles_notation.py [--dry-run] [--batch-size SIZE] [--max-molecules MAX]
```

### 3. Team Assignment

Assigns users who don't have a team to the default team.

```bash
python assign_users_to_teams.py [--dry-run] [--role ROLE]
```

## Documentation

For detailed information, see:

- [Database Quality Improvement Guide](DATABASE_QUALITY_IMPROVEMENT_GUIDE.md) - Comprehensive overview of all tools
- [Molecule Consolidation Guide](MOLECULE_CONSOLIDATION_GUIDE.md) - Details on molecule deduplication
- [Consolidated Molecules Schema](CONSOLIDATED_MOLECULES_SCHEMA.md) - Technical schema information

## Requirements

- Python 3.8+
- PostgreSQL connection to Supabase
- Environment variables for database connection:
  - `SUPABASE_DB_HOST`
  - `SUPABASE_DB_PORT`
  - `SUPABASE_DB_NAME`
  - `SUPABASE_DB_USER`
  - `SUPABASE_DB_PASSWORD`
- RDKit (for SMILES standardization only)

## Installation

1. Clone this repository
2. Install dependencies:
   ```bash
   pip install psycopg2-binary python-dotenv
   
   # For SMILES standardization
   conda install -c conda-forge rdkit
   # or
   pip install rdkit
   ```

## Audit Trail

All operations performed by these tools create entries in the `scientific_data_audit` table, which can be queried to review the changes made:

```sql
-- View recent audit records
SELECT * FROM scientific_data_audit ORDER BY timestamp DESC LIMIT 100;
```

## Monitoring and Maintenance

These tools are designed to be run periodically:

- After major data imports
- When adding new molecule sources
- After user onboarding processes

The tools are idempotent, meaning they can be run multiple times without causing problems.