# CryoProtect v2 - Database Population Guide

This guide explains how to populate the Supabase database with scientifically accurate cryoprotectant data using the provided scripts.

## Overview

The database population scripts insert real-world cryoprotectant data into the Supabase database, including:

- Cryoprotectant molecules with accurate chemical identifiers and properties
- Molecular properties (LogP, melting point, glass transition temperature, etc.)
- Cryoprotectant mixtures with detailed compositions
- Experimental data from scientific literature
- Predictions based on computational models

## Prerequisites

Before running the population scripts, ensure you have:

1. A Supabase project set up with the correct schema
2. Python 3.8+ installed
3. Required Python packages installed:
   ```
   pip install python-dotenv supabase
   ```
4. A `.env` file with Supabase credentials:
   ```
   SUPABASE_URL=your_supabase_url
   SUPABASE_KEY=your_supabase_key
   SUPABASE_USER=your_supabase_user_email
   SUPABASE_PASSWORD=your_supabase_user_password
   ```

## Available Scripts

### Main Population Script

The main script `populate_database_supabase.py` populates all tables in the correct order:

1. Property types
2. Molecules and molecular properties
3. Mixtures and mixture components
4. Calculation methods
5. Predictions
6. Experiments and experiment properties

### Usage

```bash
# Populate the database with real data
python populate_database_supabase.py

# Dry run (print actions without executing)
python populate_database_supabase.py --dry-run

# Associate data with a specific project
python populate_database_supabase.py --project-id your_project_id
```

## Data Sources

The data in these scripts is based on scientific literature and includes:

- **Molecules**: Common cryoprotectants like DMSO, glycerol, ethylene glycol, trehalose, etc.
- **Mixtures**: Standard vitrification solutions like VS55, M22, EAFS, etc.
- **Properties**: Physical and chemical properties relevant to cryopreservation
- **Experiments**: Protocols and results from published cryopreservation studies

## Verification

After running the population scripts, you can verify the data was inserted correctly by:

1. Checking the Supabase dashboard to see the number of rows in each table
2. Running SQL queries to check relationships between tables
3. Using the CryoProtect v2 application to browse the populated data

Example verification query:

```sql
-- Check molecule count
SELECT COUNT(*) FROM molecules;

-- Check mixture components
SELECT m.name AS mixture_name, mol.name AS molecule_name, mc.concentration, mc.concentration_unit
FROM mixtures m
JOIN mixture_components mc ON m.id = mc.mixture_id
JOIN molecules mol ON mc.molecule_id = mol.id
LIMIT 10;
```

## Troubleshooting

If you encounter issues:

1. **Authentication errors**: Ensure your Supabase credentials are correct in the `.env` file
2. **Permission errors**: Verify your user has the necessary permissions in Supabase
3. **Schema errors**: Make sure the database schema matches what the scripts expect
4. **Rate limiting**: The scripts include delays to avoid rate limiting, but you may need to adjust these

## Extending the Data

To add more data:

1. Edit the data structures in `populate_database_supabase.py` (CRYOPROTECTANTS, MIXTURES, etc.)
2. Add new property types to the PROPERTY_TYPES list
3. Add new prediction methods to the PREDICTION_METHODS list
4. Run the script with the `--dry-run` flag to verify your changes
5. Run the script without the flag to insert the new data

## Project ID: tsdlmynydfuypiugmkev

The Supabase project ID for CryoProtect v2 is `tsdlmynydfuypiugmkev`. You can use this ID to access the database directly through the Supabase dashboard or API.