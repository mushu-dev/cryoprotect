# CryoProtect Analyzer - Supabase Version

This repository contains a refactored version of the CryoProtect Analyzer that uses Supabase for data storage instead of CSV and JSON files.

## Overview

The CryoProtect Analyzer is a tool for analyzing potential cryoprotectant molecules from PubChem. It:

1. Retrieves molecule data from PubChem
2. Filters molecules based on predefined criteria (LogP, molecular weight, TPSA, functional groups)
3. Scores molecules based on various properties
4. Stores the results in a Supabase database

## Prerequisites

- Python 3.6+
- A Supabase project with the CryoProtect schema applied
- Required Python packages:
  - `supabase-py`
  - `requests`
  - `python-dotenv`

## Setup

1. Clone this repository
2. Create a `.env` file based on the `.env.example` template:

```
# Supabase Connection Details
SUPABASE_URL=https://your-project-ref.supabase.co
SUPABASE_KEY=your-anon-key

# Authentication (optional, for running tests)
SUPABASE_USER=user@example.com
SUPABASE_PASSWORD=your-password

# PubChem API Settings (optional)
PUBCHEM_API_DELAY=0.2  # Delay between API calls in seconds
```

3. Install required packages:

```bash
pip install supabase requests python-dotenv
```

4. Ensure your Supabase project has the correct schema applied. You can use the migration script in the `migrations` folder.

## Usage

Run the script:

```bash
python PubChem_CryoProtectants_Supabase.py
```

The script will:
1. Connect to your Supabase project
2. Authenticate (if credentials are provided)
3. Read CIDs from the local file
4. Process each molecule:
   - Fetch properties from PubChem
   - Apply filtering criteria
   - Calculate scores
   - Store data in the Supabase database

## Database Schema

The script uses the following tables:

- `molecules`: Stores basic molecule information (CID, formula, SMILES)
- `property_types`: Defines different types of properties
- `molecular_properties`: Stores properties of molecules

## Key Differences from Original Version

- Data is stored in a Supabase database instead of CSV/JSON files
- Improved error handling and logging
- Authentication support for Supabase
- Database operations follow the schema defined in the migrations

## Troubleshooting

- If you encounter authentication errors, check your Supabase credentials in the `.env` file
- If database operations fail, ensure your Supabase project has the correct schema applied
- For PubChem API rate limiting issues, increase the `PUBCHEM_API_DELAY` value in the `.env` file

## License

[MIT License](LICENSE)