# CryoProtect Analyzer

A database and analysis system for cryoprotectant molecules and mixtures.

## Overview

CryoProtect Analyzer is a comprehensive system for analyzing, storing, and predicting properties of cryoprotectant molecules and mixtures. The system integrates with PubChem to fetch molecular data, applies filtering and scoring algorithms, and provides a robust database schema for storing and analyzing results.

## Project Structure

- **migrations/** - Database migration scripts for Supabase/PostgreSQL
  - `001_initial_schema.sql` - Initial database schema creation script
  - `README.md` - Documentation for the database schema
  - `apply_migration.js` - Utility script to apply migrations to a Supabase project

- **examples/** - Example code for interacting with the database
  - `cryoprotect_db_example.js` - JavaScript/Node.js example using Supabase client
  - `cryoprotect_db_example.py` - Python example using Supabase client

- **PubChem_CryoProtectants.py** - Python script for fetching and analyzing cryoprotectant molecules from PubChem

## Database Schema

The CryoProtect Analyzer database is designed to store:

1. **Molecules** - Individual molecules with their basic information
2. **Molecular Properties** - Properties of molecules (LogP, MW, TPSA, etc.)
3. **Mixtures** - Combinations of molecules in specific concentrations
4. **Mixture Components** - Which molecules are in which mixtures
5. **Experiments** - Experimental data for mixtures
6. **Predictions** - Computational predictions of mixture properties
7. **Calculation Methods** - Different methods used for predictions
8. **Property Types** - Different properties that can be measured or predicted

For detailed information about the database schema, see [migrations/README.md](migrations/README.md).

## Getting Started

### Setting Up the Database

1. Create a Supabase project at [https://app.supabase.io/](https://app.supabase.io/)
2. Apply the database migration using one of these methods:
   - Use the `apply_migration.js` script:
     ```
     cd migrations
     node apply_migration.js
     ```
   - Manually apply the SQL from `migrations/001_initial_schema.sql` in the Supabase dashboard

### Using the API

See the example files for how to interact with the database:

- JavaScript: [examples/cryoprotect_db_example.js](examples/cryoprotect_db_example.js)
- Python: [examples/cryoprotect_db_example.py](examples/cryoprotect_db_example.py)

### Running the PubChem Analysis

The `PubChem_CryoProtectants.py` script fetches and analyzes cryoprotectant molecules from PubChem:

```bash
python PubChem_CryoProtectants.py
```

This will:
1. Fetch molecules from PubChem
2. Filter them based on cryoprotectant-relevant properties
3. Score them according to weighted criteria
4. Save the results to CSV and JSON files

## Cryoprotectant Analysis Criteria

The system uses the following criteria for analyzing cryoprotectants:

### Core Filtering Criteria
- LogP range: -1.5 to 0
- Molecular weight range: 50 to 150 g/mol
- TPSA (Topological Polar Surface Area) range: 40 to 100 Å²
- Required functional groups: OH, CONH2, S=O

### Scoring Weights (Total: 200 points)
- Hydrogen bonding: 50 points
- Solubility & polarity: 40 points
- Membrane permeability: 40 points
- Toxicity & biocompatibility: 30 points
- Protein stabilization: 20 points
- Stability & reactivity: 10 points
- Environmental safety: 10 points

## Features

- **PubChem Integration**: Fetch molecular data directly from PubChem
- **Molecular Filtering**: Filter molecules based on cryoprotectant-relevant properties
- **Scoring Algorithm**: Score molecules based on weighted criteria
- **Database Storage**: Store molecules, properties, mixtures, and experimental data
- **Prediction System**: Predict properties of mixtures based on component properties
- **Experiment Tracking**: Record and compare experimental results with predictions
- **Row-Level Security**: Secure access control for multi-user environments

## Technologies Used

- **Database**: PostgreSQL (via Supabase)
- **API**: Supabase REST API
- **Languages**: Python, JavaScript/Node.js
- **External APIs**: PubChem REST API

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.