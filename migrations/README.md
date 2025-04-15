# CryoProtect Analyzer Database Schema

This directory contains database migration scripts for the CryoProtect Analyzer project. The database is designed to store and analyze cryoprotectant molecules, their properties, mixtures, and experimental results.

## Schema Overview

The database schema includes the following tables:

1. **molecules** - Stores basic information about molecules from PubChem
   - Primary identifier: UUID
   - PubChem CID (Compound ID)
   - Molecular formula, SMILES notation
   - Auto-generated PubChem link

2. **property_types** - Defines different types of properties that can be measured
   - Examples: Molecular Weight, LogP, TPSA, H-Bond Donors, etc.
   - Data type specification (numeric, text, boolean)

3. **molecular_properties** - Stores properties of individual molecules
   - Links molecules to their properties
   - Supports different data types (numeric, text, boolean)

4. **mixtures** - Stores information about mixtures of molecules
   - Name, description
   - Created by user tracking

5. **mixture_components** - Stores which molecules are in which mixtures
   - Concentration and concentration units
   - Links mixtures to molecules

6. **calculation_methods** - Stores different methods used for predictions
   - Method name, description, version

7. **predictions** - Stores computational predictions of mixture properties
   - Links to mixtures, property types, and calculation methods
   - Confidence level for predictions

8. **experiments** - Stores experimental data
   - Links to mixtures and property types
   - Experimental conditions and date performed

## Views

The schema includes two helpful views:

1. **molecule_with_properties** - Provides molecules with all their properties in a JSON format
2. **mixture_with_components** - Provides mixtures with all their components in a JSON format

## Functions

Several utility functions are included:

1. **import_molecule_from_pubchem** - Imports a molecule from PubChem by CID
2. **calculate_mixture_score** - Calculates a score for a mixture based on its components
3. **compare_prediction_with_experiment** - Compares prediction results with experimental data

## Row-Level Security (RLS)

The schema implements Row-Level Security policies:

- All data is viewable by everyone
- Only authenticated users can insert new data
- Only creators can update or delete their own data

## How to Apply the Migration

### Using Supabase CLI

1. Install the Supabase CLI if you haven't already:
   ```
   npm install -g supabase
   ```

2. Link your project:
   ```
   supabase link --project-ref your-project-ref
   ```

3. Apply the migration:
   ```
   supabase db push
   ```

### Using Supabase Dashboard

1. Navigate to your Supabase project dashboard
2. Go to the SQL Editor
3. Copy the contents of `001_initial_schema.sql`
4. Paste into the SQL Editor and run the query

### Using psql (Direct PostgreSQL Access)

1. Connect to your database:
   ```
   psql -h your-host -U postgres -d postgres
   ```

2. Run the migration:
   ```
   \i /path/to/001_initial_schema.sql
   ```

## Initial Data

The migration script includes initial data for:

- Common property types based on the CryoProtect analysis script
- Basic calculation methods

## Notes on Usage

- When adding a new molecule, use the `import_molecule_from_pubchem` function
- When creating mixtures, add components with their concentrations
- Use the views for easy querying of molecules with properties or mixtures with components
- The `calculate_mixture_score` function can be used to evaluate mixtures
- Compare predictions with experimental results using the `compare_prediction_with_experiment` function