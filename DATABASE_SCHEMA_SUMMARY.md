# CryoProtect Database Schema Summary

## Overview

The CryoProtect database is hosted on Supabase and contains a comprehensive schema for managing cryoprotectant molecules, their properties, mixtures, and experimental data. This document provides a summary of the key tables and their relationships.

## Project Information

- **Project ID**: tsdlmynydfuypiugmkev
- **Project Name**: CryoProtect
- **Status**: ACTIVE_HEALTHY
- **Database Host**: db.tsdlmynydfuypiugmkev.supabase.co
- **Database Version**: 15.8.1.073

## Key Tables

### Molecules

The `molecules` table is the central table storing information about individual cryoprotectant molecules.

**Count**: 733 molecules

**Columns**:
- `id` (uuid) - Primary key
- `name` (varchar) - Molecule name
- `smiles` (varchar) - SMILES notation
- `inchi` (text) - InChI identifier
- `inchikey` (varchar) - InChI key
- `formula` (varchar) - Molecular formula
- `molecular_weight` (numeric) - Molecular weight
- `pubchem_cid` (integer) - PubChem Compound ID
- `chembl_id` (varchar) - ChEMBL ID
- `properties` (jsonb) - JSON properties
- Various metadata fields (created_by, created_at, updated_at, etc.)

### Molecular Properties

The `molecular_properties` table stores various properties of molecules.

**Count**: 741 properties

**Columns**:
- `id` (uuid) - Primary key
- `molecule_id` (uuid) - Foreign key to molecules
- `property_type_id` (uuid) - Foreign key to property_types
- `property_name` (varchar) - Name of the property
- `property_value` (text) - Value of the property
- `property_type` (varchar) - Type of property
- `numeric_value`, `text_value`, `boolean_value` - Typed values
- `unit` (varchar) - Unit of measurement
- Various metadata fields

### Property Types

The `property_types` table defines the types of properties that can be measured.

**Key Property Types**:
- Molecular Weight
- LogP
- TPSA
- Glass Transition Temperature
- Vitrification Concentration
- Toxicity
- Cell Permeability
- Viscosity Effect
- Ice Nucleation Inhibition
- Ice Crystal Growth Inhibition

### Mixtures

The `mixtures` table stores information about mixtures of cryoprotectant molecules.

**Count**: 8 mixtures

**Columns**:
- `id` (uuid) - Primary key
- `name` (varchar) - Mixture name
- `description` (text) - Description
- `properties` (jsonb) - JSON properties
- Various metadata fields

### Mixture Components

The `mixture_components` table links mixtures to their component molecules with concentration information.

**Count**: 3 components

**Example**: VS55 Vitrification Solution contains:
- Propylene glycol (2.2)
- Formamide (1.4)
- Dimethyl sulfoxide (8.4)

### Other Important Tables

- `experiments` - Stores experimental data
- `predictions` - Stores prediction data
- `lab_verifications` - Stores lab verification data
- `teams` and `team_members` - Manages team information
- `projects` - Stores project information
- `cryoprotection_scores` - Stores cryoprotection scoring data

## Database Relationships

- Molecules have many Molecular Properties (one-to-many)
- Mixtures have many Mixture Components (one-to-many)
- Mixture Components link to Molecules (many-to-one)
- Experiments can involve multiple Molecules (many-to-many through molecule_experiments)
- Molecules can be related to multiple Proteins (many-to-many through molecule_proteins)

## Notes

1. The database appears to be populated with real molecular data, primarily from PubChem and ChEMBL.
2. Some property values are currently null, indicating that calculation or data import may still be in progress.
3. The VS55 Vitrification Solution is a well-defined mixture with three components at specific concentrations.
4. The schema includes comprehensive support for experimental data, predictions, and lab verifications.

## Next Steps for Development

1. Complete the calculation of molecular properties where values are currently null
2. Add more mixtures and their components
3. Populate experimental data
4. Link molecules to prediction models
5. Implement advanced queries for molecular filtering based on properties