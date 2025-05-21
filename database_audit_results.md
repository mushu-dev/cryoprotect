# CryoProtect Database Audit Results

## Database Schema Overview

The CryoProtect database contains 29 tables with varying degrees of data population. The largest tables by size are:

1. **scientific_data_audit** (15 MB)
2. **molecular_properties** (3.2 MB)
3. **molecules** (2.4 MB)

## Key Tables and Their Relationships

### Core Scientific Data Tables

- **molecules**: Contains 1,480 molecules with basic information (name, formula, SMILES, etc.)
- **molecular_properties**: Contains 7,473 property values associated with molecules
- **property_types**: Defines 58 different types of properties that can be measured/calculated
- **mixtures**: Contains 10 mixtures, but only 1 has defined components
- **mixture_components**: Contains 3 components for a single mixture
- **experiments**: Contains 5 experiment entries
- **experiment_properties**: Contains 3 property measurements from experiments
- **predictions**: Contains 15 prediction records

### Support Tables

- **calculation_methods**: Defines 5 different calculation methods
- **calculation_method**: Contains 2 calculation methods (redundant with calculation_methods)
- **teams**: Defines user/researcher teams
- **projects**: Links teams to specific research projects
- **chembl_import_batches**: Tracks ChEMBL data import operations
- **chembl_import_logs**: Contains logs from ChEMBL imports

## Identified Issues

### 1. Redundant Tables

- **calculation_method** and **calculation_methods** contain similar data but with slight schema differences:
  - Both have core fields: id, name, description, method_type, reference, created_at, updated_at, parameters, version
  - **calculation_methods** adds a created_by field
  - **calculation_methods** has name field as varchar(100) while **calculation_method** uses text
  - Foreign key references exist to both tables

### 2. Complex Constraints

- **molecular_properties** has a check constraint that enforces exactly one of three value fields to be non-null:
  ```sql
  CHECK ((((numeric_value IS NOT NULL) AND (text_value IS NULL) AND (boolean_value IS NULL)) OR 
         ((numeric_value IS NULL) AND (text_value IS NOT NULL) AND (boolean_value IS NULL)) OR 
         ((numeric_value IS NULL) AND (text_value IS NULL) AND (boolean_value IS NOT NULL))))
  ```
  This enforces data quality but makes schema extension challenging.

### 3. Underutilized Tables

- **mixtures**: Only 1 out of 10 mixtures has defined components
- **experiment_properties**: Only 3 records despite having a complete experiment infrastructure
- Several junction tables appear to have minimal usage

### 4. Missing Relationships

Some tables that logically should have foreign key relationships don't have them:
- No explicit link between calculation_method and molecular_properties
- scientific_data_audit doesn't have foreign keys despite apparently tracking changes

## Data Quality Issues

1. **ChEMBL Source Identifiers**: Most ChEMBL sources are prefixed with "vUnknown ID:" suggesting import issues
2. **Placeholder Data**: Several mixtures exist without components
3. **Missing Scientific Data**: Critical properties like glass transition temperature have minimal data

## Next Steps

1. **Consolidate Redundant Tables**: Merge calculation_method into calculation_methods
2. **Clean Source Identifiers**: Remove "vUnknown ID:" prefix from ChEMBL sources
3. **Complete or Remove Placeholder Data**: Either populate mixture components or remove empty mixtures
4. **Add Missing Scientific Classifications**: Add fields to identify known cryoprotectants and their types
5. **Create Data Integrity Test Suite**: Build automated tests to verify data consistency during remediation

## Detailed Table Information

| Table Name | Column Count | Size | PK Count | FK Count | Check Count |
|------------|--------------|------|----------|----------|-------------|
| scientific_data_audit | 10 | 15 MB | 1 | 0 | 6 |
| molecular_properties | 18 | 3240 kB | 1 | 2 | 6 |
| molecules | 20 | 2408 kB | 1 | 0 | 4 |
| predictions | 18 | 216 kB | 1 | 4 | 7 |
| mixtures | 12 | 168 kB | 1 | 0 | 4 |
| experiments | 25 | 160 kB | 1 | 3 | 4 |
| mixture_components | 11 | 144 kB | 1 | 2 | 7 |
| property_types | 8 | 136 kB | 1 | 0 | 5 |
| team_members | 7 | 96 kB | 1 | 3 | 7 |
| projects | 7 | 96 kB | 1 | 1 | 4 |
| calculation_methods | 10 | 80 kB | 1 | 0 | 4 |
| user_profile | 9 | 80 kB | 1 | 1 | 6 |
| property_calculation_queue | 9 | 64 kB | 1 | 1 | 3 |
| experiment_properties | 13 | 48 kB | 1 | 2 | 4 |
| cryoprotection_scores | 9 | 48 kB | 1 | 2 | 3 |
| toxicity_data_source | 8 | 48 kB | 1 | 0 | 2 |
| calculation_method | 9 | 48 kB | 1 | 0 | 2 |
| migrations | 4 | 48 kB | 1 | 0 | 3 |
| proteins | 9 | 48 kB | 1 | 1 | 4 |
| molecule_proteins | 8 | 48 kB | 1 | 3 | 5 |
| molecule_experiments | 9 | 48 kB | 1 | 3 | 5 |
| chembl_import_batches | 13 | 40 kB | 1 | 0 | 2 |
| chembl_import_logs | 10 | 40 kB | 1 | 2 | 1 |
| teams | 5 | 32 kB | 1 | 0 | 4 |
| rls_verification_reports | 4 | 32 kB | 1 | 0 | 1 |
| optimization_reports | 5 | 32 kB | 1 | 0 | 4 |
| shared_resources | 6 | 24 kB | 1 | 2 | 6 |
| shares | 10 | 16 kB | 1 | 1 | 7 |
| lab_verifications | 8 | 16 kB | 1 | 2 | 6 |