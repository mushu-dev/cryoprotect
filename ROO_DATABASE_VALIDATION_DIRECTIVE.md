# ROO DATABASE VALIDATION DIRECTIVE

## TASK OVERVIEW
Perform a comprehensive validation of the database population and scientific data integration in CryoProtect v2. This validation is designed to confirm the quality, completeness, and scientific relevance of imported data from PubChem and ChEMBL.

## SUCCESS CRITERIA
- All reference compounds are correctly imported with complete property sets
- At least 1,000 cryoprotectant compounds are properly imported from each source
- Molecular properties contain scientifically relevant data fields (LogP, H-bond donors/acceptors, etc.)
- Molecular visualization works correctly for all imported compounds
- Database queries demonstrate acceptable performance (<100ms average response time)
- No critical data integrity issues are present in the database
- Cross-references between different identifier systems are accurately maintained

## FILE REFERENCES

### Key Files to Execute:
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/verify_imported_data.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/generate_data_quality_report.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/database/verification/runner.py`

### Key Files to Reference:
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/007_seed_scientific_data.sql`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/008_seed_scientific_data_expanded.sql`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/chembl/reference_compounds.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/cryoprotectant_identifiers.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/chembl_verification_20250429_*.json`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/pubchem_import_report_*.json`

## VALIDATION TASKS

### Task 1: Database Schema Validation
Verify that the database schema is correctly configured to support scientific data storage and retrieval.

**Specific Verification Instructions:**

1. Run the database schema verification script:
```bash
python database/verification/schema.py
```

2. Verify the following tables are present with the correct structure:
   - molecules
   - molecular_properties
   - mixture_components
   - mixtures
   - calculation_methods
   - experiments
   - predictions

3. Check the constraints and relationships between tables:
```bash
python database/verification/constraints.py
```

4. Verify that the JSON property columns in the molecules table support the expected property structure:
```sql
SELECT 
  count(*) as total_count,
  count(properties->'pubchem') as pubchem_count,
  count(properties->'chembl') as chembl_count,
  count(properties->'rdkit') as rdkit_count
FROM molecules;
```

### Task 2: Reference Compound Verification
Verify that all reference cryoprotectant compounds have been correctly imported.

**Specific Verification Instructions:**

1. Run the reference compound verification:
```bash
python verify_imported_data.py --check-reference
```

2. Examine the key reference compounds for complete data:
```sql
SELECT 
  id, 
  name, 
  properties->'pubchem'->'basic'->>'molecular_weight' as molecular_weight,
  properties->'pubchem'->'identifiers'->>'inchi_key' as inchi_key,
  properties->'pubchem'->'properties'->>'logP' as logP,
  properties->'pubchem'->'properties'->>'h_bond_donor_count' as h_donors,
  properties->'pubchem'->'properties'->>'h_bond_acceptor_count' as h_acceptors
FROM molecules
WHERE id IN ('CRYO001', 'CRYO002', 'CRYO003', 'CRYO004', 'CRYO005');
```

3. Check for missing or null values in critical properties:
```sql
SELECT 
  count(*) as molecules_with_incomplete_data
FROM molecules
WHERE 
  properties->'pubchem'->'basic'->>'molecular_weight' IS NULL OR
  properties->'pubchem'->'identifiers'->>'inchi_key' IS NULL OR
  properties->'pubchem'->'identifiers'->>'smiles' IS NULL;
```

### Task 3: PubChem Import Validation
Verify the data quality and completeness of PubChem imports.

**Specific Verification Instructions:**

1. Run the PubChem validation script:
```bash
python verify_imported_data.py --source pubchem
```

2. Check for expected property completeness:
```sql
SELECT
  count(*) as total_pubchem_compounds,
  sum(CASE WHEN properties->'pubchem'->'properties'->>'logP' IS NOT NULL THEN 1 ELSE 0 END) as logp_count,
  sum(CASE WHEN properties->'pubchem'->'properties'->>'h_bond_donor_count' IS NOT NULL THEN 1 ELSE 0 END) as hbond_donor_count,
  sum(CASE WHEN properties->'pubchem'->'properties'->>'h_bond_acceptor_count' IS NOT NULL THEN 1 ELSE 0 END) as hbond_acceptor_count
FROM molecules
WHERE properties->>'pubchem' IS NOT NULL;
```

3. Generate a visualization of property completeness:
```bash
python generate_data_quality_report.py --source pubchem --output reports/pubchem_data_quality.md
```

### Task 4: ChEMBL Import Validation
Verify the data quality and completeness of ChEMBL imports.

**Specific Verification Instructions:**

1. Run the ChEMBL validation script:
```bash
python verify_imported_data.py --source chembl
```

2. Check for property completeness:
```sql
SELECT
  count(*) as total_chembl_compounds,
  sum(CASE WHEN properties->'chembl'->'molecule_properties'->>'full_mwt' IS NOT NULL THEN 1 ELSE 0 END) as mw_count,
  sum(CASE WHEN properties->'chembl'->'molecule_properties'->>'alogp' IS NOT NULL THEN 1 ELSE 0 END) as alogp_count,
  sum(CASE WHEN properties->'chembl'->'molecule_properties'->>'hba' IS NOT NULL THEN 1 ELSE 0 END) as hba_count,
  sum(CASE WHEN properties->'chembl'->'molecule_properties'->>'hbd' IS NOT NULL THEN 1 ELSE 0 END) as hbd_count
FROM molecules
WHERE properties->>'chembl' IS NOT NULL;
```

3. Generate a visualization of property completeness:
```bash
python generate_data_quality_report.py --source chembl --output reports/chembl_data_quality.md
```

### Task 5: Cross-Reference Validation
Verify that molecule identifiers are correctly cross-referenced between different systems.

**Specific Verification Instructions:**

1. Check identifier mapping consistency:
```bash
python cryoprotectant_identifiers.py --validate
```

2. Verify molecules with multiple identifier sources:
```sql
SELECT
  count(*) as total_molecules,
  sum(CASE WHEN properties->>'pubchem' IS NOT NULL AND properties->>'chembl' IS NOT NULL THEN 1 ELSE 0 END) as molecules_with_both_sources,
  sum(CASE WHEN properties->>'pubchem' IS NOT NULL AND properties->>'chembl' IS NULL THEN 1 ELSE 0 END) as pubchem_only,
  sum(CASE WHEN properties->>'pubchem' IS NULL AND properties->>'chembl' IS NOT NULL THEN 1 ELSE 0 END) as chembl_only
FROM molecules;
```

3. Check for identifier conflicts:
```bash
python verify_imported_data.py --check-identifiers
```

### Task 6: Query Performance Validation
Verify that database queries perform within acceptable parameters.

**Specific Verification Instructions:**

1. Run the query performance test:
```bash
python test_database_performance.py
```

2. Analyze query execution plans:
```bash
python analyze_query_plans.py
```

3. Verify pagination performance for large result sets:
```bash
python test_database_performance.py --test-pagination
```

### Task 7: Visualization Integration Validation
Verify that molecular visualization correctly handles the imported molecules.

**Specific Verification Instructions:**

1. Test visualization with reference molecules:
```bash
python verify_rdkit.py --test-visualization
```

2. Verify 3D coordinate generation:
```bash
python verify_rdkit.py --test-3d-coords
```

3. Check visualization for a sample of imported compounds:
```bash
python verify_rdkit.py --test-random-samples 50
```

## COMPREHENSIVE VALIDATION

After completing individual validation tasks, perform a comprehensive validation:

```bash
python verify_database_integrity.py
```

Then generate a comprehensive report:

```bash
python generate_data_quality_report.py --full --output reports/database_quality_report.md
```

## EXPECTED OUTPUTS

The validation should produce the following outputs:

1. Validation success confirmation for each task
2. Data quality visualization showing property completeness
3. Performance metrics for database queries
4. List of any issues found, categorized by severity
5. Recommendations for data quality improvement

## TROUBLESHOOTING GUIDELINES

If you encounter validation failures:

1. For schema issues:
   - Check migration files and execution order
   - Verify that all migrations have been successfully applied

2. For data quality issues:
   - Check import logs in reports/chembl_* and reports/pubchem_*
   - Verify rate limiting and error handling in import scripts

3. For performance issues:
   - Check the query execution plans
   - Verify indexes on commonly queried columns
   - Test with smaller datasets to isolate bottlenecks

4. For visualization issues:
   - Verify SMILES and InChI strings for problematic molecules
   - Check RDKit integration and configuration

Please document any issues found during validation in the quality report, along with recommended remediation steps.