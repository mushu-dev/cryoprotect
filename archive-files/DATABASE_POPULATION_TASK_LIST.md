# Database Population Task List

This document outlines the comprehensive task list for completing the database population process in CryoProtect v2. These tasks address the current blockers and ensure a reliable database setup for the project.

## Current Status

According to `project_state.json`, the database population is currently blocked by connection issues:

- Overall Status: "Blocked - DB Authentication"
- Verification attempts have failed due to database connection issues
- Critical tasks like reference compound import, ChEMBL data import, and PubChem property enhancement are marked as "Done" but cannot be verified

## Priority Tasks

### 1. Fix Database Connection Issues

**Task ID:** `db-conn-fix-01`  
**Status:** High Priority  
**Assignee:** Roo Agent Team

**Implementation steps:**
1. Fix DNS resolution for Supabase database hostname
2. Ensure consistent environment variables (DB_* and SUPABASE_DB_*)
3. Enhance connection pool wrapper to handle DNS failures
4. Implement IP-based fallback mechanism
5. Update config.py to use both environment variable patterns

**Reference file:** `SUPABASE_CONNECTION_IMPLEMENTATION_GUIDE.md`

### 2. Verify Database Schema

**Task ID:** `db-schema-verify-01`  
**Status:** Blocked by Connection Issues  
**Assignee:** Roo Agent Team

**Implementation steps:**
1. Run the schema verification script once connection is fixed
2. Verify molecule table structure matches the normalized schema
3. Verify property_types and molecular_properties tables exist
4. Ensure all required columns are present in the schema
5. Validate foreign key relationships

**Reference files:**
- `verify_imported_data.py`
- `reports/schema_investigation_report.md`

### 3. Test Data Import Process

**Task ID:** `data-import-test-01`  
**Status:** Blocked by Connection Issues  
**Assignee:** Roo Agent Team

**Implementation steps:**
1. Test reference compound import with a small subset (1-2 compounds)
2. Verify that PropertyManager correctly handles property insertion
3. Validate that both molecule and property data is inserted correctly
4. Test ChEMBL import with a small subset
5. Test PubChem property enhancement with a small subset

**Reference files:**
- `import_reference_compounds.py`
- `import_full_chembl.py`
- `enhance_pubchem_properties.py`
- `property_utils.py`

### 4. Execute Full Data Import

**Task ID:** `data-import-full-01`  
**Status:** Blocked by Connection Issues  
**Assignee:** Roo Agent Team

**Implementation steps:**
1. Run `import_reference_compounds.py` to import reference compounds
2. Run `import_full_chembl.py` to import ChEMBL data (~500 compounds)
3. Run `enhance_pubchem_properties.py` to add properties to PubChem molecules
4. Run `reconcile_chembl_properties.py` to link identical molecules
5. Run `add_performance_indexes.py` to optimize database performance

**Reference files:**
- Tasks `task-dbpop-09` through `task-dbpop-13` in `project_state.json`

### 5. Verify Database Population

**Task ID:** `data-verify-01`  
**Status:** Blocked by Connection Issues  
**Assignee:** Roo Agent Team

**Implementation steps:**
1. Run `verify_imported_data.py --update-project-state` to validate population
2. Verify all success criteria are met:
   - ChEMBL data count (≥500 molecules)
   - Property completeness (LogP, H-bond donors/acceptors)
   - Reference compounds presence (all 9)
   - Cross-referencing (PubChem/ChEMBL IDs linked)
   - Query performance (<100ms)
3. Generate a comprehensive verification report

**Reference files:**
- `verify_imported_data.py`
- `reports/final_data_quality_report.md`

## Success Criteria

The database population process will be considered successful when:

1. **Connection issues are resolved** - Both Supabase and local connections work reliably
2. **Schema verification passes** - All required tables and columns exist
3. **Data import is complete** - All data sources are imported correctly
4. **Data verification passes** - All verification checks pass with SUCCESS status
5. **Performance is acceptable** - Queries complete within 100ms

## Dependency Map

```
┌──────────────────┐
│ Fix Connection   │
│ Issues           │
└──────────────────┘
          │
          ▼
┌──────────────────┐
│ Verify Database  │
│ Schema           │
└──────────────────┘
          │
          ▼
┌──────────────────┐
│ Test Data Import │
│ Process          │
└──────────────────┘
          │
          ▼
┌──────────────────┐
│ Execute Full     │
│ Data Import      │
└──────────────────┘
          │
          ▼
┌──────────────────┐
│ Verify Database  │
│ Population       │
└──────────────────┘
```

## Useful Commands

### Connection Testing

```bash
# Test DNS resolution
python -c "import socket; print(socket.gethostbyname('db.tsdlmynydfuypiugmkev.supabase.co'))"

# Test database connection
python -c "import psycopg2; conn = psycopg2.connect(host='db.tsdlmynydfuypiugmkev.supabase.co', port=5432, dbname='postgres', user='postgres', password='YOUR_PASSWORD'); print('Connection successful!')"

# Test with IP address
python -c "import psycopg2; conn = psycopg2.connect(host='172.64.32.247', port=5432, dbname='postgres', user='postgres', password='YOUR_PASSWORD'); print('Connection successful!')"
```

### Data Import Commands

```bash
# Import reference compounds
python import_reference_compounds.py

# Import ChEMBL data
python import_full_chembl.py

# Enhance PubChem properties
python enhance_pubchem_properties.py

# Reconcile cross-references
python reconcile_chembl_properties.py

# Add performance indexes
python add_performance_indexes.py
```

### Verification Commands

```bash
# Run verification with report
python verify_imported_data.py --report=reports/final_data_quality_report.md

# Run verification and update project state
python verify_imported_data.py --update-project-state
```

## Notes on Schema Normalization

The database schema uses a normalized approach with these main tables:

1. **molecules** - Core molecule information
   - id (PK)
   - name
   - smiles
   - formula
   - pubchem_cid
   - chembl_id
   - inchi_key
   - created_at
   - updated_at

2. **property_types** - Types of molecular properties
   - id (PK)
   - name
   - description
   - units
   - data_type (enum: 'numeric', 'text', 'boolean')
   - created_at
   - updated_at

3. **molecular_properties** - Properties associated with molecules
   - id (PK)
   - molecule_id (FK to molecules.id)
   - property_type_id (FK to property_types.id)
   - numeric_value (nullable)
   - text_value (nullable)
   - boolean_value (nullable)
   - created_at
   - updated_at

All property operations must use the `PropertyManager` class from `property_utils.py` to ensure proper handling of the normalized schema.

## References

- **Connection Fix Guide:** `SUPABASE_CONNECTION_IMPLEMENTATION_GUIDE.md`
- **Project State:** `project_state.json`
- **Schema Investigation:** `reports/schema_investigation_report.md`
- **Property Utils:** `property_utils.py`
- **Verification Script:** `verify_imported_data.py`