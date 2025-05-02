# ChEMBL Data Verification Specification

**Spec Version:** 1.0  
**Date:** 2025-04-26  
**Author:** Solution Architect  
**Related Task:** task-wv-1-1-database-verification

---

## Purpose

To verify that real ChEMBL data is present in the production database and is accessible via the website/API endpoints. This ensures the scientific integrity and usability of the CryoProtect v2 platform.

---

## Scope

- **Database Presence:** Confirm that canonical ChEMBL tables are populated with real, non-trivial data.
- **Website/API Accessibility:** Confirm that the website/API can retrieve and serve real ChEMBL data to end users.
- **(Optional) Data Integrity:** Spot-check a sample of records for expected ChEMBL values.

---

## Verification Steps

### 1. Database Presence Verification

#### 1.1 Identify Canonical ChEMBL Tables

- `molecules`
- `molecular_properties`
- `property_types`
- (Add others as relevant per schema)

#### 1.2 SQL Queries

- **Row Count Check:**
  ```sql
  SELECT COUNT(*) FROM public.molecules;
  SELECT COUNT(*) FROM public.molecular_properties;
  ```
  - Acceptance: Each table contains >1000 records.

- **Non-Null Key Fields:**
  ```sql
  SELECT COUNT(*) FROM public.molecules WHERE chembl_id IS NOT NULL AND smiles IS NOT NULL;
  ```
  - Acceptance: >95% of records have non-null `chembl_id` and `smiles`.

- **Known ChEMBL ID Check:**
  ```sql
  SELECT * FROM public.molecules WHERE chembl_id = 'CHEMBL25';
  ```
  - Acceptance: At least one known ChEMBL ID is present.

#### 1.3 Data Integrity Spot-Check (Optional)

- Manually inspect a few records for expected field values (e.g., molecular weight, LogP).

---

### 2. Website/API Accessibility Verification

#### 2.1 Identify Endpoints

- `/api/molecules`
- `/api/molecular_properties`
- (Adjust as per actual API routes)

#### 2.2 Test Requests

- **Example:**
  ```
  GET /api/molecules?limit=10
  ```
  - Expected: JSON array of molecule objects with real ChEMBL data.

- **Check for:**
  - Non-empty response
  - Presence of canonical fields (`chembl_id`, `smiles`, etc.)
  - Data matches what is in the database

#### 2.3 Authentication/Access

- Test both authenticated and (if allowed) anonymous access.
- Confirm RLS/public policies do not block required data.

---

## Acceptance Criteria

- [ ] At least one canonical ChEMBL table contains >1000 records with non-null key fields.
- [ ] Website/API endpoints return real ChEMBL data matching DB content.
- [ ] (Optional) Spot-checked records match expected ChEMBL values.
- [ ] No placeholder or test data is returned.

---

## References

- Import Script: `ChEMBL_Integrated_Import.py`
- Test Workflow: `test_real_data_workflows_final.py`
- Schema/Remediation Specs: `.specs/chembl_molecules_schema_remediation.md`, `.specs/chembl_property_types_remediation.md`
- Reports: `reports/chembl_import_report.md`, `reports/chembl_molecules_schema_validation.md`
- project_state.json

---

## Change Log

- 2025-04-26: Initial version created for task-wv-1-1-database-verification.