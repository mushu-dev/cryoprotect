# ChEMBL Molecules Table Schema Remediation Specification

**Author:** Solution Architect  
**Date:** 2025-04-26  
**Related Tasks:** task-chembl-remediate-schema  
**References:**  
- `migrations/001_initial_schema.sql` (canonical schema)
- `ChEMBL_Integrated_Import.py` (import script)
- `project_state.json` (task and log context)

---

## 1. Canonical Schema for `public.molecules`

```sql
CREATE TABLE public.molecules (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    cid INTEGER UNIQUE NOT NULL,  -- PubChem Compound ID
    name TEXT,
    molecular_formula TEXT,
    smiles TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id),
    pubchem_link TEXT GENERATED ALWAYS AS ('https://pubchem.ncbi.nlm.nih.gov/compound/' || cid) STORED
);
```

**Canonical columns:**
- id (UUID, PK)
- cid (INTEGER, UNIQUE, NOT NULL) — PubChem Compound ID
- name (TEXT)
- molecular_formula (TEXT)
- smiles (TEXT)
- created_at (TIMESTAMPTZ, NOT NULL, DEFAULT NOW())
- updated_at (TIMESTAMPTZ, NOT NULL, DEFAULT NOW())
- created_by (UUID, FK)
- pubchem_link (TEXT, generated from cid)

## 2. Schema Drift: `pubchem_cid` vs `cid`

- **Canonical:** Only `cid` is defined for PubChem Compound ID.
- **Observed Issue:** Import scripts and logs reference a missing `pubchem_cid` column, causing database errors.
- **Root Cause:** Code expects `pubchem_cid`, but schema only provides `cid`.

## 3. Remediation Scenarios

### A. Code Can Be Refactored to Use `cid`
- **Preferred:** Refactor all code (import scripts, queries, API) to use `cid` instead of `pubchem_cid`.
- **Pros:** No schema drift, aligns with canonical schema.
- **Cons:** Requires code audit and refactor; may break external integrations if they expect `pubchem_cid`.

### B. Add `pubchem_cid` as a Generated or Alias Column
- **If code refactor is not feasible:** Add `pubchem_cid` as a generated column (if supported) or as a direct copy of `cid`.
- **PostgreSQL 12+ Syntax:**
  ```sql
  ALTER TABLE public.molecules
    ADD COLUMN pubchem_cid INTEGER GENERATED ALWAYS AS (cid) STORED;
  ```
- **If generated columns are not supported:** Add as a regular column and synchronize via triggers or during import.

### C. Both Columns Exist
- **If both `cid` and `pubchem_cid` exist:** Ensure they are always equal. Add a check constraint or trigger to enforce consistency.

## 4. Remediation Plan

1. **Schema Inspection:**
   - Query the current schema for `public.molecules`.
   - Detect presence and type of `cid` and `pubchem_cid`.

2. **Decision Logic:**
   - If only `cid` exists and code can be refactored, proceed with code update.
   - If `pubchem_cid` is required, add as a generated column or regular column mirroring `cid`.

3. **Migration:**
   - Generate and apply migration SQL to add `pubchem_cid` if needed.
   - Ensure all canonical columns exist with correct types and constraints.

4. **Validation:**
   - Confirm schema matches canonical definition.
   - Run import scripts and verify no errors related to `pubchem_cid` or `cid`.

5. **Documentation:**
   - Update schema documentation and code comments to clarify the canonical approach.

## 5. Acceptance Criteria

- The `public.molecules` table contains all canonical columns as per `migrations/001_initial_schema.sql`.
- The `cid` column is present, unique, and not null.
- The `pubchem_cid` column is present if and only if required by code, and always mirrors `cid`.
- No import or application code references a missing column.
- All schema changes are idempotent and safe for production data.
- Validation queries and import scripts run without schema-related errors.

---

## 6. References

- [migrations/001_initial_schema.sql](../migrations/001_initial_schema.sql)
- [ChEMBL_Integrated_Import.py](../ChEMBL_Integrated_Import.py)
- [project_state.json](../project_state.json)