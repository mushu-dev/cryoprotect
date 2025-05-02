# ChEMBL Data Remediation Plan

## Overview

This document specifies the technical remediation plan for ChEMBL data issues in CryoProtect v2, based on CHEMBL_REMEDIATION_GUIDE.md, integrity reports, current codebase analysis, and best practices from Context7 MCP (Model Context Protocol). The plan addresses schema changes, code modifications, reconciliation scripts, and verification steps to ensure ChEMBL data is accurate, complete, queryable, and fully auditable. **All database operations must be executed via MCP tools as defined in `.roo/mcp.json` for maximum automation, robustness, and traceability.**

---

## 1. Objectives

- Ensure ≥1000 ChEMBL molecules are imported, including key reference compounds (e.g., CHEMBL25, CHEMBL1118, CHEMBL1234, CHEMBL444, CHEMBL230130, CHEMBL9335, CHEMBL15151).
- Add a dedicated `chembl_id` column to the `public.molecules` table, backfill from `data_source`, and index it.
- Store ChEMBL IDs and property provenance explicitly in the database.
- Reconcile all molecular property values with official ChEMBL data, correcting discrepancies.
- Provide robust verification and reporting for all remediation steps.
- **Mandate that all remediation steps (migrations, imports, reconciliation, verification) are performed via MCP tools for automation and auditability.**

---

## 2. Schema Changes

**Migration: `migrations/chembl_add_chembl_id.sql` (to be executed via MCP)**

```sql
-- Add chembl_id column to molecules table
ALTER TABLE public.molecules ADD COLUMN IF NOT EXISTS chembl_id VARCHAR(20);

-- Backfill chembl_id from data_source
UPDATE public.molecules
SET chembl_id = SUBSTRING(data_source, 12)
WHERE data_source LIKE 'ChEMBL ID: %' AND chembl_id IS NULL;

-- Create index for fast lookup
CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON public.molecules(chembl_id);
```

**All schema migrations must be applied using MCP tools (e.g., `use_mcp_tool` or `supabase_mcp_tools.py`), referencing `.roo/mcp.json` for configuration.**

---

## 3. Code Modifications

**Target: `ChEMBL_Integrated_Import.py`**

### 3.1 Store ChEMBL ID in dedicated column

- In `transform_chembl_to_molecule` (lines ~477-519):  
  Add `"chembl_id": compound.get('molecule_chembl_id')` to the returned dictionary.

### 3.2 Add Standard Reference Compounds

- After checkpoint directory creation (lines ~328):  
  Add logic to fetch and insert standard ChEMBL compounds by ID (see guide for list).

### 3.3 Track Property Source

- In `transform_chembl_to_properties` (lines ~607-624):  
  Modify property data to include `data_source` as  
  `f"ChEMBL: {compound.get('molecule_chembl_id')}, property: {prop_key}"`.

### 3.4 Increase Default Import Limit

- In `main` (lines ~978-983):  
  Change `parser.add_argument("--limit", type=int, default=2000, ...)`.

**All database writes and queries in import scripts must use MCP SQL execution functions for consistency and auditability.**

---

## 4. Reconciliation Script

**File: `reconcile_chembl_properties.py`**

- For each molecule with a ChEMBL ID:
  - Fetch official data from ChEMBL API.
  - Compare all key properties (LogP, Molecular Weight, HBA, HBD, PSA, RTB).
  - Update local values if they differ from ChEMBL by more than a small tolerance.
  - Log all updates for audit trail (MCP logs).
- Use the following property mapping:
  - "LogP" → ChEMBL "alogp"
  - "Molecular Weight" → "full_mwt"
  - "Hydrogen Bond Acceptor Count" → "hba"
  - "Hydrogen Bond Donor Count" → "hbd"
  - "Topological Polar Surface Area" → "psa"
  - "Rotatable Bond Count" → "rtb"

**All updates and queries must be performed via MCP tools, and all actions/errors must be logged to MCP-managed logs.**

---

## 5. Main Remediation Script

**File: `chembl_remediation_main.py`**

- Orchestrates the following steps:
  1. Applies schema changes (runs migration via MCP).
  2. Runs the full ChEMBL import with 2000+ compounds (using MCP for all DB operations).
  3. Runs the reconciliation script (using MCP).
  4. Verifies all requirements are met (see below).
  5. Generates a comprehensive remediation report (MCP log + JSON artifact).

- **All steps must be automated and idempotent, with full traceability via MCP logs.**

---

## 6. Verification Steps

- Run the following SQL queries to verify remediation (all via MCP):

```sql
-- Check molecule count
SELECT COUNT(*) FROM molecules;

-- Check ChEMBL ID presence
SELECT COUNT(*) FROM molecules WHERE chembl_id IS NOT NULL;

-- Check reference compounds
SELECT * FROM molecules WHERE chembl_id IN ('CHEMBL25', 'CHEMBL1118', 'CHEMBL1234', 'CHEMBL444', 'CHEMBL230130', 'CHEMBL9335', 'CHEMBL15151');

-- Check property counts
SELECT COUNT(*) FROM molecular_properties;

-- Check property sources
SELECT data_source, COUNT(*) FROM molecular_properties GROUP BY data_source;

-- Check LogP values for key molecules
SELECT m.chembl_id, mp.numeric_value
FROM molecules m
JOIN molecular_properties mp ON m.id = mp.molecule_id
JOIN property_types pt ON mp.property_type_id = pt.id
WHERE pt.name = 'LogP'
AND m.chembl_id IN ('CHEMBL25', 'CHEMBL1118', 'CHEMBL1234', 'CHEMBL444');
```

- Acceptance criteria:
  - ≥1000 molecules with valid ChEMBL IDs.
  - All reference compounds present.
  - Property values match ChEMBL within tolerance.
  - All changes logged and auditable via MCP.

- **Verification results and reports must be generated and stored in a standard location (e.g., `reports/chembl_remediation_report_<date>.json`), with all queries and results traceable via MCP logs.**

---

## 7. Error Handling, Resilience, and Automation

- Handle ChEMBL API rate limiting (exponential backoff, sleep).
- Use checkpointing for large imports.
- Dynamically insert missing property types as needed.
- Ensure database credentials and RLS policies are correct.
- **All remediation scripts must implement fallback and retry logic for MCP failures, with clear escalation paths (e.g., alerting, error logs).**
- **Integrate remediation scripts into CI/CD pipelines using MCP for automated, repeatable, and auditable runs.**
- **Reference `.roo/mcp.json` as the canonical configuration for all MCP operations.**

---

## 8. Artifacts

- Migration: `migrations/chembl_add_chembl_id.sql`
- Updated: `ChEMBL_Integrated_Import.py`
- New: `reconcile_chembl_properties.py`
- New: `chembl_remediation_main.py`
- Remediation report: `reports/chembl_remediation_report_<date>.json`
- **MCP logs: `logs/mcp_tool.log` and related artifacts for full auditability**

---

## 9. Diagram

```mermaid
flowchart TD
    A[Apply Schema Migration via MCP] --> B[Run ChEMBL Import (2000+) via MCP]
    B --> C[Insert Reference Compounds via MCP]
    C --> D[Reconcile Properties via MCP]
    D --> E[Verification & Reporting via MCP]
```

---

## 10. References

- [CHEMBL_REMEDIATION_GUIDE.md](CHEMBL_REMEDIATION_GUIDE.md)
- [Context7 MCP API/Tech Stack Docs] (see `.roo/mcp.json` and project documentation)
- [supabase_mcp_tools.py](supabase_mcp_tools.py)
- [use_mcp_tool.py](use_mcp_tool.py)

---

**All steps in this plan must be executed using MCP tools for maximum automation, robustness, and auditability. All actions, errors, and results must be logged and traceable via MCP.**