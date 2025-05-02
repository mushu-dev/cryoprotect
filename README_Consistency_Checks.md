# CryoProtect v2 Database Consistency Checks

This document describes the automated consistency checks for the CryoProtect v2 Supabase/Postgres database, how to run them, and how to interpret the results.

---

## Overview

The consistency checks ensure data integrity after database population. They can be run as a standalone script:

```
python check_database_consistency.py
```

The script outputs:
- A human-readable report to the terminal
- A machine-readable JSON report (`consistency_report.json`)

---

## Checks Performed

### 1. Referential Integrity

- **molecular_properties.molecule_id** must reference a valid molecule in `molecules`.
- **molecular_properties.property_type_id** must reference a valid property type in `property_types`.
- **mixture_components.mixture_id** must reference a valid mixture in `mixtures`.
- **mixture_components.molecule_id** must reference a valid molecule in `molecules`.

**Interpretation:**  
PASS means all foreign keys are valid.  
FAIL means there are orphaned records (e.g., a property references a non-existent molecule).

---

### 2. No Duplicate Molecules

- No two molecules should have the same `inchikey`.

**Interpretation:**  
PASS means all molecules are unique by inchikey.  
FAIL means there are duplicate molecules (potential data quality issue).

---

### 3. Required Fields

- **Molecules:** `inchikey`, `name`, and `cid` must be non-null.
- **MolecularProperty:** `molecule_id`, `property_type_id`, and at least one value field (`numeric_value`, `text_value`, `boolean_value`) must be non-null.
- **MixtureComponent:** `mixture_id`, `molecule_id`, `concentration`, and `concentration_unit` must be non-null.

**Interpretation:**  
PASS means all required fields are present.  
FAIL means some records are missing required data.

---

### 4. Logical Consistency

- **MixtureComponent totals:** For mixtures where `concentration_unit` is percent, the sum of concentrations should be close to 100%. For other units, totals should be positive and reasonable.
- **No negative concentrations** in mixture components.

**Interpretation:**  
PASS means all mixtures/components are logically consistent.  
FAIL means some mixtures have inconsistent totals or negative concentrations.

---

## Example Output

### Human-Readable Report

```
============================================================
CryoProtect v2 Database Consistency Check Report
============================================================

--- Referential Integrity ---
[PASS] molecular_properties.molecule_id references molecules.id: All OK
[PASS] molecular_properties.property_type_id references property_types.id: All OK
[PASS] mixture_components.mixture_id references mixtures.id: All OK
[FAIL] mixture_components.molecule_id references molecules.id: 2 orphaned mixture_components

--- No Duplicate Molecules ---
[PASS] No duplicate molecules by inchikey: All OK

--- Required Fields ---
[PASS] Molecules required fields (inchikey, name, cid) non-null: All OK
[FAIL] MolecularProperty required fields (molecule_id, property_type_id, value) non-null: 1 molecular_properties missing required fields
[PASS] MixtureComponent required fields (mixture_id, molecule_id, concentration, concentration_unit) non-null: All OK

--- Logical Consistency ---
[PASS] Mixture component totals (percent mixtures sum to ~100%): All OK
[PASS] No negative concentrations in mixture_components: All OK

See 'consistency_report.json' for full details.
```

### JSON Report (consistency_report.json)

```json
{
  "referential_integrity": [
    {
      "check": "molecular_properties.molecule_id references molecules.id",
      "status": "PASS",
      "details": "All OK",
      "orphans": []
    },
    {
      "check": "molecular_properties.property_type_id references property_types.id",
      "status": "PASS",
      "details": "All OK",
      "orphans": []
    },
    {
      "check": "mixture_components.mixture_id references mixtures.id",
      "status": "PASS",
      "details": "All OK",
      "orphans": []
    },
    {
      "check": "mixture_components.molecule_id references molecules.id",
      "status": "FAIL",
      "details": "2 orphaned mixture_components",
      "orphans": [123, 456]
    }
  ],
  "no_duplicate_molecules": {
    "check": "No duplicate molecules by inchikey",
    "status": "PASS",
    "details": "All OK",
    "duplicates": []
  },
  "required_fields": [
    {
      "check": "Molecules required fields (inchikey, name, cid) non-null",
      "status": "PASS",
      "details": "All OK",
      "molecule_ids": []
    },
    {
      "check": "MolecularProperty required fields (molecule_id, property_type_id, value) non-null",
      "status": "FAIL",
      "details": "1 molecular_properties missing required fields",
      "property_ids": [789]
    },
    {
      "check": "MixtureComponent required fields (mixture_id, molecule_id, concentration, concentration_unit) non-null",
      "status": "PASS",
      "details": "All OK",
      "component_ids": []
    }
  ],
  "logical_consistency": [
    {
      "check": "Mixture component totals (percent mixtures sum to ~100%)",
      "status": "PASS",
      "details": "All OK",
      "mixtures": []
    },
    {
      "check": "No negative concentrations in mixture_components",
      "status": "PASS",
      "details": "All OK",
      "components": []
    }
  ]
}
```

---

## How to Interpret Results

- **PASS:** No issues found for this check.
- **FAIL:** Issues found; see details and IDs for affected records.
- **Orphans:** IDs listed should be investigated and resolved.
- **Duplicates:** Molecules with the same inchikey should be deduplicated.
- **Missing Required Fields:** Records listed should be completed or removed.
- **Logical Consistency:** Mixtures/components with issues should be reviewed for data entry or calculation errors.

---

## Issues & Considerations

- The script requires direct database access; ensure credentials are set in your environment or `.env` file.
- If the schema changes, update the script accordingly.
- Some checks (e.g., logical totals) may need adjustment for new units or business rules.
- The script does not attempt to fix issues, only report them.

---

For questions or improvements, see the project documentation or contact the CryoProtect v2 team.