# ChEMBL Property Types Remediation Specification

**Author:** Solution Architect  
**Date:** 2025-04-26  
**Related Tasks:** task-chembl-remediate-property-types  
**References:**  
- migrations/001_initial_schema.sql  
- ChEMBL_Integrated_Import.py  
- project_state.json

---

## 1. Objective

Resolve property type retrieval failures in ChEMBL integration by:
- Verifying and remediating the `property_types` table.
- Ensuring all required property types (LogP, Molecular Weight, HBA, HBD, etc.) are present.
- Implementing a robust fallback for unknown property types.

---

## 2. Canonical Property Types

The following property types **must** exist in `public.property_types`:

| Name                          | Description                                 | Units    | Data Type |
|-------------------------------|---------------------------------------------|----------|-----------|
| Molecular Weight              | The mass of a molecule                      | g/mol    | numeric   |
| LogP                          | Octanol-water partition coefficient         |          | numeric   |
| LogD                          | Distribution coefficient                    |          | numeric   |
| TPSA                          | Topological Polar Surface Area              | Å²       | numeric   |
| Hydrogen Bond Donor Count     | Number of hydrogen bond donor groups        |          | numeric   |
| Hydrogen Bond Acceptor Count  | Number of hydrogen bond acceptor groups     |          | numeric   |
| Rotatable Bond Count          | Number of rotatable bonds                   |          | numeric   |
| Toxicity                      | Toxicity information                        |          | text      |
| Stability                     | Chemical stability information              |          | text      |
| Environmental Safety          | Environmental impact information            |          | text      |
| Hydrogen Bonding Score        | Score for hydrogen bonding capability       |          | numeric   |
| Solubility Score              | Score for solubility and polarity           |          | numeric   |
| Membrane Permeability Score   | Score for membrane permeability             |          | numeric   |
| Toxicity Score                | Score for toxicity and biocompatibility     |          | numeric   |
| Stability Score               | Score for stability and reactivity          |          | numeric   |
| Environmental Score           | Score for environmental safety              |          | numeric   |
| Total Score                   | Overall cryoprotectant score                |          | numeric   |

**Note:** Names must match those used in ChEMBL_Integrated_Import.py property mappings (case-insensitive).

---

## 3. Fallback Logic for Unknown Property Types

- If a property type is encountered during import that does **not** exist in `property_types`:
  1. Attempt to insert a new row into `property_types` with:
     - `name`: The new property type name (as encountered)
     - `description`: "Auto-added by ChEMBL import"
     - `units`: NULL (unless known)
     - `data_type`: Infer from value type (`numeric`, `text`, `boolean`)
  2. If insertion fails (e.g., due to constraint), log the error and skip the property.
  3. All fallback insertions must be logged for audit.

- If dynamic insertion is not possible (e.g., due to permissions), assign the property to a generic "Unknown Property Type" entry (which must exist in the table).

---

## 4. Implementation Tasks

1. **Verify property_types table**
   - Query all rows; check for presence of all canonical property types.
   - Add any missing types with correct schema.

2. **Update ChEMBL_Integrated_Import.py**
   - Before property insertion, check if property type exists.
   - If not, attempt to insert it dynamically (see fallback logic).
   - If dynamic insertion fails, use "Unknown Property Type" as fallback.

3. **Testing & Validation**
   - Simulate import with known and unknown property types.
   - Confirm all properties are inserted or logged as skipped with clear error.

---

## 5. Acceptance Criteria

- All canonical property types are present in `property_types`.
- ChEMBL import does not fail due to missing property types.
- Unknown property types are handled per fallback logic.
- All fallback insertions and skips are logged for traceability.

---

## 6. References

- [migrations/001_initial_schema.sql](../migrations/001_initial_schema.sql)
- [ChEMBL_Integrated_Import.py](../ChEMBL_Integrated_Import.py)
- [project_state.json](../project_state.json)