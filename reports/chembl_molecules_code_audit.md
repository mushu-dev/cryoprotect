# ChEMBL Molecules Code Audit Report: `cid` vs `pubchem_cid`

**Date:** 2025-04-26  
**Task:** task-imp-chembl-schema-003  
**Author:** Apex Implementer  
**References:**
- `.specs/chembl_molecules_schema_remediation.md`
- `reports/chembl_molecules_schema_inspection.md`
- `reports/chembl_molecules_schema_migration_log.md`
- `ChEMBL_Integrated_Import.py`

## 1. Executive Summary

This audit examines the usage of `cid` and `pubchem_cid` columns across the codebase following the schema remediation that added both columns to the `public.molecules` table. The audit identifies all code references to these columns and provides recommendations for maintaining consistency.

**Key Findings:**
- The canonical schema uses `cid` for PubChem Compound ID
- The migration added `pubchem_cid` as the primary storage column and `cid` as a generated column
- Multiple files reference both column names, with newer ChEMBL integration code primarily using `pubchem_cid`
- The current database schema supports both column names through the generated column approach

## 2. Database Schema Current State

The current database schema for `public.molecules` includes both columns:

```sql
-- Primary storage column
pubchem_cid INTEGER UNIQUE

-- Generated column based on pubchem_cid
cid INTEGER GENERATED ALWAYS AS (pubchem_cid) STORED

-- Generated link using pubchem_cid
pubchem_link TEXT GENERATED ALWAYS AS ('https://pubchem.ncbi.nlm.nih.gov/compound/' || pubchem_cid) STORED
```

This approach allows code to reference either column name while maintaining data consistency.

## 3. Code References Analysis

### 3.1 Import Scripts

#### ChEMBL_Integrated_Import.py
This script does not currently reference either `cid` or `pubchem_cid` in its molecule transformation function:

```python
def transform_chembl_to_molecule(compound: Dict[str, Any], user_profile_id: str) -> Dict[str, Any]:
    # ...
    return {
        "name": compound.get('pref_name') or compound.get('molecule_chembl_id'),
        "smiles": structures.get('canonical_smiles'),
        "inchi": structures.get('standard_inchi'),
        "inchikey": structures.get('standard_inchi_key'),
        "formula": compound.get('molecule_properties', {}).get('full_molformula'),
        "molecular_weight": compound.get('molecule_properties', {}).get('full_mwt'),
        "created_by": user_profile_id,
        "data_source": f"ChEMBL ID: {compound.get('molecule_chembl_id')}",
        # No pubchem_cid or cid field
        # ...
    }
```

**Finding:** The script needs to be updated to include the `pubchem_cid` field in the transformed data.

#### ChEMBL_CryoProtectants_Supabase.py
This script references `pubchem_cid` in its molecule transformation:

```python
molecules_to_insert.append({
    "name": molecule.get("Name") or f"ChEMBL: {molecule['ChEMBL ID']}",
    "smiles": molecule.get("SMILES"),
    "inchi": molecule.get("InChI"),
    "inchikey": molecule.get("InChIKey"),
    "formula": molecule.get("Molecular Formula"),
    "molecular_weight": float(molecule.get("Molecular Weight")) if molecule.get("Molecular Weight") else None,
    "created_by": user_id,
    "data_source": "ChEMBL",
    # No pubchem_cid or cid field
    # ...
})
```

**Finding:** The script needs to be updated to include the `pubchem_cid` field.

#### populate_database_supabase.py
This script explicitly uses `pubchem_cid`:

```python
"formula": "C2H6OS",
"pubchem_cid": 679,
"molecular_weight": 78.13,
```

**Finding:** This script is already using the correct column name.

### 3.2 API and Models

#### api/models.py
The models use `cid` for the PubChem Compound ID:

```python
'id': fields.String,
'cid': fields.Integer,
'name': fields.String,
```

And in the schema:

```python
cid = ma_fields.Integer(
    description="PubChem Compound ID",
)
```

**Finding:** The API models use `cid` instead of `pubchem_cid`.

#### api/resources.py
The API resources use `cid` consistently:

```python
molecule_parser.add_argument('cid', type=int, required=True,
                           help='PubChem Compound ID is required')
```

And for importing from PubChem:

```python
"import_molecule_from_pubchem",
{"p_cid": cid, "p_user_id": user_id}
```

**Finding:** The API resources use `cid` consistently.

### 3.3 Database Functions

The database function for importing from PubChem uses `cid`:

```sql
-- Check if molecule already exists
SELECT id INTO v_molecule_id FROM public.molecules WHERE cid = p_cid;

-- Insert new molecule
INSERT INTO public.molecules (cid, created_by)
VALUES (p_cid, p_user_id)
```

**Finding:** The database functions use `cid` instead of `pubchem_cid`.

### 3.4 Frontend Code

The frontend JavaScript code uses `cid` consistently:

```javascript
document.getElementById('molecule-name')?.textContent = molecule.name || `CID: ${molecule.cid}`;
```

And for PubChem links:

```javascript
if (pubchemLink && molecule.cid) {
  pubchemLink.href = `https://pubchem.ncbi.nlm.nih.gov/compound/${molecule.cid}`;
}
```

**Finding:** The frontend code uses `cid` consistently.

### 3.5 Tests

The test files use `cid` consistently:

```python
'id': self.sample_molecule_id,
'cid': 123456,
'name': 'Glycerol',
```

**Finding:** The test files use `cid` consistently.

## 4. Error Logs Analysis

The error logs show issues with missing `pubchem_cid` column:

```
2025-04-26 21:09:29,624 [ERROR] chembl_import - [1b2f7d0e-814b-4cd9-9d88-dcf3a8feef0b] Database error: Error inserting molecule: {'code': 'PGRST204', 'details': None, 'hint': None, 'message': "Could not find the 'pubchem_cid' column of 'molecules' in the schema cache"}
```

**Finding:** The ChEMBL import script is attempting to use `pubchem_cid` but encountering errors.

## 5. Recommendations

Based on the audit findings, the following recommendations are made:

### 5.1 Maintain Dual Column Approach

The current approach of having both `cid` and `pubchem_cid` columns with `cid` generated from `pubchem_cid` is appropriate given the mixed usage across the codebase. This approach allows for a gradual transition while maintaining compatibility.

### 5.2 Update Import Scripts

1. **ChEMBL_Integrated_Import.py**: Update the `transform_chembl_to_molecule` function to include the `pubchem_cid` field:

```python
def transform_chembl_to_molecule(compound: Dict[str, Any], user_profile_id: str) -> Dict[str, Any]:
    # ...
    return {
        # ... existing fields ...
        "pubchem_cid": None,  # Add this field, can be populated later if available
        # ...
    }
```

2. **ChEMBL_CryoProtectants_Supabase.py**: Update the molecule transformation to include the `pubchem_cid` field:

```python
molecules_to_insert.append({
    # ... existing fields ...
    "pubchem_cid": None,  # Add this field, can be populated later if available
    # ...
})
```

### 5.3 Documentation Updates

1. Add comments in key files explaining the relationship between `cid` and `pubchem_cid`:

```python
# Note: The database schema supports both 'cid' and 'pubchem_cid' columns.
# 'pubchem_cid' is the primary storage column, and 'cid' is a generated column based on 'pubchem_cid'.
# Either column name can be used in queries, but new code should prefer 'pubchem_cid'.
```

2. Update API documentation to clarify that both column names are supported.

### 5.4 Long-term Strategy

For long-term consistency, consider one of the following approaches:

1. **Gradual Migration to `pubchem_cid`**: Gradually update all code to use `pubchem_cid` as the canonical column name.
2. **Standardize on `cid`**: Update ChEMBL integration code to use `cid` to match the rest of the codebase.

The recommended approach is to standardize on `pubchem_cid` as it more clearly indicates the source of the identifier and aligns with the current database schema design where `pubchem_cid` is the primary storage column.

## 6. Implementation Plan

1. **Update Import Scripts**: Modify the ChEMBL import scripts to include the `pubchem_cid` field.
2. **Add Documentation**: Add comments explaining the dual column approach in key files.
3. **Update Error Handling**: Ensure error handling in import scripts can handle both column names.
4. **Test Changes**: Verify that the updated scripts work correctly with the current database schema.
5. **Monitor Usage**: Continue monitoring usage of both column names in new code.

## 7. Conclusion

The audit has identified mixed usage of `cid` and `pubchem_cid` across the codebase. The current database schema supports both column names through the generated column approach, which provides a good balance between compatibility and data consistency.

The recommended approach is to maintain the dual column approach in the short term while gradually standardizing on `pubchem_cid` as the canonical column name for new code. This approach minimizes disruption while moving toward a more consistent naming convention.