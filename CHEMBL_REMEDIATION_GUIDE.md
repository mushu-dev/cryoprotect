# ChEMBL Data Remediation Technical Guide

This guide provides comprehensive technical details for remediating the ChEMBL data issues in CryoProtect v2.

## 1. Key Files

| File | Path | Purpose |
|------|------|---------|
| ChEMBL Import | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ChEMBL_Integrated_Import.py` | Main import script |
| Schema | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/001_initial_schema.sql` | DB schema reference |
| Verification Report | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/chembl_integrity_check_20250427_001142.json` | Current issues |

## 2. Issue Details

1. **Insufficient Data Volume**: 
   - Current: 30 molecules 
   - Required: 1000+ molecules
   - The import was likely run with a small limit setting

2. **Missing ChEMBL ID Column**:
   - Current: ChEMBL IDs stored in `data_source` as "ChEMBL ID: CHEMBLXXX"
   - Required: Dedicated `chembl_id` column for direct querying

3. **Property Value Discrepancies**:
   - 0% match between our values and official ChEMBL
   - LogP values show consistent differences (see integrity report)
   - Current code doesn't track which ChEMBL property was used (alogp vs cx_logp)

4. **Missing Reference Compounds**:
   - Standard compounds like CHEMBL25 (Aspirin) are missing
   - These are essential for validation and comparison

## 3. Schema Enhancement

```sql
-- Add dedicated chembl_id column to molecules table
ALTER TABLE molecules 
ADD COLUMN chembl_id VARCHAR(20);

-- Backfill chembl_id from data_source field
UPDATE molecules
SET chembl_id = SUBSTRING(data_source, 12) 
WHERE data_source LIKE 'ChEMBL ID: %';

-- Add index on chembl_id for performance
CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON molecules(chembl_id);
```

## 4. Code Modifications

### 4.1 Store ChEMBL ID in dedicated column

**File**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ChEMBL_Integrated_Import.py`  
**Lines**: 477-519 (transform_chembl_to_molecule function)

```python
# Modify return dictionary to include chembl_id
return {
    "name": compound.get('pref_name') or compound.get('molecule_chembl_id'),
    "smiles": structures.get('canonical_smiles'),
    "inchi": structures.get('standard_inchi'),
    "inchikey": structures.get('standard_inchi_key'),
    "formula": compound.get('molecule_properties', {}).get('full_molformula'),
    "molecular_weight": compound.get('molecule_properties', {}).get('full_mwt'),
    "pubchem_cid": None,
    "chembl_id": compound.get('molecule_chembl_id'),  # Add this line
    "created_by": user_profile_id,
    "data_source": f"ChEMBL ID: {compound.get('molecule_chembl_id')}",
    "version": 1,
    "modification_history": json.dumps([{
        "timestamp": datetime.now().isoformat(),
        "action": "created",
        "user_id": user_profile_id
    }])
}
```

### 4.2 Add Standard Reference Compounds

**File**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ChEMBL_Integrated_Import.py`  
**Lines**: Insert after line 328 (after creating checkpoint directory)

```python
# Add specific known ChEMBL IDs for standard compounds
standard_chembl_ids = [
    "CHEMBL25",     # Aspirin - common reference
    "CHEMBL1118",   # Glycerol - important cryoprotectant
    "CHEMBL1234",   # DMSO - important cryoprotectant
    "CHEMBL444",    # Ethylene glycol
    "CHEMBL230130", # Propylene glycol
    "CHEMBL9335",   # Trehalose
    "CHEMBL15151"   # Sucrose
]

# Fetch standard compounds by ID first
for chembl_id in standard_chembl_ids:
    try:
        compound = molecule.get(chembl_id)
        if compound:
            # Add properties list if needed
            if 'molecule_properties' in compound:
                compound['properties'] = []
                mol_props = compound.get('molecule_properties', {})
                for prop_key, prop_value in mol_props.items():
                    if prop_value is not None:
                        compound['properties'].append({
                            'property_name': prop_key,
                            'value': prop_value
                        })
            
            all_compounds.append(compound)
            logger.info(f"Fetched standard compound: {chembl_id}")
    except Exception as e:
        log_error(
            error_type="API",
            message=f"Error fetching standard compound {chembl_id}: {str(e)}",
            context={
                "exception": e,
                "compound_id": chembl_id,
                "source": "fetch_cryoprotectant_compounds.fetch_standard"
            }
        )
```

### 4.3 Track Property Source 

**File**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ChEMBL_Integrated_Import.py`  
**Lines**: 607-624 (inside transform_chembl_to_properties function)

```python
# Modify property data dictionary to include source property
properties.append({
    "id": str(uuid.uuid4()),
    "molecule_id": molecule_id,
    "property_type_id": property_type_id,
    "numeric_value": float(mol_props[prop_key]),
    "text_value": None,
    "boolean_value": None,
    "created_by": user_profile_id,
    "data_source": f"ChEMBL: {compound.get('molecule_chembl_id')}, property: {prop_key}",  # Modified line
    "version": 1,
    "modification_history": json.dumps([{
        "timestamp": datetime.now().isoformat(),
        "action": "created",
        "user_id": user_profile_id
    }])
})
```

### 4.4 Increase Default Import Limit

**File**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ChEMBL_Integrated_Import.py`  
**Lines**: 978-983 (in main function)

```python
parser.add_argument("--limit", type=int, default=2000, help="Maximum number of compounds to import")
```

## 5. Property Reconciliation Implementation

Create a new file: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reconcile_chembl_properties.py`

This script should:
1. Get all molecules with ChEMBL IDs from our database
2. For each molecule, fetch official data from ChEMBL API
3. Compare our property values with official values
4. Update our values if they differ by more than a small tolerance
5. Log all updates for audit trail

Key conversion map to use:
- Our "LogP" → ChEMBL "alogp"
- Our "Molecular Weight" → ChEMBL "full_mwt"
- Our "Hydrogen Bond Acceptor Count" → ChEMBL "hba"
- Our "Hydrogen Bond Donor Count" → ChEMBL "hbd"
- Our "Topological Polar Surface Area" → ChEMBL "psa"
- Our "Rotatable Bond Count" → ChEMBL "rtb"

## 6. Full Remediation Script

Create a main remediation script that:
1. Applies schema changes (add chembl_id column)
2. Runs the full ChEMBL import with 2000+ compounds
3. Reconciles property values with official API
4. Verifies all requirements are met
5. Generates a comprehensive report

## 7. Execution Order

1. **Schema Enhancement**:
   - Add chembl_id column
   - Backfill from data_source
   - Create index

2. **Code Modification**:
   - Update transform_chembl_to_molecule to store chembl_id
   - Add standard compound fetching
   - Improve property source tracking
   - Increase default import limit

3. **Run Import**:
   - Execute with limit=2000
   - Ensure checkpoint capability for resilience

4. **Reconcile Properties**:
   - Run reconciliation script
   - Generate reconciliation report

5. **Verification**:
   - Check molecule count (≥1000)
   - Verify reference compounds presence
   - Confirm property value matches

## 8. Verification Queries

```sql
-- Check molecule count
SELECT COUNT(*) FROM molecules;

-- Check ChEMBL ID presence
SELECT COUNT(*) FROM molecules WHERE chembl_id IS NOT NULL;

-- Check reference compounds
SELECT * FROM molecules 
WHERE chembl_id IN ('CHEMBL25', 'CHEMBL1118', 'CHEMBL1234', 'CHEMBL444');

-- Check property counts
SELECT COUNT(*) FROM molecular_properties;

-- Check property sources
SELECT data_source, COUNT(*) 
FROM molecular_properties 
GROUP BY data_source;

-- Check LogP values for key molecules
SELECT m.chembl_id, mp.numeric_value
FROM molecules m
JOIN molecular_properties mp ON m.id = mp.molecule_id
JOIN property_types pt ON mp.property_type_id = pt.id
WHERE pt.name = 'LogP'
AND m.chembl_id IN ('CHEMBL25', 'CHEMBL1118', 'CHEMBL1234', 'CHEMBL444');
```

## 9. Common Issues & Solutions

1. **ChEMBL API Rate Limiting**:
   - If encountering API errors, increase sleep time (line 432)
   - Consider implementing exponential backoff

2. **Large Data Volumes**:
   - Process in smaller batches if memory issues occur
   - Use checkpointing to resume interrupted imports

3. **Property Type Mismatches**:
   - Ensure property_types table has all required types
   - Add dynamic property type insertion for missing types

4. **Database Connection Issues**:
   - Check .env file for correct credentials
   - Verify RLS/security settings don't block operations