# Molecule Consolidation Process

This document describes the process of identifying and consolidating duplicate molecules in the CryoProtect database.

## Background

During database population from multiple sources (PubChem, ChEMBL, etc.), duplicate molecules with the same InChIKey were inadvertently created. These duplicates can cause confusion in the UI, inconsistent property assignments, and inefficient storage.

The consolidation process addresses this by:
1. Identifying molecules with duplicate InChIKeys
2. Selecting a primary molecule from each group (prioritizing those with PubChem CIDs)
3. Updating all dependent tables to reference only the primary molecule
4. Maintaining a consolidated_molecules table to track the relationships

## Implementation Details

The `consolidate_duplicate_molecules.py` script implements the consolidation process with these key components:

### 1. Identification

```sql
SELECT inchikey, COUNT(*) as count, array_agg(id) as molecule_ids, 
       array_agg(name) as names, array_agg(pubchem_cid) as pubchem_cids
FROM molecules
WHERE inchikey IS NOT NULL
GROUP BY inchikey
HAVING COUNT(*) > 1
```

### 2. Primary Molecule Selection

The primary molecule is selected based on the following criteria:
1. Presence of a PubChem CID (molecules with CIDs are prioritized)
2. Name quality (test molecules are deprioritized)
3. Creation date (older molecules are preferred)

### 3. Consolidated Molecules Table

The `consolidated_molecules` table maintains a record of the consolidation:
- Each molecule gets a record
- Primary molecules have `molecule_status = 'primary'`
- Duplicates have `molecule_status = 'duplicate'` and reference their primary molecule

### 4. Dependent Table Updates

The following tables are updated to reference only the primary molecule:
- `molecular_properties`
- `mixture_components`
- `predictions`

All updates include audit trail records in the `scientific_data_audit` table.

## Current Status

Our database has 13 InChIKeys with duplicate molecules, affecting a total of 31 molecules. The breakdown is:
- 7 InChIKeys with 2 molecules each
- 5 InChIKeys with 3 molecules each
- 1 InChIKey with 4 molecules

## Running the Consolidation Process

### Dry Run Mode

To preview what changes would be made without actually modifying the database:

```bash
python run_duplicate_molecule_consolidation.py --dry-run
```

### Execution Mode

To perform the actual consolidation:

```bash
python run_duplicate_molecule_consolidation.py
```

## Verification

After running the consolidation process, verify the results:

1. Check the `consolidation_summary.json` file for a detailed report
2. Query the `consolidated_molecules` table to view consolidated relationships:

```sql
SELECT primary_molecule_id, primary_molecule_name, COUNT(*) as duplicate_count
FROM consolidated_molecules
WHERE molecule_status = 'duplicate'
GROUP BY primary_molecule_id, primary_molecule_name
ORDER BY duplicate_count DESC;
```

3. Verify that dependent tables now reference only primary molecules:

```sql
SELECT m.id, m.name, COUNT(mp.id) as property_count
FROM molecules m
JOIN molecular_properties mp ON m.id = mp.molecule_id
WHERE m.id IN (SELECT id FROM consolidated_molecules WHERE molecule_status = 'primary')
GROUP BY m.id, m.name;
```

## Next Steps

- Add a UI layer to show all names for consolidated molecules in the interface
- Consider implementing molecule name standardization to further improve data quality
- Implement a periodic cleanup process to identify and consolidate new duplicates