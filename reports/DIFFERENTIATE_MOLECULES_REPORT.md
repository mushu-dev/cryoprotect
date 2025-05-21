# Molecule Differentiation Report

## Overview

This report documents the implementation of the "differentiate" strategy for molecules that were identified as similar but having different chemical structures. This completes the consolidation phase of the Molecule Data Quality Enhancement project.

## Differentiation Process

The differentiation process was implemented to handle molecules with similar names or identifiers that have distinct chemical structures. Unlike the consolidation process for true duplicates, differentiation preserves each molecule as separate but adds properties to clarify their differences.

### Implementation Approach

1. **Analysis**: Molecules with matching names or identifiers but different SMILES strings or molecular formulas were classified into differentiation groups.

2. **Property Addition**: Each molecule was enhanced with differentiation properties in a structured JSON format:
   ```json
   {
     "differentiation_group": "group_uuid",
     "has_structural_differences": true,
     "smiles_structure": "C1CC(=O)OCC1",
     "formula": "C6H10O2"
   }
   ```

3. **Batch Processing**: For large groups, we implemented a batched approach to handle database load efficiently.

## Results

The differentiation process successfully processed:

- **4 differentiation groups**
- **248 total molecules** that needed differentiation
- **100% completion** with no errors or issues

### Group Details

1. **Group 1d968532-913b-4744-b67a-65ee497767a8**:
   - 5 molecules differentiated
   - Included variants of Sucrose and Trehalose with distinct chemical structures
   - Example molecules:
     - `(2R,3R,4S,5S,6R)-2-[(2S,3S,4S,5R)-3,4-dihydroxy-2,5-bis(hydroxymethyl)oxolan-2-yl]oxy-6-(hydroxymethyl)oxane-3,4,5-triol`
     - `Sucrose`
     - `Trehalose`

2. **Group 476326ba-5597-442b-b21b-b5a92b82ec65**:
   - 2 molecules differentiated
   - Different IUPAC representations of similar complex sugar molecules

3. **Group 76397481-971c-41b5-a695-82072235de5b**:
   - 3 molecules differentiated
   - Variants of amino acids with different stereochemistry
   - Example molecules:
     - `2-azaniumylpentanedioate`
     - `(2R)-2-amino-5-hydroxy-5-oxopentanoate`
     - `(2S)-2-amino-5-hydroxy-5-oxopentanoate`

4. **Group b278ce51-6575-46f2-ae4e-79369475cc8d**:
   - 238 molecules differentiated
   - Large group with diverse chemical structures
   - Required optimized batch processing
   - Completed successfully with all molecules properly tagged

## Database Impact

The differentiation process made the following changes to the database:

1. **No Relationships Changed**: Unlike consolidation, differentiation doesn't modify or migrate relationships between molecules.

2. **Added Properties**: Each molecule received a new `differentiation` property object in its JSONB properties field.

3. **Search Improvements**: The added properties enable searching for molecules within the same differentiation group.

## Application Integration

Applications can now use the differentiation properties to:

1. **Display Structural Differences**: Show users that similarly named molecules have different chemical structures.

2. **Group Related Variants**: Present related chemical variants together, while maintaining their distinct identities.

3. **Filter by Structural Properties**: Allow searching by specific structural variations within similarly named compounds.

## Code Example

```python
def get_differentiated_variants(conn, molecule_id):
    """
    Get all differentiated variants for a molecule.
    """
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # First, check if this molecule is part of a differentiation group
        cursor.execute("""
            SELECT properties->'differentiation'->>'differentiation_group' AS group_id
            FROM molecules
            WHERE id = %s AND properties->'differentiation' IS NOT NULL
        """, (molecule_id,))
        
        result = cursor.fetchone()
        if not result or not result['group_id']:
            return []
            
        group_id = result['group_id']
        
        # Find all molecules in this differentiation group
        cursor.execute("""
            SELECT id, name, molecular_formula, smiles, 
                   properties->'differentiation' AS differentiation_info
            FROM molecules
            WHERE properties->'differentiation'->>'differentiation_group' = %s
              AND id != %s
        """, (group_id, molecule_id))
        
        return cursor.fetchall()
```

## Conclusion

The differentiation implementation successfully completes the molecule consolidation and differentiation phase of the Data Quality Enhancement project. With this milestone, the database now provides:

1. **Consolidated true duplicates**: Molecules that represent the same chemical entity are now properly consolidated with a primary/secondary relationship.

2. **Differentiated similar molecules**: Molecules with similar names but distinct structures are properly marked and can be displayed appropriately.

3. **Improved data quality**: The overall data quality is enhanced by providing clear relationships between similar or identical molecules.

4. **Better user experience**: Applications can now provide clearer information to users about the relationships between molecules.

This work forms a solid foundation for further data quality enhancements in future phases of the project.