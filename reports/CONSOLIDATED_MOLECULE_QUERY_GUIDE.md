# Consolidated Molecule Query Guide

This guide explains how to query the database after the molecule consolidation process.

## Overview

Duplicate molecules in the database have been consolidated using a "primary/secondary" approach:

- **Primary molecules** are the canonical representation of each unique chemical entity
- **Secondary molecules** point to their primary counterpart using the `consolidated_to` property

This approach preserves all molecule IDs while reducing redundancy.

## Database Schema Changes

The following database objects have been added to support molecule consolidation:

### Helper Functions

```sql
-- Get the primary molecule ID for any molecule
CREATE OR REPLACE FUNCTION get_primary_molecule_id(molecule_id UUID)
RETURNS UUID AS $$
DECLARE
    primary_id UUID;
BEGIN
    -- Check if the molecule has a consolidated_to property
    SELECT (properties->>'consolidated_to')::UUID INTO primary_id
    FROM molecules
    WHERE id = molecule_id;
    
    -- If not, it's already a primary molecule
    IF primary_id IS NULL THEN
        RETURN molecule_id;
    ELSE
        RETURN primary_id;
    END IF;
END;
$$ LANGUAGE plpgsql;

-- Check if a molecule is a primary molecule
CREATE OR REPLACE FUNCTION is_primary_molecule(molecule_id UUID)
RETURNS BOOLEAN AS $$
BEGIN
    -- Check if the molecule has no consolidated_to property
    RETURN NOT EXISTS (
        SELECT 1
        FROM molecules
        WHERE id = molecule_id
        AND properties ? 'consolidated_to'
    );
END;
$$ LANGUAGE plpgsql;
```

### Database Views

```sql
-- View that handles consolidated molecules transparently
CREATE OR REPLACE VIEW consolidated_molecules AS
WITH consolidated_mapping AS (
    SELECT 
        m.id,
        CASE 
            WHEN m.properties->>'consolidated_to' IS NOT NULL THEN 
                m.properties->>'consolidated_to'
            ELSE 
                m.id::text
        END AS primary_id,
        m.properties->>'consolidated_to' IS NOT NULL AS is_consolidated
    FROM 
        molecules m
)
SELECT 
    m.id,
    m.name,
    m.pubchem_cid,
    m.smiles,
    m.molecular_formula,
    m.data_source,
    m.created_at,
    m.updated_at,
    m.properties,
    cm.is_consolidated,
    CASE WHEN cm.is_consolidated THEN pm.name ELSE m.name END AS consolidated_name,
    CASE WHEN cm.is_consolidated THEN pm.id ELSE m.id END AS primary_molecule_id
FROM 
    molecules m
JOIN 
    consolidated_mapping cm ON m.id = cm.id
LEFT JOIN 
    molecules pm ON pm.id = cm.primary_id::uuid;

-- View that returns only primary molecules
CREATE OR REPLACE VIEW primary_molecules AS
SELECT 
    m.*
FROM 
    molecules m
WHERE 
    NOT (m.properties ? 'consolidated_to');
```

## Query Patterns

### 1. Get Primary Molecule for Any ID

```sql
-- Using the helper function
SELECT get_primary_molecule_id('c49fd8f6-2e3b-5c4d-9f4a-d0a1e6f7a8b9');

-- Or with a manual query
SELECT 
    CASE 
        WHEN properties->>'consolidated_to' IS NOT NULL THEN 
            (properties->>'consolidated_to')::uuid
        ELSE 
            id
    END AS primary_id
FROM molecules
WHERE id = 'c49fd8f6-2e3b-5c4d-9f4a-d0a1e6f7a8b9';
```

### 2. Check if a Molecule is Primary

```sql
-- Using the helper function
SELECT is_primary_molecule('c49fd8f6-2e3b-5c4d-9f4a-d0a1e6f7a8b9');

-- Or with a manual query
SELECT 
    properties->>'consolidated_to' IS NULL AS is_primary
FROM molecules
WHERE id = 'c49fd8f6-2e3b-5c4d-9f4a-d0a1e6f7a8b9';
```

### 3. Find Only Primary Molecules

```sql
-- Using the view
SELECT * FROM primary_molecules
WHERE name ILIKE '%glycol%';

-- Or with a manual query
SELECT * FROM molecules
WHERE name ILIKE '%glycol%'
AND (properties->>'consolidated_to' IS NULL);
```

### 4. Find Secondary Molecules for a Primary

```sql
SELECT * FROM molecules
WHERE properties->>'consolidated_to' = '9b5b9bbd-ad2c-46ec-a1d9-279fe46a0549';
```

### 5. Query with Transparent Handling of Consolidation

```sql
-- Using the consolidated_molecules view
SELECT * FROM consolidated_molecules
WHERE name ILIKE '%glycol%';
```

## Python Example

```python
def get_molecule_with_primary_info(conn, molecule_id):
    """
    Get a molecule with information about its primary molecule if it's consolidated.
    """
    # Get the primary ID
    with conn.cursor() as cursor:
        cursor.execute("""
            SELECT get_primary_molecule_id(%s) as primary_id
        """, (molecule_id,))
        primary_id = cursor.fetchone()[0]
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Get the original molecule
        cursor.execute("""
            SELECT * FROM molecules WHERE id = %s
        """, (molecule_id,))
        
        molecule = cursor.fetchone()
        if not molecule:
            return None
        
        # Add consolidated info
        molecule['is_consolidated'] = primary_id != molecule_id
        
        # If consolidated, get primary molecule info
        if molecule['is_consolidated']:
            cursor.execute("""
                SELECT * FROM molecules WHERE id = %s
            """, (primary_id,))
            
            primary_molecule = cursor.fetchone()
            if primary_molecule:
                molecule['primary_molecule'] = primary_molecule
                molecule['primary_id'] = primary_id
        
        return molecule
```

## Best Practices

1. **Use Primary Molecules for New Relationships**
   - When creating new mixture components, predictions, or other relationships, use the primary molecule ID

2. **Transparently Handle Both Primary and Secondary**
   - Applications should check if a molecule is consolidated and use its primary when needed
   - Using the provided views and functions simplifies this process

3. **Search Both Primary and Secondary Names**
   - When searching by name, consider searching both the original and primary molecule names

4. **Consider Migration for Legacy Data**
   - For existing data with references to secondary molecules, consider migrating those references to primary molecules

## Conclusion

This consolidated molecule approach preserves all existing molecule IDs while reducing redundancy in the database. The provided database functions, views, and query patterns make it easy to work with the consolidated data model.

For further examples, see the `query_consolidated_molecules_fixed.py` script in the project repository.