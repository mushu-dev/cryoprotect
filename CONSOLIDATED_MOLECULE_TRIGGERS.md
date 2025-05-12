# Consolidated Molecule Triggers

This document explains the database trigger system implemented to maintain the integrity of consolidated molecules in the CryoProtect database.

## Background

As part of the Molecule Data Quality Enhancement initiative (Phase 2), we implemented two strategies for handling duplicate molecules:

1. **Consolidate**: For true duplicates with identical structures, we established a primary/secondary relationship, migrating all data to the primary molecule.
2. **Differentiate**: For molecules with similar names but different structures, we added properties to clearly distinguish them.

To ensure data integrity is maintained in ongoing operations, we've implemented database triggers that:

- Automatically redirect operations to primary molecules
- Prevent modifications to secondary molecules
- Enforce data integrity rules for consolidated molecules
- Prevent deletion of primary molecules that have secondaries 

## Implemented Triggers

### 1. Primary Molecule Redirection

When inserting or updating records that reference molecules, the `redirect_to_primary_molecule` trigger automatically ensures that references to secondary molecules are redirected to their primary molecules:

```sql
CREATE OR REPLACE FUNCTION redirect_to_primary_molecule()
RETURNS TRIGGER AS $$
DECLARE
    primary_id UUID;
BEGIN
    -- Only process if we have a molecule_id in the record
    IF NEW.molecule_id IS NOT NULL THEN
        -- Get the primary molecule ID
        SELECT get_primary_molecule_id(NEW.molecule_id) INTO primary_id;
        
        -- If it's different from the provided ID, use the primary
        IF primary_id != NEW.molecule_id THEN
            NEW.molecule_id := primary_id;
        END IF;
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;
```

This trigger is applied to:
- `mixture_components`
- `molecular_properties`
- `predictions`
- `experiment_properties`

### 2. Consolidated Molecule Validation

The `validate_consolidated_molecule_operation` trigger enforces rules for operations on consolidated molecules:

- Prevents modification of secondary molecules
- Prevents circular references (a molecule can't be its own primary)
- Prevents setting a secondary molecule as a primary

```sql
CREATE OR REPLACE FUNCTION validate_consolidated_molecule_operation()
RETURNS TRIGGER AS $$
BEGIN
    -- Prevent modification of secondary molecules
    IF TG_OP IN ('UPDATE', 'DELETE') AND EXISTS (
        SELECT 1 FROM molecules 
        WHERE id = OLD.id AND primary_molecule_id IS NOT NULL
    ) THEN
        RAISE EXCEPTION 'Cannot modify or delete a secondary molecule (ID: %). Use the primary molecule instead.', OLD.id;
    END IF;
    
    -- For insert or update operations, ensure we're not creating circular references
    IF TG_OP IN ('INSERT', 'UPDATE') AND NEW.primary_molecule_id IS NOT NULL THEN
        -- Check that primary_molecule_id doesn't point to a secondary molecule
        IF EXISTS (
            SELECT 1 FROM molecules 
            WHERE id = NEW.primary_molecule_id AND primary_molecule_id IS NOT NULL
        ) THEN
            RAISE EXCEPTION 'Cannot set primary_molecule_id to a secondary molecule (ID: %).', NEW.primary_molecule_id;
        END IF;
        
        -- Check for circular references
        IF NEW.id = NEW.primary_molecule_id THEN
            RAISE EXCEPTION 'Cannot set primary_molecule_id to itself (ID: %).', NEW.id;
        END IF;
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;
```

### 3. Primary Molecule Deletion Protection

The `prevent_primary_molecule_deletion` trigger prevents the deletion of primary molecules that have secondary molecules referencing them:

```sql
CREATE OR REPLACE FUNCTION prevent_primary_molecule_deletion()
RETURNS TRIGGER AS $$
BEGIN
    -- Check if this molecule is referenced as primary by others
    IF EXISTS (
        SELECT 1 FROM molecules 
        WHERE primary_molecule_id = OLD.id
    ) THEN
        RAISE EXCEPTION 'Cannot delete a primary molecule (ID: %) that has secondary molecules referencing it.', OLD.id;
    END IF;
    
    RETURN OLD;
END;
$$ LANGUAGE plpgsql;
```

### 4. Consolidation Handler

The `handle_molecule_consolidation` trigger ensures proper flags are set when establishing a consolidation relationship:

```sql
CREATE OR REPLACE FUNCTION handle_molecule_consolidation()
RETURNS TRIGGER AS $$
BEGIN
    -- Only run when setting or updating primary_molecule_id
    IF (TG_OP = 'UPDATE' AND NEW.primary_molecule_id IS DISTINCT FROM OLD.primary_molecule_id) OR
       (TG_OP = 'INSERT' AND NEW.primary_molecule_id IS NOT NULL) THEN
        -- Mark the molecule as secondary
        NEW.is_secondary := TRUE;
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;
```

## Helper Functions

### get_primary_molecule_id

This function returns the primary molecule ID for any molecule. If the molecule is a secondary, it returns its primary. Otherwise, it returns the original ID:

```sql
CREATE OR REPLACE FUNCTION get_primary_molecule_id(molecule_uuid UUID)
RETURNS UUID AS $$
DECLARE
    primary_id UUID;
BEGIN
    -- Check if this molecule is consolidated (has a primary_molecule_id)
    SELECT primary_molecule_id INTO primary_id
    FROM molecules
    WHERE id = molecule_uuid AND primary_molecule_id IS NOT NULL;
    
    -- If it has a primary, return that
    IF primary_id IS NOT NULL THEN
        RETURN primary_id;
    END IF;
    
    -- Otherwise, return the original ID (it's already a primary or not consolidated)
    RETURN molecule_uuid;
END;
$$ LANGUAGE plpgsql;
```

## Using Consolidated Molecules in Application Code

### Database Access

When writing SQL queries, you can use the `get_primary_molecule_id` function to ensure you're always working with the primary molecule:

```sql
-- Example: Get molecular properties for a molecule (or its primary)
SELECT * FROM molecular_properties
WHERE molecule_id = get_primary_molecule_id('molecule-uuid-here');
```

### API Endpoints

For API endpoints, you should implement a middleware function that:

1. Resolves molecule IDs to their primary molecules
2. Returns appropriate information about consolidation status

Example middleware function:

```python
def resolve_molecule_id(molecule_id):
    """Resolve a molecule ID to its primary, if it's consolidated."""
    with get_db_connection() as conn:
        with conn.cursor() as cursor:
            cursor.execute("""
                SELECT 
                    id, 
                    primary_molecule_id,
                    CASE WHEN primary_molecule_id IS NOT NULL THEN true ELSE false END as is_consolidated
                FROM molecules
                WHERE id = %s
            """, (molecule_id,))
            result = cursor.fetchone()
            
            if not result:
                return None, False
                
            if result[1]:  # Has a primary_molecule_id
                return result[1], True
            
            return result[0], False
```

## Testing

A comprehensive test script (`test_consolidated_molecule_triggers.py`) is provided to verify that all triggers are working correctly. It tests:

1. Automatic redirection to primary molecules
2. Protection against modifying secondary molecules
3. Prevention of circular references
4. Protection against deleting primary molecules with secondaries

Run the test script after applying the migration to ensure everything is working correctly:

```bash
./apply_consolidated_molecule_triggers.sh
```

## Troubleshooting

If you encounter issues with the triggers:

1. Check the PostgreSQL logs for detailed error messages
2. Verify that all triggers were created successfully by querying the `pg_trigger` system catalog:

```sql
SELECT tgname, tgrelid::regclass
FROM pg_trigger
WHERE tgrelid = 'molecules'::regclass OR tgrelid = 'molecular_properties'::regclass;
```

3. Run the test script to identify specific issues

## Next Steps

With these triggers in place, the database will maintain the integrity of consolidated molecules. However, you should also:

1. Update API endpoints to be aware of consolidated molecules
2. Update the user interface to display consolidation status
3. Add logging for consolidation operations to track changes
4. Consider adding a way to "unconsolidate" molecules if needed in the future