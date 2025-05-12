-- 023_create_consolidated_molecule_triggers.sql
-- Migration to add database triggers that maintain the integrity of consolidated molecules

-- First, create the helper function to get primary molecule ID if it exists
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

-- Create function that redirects references to use primary molecule IDs
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

-- Create triggers for each table that references molecules

-- Trigger for mixture_components
DROP TRIGGER IF EXISTS use_primary_molecule_in_mixture_components ON mixture_components;
CREATE TRIGGER use_primary_molecule_in_mixture_components
BEFORE INSERT OR UPDATE ON mixture_components
FOR EACH ROW EXECUTE FUNCTION redirect_to_primary_molecule();

-- Trigger for molecular_properties
DROP TRIGGER IF EXISTS use_primary_molecule_in_molecular_properties ON molecular_properties;
CREATE TRIGGER use_primary_molecule_in_molecular_properties
BEFORE INSERT OR UPDATE ON molecular_properties
FOR EACH ROW EXECUTE FUNCTION redirect_to_primary_molecule();

-- Trigger for predictions
DROP TRIGGER IF EXISTS use_primary_molecule_in_predictions ON predictions;
CREATE TRIGGER use_primary_molecule_in_predictions
BEFORE INSERT OR UPDATE ON predictions
FOR EACH ROW EXECUTE FUNCTION redirect_to_primary_molecule();

-- Trigger for experiment_properties that reference molecules
DROP TRIGGER IF EXISTS use_primary_molecule_in_experiment_properties ON experiment_properties;
CREATE TRIGGER use_primary_molecule_in_experiment_properties
BEFORE INSERT OR UPDATE ON experiment_properties
FOR EACH ROW
WHEN (NEW.molecule_id IS NOT NULL)
EXECUTE FUNCTION redirect_to_primary_molecule();

-- Create function to validate operations on consolidated molecules
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

-- Apply validation trigger to molecules table
DROP TRIGGER IF EXISTS validate_consolidated_molecule_operations ON molecules;
CREATE TRIGGER validate_consolidated_molecule_operations
BEFORE INSERT OR UPDATE OR DELETE ON molecules
FOR EACH ROW EXECUTE FUNCTION validate_consolidated_molecule_operation();

-- Create a trigger function to prevent deletion of primary molecules that have secondaries
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

-- Apply the deletion protection trigger
DROP TRIGGER IF EXISTS protect_primary_molecule_deletion ON molecules;
CREATE TRIGGER protect_primary_molecule_deletion
BEFORE DELETE ON molecules
FOR EACH ROW EXECUTE FUNCTION prevent_primary_molecule_deletion();

-- Create a cascade function to handle property rerouting when establishing consolidation
CREATE OR REPLACE FUNCTION handle_molecule_consolidation()
RETURNS TRIGGER AS $$
BEGIN
    -- Only run when setting or updating primary_molecule_id
    IF (TG_OP = 'UPDATE' AND NEW.primary_molecule_id IS DISTINCT FROM OLD.primary_molecule_id) OR
       (TG_OP = 'INSERT' AND NEW.primary_molecule_id IS NOT NULL) THEN
        -- Mark the molecule as secondary
        NEW.is_secondary := TRUE;
        
        -- We don't actually migrate properties here - that should be done explicitly
        -- This trigger merely ensures proper database consistency
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Apply the consolidation handler trigger
DROP TRIGGER IF EXISTS handle_molecule_consolidation_trigger ON molecules;
CREATE TRIGGER handle_molecule_consolidation_trigger
BEFORE INSERT OR UPDATE ON molecules
FOR EACH ROW
WHEN (NEW.primary_molecule_id IS NOT NULL)
EXECUTE FUNCTION handle_molecule_consolidation();

-- Add documentation comment
COMMENT ON FUNCTION get_primary_molecule_id(UUID) IS 'Gets the primary molecule ID for a molecule. If the molecule is secondary, returns its primary. Otherwise returns the original molecule ID.';
COMMENT ON FUNCTION redirect_to_primary_molecule() IS 'Trigger function to automatically use primary molecule IDs when referencing consolidated molecules.';
COMMENT ON FUNCTION validate_consolidated_molecule_operation() IS 'Ensures operations on consolidated molecules follow data integrity rules.';
COMMENT ON FUNCTION prevent_primary_molecule_deletion() IS 'Prevents deletion of primary molecules that have secondary molecules referring to them.';
COMMENT ON FUNCTION handle_molecule_consolidation() IS 'Sets proper flags when establishing a consolidation relationship.';