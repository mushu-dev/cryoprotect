-- Migration: 026_restore_foreign_keys.sql
-- Purpose: Restore foreign key constraints that were previously dropped

-- First, create a function to help with safely adding foreign keys
-- This function checks if the foreign key already exists before adding it
CREATE OR REPLACE FUNCTION safe_add_foreign_key(
    source_table TEXT,
    source_column TEXT, 
    target_table TEXT, 
    target_column TEXT,
    constraint_name TEXT DEFAULT NULL
) RETURNS VOID AS $$
DECLARE
    actual_constraint_name TEXT;
    source_data_type TEXT;
    target_data_type TEXT;
    matching_count INTEGER;
    orphaned_count INTEGER;
BEGIN
    -- Generate a constraint name if not provided
    IF constraint_name IS NULL THEN
        actual_constraint_name := source_table || '_' || source_column || '_fkey';
    ELSE
        actual_constraint_name := constraint_name;
    END IF;
    
    -- Check if constraint already exists
    PERFORM 1
    FROM information_schema.table_constraints
    WHERE constraint_name = actual_constraint_name
    AND table_name = source_table
    AND constraint_type = 'FOREIGN KEY';
    
    IF FOUND THEN
        RAISE NOTICE 'Foreign key constraint % already exists', actual_constraint_name;
        RETURN;
    END IF;
    
    -- Get column data types to ensure they match
    EXECUTE format('SELECT data_type FROM information_schema.columns 
                   WHERE table_name = %L AND column_name = %L', 
                   source_table, source_column)
    INTO source_data_type;
    
    EXECUTE format('SELECT data_type FROM information_schema.columns 
                   WHERE table_name = %L AND column_name = %L', 
                   target_table, target_column)
    INTO target_data_type;
    
    IF source_data_type IS NULL THEN
        RAISE EXCEPTION 'Source column %.% does not exist', source_table, source_column;
    END IF;
    
    IF target_data_type IS NULL THEN
        RAISE EXCEPTION 'Target column %.% does not exist', target_table, target_column;
    END IF;
    
    IF source_data_type <> target_data_type THEN
        RAISE EXCEPTION 'Data type mismatch: %.% is %, but %.% is %', 
                        source_table, source_column, source_data_type,
                        target_table, target_column, target_data_type;
    END IF;
    
    -- Check for orphaned records
    EXECUTE format('
        SELECT COUNT(*) 
        FROM %I 
        WHERE %I IS NOT NULL 
        AND NOT EXISTS (
            SELECT 1 
            FROM %I 
            WHERE %I = %I.%I
        )', 
        source_table, source_column, 
        target_table, target_column, source_table, source_column
    ) INTO orphaned_count;
    
    IF orphaned_count > 0 THEN
        RAISE WARNING 'Found % orphaned records in %.% that do not exist in %.%',
                     orphaned_count, source_table, source_column, target_table, target_column;
        RETURN;
    END IF;
    
    -- Add the foreign key constraint
    EXECUTE format('ALTER TABLE %I ADD CONSTRAINT %I 
                   FOREIGN KEY (%I) REFERENCES %I(%I) ON DELETE CASCADE',
                   source_table, actual_constraint_name,
                   source_column, target_table, target_column);
    
    RAISE NOTICE 'Added foreign key constraint % on %.% referencing %.%',
                 actual_constraint_name, source_table, source_column, target_table, target_column;
END;
$$ LANGUAGE plpgsql;

-- Now, restore the foreign key constraints
-- First, update RBAC tables
SELECT safe_add_foreign_key('team_members', 'team_id', 'teams', 'id');
SELECT safe_add_foreign_key('user_profile', 'team_id', 'teams', 'id');
SELECT safe_add_foreign_key('projects', 'team_id', 'teams', 'id');

-- Then, update the core molecule tables
SELECT safe_add_foreign_key('molecular_properties', 'molecule_id', 'molecules', 'id');
SELECT safe_add_foreign_key('molecular_properties', 'property_type_id', 'property_types', 'id');

SELECT safe_add_foreign_key('mixtures', 'project_id', 'projects', 'id');
SELECT safe_add_foreign_key('mixture_components', 'mixture_id', 'mixtures', 'id');
SELECT safe_add_foreign_key('mixture_components', 'molecule_id', 'molecules', 'id');

-- Update experimental data tables
SELECT safe_add_foreign_key('experiment_properties', 'experiment_id', 'molecule_experiments', 'id');
SELECT safe_add_foreign_key('molecule_experiments', 'molecule_id', 'molecules', 'id');
SELECT safe_add_foreign_key('molecule_proteins', 'molecule_id', 'molecules', 'id');
SELECT safe_add_foreign_key('molecule_proteins', 'protein_id', 'proteins', 'id');
SELECT safe_add_foreign_key('predictions', 'molecule_id', 'molecules', 'id');
SELECT safe_add_foreign_key('predictions', 'mixture_id', 'mixtures', 'id');

-- Connect consolidated_molecules to molecules
SELECT safe_add_foreign_key('consolidated_molecules', 'primary_molecule_id', 'molecules', 'id');

-- Add any missing indexes on the foreign key columns
CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id 
ON molecular_properties(molecule_id);

CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_type_id 
ON molecular_properties(property_type_id);

CREATE INDEX IF NOT EXISTS idx_mixture_components_mixture_id 
ON mixture_components(mixture_id);

CREATE INDEX IF NOT EXISTS idx_mixture_components_molecule_id 
ON mixture_components(molecule_id);

CREATE INDEX IF NOT EXISTS idx_experiment_properties_experiment_id 
ON experiment_properties(experiment_id);

CREATE INDEX IF NOT EXISTS idx_molecule_experiments_molecule_id 
ON molecule_experiments(molecule_id);

CREATE INDEX IF NOT EXISTS idx_molecule_proteins_molecule_id 
ON molecule_proteins(molecule_id);

CREATE INDEX IF NOT EXISTS idx_molecule_proteins_protein_id 
ON molecule_proteins(protein_id);

CREATE INDEX IF NOT EXISTS idx_predictions_molecule_id 
ON predictions(molecule_id);

CREATE INDEX IF NOT EXISTS idx_predictions_mixture_id 
ON predictions(mixture_id);

CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_primary_molecule_id 
ON consolidated_molecules(primary_molecule_id);

-- Drop the helper function when we're done
DROP FUNCTION IF EXISTS safe_add_foreign_key(TEXT, TEXT, TEXT, TEXT, TEXT);