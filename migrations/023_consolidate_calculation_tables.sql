-- Migration: Consolidate calculation_method into calculation_methods
-- This migration merges data from calculation_method into calculation_methods
-- and updates foreign key references to point to calculation_methods

-- Step 1: Create a savepoint to make this transaction atomic
SAVEPOINT consolidate_calculation_tables;

-- Step 2: Insert data from calculation_method into calculation_methods if not already there
DO $$
DECLARE
    method_record RECORD;
BEGIN
    FOR method_record IN (SELECT * FROM calculation_method) LOOP
        -- Check if method already exists in calculation_methods
        IF NOT EXISTS (SELECT 1 FROM calculation_methods WHERE id = method_record.id) THEN
            INSERT INTO calculation_methods (
                id, 
                name, 
                description, 
                version, 
                created_at, 
                updated_at, 
                created_by, 
                method_type, 
                reference, 
                parameters
            ) VALUES (
                method_record.id,
                method_record.name,
                method_record.description,
                method_record.version,
                COALESCE(method_record.created_at, NOW()),
                COALESCE(method_record.updated_at, NOW()),
                NULL, -- created_by is new in calculation_methods
                method_record.method_type,
                method_record.reference,
                method_record.parameters
            );
            
            RAISE NOTICE 'Migrated calculation method: %', method_record.name;
        ELSE
            RAISE NOTICE 'Calculation method already exists: %', method_record.name;
        END IF;
    END LOOP;
END $$;

-- Step 3: Update foreign key references in predictions table (if they reference calculation_method)
DO $$
DECLARE
    fk_constraint_name TEXT;
    fk_def TEXT;
BEGIN
    -- Find foreign key constraints referencing calculation_method
    SELECT 
        tc.constraint_name,
        pg_get_constraintdef(tc.oid) AS constraint_def
    INTO
        fk_constraint_name,
        fk_def
    FROM 
        pg_constraint tc
    JOIN 
        pg_class rc ON tc.conrelid = rc.oid
    JOIN 
        pg_namespace rn ON rc.relnamespace = rn.oid
    JOIN 
        pg_class fc ON tc.confrelid = fc.oid
    JOIN 
        pg_namespace fn ON fc.relnamespace = fn.oid
    WHERE 
        tc.contype = 'f'
        AND rn.nspname = 'public'
        AND rc.relname = 'predictions'
        AND fn.nspname = 'public'
        AND fc.relname = 'calculation_method';
    
    -- If a foreign key constraint was found, drop it and create a new one referencing calculation_methods
    IF fk_constraint_name IS NOT NULL THEN
        -- Extract column name from constraint definition
        -- The constraint def will look like: FOREIGN KEY (calculation_method_id) REFERENCES public.calculation_method(id)
        DECLARE
            column_name TEXT;
        BEGIN
            column_name := substring(fk_def FROM 'FOREIGN KEY \(([^\)]+)\)');
            
            -- Drop the old constraint
            EXECUTE 'ALTER TABLE predictions DROP CONSTRAINT ' || fk_constraint_name;
            RAISE NOTICE 'Dropped foreign key constraint: %', fk_constraint_name;
            
            -- Create new constraint referencing calculation_methods
            EXECUTE 'ALTER TABLE predictions ADD CONSTRAINT predictions_' || column_name || '_fkey ' ||
                    'FOREIGN KEY (' || column_name || ') REFERENCES calculation_methods(id)';
            RAISE NOTICE 'Created new foreign key constraint referencing calculation_methods';
        END;
    ELSE
        RAISE NOTICE 'No foreign key constraint found referencing calculation_method';
    END IF;
END $$;

-- Step 4: Verify that all data has been migrated correctly
DO $$
DECLARE
    methods_count INT;
    missing_count INT;
BEGIN
    -- Check count of methods in calculation_method
    SELECT COUNT(*) INTO methods_count FROM calculation_method;
    
    -- Check count of methods from calculation_method not in calculation_methods
    SELECT COUNT(*) 
    INTO missing_count 
    FROM calculation_method cm
    WHERE NOT EXISTS (
        SELECT 1 FROM calculation_methods cms WHERE cms.id = cm.id
    );
    
    -- If all data has been migrated, we can drop the original table
    IF missing_count = 0 THEN
        -- Only drop if there's data to migrate (to avoid dropping an empty table unnecessarily)
        IF methods_count > 0 THEN
            -- Create a backup of the table data just in case
            CREATE TABLE IF NOT EXISTS calculation_method_backup AS SELECT * FROM calculation_method;
            
            -- Drop the original table
            DROP TABLE calculation_method;
            
            RAISE NOTICE 'Successfully dropped calculation_method table after migrating % records', methods_count;
        ELSE
            RAISE NOTICE 'calculation_method table is empty, safe to drop';
            DROP TABLE calculation_method;
        END IF;
    ELSE
        -- If migration is not complete, raise an exception to roll back
        RAISE EXCEPTION 'Migration incomplete: % records not migrated from calculation_method', missing_count;
    END IF;
END $$;

-- Step 5: Commit the transaction if all steps succeed (implicit, as we're not in an explicit transaction)
-- If any step fails, PostgreSQL will roll back to the savepoint automatically
RELEASE SAVEPOINT consolidate_calculation_tables;

-- Step 6: Record the migration in the migrations table
INSERT INTO migrations (name, created_at)
VALUES ('023_consolidate_calculation_tables', NOW());