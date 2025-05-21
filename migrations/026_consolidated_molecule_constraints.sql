-- Check if table exists
DO $$ 
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'consolidated_molecules') THEN
        -- First validate the current data
        IF EXISTS (
            SELECT 1
            FROM consolidated_molecules 
            WHERE molecule_status NOT IN ('original', 'primary', 'duplicate')
        ) THEN
            RAISE EXCEPTION 'Invalid molecule_status values found in consolidated_molecules table. Please fix before applying constraint.';
        END IF;

        -- Add CHECK constraint for molecule_status
        IF NOT EXISTS (
            SELECT 1
            FROM information_schema.table_constraints 
            WHERE table_name = 'consolidated_molecules' 
            AND constraint_type = 'CHECK'
            AND constraint_name LIKE '%valid_status%'
        ) THEN
            ALTER TABLE consolidated_molecules 
            ADD CONSTRAINT valid_status 
            CHECK (molecule_status IN ('original', 'primary', 'duplicate'));
            
            RAISE NOTICE 'Added CHECK constraint for molecule_status column.';
        ELSE
            RAISE NOTICE 'CHECK constraint for molecule_status already exists.';
        END IF;

        -- Create index for inchikey if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_inchikey%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_inchikey 
            ON public.consolidated_molecules(inchikey);
            
            RAISE NOTICE 'Created index for inchikey column.';
        ELSE
            RAISE NOTICE 'Index for inchikey column already exists.';
        END IF;

        -- Create index for primary_molecule_id if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_primary_id%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_primary_id 
            ON public.consolidated_molecules(primary_molecule_id);
            
            RAISE NOTICE 'Created index for primary_molecule_id column.';
        ELSE
            RAISE NOTICE 'Index for primary_molecule_id column already exists.';
        END IF;

        -- Create index for molecule_status if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_status%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_status 
            ON public.consolidated_molecules(molecule_status);
            
            RAISE NOTICE 'Created index for molecule_status column.';
        ELSE
            RAISE NOTICE 'Index for molecule_status column already exists.';
        END IF;
    ELSE
        RAISE NOTICE 'consolidated_molecules table does not exist. Skipping constraint and index creation.';
    END IF;
END $$;
