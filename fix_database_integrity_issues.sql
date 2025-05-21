-- fix_database_integrity_issues.sql
--
-- This script addresses common data integrity issues identified by the database verification tool.
-- Run this script to fix data integrity problems in the CryoProtect database.
--
-- Usage: 
--   export $(grep -v '^#' .env | xargs) && PGPASSWORD=$SUPABASE_DB_PASSWORD psql -h $SUPABASE_DB_HOST -p $SUPABASE_DB_PORT -U $SUPABASE_DB_USER -d $SUPABASE_DB_NAME -f fix_database_integrity_issues.sql

-- Begin transaction
BEGIN;

-- Set search path
SET search_path TO public;

-- Create temporary function to log repair operations
CREATE OR REPLACE FUNCTION temp_log_repair(operation text, details text) RETURNS void AS $$
BEGIN
    RAISE NOTICE '% - %', operation, details;
END;
$$ LANGUAGE plpgsql;

-- 1. Fix NULL values in required fields
DO $$
DECLARE
    tbl RECORD;
    col RECORD;
    null_count INT;
    update_sql TEXT;
BEGIN
    RAISE NOTICE 'Checking for NULL values in required fields...';
    
    FOR tbl IN 
        SELECT table_name FROM information_schema.tables 
        WHERE table_schema = 'public' 
        AND table_type = 'BASE TABLE'
    LOOP
        FOR col IN 
            SELECT column_name, data_type 
            FROM information_schema.columns 
            WHERE table_schema = 'public' 
            AND table_name = tbl.table_name 
            AND is_nullable = 'NO'
            AND column_default IS NULL
        LOOP
            EXECUTE format('SELECT COUNT(*) FROM %I WHERE %I IS NULL', 
                          tbl.table_name, col.column_name) 
            INTO null_count;
            
            IF null_count > 0 THEN
                PERFORM temp_log_repair(
                    'REQUIRED FIELD REPAIR', 
                    format('Table %s has %s NULL values in required column %s', 
                          tbl.table_name, null_count, col.column_name)
                );
                
                -- Generate appropriate default value based on data type
                CASE 
                    WHEN col.data_type LIKE '%char%' OR col.data_type = 'text' THEN
                        update_sql := format('UPDATE %I SET %I = ''[Fixed - was NULL]'' WHERE %I IS NULL', 
                                          tbl.table_name, col.column_name, col.column_name);
                    WHEN col.data_type = 'uuid' THEN
                        update_sql := format('UPDATE %I SET %I = uuid_generate_v4() WHERE %I IS NULL', 
                                          tbl.table_name, col.column_name, col.column_name);
                    WHEN col.data_type IN ('integer', 'bigint', 'smallint', 'numeric', 'decimal') THEN
                        update_sql := format('UPDATE %I SET %I = 0 WHERE %I IS NULL', 
                                          tbl.table_name, col.column_name, col.column_name);
                    WHEN col.data_type = 'boolean' THEN
                        update_sql := format('UPDATE %I SET %I = FALSE WHERE %I IS NULL', 
                                          tbl.table_name, col.column_name, col.column_name);
                    WHEN col.data_type LIKE '%timestamp%' THEN
                        update_sql := format('UPDATE %I SET %I = now() WHERE %I IS NULL', 
                                          tbl.table_name, col.column_name, col.column_name);
                    ELSE
                        -- Skip unknown types
                        PERFORM temp_log_repair(
                            'SKIPPED', 
                            format('Unknown data type %s for column %s.%s', 
                                  col.data_type, tbl.table_name, col.column_name)
                        );
                        CONTINUE;
                END CASE;
                
                EXECUTE update_sql;
                PERFORM temp_log_repair(
                    'FIXED', 
                    format('Updated %s NULL values in %s.%s', 
                          null_count, tbl.table_name, col.column_name)
                );
            END IF;
        END LOOP;
    END LOOP;
END;
$$;

-- 2. Address foreign key violations
DO $$
DECLARE
    fk_record RECORD;
    ref_record RECORD;
    violation_count INT;
    fix_sql TEXT;
BEGIN
    RAISE NOTICE 'Checking for foreign key violations...';
    
    FOR fk_record IN 
        SELECT
            tc.table_name,
            kcu.column_name,
            ccu.table_name AS ref_table_name,
            ccu.column_name AS ref_column_name,
            tc.constraint_name
        FROM
            information_schema.table_constraints AS tc
            JOIN information_schema.key_column_usage AS kcu
              ON tc.constraint_name = kcu.constraint_name
            JOIN information_schema.constraint_column_usage AS ccu
              ON ccu.constraint_name = tc.constraint_name
        WHERE
            tc.constraint_type = 'FOREIGN KEY'
            AND tc.table_schema = 'public'
    LOOP
        -- Check for existing violations
        EXECUTE format(
            'SELECT COUNT(*) FROM %I t1 
             WHERE t1.%I IS NOT NULL AND 
             NOT EXISTS (SELECT 1 FROM %I t2 WHERE t2.%I = t1.%I)',
            fk_record.table_name, 
            fk_record.column_name,
            fk_record.ref_table_name,
            fk_record.ref_column_name,
            fk_record.column_name
        ) INTO violation_count;
        
        IF violation_count > 0 THEN
            PERFORM temp_log_repair(
                'FOREIGN KEY VIOLATION', 
                format('Table %s has %s rows with invalid foreign keys in column %s referencing %s.%s',
                      fk_record.table_name, violation_count, fk_record.column_name, 
                      fk_record.ref_table_name, fk_record.ref_column_name)
            );
            
            -- Try to find a valid reference value to use as replacement
            EXECUTE format(
                'SELECT %I FROM %I LIMIT 1',
                fk_record.ref_column_name,
                fk_record.ref_table_name
            ) INTO ref_record;
            
            IF ref_record IS NOT NULL THEN
                -- Option 1: Set to NULL if column allows nulls
                EXECUTE format(
                    'SELECT COUNT(*) FROM information_schema.columns 
                     WHERE table_schema = ''public'' 
                     AND table_name = %L 
                     AND column_name = %L 
                     AND is_nullable = ''YES''',
                    fk_record.table_name, 
                    fk_record.column_name
                ) INTO violation_count;
                
                IF violation_count > 0 THEN
                    -- Set invalid references to NULL
                    fix_sql := format(
                        'UPDATE %I t1 SET %I = NULL 
                         WHERE t1.%I IS NOT NULL AND 
                         NOT EXISTS (SELECT 1 FROM %I t2 WHERE t2.%I = t1.%I)',
                        fk_record.table_name, 
                        fk_record.column_name,
                        fk_record.column_name,
                        fk_record.ref_table_name,
                        fk_record.ref_column_name,
                        fk_record.column_name
                    );
                    
                    EXECUTE fix_sql;
                    PERFORM temp_log_repair(
                        'FIXED', 
                        format('Set invalid foreign keys to NULL in %s.%s',
                              fk_record.table_name, fk_record.column_name)
                    );
                ELSE
                    -- Option 2: For required FK, update to a valid value
                    EXECUTE format(
                        'UPDATE %I t1 SET %I = %L
                         WHERE t1.%I IS NOT NULL AND 
                         NOT EXISTS (SELECT 1 FROM %I t2 WHERE t2.%I = t1.%I)',
                        fk_record.table_name, 
                        fk_record.column_name,
                        ref_record,
                        fk_record.column_name,
                        fk_record.ref_table_name,
                        fk_record.ref_column_name,
                        fk_record.column_name
                    );
                    
                    PERFORM temp_log_repair(
                        'FIXED', 
                        format('Updated invalid foreign keys in %s.%s to valid reference %s',
                              fk_record.table_name, fk_record.column_name, ref_record)
                    );
                END IF;
            ELSE
                -- Option 3: No valid references found, handle case by case
                PERFORM temp_log_repair(
                    'MANUAL FIX NEEDED', 
                    format('No valid references found in %s.%s to fix violations in %s.%s',
                          fk_record.ref_table_name, fk_record.ref_column_name,
                          fk_record.table_name, fk_record.column_name)
                );
            END IF;
        END IF;
    END LOOP;
END;
$$;

-- 3. Fix duplicate records
DO $$
DECLARE
    tbl RECORD;
    duplicate_count INT;
BEGIN
    RAISE NOTICE 'Checking for duplicate molecules by InChIKey...';
    
    -- Check for duplicate molecules by InChIKey
    SELECT COUNT(*) INTO duplicate_count
    FROM (
        SELECT inchikey, COUNT(*) 
        FROM molecules 
        WHERE inchikey IS NOT NULL 
        GROUP BY inchikey 
        HAVING COUNT(*) > 1
    ) AS duplicates;
    
    IF duplicate_count > 0 THEN
        PERFORM temp_log_repair(
            'DUPLICATE MOLECULES', 
            format('Found %s InChIKeys with multiple molecule entries',
                  duplicate_count)
        );
        
        -- Create a temporary table to track duplicates
        CREATE TEMP TABLE duplicate_molecules AS
        SELECT array_agg(id) AS duplicate_ids, inchikey, MIN(id) AS keep_id
        FROM molecules
        WHERE inchikey IS NOT NULL
        GROUP BY inchikey
        HAVING COUNT(*) > 1;
        
        -- Update references to use the keeper ID
        FOR tbl IN 
            SELECT
                tc.table_name,
                kcu.column_name
            FROM
                information_schema.table_constraints AS tc
                JOIN information_schema.key_column_usage AS kcu
                  ON tc.constraint_name = kcu.constraint_name
                JOIN information_schema.constraint_column_usage AS ccu
                  ON ccu.constraint_name = tc.constraint_name
            WHERE
                tc.constraint_type = 'FOREIGN KEY'
                AND tc.table_schema = 'public'
                AND ccu.table_name = 'molecules'
                AND ccu.column_name = 'id'
        LOOP
            EXECUTE format(
                'WITH updates AS (
                    SELECT dm.keep_id, unnest(dm.duplicate_ids) AS dup_id
                    FROM duplicate_molecules dm
                )
                UPDATE %I t
                SET %I = u.keep_id
                FROM updates u
                WHERE t.%I = u.dup_id
                AND t.%I != u.keep_id',
                tbl.table_name,
                tbl.column_name,
                tbl.column_name,
                tbl.column_name
            );
            
            PERFORM temp_log_repair(
                'FIXED REFERENCES', 
                format('Updated references in %s.%s to use primary molecule IDs',
                      tbl.table_name, tbl.column_name)
            );
        END LOOP;
        
        -- Delete duplicate molecules (this will only work if all references are updated)
        BEGIN
            EXECUTE 
                'DELETE FROM molecules m
                 USING duplicate_molecules dm
                 WHERE m.id = ANY(dm.duplicate_ids)
                 AND m.id != dm.keep_id';
                 
            PERFORM temp_log_repair(
                'REMOVED DUPLICATES', 
                'Successfully deleted duplicate molecule records'
            );
        EXCEPTION WHEN OTHERS THEN
            PERFORM temp_log_repair(
                'FAILED TO REMOVE DUPLICATES', 
                'Error removing duplicate molecules - foreign key constraints still exist'
            );
        END;
        
        -- Clean up
        DROP TABLE duplicate_molecules;
    END IF;
END;
$$;

-- 4. Fix chemical data issues
DO $$
BEGIN
    RAISE NOTICE 'Fixing chemical data issues...';
    
    -- Fix molecules with invalid SMILES
    UPDATE molecules
    SET smiles = NULL
    WHERE smiles IS NOT NULL 
    AND smiles NOT SIMILAR TO '[A-Za-z0-9@+\\-\\[\\]\\(\\)\\\\\\/%=#$]+';
    
    PERFORM temp_log_repair(
        'FIXED SMILES DATA', 
        'Set invalid SMILES strings to NULL'
    );
    
    -- Fix molecules with invalid InChI
    UPDATE molecules
    SET inchi = NULL
    WHERE inchi IS NOT NULL 
    AND inchi NOT LIKE 'InChI=1%';
    
    PERFORM temp_log_repair(
        'FIXED INCHI DATA', 
        'Set invalid InChI strings to NULL'
    );
    
    -- Fix molecules with invalid InChIKey
    UPDATE molecules
    SET inchikey = NULL
    WHERE inchikey IS NOT NULL 
    AND inchikey NOT SIMILAR TO '[A-Z]{14}-[A-Z]{10}-[A-Z]';
    
    PERFORM temp_log_repair(
        'FIXED INCHIKEY DATA', 
        'Set invalid InChIKey strings to NULL'
    );
    
    -- Fix molecules with suspicious molecular weight
    UPDATE molecules
    SET molecular_weight = NULL
    WHERE molecular_weight IS NOT NULL 
    AND (molecular_weight <= 0 OR molecular_weight > 5000);
    
    PERFORM temp_log_repair(
        'FIXED MOLECULAR WEIGHT', 
        'Set suspicious molecular weight values to NULL'
    );
end;
$$;

-- 5. Fix mixture issues
DO $$
BEGIN
    RAISE NOTICE 'Fixing mixture issues...';
    
    -- Standardize concentration units
    UPDATE mixture_components
    SET concentration_unit = '%'
    WHERE concentration_unit IN ('percent', 'pct');
    
    PERFORM temp_log_repair(
        'FIXED CONCENTRATION UNITS', 
        'Standardized concentration units to use %'
    );
    
    -- Fix percentage mixtures that don't sum to 100%
    WITH mixture_totals AS (
        SELECT 
            mixture_id,
            SUM(concentration) AS total
        FROM mixture_components
        WHERE concentration_unit = '%'
        GROUP BY mixture_id
        HAVING SUM(concentration) != 100
    )
    UPDATE mixture_components mc
    SET concentration = concentration * (100.0 / mt.total)
    FROM mixture_totals mt
    WHERE mc.mixture_id = mt.mixture_id
    AND mc.concentration_unit = '%';
    
    PERFORM temp_log_repair(
        'FIXED PERCENTAGE MIXTURES', 
        'Adjusted component concentrations to sum to 100%'
    );
end;
$$;

-- 6. Fix user/team relationships
DO $$
DECLARE
    default_team_id UUID;
BEGIN
    RAISE NOTICE 'Fixing user and team relationships...';
    
    -- Find or create a default team
    SELECT id INTO default_team_id FROM teams LIMIT 1;
    
    IF default_team_id IS NULL THEN
        INSERT INTO teams (name, description)
        VALUES ('Default Team', 'Automatically created default team')
        RETURNING id INTO default_team_id;
        
        PERFORM temp_log_repair(
            'CREATED DEFAULT TEAM', 
            format('Created a default team with ID %s', default_team_id)
        );
    END IF;
    
    -- Assign users without team to the default team
    UPDATE user_profile
    SET team_id = default_team_id
    WHERE team_id IS NULL;
    
    PERFORM temp_log_repair(
        'FIXED USER TEAMS', 
        'Assigned users without a team to the default team'
    );
    
    -- Assign projects without team to the default team
    UPDATE projects
    SET team_id = default_team_id
    WHERE team_id IS NULL;
    
    PERFORM temp_log_repair(
        'FIXED PROJECT TEAMS', 
        'Assigned projects without a team to the default team'
    );
    
    -- Ensure each team has at least one member
    INSERT INTO team_members (team_id, user_id, role)
    SELECT 
        t.id,
        (SELECT id FROM user_profile LIMIT 1),
        'member'
    FROM teams t
    LEFT JOIN team_members tm ON t.id = tm.team_id
    WHERE tm.id IS NULL;
    
    PERFORM temp_log_repair(
        'FIXED TEAM MEMBERS', 
        'Added at least one member to each team'
    );
END;
$$;

-- Drop the temporary function
DROP FUNCTION IF EXISTS temp_log_repair(text, text);

-- Commit all changes
COMMIT;

-- Final verification
DO $$
BEGIN
    RAISE NOTICE 'Database integrity fixes applied successfully.';
    RAISE NOTICE 'Run the verification script again to check the results.';
END;
$$;