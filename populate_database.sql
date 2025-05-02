-- CryoProtect v2 - Database Population Script
-- This script populates the CryoProtect database with initial data
-- Execute this script to populate all tables in the correct order

-- Begin transaction
BEGIN;

-- 1. Populate molecules table
\echo 'Populating molecules table...'
\i populate_molecules.sql

-- 2. Populate molecular_properties table
\echo 'Populating molecular_properties table...'
\i populate_molecular_properties.sql

-- 3. Populate mixtures table
\echo 'Populating mixtures table...'
\i populate_mixtures.sql

-- 4. Populate mixture_components table
\echo 'Populating mixture_components table...'
\i populate_mixture_components.sql

-- 5. Populate experiments table
\echo 'Populating experiments table...'
\i populate_experiments.sql

-- 6. Populate experiment_properties table
\echo 'Populating experiment_properties table...'
\i populate_experiment_properties.sql

-- 7. Populate predictions table
\echo 'Populating predictions table...'
\i populate_predictions.sql

-- Commit transaction
COMMIT;

-- Verify data was inserted correctly
\echo 'Verifying data was inserted correctly...'

-- Count molecules
SELECT COUNT(*) AS molecule_count FROM molecules;

-- Count molecular properties
SELECT COUNT(*) AS molecular_property_count FROM molecular_properties;

-- Count mixtures
SELECT COUNT(*) AS mixture_count FROM mixtures;

-- Count mixture components
SELECT COUNT(*) AS mixture_component_count FROM mixture_components;

-- Count experiments
SELECT COUNT(*) AS experiment_count FROM experiments;

-- Count predictions
SELECT COUNT(*) AS prediction_count FROM predictions;

\echo 'Database population complete!'