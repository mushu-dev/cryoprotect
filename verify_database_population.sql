-- CryoProtect v2 - Database Population Verification Script
-- This script verifies that the database was populated correctly

-- Count records in each table
\echo 'Counting records in each table...'
SELECT 'molecules' AS table_name, COUNT(*) AS record_count FROM molecules
UNION ALL
SELECT 'molecular_properties' AS table_name, COUNT(*) AS record_count FROM molecular_properties
UNION ALL
SELECT 'mixtures' AS table_name, COUNT(*) AS record_count FROM mixtures
UNION ALL
SELECT 'mixture_components' AS table_name, COUNT(*) AS record_count FROM mixture_components
UNION ALL
SELECT 'experiments' AS table_name, COUNT(*) AS record_count FROM experiments
UNION ALL
SELECT 'predictions' AS table_name, COUNT(*) AS record_count FROM predictions
ORDER BY table_name;

-- Verify molecules
\echo '\nVerifying molecules...'
SELECT name, formula, molecular_weight FROM molecules ORDER BY name;

-- Verify molecular properties
\echo '\nVerifying molecular properties...'
SELECT m.name AS molecule_name, pt.name AS property_name, mp.numeric_value, mp.unit
FROM molecular_properties mp
JOIN molecules m ON mp.molecule_id = m.id
JOIN property_types pt ON mp.property_type_id = pt.id
ORDER BY m.name, pt.name
LIMIT 10;

-- Verify mixtures
\echo '\nVerifying mixtures...'
SELECT name, description FROM mixtures ORDER BY name;

-- Verify mixture components
\echo '\nVerifying mixture components...'
SELECT mix.name AS mixture_name, mol.name AS molecule_name, mc.concentration, mc.concentration_unit, mc.role
FROM mixture_components mc
JOIN mixtures mix ON mc.mixture_id = mix.id
JOIN molecules mol ON mc.molecule_id = mol.id
ORDER BY mix.name, mol.name;

-- Verify experiments
\echo '\nVerifying experiments...'
SELECT e.id, mix.name AS mixture_name, pt.name AS property_type, e.experimental_conditions
FROM experiments e
JOIN mixtures mix ON e.mixture_id = mix.id
JOIN property_types pt ON e.property_type_id = pt.id
ORDER BY mix.name;

-- Verify predictions
\echo '\nVerifying predictions...'
SELECT m.name AS molecule_name, pt.name AS property_name, cm.name AS method_name, 
       p.numeric_value, p.unit, p.confidence
FROM predictions p
JOIN molecules m ON p.molecule_id = m.id
JOIN property_types pt ON p.property_type_id = pt.id
JOIN calculation_methods cm ON p.calculation_method_id = cm.id
ORDER BY m.name, pt.name
LIMIT 10;

-- Verify relationships
\echo '\nVerifying relationships...'
-- Check that all mixture components reference valid mixtures and molecules
SELECT COUNT(*) AS invalid_mixture_components
FROM mixture_components mc
WHERE NOT EXISTS (SELECT 1 FROM mixtures m WHERE m.id = mc.mixture_id)
   OR NOT EXISTS (SELECT 1 FROM molecules m WHERE m.id = mc.molecule_id);

-- Check that all molecular properties reference valid molecules
SELECT COUNT(*) AS invalid_molecular_properties
FROM molecular_properties mp
WHERE NOT EXISTS (SELECT 1 FROM molecules m WHERE m.id = mp.molecule_id);

-- Check that all experiments reference valid mixtures
SELECT COUNT(*) AS invalid_experiments
FROM experiments e
WHERE NOT EXISTS (SELECT 1 FROM mixtures m WHERE m.id = e.mixture_id);

-- Check that all predictions reference valid molecules
SELECT COUNT(*) AS invalid_predictions
FROM predictions p
WHERE NOT EXISTS (SELECT 1 FROM molecules m WHERE m.id = p.molecule_id);

\echo '\nVerification complete!'