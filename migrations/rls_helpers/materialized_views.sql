-- Materialized views for frequently accessed data
-- These views cache common query results and improve performance
-- They need to be refreshed periodically

-- 1. Public Molecules Summary View
CREATE MATERIALIZED VIEW IF NOT EXISTS public_molecules_summary AS
SELECT 
    m.id, 
    m.name, 
    m.molecular_formula, 
    m.smiles, 
    m.cid, 
    m.pubchem_link,
    m.is_public,
    COUNT(mp.id) AS property_count
FROM 
    molecules m
LEFT JOIN 
    molecular_properties mp ON m.id = mp.molecule_id
WHERE 
    m.is_public = true
GROUP BY 
    m.id;

-- Create index on the materialized view
CREATE UNIQUE INDEX IF NOT EXISTS idx_public_molecules_summary_id 
ON public_molecules_summary(id);

-- Create indexes for common access patterns
CREATE INDEX IF NOT EXISTS idx_public_molecules_summary_name
ON public_molecules_summary(name);

CREATE INDEX IF NOT EXISTS idx_public_molecules_summary_cid
ON public_molecules_summary(cid);

-- 2. Public Molecular Properties View
CREATE MATERIALIZED VIEW IF NOT EXISTS public_molecular_properties AS
SELECT 
    mp.id,
    mp.molecule_id,
    mp.property_name,
    mp.property_value,
    mp.property_unit,
    mp.property_type,
    mp.calculation_method,
    mp.reference_doi
FROM 
    molecular_properties mp
JOIN 
    molecules m ON mp.molecule_id = m.id
WHERE 
    m.is_public = true
    AND (mp.sensitivity_level IS NULL OR mp.sensitivity_level <= 'low');

-- Create indexes on the materialized view
CREATE INDEX IF NOT EXISTS idx_public_molecular_properties_molecule_id 
ON public_molecular_properties(molecule_id);

CREATE INDEX IF NOT EXISTS idx_public_molecular_properties_property_name 
ON public_molecular_properties(property_name);

-- 3. Public Mixtures Summary View
CREATE MATERIALIZED VIEW IF NOT EXISTS public_mixtures_summary AS
SELECT 
    m.id, 
    m.name, 
    m.description,
    m.created_at,
    m.is_public,
    COUNT(mc.id) AS component_count
FROM 
    mixtures m
LEFT JOIN 
    mixture_components mc ON m.id = mc.mixture_id
WHERE 
    m.is_public = true
GROUP BY 
    m.id;

-- Create indexes on the materialized view
CREATE UNIQUE INDEX IF NOT EXISTS idx_public_mixtures_summary_id 
ON public_mixtures_summary(id);

CREATE INDEX IF NOT EXISTS idx_public_mixtures_summary_name
ON public_mixtures_summary(name);

-- 4. Popular Molecules View (frequently accessed molecules)
CREATE MATERIALIZED VIEW IF NOT EXISTS popular_molecules AS
SELECT 
    m.id,
    m.name,
    m.molecular_formula,
    m.smiles,
    m.cid,
    m.pubchem_link,
    m.is_public,
    COUNT(DISTINCT e.id) as experiment_count,
    COUNT(DISTINCT mx.id) as mixture_count,
    COUNT(DISTINCT mp.id) as property_count
FROM 
    molecules m
LEFT JOIN 
    experiment_molecules em ON m.id = em.molecule_id
LEFT JOIN 
    experiments e ON em.experiment_id = e.id
LEFT JOIN 
    mixture_components mc ON m.id = mc.molecule_id
LEFT JOIN 
    mixtures mx ON mc.mixture_id = mx.id
LEFT JOIN 
    molecular_properties mp ON m.id = mp.molecule_id
WHERE 
    m.is_public = true
GROUP BY 
    m.id
HAVING 
    COUNT(DISTINCT e.id) + COUNT(DISTINCT mx.id) > 5
ORDER BY 
    COUNT(DISTINCT e.id) + COUNT(DISTINCT mx.id) DESC
LIMIT 500;

-- Create indexes on the materialized view
CREATE UNIQUE INDEX IF NOT EXISTS idx_popular_molecules_id
ON popular_molecules(id);

-- 5. Property Statistics View (aggregate statistics on molecular properties)
CREATE MATERIALIZED VIEW IF NOT EXISTS property_statistics AS
SELECT 
    property_name,
    property_type,
    COUNT(*) as count,
    MIN(CASE WHEN property_value ~ '^-?[0-9]+(\.[0-9]+)?$' THEN property_value::numeric ELSE NULL END) as min_numeric_value,
    MAX(CASE WHEN property_value ~ '^-?[0-9]+(\.[0-9]+)?$' THEN property_value::numeric ELSE NULL END) as max_numeric_value,
    AVG(CASE WHEN property_value ~ '^-?[0-9]+(\.[0-9]+)?$' THEN property_value::numeric ELSE NULL END) as avg_numeric_value,
    array_agg(DISTINCT property_unit) FILTER (WHERE property_unit IS NOT NULL) as units,
    COUNT(DISTINCT molecule_id) as molecule_count
FROM 
    molecular_properties
GROUP BY 
    property_name, property_type;

-- Create indexes on the materialized view
CREATE UNIQUE INDEX IF NOT EXISTS idx_property_statistics_name_type
ON property_statistics(property_name, property_type);

-- 6. User Favorites Statistics (helps with personalized recommendations)
CREATE MATERIALIZED VIEW IF NOT EXISTS user_favorites_statistics AS
SELECT
    up.auth_user_id,
    COUNT(DISTINCT m.id) as saved_molecules,
    COUNT(DISTINCT mx.id) as saved_mixtures,
    COUNT(DISTINCT e.id) as saved_experiments,
    COUNT(DISTINCT mc.molecule_id) as component_molecules
FROM
    user_profile up
LEFT JOIN
    user_saved_molecules usm ON up.id = usm.user_id
LEFT JOIN
    molecules m ON usm.molecule_id = m.id
LEFT JOIN
    user_saved_mixtures uxm ON up.id = uxm.user_id
LEFT JOIN
    mixtures mx ON uxm.mixture_id = mx.id
LEFT JOIN
    user_saved_experiments use ON up.id = use.user_id
LEFT JOIN
    experiments e ON use.experiment_id = e.id
LEFT JOIN
    mixture_components mc ON mx.id = mc.mixture_id
GROUP BY
    up.auth_user_id;

-- Create refresh function to update all materialized views
CREATE OR REPLACE FUNCTION refresh_materialized_views()
RETURNS void AS $$
BEGIN
    REFRESH MATERIALIZED VIEW CONCURRENTLY public_molecules_summary;
    REFRESH MATERIALIZED VIEW CONCURRENTLY public_molecular_properties;
    REFRESH MATERIALIZED VIEW CONCURRENTLY public_mixtures_summary;
    REFRESH MATERIALIZED VIEW CONCURRENTLY popular_molecules;
    REFRESH MATERIALIZED VIEW CONCURRENTLY property_statistics;
    REFRESH MATERIALIZED VIEW CONCURRENTLY user_favorites_statistics;
    
    -- Log the refresh
    INSERT INTO materialized_view_refresh_log (refresh_time)
    VALUES (NOW());
END;
$$ LANGUAGE plpgsql;

-- Create a table to track refresh history
CREATE TABLE IF NOT EXISTS materialized_view_refresh_log (
    id SERIAL PRIMARY KEY,
    refresh_time TIMESTAMP WITH TIME ZONE NOT NULL
);

-- To schedule refresh:
-- COMMENT OUT FOR NOW: SELECT cron.schedule('refresh_materialized_views', '0 * * * *', 'SELECT refresh_materialized_views()');
-- This requires pg_cron extension to be activated in the database