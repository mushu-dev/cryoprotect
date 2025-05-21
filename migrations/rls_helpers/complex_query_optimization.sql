-- Complex Query Optimization SQL
-- This adds specialized functions and indexes for complex query patterns

-- 1. Molecular Property Range Query Optimization
CREATE OR REPLACE FUNCTION find_molecules_by_property_range(
    property_name_param TEXT,
    min_value NUMERIC,
    max_value NUMERIC
) RETURNS SETOF uuid AS $$
BEGIN
    RETURN QUERY
    SELECT DISTINCT m.id
    FROM molecules m
    JOIN molecular_properties mp ON m.id = mp.molecule_id
    WHERE mp.property_name = property_name_param
    AND mp.property_value::numeric >= min_value
    AND mp.property_value::numeric <= max_value
    AND (m.is_public = true OR m.created_by = auth.uid() OR 
         EXISTS (
             SELECT 1 FROM project_molecules pm
             JOIN team_projects tp ON pm.project_id = tp.project_id
             JOIN user_profile up ON tp.team_id = up.team_id
             WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
         ));
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 2. Full Text Search Optimization for molecules
CREATE INDEX IF NOT EXISTS idx_molecules_name_description_tsvector 
ON molecules USING gin(to_tsvector('english', COALESCE(name, '') || ' ' || COALESCE(description, '')));

CREATE OR REPLACE FUNCTION search_molecules_text(search_term TEXT)
RETURNS SETOF uuid AS $$
BEGIN
    RETURN QUERY
    SELECT id FROM molecules
    WHERE to_tsvector('english', COALESCE(name, '') || ' ' || COALESCE(description, '')) @@ plainto_tsquery('english', search_term)
    AND (is_public = true OR created_by = auth.uid() OR
         EXISTS (
             SELECT 1 FROM project_molecules pm
             JOIN team_projects tp ON pm.project_id = tp.project_id
             JOIN user_profile up ON tp.team_id = up.team_id
             WHERE pm.molecule_id = molecules.id AND up.auth_user_id = auth.uid()
         ));
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 3. Bulk Molecule Access Check
CREATE OR REPLACE FUNCTION filter_accessible_molecules(molecule_ids uuid[])
RETURNS SETOF uuid AS $$
BEGIN
    RETURN QUERY
    SELECT id FROM molecules
    WHERE id = ANY(molecule_ids)
    AND (is_public = true OR created_by = auth.uid() OR
         EXISTS (
             SELECT 1 FROM project_molecules pm
             JOIN team_projects tp ON pm.project_id = tp.project_id
             JOIN user_profile up ON tp.team_id = up.team_id
             WHERE pm.molecule_id = molecules.id AND up.auth_user_id = auth.uid()
         ));
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 4. Optimization for Common Query Patterns
-- 4.1 Molecule-with-properties query
CREATE OR REPLACE FUNCTION get_molecules_with_properties(
    limit_count INTEGER DEFAULT 100,
    offset_val INTEGER DEFAULT 0
) RETURNS TABLE (
    id uuid,
    name TEXT, 
    smiles TEXT, 
    molecular_formula TEXT,
    property_count INTEGER
) AS $$
DECLARE
    user_id uuid;
BEGIN
    user_id := auth.uid();
    
    RETURN QUERY
    WITH accessible_molecules AS (
        SELECT m.id
        FROM molecules m
        WHERE m.is_public = true 
           OR m.created_by = user_id
           OR EXISTS (
               SELECT 1 FROM project_molecules pm
               JOIN team_projects tp ON pm.project_id = tp.project_id
               JOIN user_profile up ON tp.team_id = up.team_id
               WHERE pm.molecule_id = m.id AND up.auth_user_id = user_id
           )
    )
    SELECT 
        m.id,
        m.name,
        m.smiles,
        m.molecular_formula,
        COUNT(mp.id)::INTEGER AS property_count
    FROM 
        molecules m
    JOIN
        accessible_molecules am ON m.id = am.id
    LEFT JOIN
        molecular_properties mp ON m.id = mp.molecule_id
    GROUP BY
        m.id, m.name, m.smiles, m.molecular_formula
    ORDER BY
        m.name
    LIMIT limit_count
    OFFSET offset_val;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 4.2 Mixture-with-components query
CREATE OR REPLACE FUNCTION get_mixtures_with_components(
    limit_count INTEGER DEFAULT 100,
    offset_val INTEGER DEFAULT 0
) RETURNS TABLE (
    id uuid,
    name TEXT, 
    description TEXT,
    component_count INTEGER
) AS $$
DECLARE
    user_id uuid;
BEGIN
    user_id := auth.uid();
    
    RETURN QUERY
    WITH accessible_mixtures AS (
        SELECT m.id
        FROM mixtures m
        WHERE m.is_public = true 
           OR m.created_by = user_id
           OR EXISTS (
               SELECT 1 FROM project_mixtures pm
               JOIN team_projects tp ON pm.project_id = tp.project_id
               JOIN user_profile up ON tp.team_id = up.team_id
               WHERE pm.mixture_id = m.id AND up.auth_user_id = user_id
           )
    )
    SELECT 
        m.id,
        m.name,
        m.description,
        COUNT(mc.id)::INTEGER AS component_count
    FROM 
        mixtures m
    JOIN
        accessible_mixtures am ON m.id = am.id
    LEFT JOIN
        mixture_components mc ON m.id = mc.mixture_id
    GROUP BY
        m.id, m.name, m.description
    ORDER BY
        m.name
    LIMIT limit_count
    OFFSET offset_val;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 5. Experiment access optimization
CREATE OR REPLACE FUNCTION has_experiment_access(experiment_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM experiments e
    WHERE e.id = experiment_id AND (
      e.is_public = true OR
      e.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_experiments pe
        JOIN team_projects tp ON pe.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pe.experiment_id = experiment_id AND up.auth_user_id = auth.uid()
      )
    )
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 6. Project access optimization
CREATE OR REPLACE FUNCTION has_project_access(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM team_projects tp
    JOIN user_profile up ON tp.team_id = up.team_id
    WHERE tp.project_id = project_id AND up.auth_user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 7. Project admin check optimization
CREATE OR REPLACE FUNCTION is_project_admin(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM project_members pm
    WHERE pm.project_id = project_id AND pm.user_id = auth.uid() AND pm.role = 'admin'
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Create hypertable indexes for efficient time-series queries (if timescaledb is available)
DO $$
BEGIN
    IF EXISTS (
        SELECT 1 FROM pg_extension WHERE extname = 'timescaledb'
    ) THEN
        -- Add timescale-specific optimizations for time-series data
        CREATE INDEX IF NOT EXISTS idx_experiments_time 
        ON experiments(created_at DESC, id);
        
        CREATE INDEX IF NOT EXISTS idx_predictions_time 
        ON predictions(created_at DESC, id);
    END IF;
END $$;

-- Additional composite indexes for complex query patterns
CREATE INDEX IF NOT EXISTS idx_molecular_properties_name_numeric 
ON molecular_properties(property_name, (property_value::numeric))
WHERE property_value ~ '^-?[0-9]+(\.[0-9]+)?$';

CREATE INDEX IF NOT EXISTS idx_molecular_properties_name_text
ON molecular_properties(property_name, property_value)
WHERE property_type = 'text';

-- Add covering indexes for common queries that involve multiple columns
CREATE INDEX IF NOT EXISTS idx_molecules_cover_common
ON molecules(id, name, smiles, molecular_formula, is_public, created_by);

CREATE INDEX IF NOT EXISTS idx_mixtures_cover_common
ON mixtures(id, name, description, is_public, created_by);

-- Create index to optimize property-based molecule filtering
CREATE INDEX IF NOT EXISTS idx_molecule_properties_filtered
ON molecular_properties (molecule_id, property_name, property_type, property_value);

-- Create index for molecular formula search
CREATE INDEX IF NOT EXISTS idx_molecules_formula_trgm
ON molecules USING gin (molecular_formula gin_trgm_ops);

-- Create indexes for batch operations
CREATE INDEX IF NOT EXISTS idx_molecules_cid
ON molecules(cid);

CREATE INDEX IF NOT EXISTS idx_molecule_batch_access
ON molecules(is_public, created_by);

-- Create indexes for similar structure searches
CREATE INDEX IF NOT EXISTS idx_molecules_smiles_trgm
ON molecules USING gin (smiles gin_trgm_ops);

-- Create expression index for quick filtering of numeric properties
CREATE INDEX IF NOT EXISTS idx_molecular_properties_numeric_value
ON molecular_properties ((property_value::numeric))
WHERE property_value ~ '^-?[0-9]+(\.[0-9]+)?$';

-- Create a cache table for storing temporary query results
CREATE TABLE IF NOT EXISTS query_cache (
    cache_key TEXT PRIMARY KEY,
    result JSONB NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    expires_at TIMESTAMP WITH TIME ZONE NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_query_cache_expiry
ON query_cache(expires_at);

-- Create function to automatically clean expired cache entries
CREATE OR REPLACE FUNCTION clean_expired_query_cache() RETURNS void AS $$
BEGIN
    DELETE FROM query_cache WHERE expires_at < CURRENT_TIMESTAMP;
END;
$$ LANGUAGE plpgsql;

-- Set up cache maintenance (will only work if pg_cron extension is available)
DO $$
BEGIN
    IF EXISTS (
        SELECT 1 FROM pg_extension WHERE extname = 'pg_cron'
    ) THEN
        -- Schedule cache cleanup every 10 minutes
        PERFORM cron.schedule('clean_query_cache', '*/10 * * * *', 'SELECT clean_expired_query_cache()');
    END IF;
END $$;

-- Create utility function to store query results in cache
CREATE OR REPLACE FUNCTION cache_query_result(
    query_id TEXT,
    result JSONB,
    ttl_minutes INTEGER DEFAULT 10
) RETURNS void AS $$
BEGIN
    -- Delete existing cached result if any
    DELETE FROM query_cache WHERE cache_key = query_id;
    
    -- Store new result
    INSERT INTO query_cache (cache_key, result, expires_at)
    VALUES (query_id, result, CURRENT_TIMESTAMP + (ttl_minutes || ' minutes')::INTERVAL);
END;
$$ LANGUAGE plpgsql;

-- Create utility function to retrieve cached query results
CREATE OR REPLACE FUNCTION get_cached_query_result(query_id TEXT)
RETURNS JSONB AS $$
DECLARE
    cached_result JSONB;
BEGIN
    -- Get cached result if not expired
    SELECT result INTO cached_result
    FROM query_cache
    WHERE cache_key = query_id AND expires_at > CURRENT_TIMESTAMP;
    
    RETURN cached_result;
END;
$$ LANGUAGE plpgsql;

-- Create materialized view for popular property types (for faster searching)
CREATE MATERIALIZED VIEW IF NOT EXISTS property_types_summary AS
SELECT 
    property_name,
    property_type,
    COUNT(*) as usage_count,
    COUNT(DISTINCT molecule_id) as molecule_count,
    MIN(CASE WHEN property_value ~ '^-?[0-9]+(\.[0-9]+)?$' THEN property_value::numeric ELSE NULL END) as min_value,
    MAX(CASE WHEN property_value ~ '^-?[0-9]+(\.[0-9]+)?$' THEN property_value::numeric ELSE NULL END) as max_value,
    array_agg(DISTINCT property_unit) FILTER (WHERE property_unit IS NOT NULL) as units
FROM 
    molecular_properties
GROUP BY 
    property_name, property_type;

CREATE UNIQUE INDEX IF NOT EXISTS idx_property_types_summary_key
ON property_types_summary(property_name, property_type);

-- Create function to get molecules by similar properties
CREATE OR REPLACE FUNCTION find_molecules_by_similar_properties(
    reference_molecule_id uuid,
    similarity_threshold FLOAT DEFAULT 0.7,
    limit_count INTEGER DEFAULT 20
) RETURNS TABLE (
    id uuid,
    name TEXT,
    similarity FLOAT
) AS $$
DECLARE
    user_id uuid;
BEGIN
    user_id := auth.uid();
    
    RETURN QUERY
    WITH reference_properties AS (
        SELECT property_name, property_value::numeric as numeric_value
        FROM molecular_properties
        WHERE molecule_id = reference_molecule_id
        AND property_value ~ '^-?[0-9]+(\.[0-9]+)?$'
    ),
    molecule_similarities AS (
        SELECT 
            m.id,
            m.name,
            COUNT(rp.property_name) FILTER (
                WHERE ABS(mp.property_value::numeric - rp.numeric_value) / 
                     (CASE WHEN rp.numeric_value = 0 THEN 1 ELSE ABS(rp.numeric_value) END) <= (1 - similarity_threshold)
            )::float / 
            NULLIF(COUNT(rp.property_name), 0) as similarity_score
        FROM 
            molecules m
        JOIN 
            molecular_properties mp ON m.id = mp.molecule_id
        JOIN 
            reference_properties rp ON mp.property_name = rp.property_name
        WHERE 
            mp.property_value ~ '^-?[0-9]+(\.[0-9]+)?$'
            AND (m.is_public = true OR m.created_by = user_id OR
                 EXISTS (
                     SELECT 1 FROM project_molecules pm
                     JOIN team_projects tp ON pm.project_id = tp.project_id
                     JOIN user_profile up ON tp.team_id = up.team_id
                     WHERE pm.molecule_id = m.id AND up.auth_user_id = user_id
                 ))
            AND m.id != reference_molecule_id
        GROUP BY 
            m.id, m.name
        HAVING 
            COUNT(rp.property_name) > 0
    )
    SELECT 
        id, 
        name, 
        similarity_score as similarity
    FROM 
        molecule_similarities
    WHERE 
        similarity_score >= similarity_threshold
    ORDER BY 
        similarity_score DESC
    LIMIT 
        limit_count;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;