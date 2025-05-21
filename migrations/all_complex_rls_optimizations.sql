-- All Complex RLS Query Optimizations for CryoProtect
-- This file contains all the optimizations needed to improve complex query performance with RLS

BEGIN;

-- 1. Create Security Definer Functions for Access Control
-- ======================================================

-- Check if user is a team member
CREATE OR REPLACE FUNCTION is_team_member(team_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE auth_user_id = auth.uid() AND user_profile.team_id = $1
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Check if user is a team admin
CREATE OR REPLACE FUNCTION is_team_admin(team_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM team_members
    WHERE team_id = $1 AND user_id = auth.uid() AND role = 'admin'
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Get all teams for the current user
CREATE OR REPLACE FUNCTION get_user_teams()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT team_id FROM user_profile
  WHERE auth_user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Check if user has access to a molecule
CREATE OR REPLACE FUNCTION has_molecule_access(molecule_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM molecules m
    WHERE m.id = molecule_id AND (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_molecules pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.molecule_id = molecule_id AND up.auth_user_id = auth.uid()
      )
    )
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Check if user has access to a mixture
CREATE OR REPLACE FUNCTION has_mixture_access(mixture_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM mixtures m
    WHERE m.id = mixture_id AND (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_mixtures pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.mixture_id = mixture_id AND up.auth_user_id = auth.uid()
      )
    )
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Check if user has required clearance level
CREATE OR REPLACE FUNCTION user_has_clearance(required_level text)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE auth_user_id = auth.uid() 
    AND (clearance_level >= required_level OR required_level IS NULL)
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Check if user has access to an experiment
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

-- Check if user has access to a project
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

-- Check if user is project admin
CREATE OR REPLACE FUNCTION is_project_admin(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM project_members pm
    WHERE pm.project_id = project_id AND pm.user_id = auth.uid() AND pm.role = 'admin'
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 2. Create Complex Query Optimization Functions
-- =============================================

-- Find molecules by property range query
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

-- Full-text search optimization function
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

-- Bulk molecule access check
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

-- Get molecules with properties optimization function
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

-- Get mixtures with components optimization function
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

-- Find molecules with similar properties
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

-- 3. Create Performance Indexes
-- ============================

-- User profile indexes
CREATE INDEX IF NOT EXISTS idx_user_profile_auth_user_id ON user_profile(auth_user_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_team_id ON user_profile(team_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_clearance_level ON user_profile(clearance_level);
CREATE INDEX IF NOT EXISTS idx_user_profile_role ON user_profile(role);

-- Team projects indexes
CREATE INDEX IF NOT EXISTS idx_team_projects_team_id ON team_projects(team_id);
CREATE INDEX IF NOT EXISTS idx_team_projects_project_id ON team_projects(project_id);

-- Project molecule indexes
CREATE INDEX IF NOT EXISTS idx_project_molecules_project_id ON project_molecules(project_id);
CREATE INDEX IF NOT EXISTS idx_project_molecules_molecule_id ON project_molecules(molecule_id);

-- Project mixture indexes
CREATE INDEX IF NOT EXISTS idx_project_mixtures_project_id ON project_mixtures(project_id);
CREATE INDEX IF NOT EXISTS idx_project_mixtures_mixture_id ON project_mixtures(mixture_id);

-- Project experiment indexes
CREATE INDEX IF NOT EXISTS idx_project_experiments_project_id ON project_experiments(project_id);
CREATE INDEX IF NOT EXISTS idx_project_experiments_experiment_id ON project_experiments(experiment_id);

-- Molecule ownership and visibility indexes
CREATE INDEX IF NOT EXISTS idx_molecules_created_by ON molecules(created_by);
CREATE INDEX IF NOT EXISTS idx_molecules_is_public ON molecules(is_public);

-- Mixture ownership and visibility indexes
CREATE INDEX IF NOT EXISTS idx_mixtures_created_by ON mixtures(created_by);
CREATE INDEX IF NOT EXISTS idx_mixtures_is_public ON mixtures(is_public);

-- Experiment ownership and visibility indexes
CREATE INDEX IF NOT EXISTS idx_experiments_created_by ON experiments(created_by);
CREATE INDEX IF NOT EXISTS idx_experiments_is_public ON experiments(is_public);

-- Molecular Properties indexes
CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties(molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_name ON molecular_properties(property_name);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_type ON molecular_properties(property_type);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_sensitivity_level ON molecular_properties(sensitivity_level);

-- Team members indexes
CREATE INDEX IF NOT EXISTS idx_team_members_team_id ON team_members(team_id);
CREATE INDEX IF NOT EXISTS idx_team_members_user_id ON team_members(user_id);
CREATE INDEX IF NOT EXISTS idx_team_members_role ON team_members(role);

-- Project members indexes
CREATE INDEX IF NOT EXISTS idx_project_members_project_id ON project_members(project_id);
CREATE INDEX IF NOT EXISTS idx_project_members_user_id ON project_members(user_id);
CREATE INDEX IF NOT EXISTS idx_project_members_role ON project_members(role);

-- Partial indexes for common access patterns
CREATE INDEX IF NOT EXISTS idx_molecules_public ON molecules(id) WHERE is_public = true;
CREATE INDEX IF NOT EXISTS idx_mixtures_public ON mixtures(id) WHERE is_public = true;
CREATE INDEX IF NOT EXISTS idx_experiments_public ON experiments(id) WHERE is_public = true;
CREATE INDEX IF NOT EXISTS idx_molecular_properties_sensitive 
ON molecular_properties(molecule_id, sensitivity_level)
WHERE sensitivity_level > 'low';

-- Composite indexes for common query patterns
CREATE INDEX IF NOT EXISTS idx_user_profile_auth_team 
ON user_profile(auth_user_id, team_id);

CREATE INDEX IF NOT EXISTS idx_project_mixtures_composite
ON project_mixtures(mixture_id, project_id);

CREATE INDEX IF NOT EXISTS idx_project_molecules_composite
ON project_molecules(molecule_id, project_id);

CREATE INDEX IF NOT EXISTS idx_project_experiments_composite
ON project_experiments(experiment_id, project_id);

-- Create full-text search index for molecules
CREATE INDEX IF NOT EXISTS idx_molecules_name_description_tsvector 
ON molecules USING gin(to_tsvector('english', COALESCE(name, '') || ' ' || COALESCE(description, '')));

-- Create indexes for property searches
CREATE INDEX IF NOT EXISTS idx_molecular_properties_name_numeric 
ON molecular_properties(property_name, (property_value::numeric))
WHERE property_value ~ '^-?[0-9]+(\.[0-9]+)?$';

CREATE INDEX IF NOT EXISTS idx_molecular_properties_name_text
ON molecular_properties(property_name, property_value)
WHERE property_type = 'text';

-- Create covering indexes for common queries
CREATE INDEX IF NOT EXISTS idx_molecules_cover_common
ON molecules(id, name, smiles, molecular_formula, is_public, created_by);

CREATE INDEX IF NOT EXISTS idx_mixtures_cover_common
ON mixtures(id, name, description, is_public, created_by);

-- 4. Create Materialized Views
-- ===========================

-- Public molecules summary
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

-- Create search indexes on the materialized view
CREATE INDEX IF NOT EXISTS idx_public_molecules_summary_name
ON public_molecules_summary(name);

CREATE INDEX IF NOT EXISTS idx_public_molecules_summary_cid
ON public_molecules_summary(cid);

-- Public molecular properties view
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

-- Popular property types summary
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

-- Create refresh materialized views function
CREATE OR REPLACE FUNCTION refresh_materialized_views()
RETURNS void AS $$
BEGIN
    REFRESH MATERIALIZED VIEW CONCURRENTLY public_molecules_summary;
    REFRESH MATERIALIZED VIEW CONCURRENTLY public_molecular_properties;
    REFRESH MATERIALIZED VIEW CONCURRENTLY property_types_summary;
    
    -- Log the refresh if the log table exists
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'materialized_view_refresh_log') THEN
        INSERT INTO materialized_view_refresh_log (refresh_time)
        VALUES (NOW());
    END IF;
END;
$$ LANGUAGE plpgsql;

-- Create refresh log table if it doesn't exist
CREATE TABLE IF NOT EXISTS materialized_view_refresh_log (
    id SERIAL PRIMARY KEY,
    refresh_time TIMESTAMP WITH TIME ZONE NOT NULL
);

-- 5. Create Query Cache System
-- ==========================

-- Create cache table for storing temporary query results
CREATE TABLE IF NOT EXISTS query_cache (
    cache_key TEXT PRIMARY KEY,
    result JSONB NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    expires_at TIMESTAMP WITH TIME ZONE NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_query_cache_expiry
ON query_cache(expires_at);

-- Create function to clean expired cache entries
CREATE OR REPLACE FUNCTION clean_expired_query_cache() RETURNS void AS $$
BEGIN
    DELETE FROM query_cache WHERE expires_at < CURRENT_TIMESTAMP;
END;
$$ LANGUAGE plpgsql;

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

-- 6. Update RLS Policies to Use Security Definer Functions
-- ======================================================

-- Molecules table policies
DROP POLICY IF EXISTS "users_can_view_molecules" ON molecules;
CREATE POLICY "users_can_view_molecules" ON molecules
FOR SELECT
TO authenticated
USING (
  has_molecule_access(id)
);

-- Mixtures table policies
DROP POLICY IF EXISTS "users_can_view_mixtures" ON mixtures;
CREATE POLICY "users_can_view_mixtures" ON mixtures
FOR SELECT
TO authenticated
USING (
  has_mixture_access(id)
);

-- Experiments table policies 
DROP POLICY IF EXISTS "users_can_view_experiments" ON experiments;
CREATE POLICY "users_can_view_experiments" ON experiments
FOR SELECT
TO authenticated
USING (
  has_experiment_access(id)
);

-- Molecular properties table policies
DROP POLICY IF EXISTS "users_can_view_properties" ON molecular_properties;
CREATE POLICY "users_can_view_properties" ON molecular_properties
FOR SELECT
TO authenticated
USING (
  has_molecule_access(molecule_id) AND
  (sensitivity_level IS NULL OR 
   sensitivity_level <= 'medium' OR 
   user_has_clearance(sensitivity_level))
);

-- Team projects table policies
DROP POLICY IF EXISTS "users_can_view_team_projects" ON team_projects;
CREATE POLICY "users_can_view_team_projects" ON team_projects
FOR SELECT
TO authenticated
USING (
  is_team_member(team_id)
);

-- Project members table policies
DROP POLICY IF EXISTS "users_can_view_project_members" ON project_members;
CREATE POLICY "users_can_view_project_members" ON project_members
FOR SELECT
TO authenticated
USING (
  has_project_access(project_id)
);

-- Team members table policies
DROP POLICY IF EXISTS "users_can_view_team_members" ON team_members;
CREATE POLICY "users_can_view_team_members" ON team_members
FOR SELECT
TO authenticated
USING (
  is_team_member(team_id)
);

COMMIT;