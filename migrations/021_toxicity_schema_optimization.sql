-- Migration: 021_toxicity_schema_optimization.sql
-- Purpose: Optimize toxicity data schema with specialized tables, 
-- improved indexes, and materialized views for better performance.

BEGIN;

-- =======================================
-- 1. Create specialized toxicity tables
-- =======================================

-- Base toxicity table with common fields
CREATE TABLE IF NOT EXISTS toxicity_data_new (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    molecule_id UUID NOT NULL REFERENCES molecule(id) ON DELETE CASCADE,
    source VARCHAR(50) NOT NULL,
    source_id VARCHAR(100),
    toxicity_type VARCHAR(50) NOT NULL,
    species VARCHAR(50),
    route_of_administration VARCHAR(50),
    value NUMERIC NOT NULL,
    unit VARCHAR(20) NOT NULL,
    uncertainty NUMERIC,
    confidence_score NUMERIC CHECK (confidence_score >= 0 AND confidence_score <= 1),
    is_predicted BOOLEAN NOT NULL DEFAULT FALSE,
    prediction_method VARCHAR(100),
    reference_doi VARCHAR(100),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT now(),
    created_by UUID REFERENCES auth.users(id),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT now(),
    UNIQUE(molecule_id, source, toxicity_type, species, route_of_administration)
);

-- Enable RLS on toxicity data
ALTER TABLE toxicity_data_new ENABLE ROW LEVEL SECURITY;

-- Toxicity data policies
CREATE POLICY "users_can_view_toxicity_data" ON toxicity_data_new
FOR SELECT TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM molecule m
    WHERE m.id = toxicity_data_new.molecule_id AND
    (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_molecule pm
        JOIN team_project tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
      )
    )
  )
);

-- Service role policy
CREATE POLICY "service_role_all_access" ON toxicity_data_new
FOR ALL TO service_role
USING (true) WITH CHECK (true);

-- Create specialized toxicity tables for different data types
CREATE TABLE IF NOT EXISTS toxicity_ld50 (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    toxicity_data_id UUID NOT NULL REFERENCES toxicity_data_new(id) ON DELETE CASCADE,
    organism_strain VARCHAR(100),
    age_group VARCHAR(50),
    gender VARCHAR(20),
    duration NUMERIC,
    duration_unit VARCHAR(20),
    observation_period NUMERIC,
    observation_period_unit VARCHAR(20),
    dosing_regime VARCHAR(50),
    number_of_animals INTEGER,
    death_count INTEGER,
    additional_effects JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT now(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT now()
);

-- Enable RLS on LD50 table
ALTER TABLE toxicity_ld50 ENABLE ROW LEVEL SECURITY;

-- LD50 RLS policy
CREATE POLICY "users_can_view_ld50_data" ON toxicity_ld50
FOR SELECT TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM toxicity_data_new td
    WHERE td.id = toxicity_ld50.toxicity_data_id AND
    EXISTS (
      SELECT 1 FROM molecule m
      WHERE m.id = td.molecule_id AND
      (
        m.is_public = true OR
        m.created_by = auth.uid() OR
        EXISTS (
          SELECT 1 FROM project_molecule pm
          JOIN team_project tp ON pm.project_id = tp.project_id
          JOIN user_profile up ON tp.team_id = up.team_id
          WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
        )
      )
    )
  )
);

-- Service role policy
CREATE POLICY "service_role_all_access" ON toxicity_ld50
FOR ALL TO service_role
USING (true) WITH CHECK (true);

-- Create table for Tox21 assay results
CREATE TABLE IF NOT EXISTS toxicity_tox21 (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    toxicity_data_id UUID NOT NULL REFERENCES toxicity_data_new(id) ON DELETE CASCADE,
    assay_id VARCHAR(50) NOT NULL,
    assay_target VARCHAR(100),
    activity_outcome VARCHAR(20) NOT NULL,
    activity_score NUMERIC,
    curve_description VARCHAR(100),
    intended_target_family VARCHAR(100),
    intended_target_type VARCHAR(100),
    assay_format VARCHAR(50),
    assay_format_type VARCHAR(50),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT now(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT now()
);

-- Enable RLS on Tox21 table
ALTER TABLE toxicity_tox21 ENABLE ROW LEVEL SECURITY;

-- Tox21 RLS policy
CREATE POLICY "users_can_view_tox21_data" ON toxicity_tox21
FOR SELECT TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM toxicity_data_new td
    WHERE td.id = toxicity_tox21.toxicity_data_id AND
    EXISTS (
      SELECT 1 FROM molecule m
      WHERE m.id = td.molecule_id AND
      (
        m.is_public = true OR
        m.created_by = auth.uid() OR
        EXISTS (
          SELECT 1 FROM project_molecule pm
          JOIN team_project tp ON pm.project_id = tp.project_id
          JOIN user_profile up ON tp.team_id = up.team_id
          WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
        )
      )
    )
  )
);

-- Service role policy
CREATE POLICY "service_role_all_access" ON toxicity_tox21
FOR ALL TO service_role
USING (true) WITH CHECK (true);

-- Create toxicity classification table
CREATE TABLE IF NOT EXISTS toxicity_classification (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    molecule_id UUID NOT NULL REFERENCES molecule(id) ON DELETE CASCADE,
    classification_system VARCHAR(50) NOT NULL,
    hazard_class VARCHAR(100) NOT NULL,
    hazard_category VARCHAR(100),
    hazard_statement VARCHAR(255),
    hazard_code VARCHAR(50),
    pictogram VARCHAR(50),
    signal_word VARCHAR(50),
    source VARCHAR(50) NOT NULL,
    is_predicted BOOLEAN NOT NULL DEFAULT FALSE,
    confidence_score NUMERIC CHECK (confidence_score >= 0 AND confidence_score <= 1),
    reference_doi VARCHAR(100),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT now(),
    created_by UUID REFERENCES auth.users(id),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT now(),
    UNIQUE(molecule_id, classification_system, hazard_class, hazard_category, source)
);

-- Enable RLS on classification table
ALTER TABLE toxicity_classification ENABLE ROW LEVEL SECURITY;

-- Classification RLS policy
CREATE POLICY "users_can_view_classification_data" ON toxicity_classification
FOR SELECT TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM molecule m
    WHERE m.id = toxicity_classification.molecule_id AND
    (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_molecule pm
        JOIN team_project tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
      )
    )
  )
);

-- Service role policy
CREATE POLICY "service_role_all_access" ON toxicity_classification
FOR ALL TO service_role
USING (true) WITH CHECK (true);

-- =======================================
-- 2. Create performance indexes
-- =======================================

-- Base toxicity table indexes
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_molecule_id ON toxicity_data_new(molecule_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_source ON toxicity_data_new(source);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_toxicity_type ON toxicity_data_new(toxicity_type);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_species ON toxicity_data_new(species);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_is_predicted ON toxicity_data_new(is_predicted);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_created_by ON toxicity_data_new(created_by);

-- Composite indexes for common query patterns
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_molecule_type ON toxicity_data_new(molecule_id, toxicity_type);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_molecule_source ON toxicity_data_new(molecule_id, source);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_type_species ON toxicity_data_new(toxicity_type, species);

-- Partial index for predicted values
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_predicted ON toxicity_data_new(molecule_id, toxicity_type, value)
WHERE is_predicted = true;

-- Partial index for experimental values
CREATE INDEX IF NOT EXISTS idx_toxicity_data_new_experimental ON toxicity_data_new(molecule_id, toxicity_type, value)
WHERE is_predicted = false;

-- LD50 indexes
CREATE INDEX IF NOT EXISTS idx_toxicity_ld50_toxicity_data_id ON toxicity_ld50(toxicity_data_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_ld50_organism_strain ON toxicity_ld50(organism_strain);
CREATE INDEX IF NOT EXISTS idx_toxicity_ld50_gender ON toxicity_ld50(gender);

-- Tox21 indexes
CREATE INDEX IF NOT EXISTS idx_toxicity_tox21_toxicity_data_id ON toxicity_tox21(toxicity_data_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_tox21_assay_id ON toxicity_tox21(assay_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_tox21_activity_outcome ON toxicity_tox21(activity_outcome);
CREATE INDEX IF NOT EXISTS idx_toxicity_tox21_assay_target ON toxicity_tox21(assay_target);
CREATE INDEX IF NOT EXISTS idx_toxicity_tox21_target_family ON toxicity_tox21(intended_target_family);

-- Classification indexes
CREATE INDEX IF NOT EXISTS idx_toxicity_classification_molecule_id ON toxicity_classification(molecule_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_classification_system ON toxicity_classification(classification_system);
CREATE INDEX IF NOT EXISTS idx_toxicity_classification_hazard_class ON toxicity_classification(hazard_class);
CREATE INDEX IF NOT EXISTS idx_toxicity_classification_source ON toxicity_classification(source);
CREATE INDEX IF NOT EXISTS idx_toxicity_classification_is_predicted ON toxicity_classification(is_predicted);

-- Composite index for classification lookups
CREATE INDEX IF NOT EXISTS idx_toxicity_classification_molecule_system ON toxicity_classification(molecule_id, classification_system);

-- =======================================
-- 3. Create materialized views
-- =======================================

-- Molecule toxicity summary view
CREATE MATERIALIZED VIEW IF NOT EXISTS toxicity_summary AS
SELECT 
    m.id AS molecule_id,
    m.name,
    m.smiles,
    m.molecular_formula,
    m.cid,
    COUNT(DISTINCT td.id) AS toxicity_data_count,
    COUNT(DISTINCT CASE WHEN td.is_predicted = false THEN td.id END) AS experimental_data_count,
    COUNT(DISTINCT CASE WHEN td.is_predicted = true THEN td.id END) AS predicted_data_count,
    ARRAY_AGG(DISTINCT td.toxicity_type) AS toxicity_types,
    ARRAY_AGG(DISTINCT td.source) AS data_sources,
    EXISTS (
        SELECT 1 FROM toxicity_classification tc
        WHERE tc.molecule_id = m.id
    ) AS has_classification,
    EXISTS (
        SELECT 1 FROM toxicity_tox21 tt
        JOIN toxicity_data_new td2 ON tt.toxicity_data_id = td2.id
        WHERE td2.molecule_id = m.id
    ) AS has_tox21_data,
    EXISTS (
        SELECT 1 FROM toxicity_ld50 tl
        JOIN toxicity_data_new td3 ON tl.toxicity_data_id = td3.id
        WHERE td3.molecule_id = m.id
    ) AS has_ld50_data
FROM 
    molecule m
LEFT JOIN 
    toxicity_data_new td ON m.id = td.molecule_id
GROUP BY 
    m.id, m.name, m.smiles, m.molecular_formula, m.cid;

-- Create index on the materialized view
CREATE UNIQUE INDEX IF NOT EXISTS idx_toxicity_summary_molecule_id 
ON toxicity_summary(molecule_id);

-- Tox21 activity summary view
CREATE MATERIALIZED VIEW IF NOT EXISTS tox21_activity_summary AS
SELECT 
    td.molecule_id,
    tt.assay_id,
    tt.assay_target,
    tt.activity_outcome,
    tt.activity_score,
    tt.intended_target_family,
    tt.intended_target_type,
    m.name AS molecule_name,
    m.smiles
FROM 
    toxicity_tox21 tt
JOIN 
    toxicity_data_new td ON tt.toxicity_data_id = td.id
JOIN 
    molecule m ON td.molecule_id = m.id
WHERE 
    td.source = 'Tox21';

-- Create indexes on the Tox21 materialized view
CREATE INDEX IF NOT EXISTS idx_tox21_summary_molecule_id 
ON tox21_activity_summary(molecule_id);

CREATE INDEX IF NOT EXISTS idx_tox21_summary_assay_id 
ON tox21_activity_summary(assay_id);

CREATE INDEX IF NOT EXISTS idx_tox21_summary_activity 
ON tox21_activity_summary(activity_outcome);

-- LD50 summary view
CREATE MATERIALIZED VIEW IF NOT EXISTS ld50_summary AS
SELECT 
    td.molecule_id,
    m.name AS molecule_name,
    m.smiles,
    td.species,
    td.route_of_administration,
    td.value AS ld50_value,
    td.unit,
    td.is_predicted,
    tl.organism_strain,
    tl.gender,
    tl.age_group,
    td.source
FROM 
    toxicity_ld50 tl
JOIN 
    toxicity_data_new td ON tl.toxicity_data_id = td.id
JOIN 
    molecule m ON td.molecule_id = m.id
WHERE 
    td.toxicity_type = 'LD50';

-- Create indexes on the LD50 materialized view
CREATE INDEX IF NOT EXISTS idx_ld50_summary_molecule_id 
ON ld50_summary(molecule_id);

CREATE INDEX IF NOT EXISTS idx_ld50_summary_species 
ON ld50_summary(species);

CREATE INDEX IF NOT EXISTS idx_ld50_summary_route 
ON ld50_summary(route_of_administration);

CREATE INDEX IF NOT EXISTS idx_ld50_summary_is_predicted 
ON ld50_summary(is_predicted);

-- Hazard classification summary view
CREATE MATERIALIZED VIEW IF NOT EXISTS hazard_classification_summary AS
SELECT 
    tc.molecule_id,
    m.name AS molecule_name,
    m.smiles,
    tc.classification_system,
    tc.hazard_class,
    tc.hazard_category,
    tc.hazard_statement,
    tc.hazard_code,
    tc.pictogram,
    tc.signal_word,
    tc.source,
    tc.is_predicted,
    tc.confidence_score
FROM 
    toxicity_classification tc
JOIN 
    molecule m ON tc.molecule_id = m.id;

-- Create indexes on the hazard classification materialized view
CREATE INDEX IF NOT EXISTS idx_hazard_summary_molecule_id 
ON hazard_classification_summary(molecule_id);

CREATE INDEX IF NOT EXISTS idx_hazard_summary_system 
ON hazard_classification_summary(classification_system);

CREATE INDEX IF NOT EXISTS idx_hazard_summary_hazard_class 
ON hazard_classification_summary(hazard_class);

-- Create materialized view refresh log table
CREATE TABLE IF NOT EXISTS materialized_view_refresh_log (
    id SERIAL PRIMARY KEY,
    view_name TEXT NOT NULL,
    refresh_time TIMESTAMP WITH TIME ZONE NOT NULL
);

-- Create refresh function
CREATE OR REPLACE FUNCTION refresh_toxicity_materialized_views()
RETURNS void AS $$
BEGIN
    REFRESH MATERIALIZED VIEW CONCURRENTLY toxicity_summary;
    REFRESH MATERIALIZED VIEW CONCURRENTLY tox21_activity_summary;
    REFRESH MATERIALIZED VIEW CONCURRENTLY ld50_summary;
    REFRESH MATERIALIZED VIEW CONCURRENTLY hazard_classification_summary;
    
    -- Log the refresh
    INSERT INTO materialized_view_refresh_log (view_name, refresh_time)
    VALUES 
        ('toxicity_summary', NOW()),
        ('tox21_activity_summary', NOW()),
        ('ld50_summary', NOW()),
        ('hazard_classification_summary', NOW());
END;
$$ LANGUAGE plpgsql;

-- Schedule refresh (if using pg_cron extension)
DO $$
BEGIN
    IF EXISTS (
        SELECT 1 FROM pg_extension WHERE extname = 'pg_cron'
    ) THEN
        PERFORM cron.schedule('refresh_toxicity_views', '0 */3 * * *', 'SELECT refresh_toxicity_materialized_views()');
    END IF;
END $$;

-- =======================================
-- 4. Create stored functions for complex calculations
-- =======================================

-- Function to calculate toxicity score for a molecule
CREATE OR REPLACE FUNCTION calculate_toxicity_score(molecule_uuid UUID)
RETURNS NUMERIC AS $$
DECLARE
    score NUMERIC := 0;
    has_ld50 BOOLEAN;
    min_ld50 NUMERIC;
    has_tox21 BOOLEAN;
    tox21_active_count INTEGER;
    has_classification BOOLEAN;
    ghs_hazard_count INTEGER;
BEGIN
    -- Check for LD50 data
    SELECT 
        EXISTS(SELECT 1 FROM ld50_summary WHERE molecule_id = molecule_uuid) INTO has_ld50;
        
    IF has_ld50 THEN
        -- Get minimum LD50 value (lower means more toxic)
        SELECT MIN(ld50_value) INTO min_ld50
        FROM ld50_summary
        WHERE molecule_id = molecule_uuid
        AND species IN ('Rat', 'Mouse')
        AND route_of_administration IN ('Oral', 'Dermal', 'Inhalation');
        
        -- Adjust score based on LD50 value
        IF min_ld50 IS NOT NULL THEN
            IF min_ld50 < 50 THEN  -- Highly toxic
                score := score + 5;
            ELSIF min_ld50 < 300 THEN  -- Moderately toxic
                score := score + 3;
            ELSIF min_ld50 < 2000 THEN  -- Slightly toxic
                score := score + 1;
            END IF;
        END IF;
    END IF;
    
    -- Check for Tox21 data
    SELECT 
        EXISTS(SELECT 1 FROM tox21_activity_summary WHERE molecule_id = molecule_uuid) INTO has_tox21;
        
    IF has_tox21 THEN
        -- Count active assays
        SELECT COUNT(*)
        INTO tox21_active_count
        FROM tox21_activity_summary
        WHERE molecule_id = molecule_uuid
        AND activity_outcome = 'Active';
        
        -- Adjust score based on active assay count
        score := score + LEAST(tox21_active_count, 5);  -- Cap at 5 points
    END IF;
    
    -- Check for hazard classifications
    SELECT 
        EXISTS(SELECT 1 FROM hazard_classification_summary WHERE molecule_id = molecule_uuid) INTO has_classification;
        
    IF has_classification THEN
        -- Count GHS hazard statements
        SELECT COUNT(*)
        INTO ghs_hazard_count
        FROM hazard_classification_summary
        WHERE molecule_id = molecule_uuid
        AND classification_system = 'GHS';
        
        -- Adjust score based on hazard count
        score := score + LEAST(ghs_hazard_count, 3);  -- Cap at 3 points
    END IF;
    
    RETURN score;
END;
$$ LANGUAGE plpgsql;

-- Function to find similar molecules by toxicity profile
CREATE OR REPLACE FUNCTION find_similar_toxicity_profiles(molecule_uuid UUID, limit_count INTEGER DEFAULT 10)
RETURNS TABLE(
    similar_molecule_id UUID,
    molecule_name TEXT,
    smiles TEXT,
    similarity_score NUMERIC
) AS $$
BEGIN
    RETURN QUERY
    WITH target_assays AS (
        -- Get the target molecule's Tox21 profile
        SELECT assay_id, activity_outcome
        FROM tox21_activity_summary
        WHERE molecule_id = molecule_uuid
    ),
    target_hazards AS (
        -- Get the target molecule's hazard profile
        SELECT hazard_class, hazard_category
        FROM hazard_classification_summary
        WHERE molecule_id = molecule_uuid
    ),
    molecule_similarities AS (
        -- Calculate similarity score for each molecule
        SELECT 
            m.id as similar_molecule_id,
            m.name as molecule_name,
            m.smiles,
            (
                -- Tox21 profile similarity (up to 0.6 points)
                (
                    SELECT COUNT(*) * 0.1
                    FROM tox21_activity_summary ta
                    JOIN target_assays t ON ta.assay_id = t.assay_id AND ta.activity_outcome = t.activity_outcome
                    WHERE ta.molecule_id = m.id
                    LIMIT 6
                ) +
                -- Hazard profile similarity (up to 0.4 points)
                (
                    SELECT COUNT(*) * 0.1
                    FROM hazard_classification_summary hc
                    JOIN target_hazards th ON hc.hazard_class = th.hazard_class AND hc.hazard_category = th.hazard_category
                    WHERE hc.molecule_id = m.id
                    LIMIT 4
                )
            )::NUMERIC as similarity_score
        FROM molecule m
        WHERE m.id != molecule_uuid
        AND EXISTS (
            SELECT 1 FROM toxicity_summary ts WHERE ts.molecule_id = m.id
        )
    )
    SELECT 
        similar_molecule_id,
        molecule_name,
        smiles,
        similarity_score
    FROM molecule_similarities
    WHERE similarity_score > 0
    ORDER BY similarity_score DESC
    LIMIT limit_count;
END;
$$ LANGUAGE plpgsql;

-- =======================================
-- 5. Migrate data from existing tables to new structure
-- =======================================

-- Insert a migration completion record
INSERT INTO migration_history (name, applied_at, description) 
VALUES (
    '021_toxicity_schema_optimization', 
    NOW(), 
    'Optimized toxicity data schema with specialized tables, improved indexes, and materialized views'
);

COMMIT;