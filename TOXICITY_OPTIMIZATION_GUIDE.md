# Toxicity Data Optimization Guide

This guide provides detailed implementation steps for optimizing toxicity data handling and querying in the CryoProtect system. It serves as a companion document to the main CRYOPROTECT_OPTIMIZATION_PLAN.md and focuses specifically on toxicity data performance improvements.

## 1. Current Toxicity Data Implementation Analysis

### Current Limitations
1. **Inefficient schema design**: Current schema doesn't efficiently support common query patterns
2. **Missing indexes**: Key columns used in filtering lack appropriate indexes
3. **Inefficient joins**: Complex joins required for toxicity data retrieval
4. **No data partitioning**: All toxicity data stored in a single table regardless of source
5. **Limited caching**: Frequently accessed toxicity data is not cached
6. **Poor query optimization**: Queries not optimized for performance

### Performance Impact
- Slow query execution for toxicity-related queries
- High database load during toxicity data retrieval
- Poor response times for toxicity screening
- Inefficient data access patterns
- Limited scalability as toxicity dataset grows

## 2. Optimization Steps

### 2.1 Schema Optimization

#### Problem
The current schema design requires complex joins and doesn't efficiently support common query patterns for toxicity data.

#### Solution
Optimize the toxicity data schema:
- Create specialized tables for different toxicity data types
- Implement appropriate indexes for common query patterns
- Use materialized views for frequently accessed data
- Implement efficient data partitioning

#### Implementation Examples

**1. Updated Toxicity Schema**

```sql
-- toxicity_schema.sql

-- Base toxicity table with common fields
CREATE TABLE IF NOT EXISTS toxicity_data (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES molecules(id) ON DELETE CASCADE,
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
ALTER TABLE toxicity_data ENABLE ROW LEVEL SECURITY;

-- Toxicity data policies
CREATE POLICY "users_can_view_toxicity_data" ON toxicity_data
FOR SELECT TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM molecules m
    WHERE m.id = toxicity_data.molecule_id AND
    (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_molecules pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
      )
    )
  )
);

-- Service role policy
CREATE POLICY "service_role_all_access" ON toxicity_data
FOR ALL TO service_role
USING (true) WITH CHECK (true);

-- Create specialized toxicity tables for different data types
CREATE TABLE IF NOT EXISTS toxicity_ld50 (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    toxicity_data_id UUID NOT NULL REFERENCES toxicity_data(id) ON DELETE CASCADE,
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
    SELECT 1 FROM toxicity_data td
    WHERE td.id = toxicity_ld50.toxicity_data_id AND
    EXISTS (
      SELECT 1 FROM molecules m
      WHERE m.id = td.molecule_id AND
      (
        m.is_public = true OR
        m.created_by = auth.uid() OR
        EXISTS (
          SELECT 1 FROM project_molecules pm
          JOIN team_projects tp ON pm.project_id = tp.project_id
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
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    toxicity_data_id UUID NOT NULL REFERENCES toxicity_data(id) ON DELETE CASCADE,
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
    SELECT 1 FROM toxicity_data td
    WHERE td.id = toxicity_tox21.toxicity_data_id AND
    EXISTS (
      SELECT 1 FROM molecules m
      WHERE m.id = td.molecule_id AND
      (
        m.is_public = true OR
        m.created_by = auth.uid() OR
        EXISTS (
          SELECT 1 FROM project_molecules pm
          JOIN team_projects tp ON pm.project_id = tp.project_id
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
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES molecules(id) ON DELETE CASCADE,
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
    SELECT 1 FROM molecules m
    WHERE m.id = toxicity_classification.molecule_id AND
    (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_molecules pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
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
```

### 2.2 Performance Indexes

#### Problem
Missing indexes on key columns used in filtering and joining toxicity data tables.

#### Solution
Create strategic indexes on columns used in toxicity data queries:
- Add indexes for foreign keys
- Add indexes for frequently filtered columns
- Create composite indexes for common query patterns
- Implement partial indexes for selective queries

#### Implementation Examples

```sql
-- toxicity_indexes.sql

-- Base toxicity table indexes
CREATE INDEX IF NOT EXISTS idx_toxicity_data_molecule_id ON toxicity_data(molecule_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_source ON toxicity_data(source);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_toxicity_type ON toxicity_data(toxicity_type);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_species ON toxicity_data(species);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_is_predicted ON toxicity_data(is_predicted);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_created_by ON toxicity_data(created_by);

-- Composite indexes for common query patterns
CREATE INDEX IF NOT EXISTS idx_toxicity_data_molecule_type ON toxicity_data(molecule_id, toxicity_type);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_molecule_source ON toxicity_data(molecule_id, source);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_type_species ON toxicity_data(toxicity_type, species);

-- Partial index for predicted values
CREATE INDEX IF NOT EXISTS idx_toxicity_data_predicted ON toxicity_data(molecule_id, toxicity_type, value)
WHERE is_predicted = true;

-- Partial index for experimental values
CREATE INDEX IF NOT EXISTS idx_toxicity_data_experimental ON toxicity_data(molecule_id, toxicity_type, value)
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
```

### 2.3 Materialized Views

#### Problem
Frequently accessed toxicity data requires complex joins and calculations, impacting performance.

#### Solution
Create materialized views for commonly accessed toxicity data:
- Create summaries of toxicity data by molecule
- Implement pre-calculated views for toxicity screening
- Set up refresh schedules for views
- Add indexes to materialized views

#### Implementation Examples

```sql
-- toxicity_materialized_views.sql

-- Molecule toxicity summary view
CREATE MATERIALIZED VIEW IF NOT EXISTS toxicity_summary AS
SELECT 
    m.id AS molecule_id,
    m.name,
    m.smiles,
    m.molecular_formula,
    m.cid,
    m.pubchem_link,
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
        JOIN toxicity_data td2 ON tt.toxicity_data_id = td2.id
        WHERE td2.molecule_id = m.id
    ) AS has_tox21_data,
    EXISTS (
        SELECT 1 FROM toxicity_ld50 tl
        JOIN toxicity_data td3 ON tl.toxicity_data_id = td3.id
        WHERE td3.molecule_id = m.id
    ) AS has_ld50_data
FROM 
    molecules m
LEFT JOIN 
    toxicity_data td ON m.id = td.molecule_id
GROUP BY 
    m.id;

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
    toxicity_data td ON tt.toxicity_data_id = td.id
JOIN 
    molecules m ON td.molecule_id = m.id
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
    toxicity_data td ON tl.toxicity_data_id = td.id
JOIN 
    molecules m ON td.molecule_id = m.id
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
    molecules m ON tc.molecule_id = m.id;

-- Create indexes on the hazard classification materialized view
CREATE INDEX IF NOT EXISTS idx_hazard_summary_molecule_id 
ON hazard_classification_summary(molecule_id);

CREATE INDEX IF NOT EXISTS idx_hazard_summary_system 
ON hazard_classification_summary(classification_system);

CREATE INDEX IF NOT EXISTS idx_hazard_summary_hazard_class 
ON hazard_classification_summary(hazard_class);

-- Create refresh function and schedule
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

-- Create table to track refresh history if it doesn't exist
CREATE TABLE IF NOT EXISTS materialized_view_refresh_log (
    id SERIAL PRIMARY KEY,
    view_name TEXT NOT NULL,
    refresh_time TIMESTAMP WITH TIME ZONE NOT NULL
);

-- Schedule refresh (if using pg_cron extension)
-- This requires the pg_cron extension to be installed
DO $$
BEGIN
    IF EXISTS (
        SELECT 1 FROM pg_extension WHERE extname = 'pg_cron'
    ) THEN
        PERFORM cron.schedule('refresh_toxicity_views', '0 */3 * * *', 'SELECT refresh_toxicity_materialized_views()');
    END IF;
END $$;
```

### 2.4 Query Optimization

#### Problem
Current queries for toxicity data are not optimized for performance, resulting in slow response times.

#### Solution
Optimize common query patterns:
- Use the materialized views for common queries
- Optimize complex queries with CTEs and subqueries
- Implement efficient filtering strategies
- Leverage database functions for complex calculations

#### Implementation Examples

**1. Optimized Queries**

```sql
-- Get toxicity summary for a molecule
-- Before: Complex JOIN across multiple tables
-- After: Simple query against materialized view
SELECT * FROM toxicity_summary
WHERE molecule_id = '123e4567-e89b-12d3-a456-426614174000';

-- Get all active compounds in Tox21 assays
-- Before: Multiple JOINs with filtering
-- After: Simple query against materialized view with efficient filtering
SELECT molecule_id, molecule_name, smiles, assay_id, assay_target
FROM tox21_activity_summary
WHERE activity_outcome = 'Active'
ORDER BY assay_target, molecule_name;

-- Get LD50 values for rat, oral route
-- Before: Complex JOINs with multiple filters
-- After: Simple query against materialized view
SELECT molecule_id, molecule_name, smiles, ld50_value, unit, is_predicted, source
FROM ld50_summary
WHERE species = 'Rat' AND route_of_administration = 'Oral'
ORDER BY ld50_value ASC;

-- Get GHS classifications for a molecule
-- Before: JOIN with multiple filters
-- After: Simple query against materialized view
SELECT *
FROM hazard_classification_summary
WHERE molecule_id = '123e4567-e89b-12d3-a456-426614174000'
AND classification_system = 'GHS';
```

**2. Stored Functions for Complex Calculations**

```sql
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
        FROM molecules m
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
```

### 2.5 Client-Side Caching & API Optimizations

#### Problem
Client applications make repeated requests for the same toxicity data, increasing server load.

#### Solution
Implement efficient client-side caching and API optimizations:
- Add proper caching headers to API responses
- Implement ETags for efficient caching
- Create specialized endpoints for common toxicity data requests
- Implement bulk endpoints for retrieving multiple molecules' toxicity data

#### Implementation Examples

**1. REST API Endpoints**

```python
# toxicity_resources.py
from flask import Blueprint, jsonify, request, make_response
from flask_restful import Resource, Api
from .enhanced_jwt_auth import jwt_required
import hashlib
import time
import functools

toxicity_bp = Blueprint('toxicity', __name__)
api = Api(toxicity_bp)

# Cache control decorator
def cache_control(max_age=3600):
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            response = make_response(f(*args, **kwargs))
            response.headers['Cache-Control'] = f'public, max-age={max_age}'
            return response
        return wrapper
    return decorator

# ETag generator
def generate_etag(data):
    return hashlib.md5(str(data).encode('utf-8')).hexdigest()

# Optimized toxicity summary endpoint
class ToxicitySummaryResource(Resource):
    @jwt_required()
    @cache_control(max_age=3600)  # 1 hour cache
    def get(self, molecule_id):
        # Generate ETag based on last update time of toxicity data
        last_update = db.execute_query(
            "SELECT MAX(updated_at) FROM toxicity_data WHERE molecule_id = %s",
            [molecule_id]
        )
        
        etag = generate_etag(f"{molecule_id}:{last_update}")
        if request.headers.get('If-None-Match') == etag:
            return '', 304  # Not Modified
        
        # Use materialized view for efficient query
        summary = db.execute_query(
            "SELECT * FROM toxicity_summary WHERE molecule_id = %s",
            [molecule_id]
        )
        
        response = make_response(jsonify(summary))
        response.headers['ETag'] = etag
        return response

# Bulk toxicity data endpoint
class BulkToxicityResource(Resource):
    @jwt_required()
    def post(self):
        data = request.get_json()
        molecule_ids = data.get('molecule_ids', [])
        
        if not molecule_ids:
            return {'error': 'No molecule IDs provided'}, 400
            
        if len(molecule_ids) > 100:
            return {'error': 'Maximum of 100 molecule IDs allowed'}, 400
        
        # Use optimized query with IN clause
        results = db.execute_query(
            "SELECT * FROM toxicity_summary WHERE molecule_id = ANY(%s)",
            [molecule_ids]
        )
        
        return jsonify(results)

# Toxicity score endpoint
class ToxicityScoreResource(Resource):
    @jwt_required()
    @cache_control(max_age=7200)  # 2 hour cache
    def get(self, molecule_id):
        # Use database function for efficient calculation
        score = db.execute_query(
            "SELECT calculate_toxicity_score(%s) AS toxicity_score",
            [molecule_id]
        )
        
        etag = generate_etag(f"{molecule_id}:score:{score}")
        if request.headers.get('If-None-Match') == etag:
            return '', 304  # Not Modified
        
        response = make_response(jsonify({'toxicity_score': score}))
        response.headers['ETag'] = etag
        return response

# Similar toxicity profiles endpoint
class SimilarToxicityResource(Resource):
    @jwt_required()
    def get(self, molecule_id):
        limit = int(request.args.get('limit', 10))
        
        # Use database function for efficient similarity calculation
        similar = db.execute_query(
            "SELECT * FROM find_similar_toxicity_profiles(%s, %s)",
            [molecule_id, limit]
        )
        
        return jsonify(similar)

# Register resources
api.add_resource(ToxicitySummaryResource, '/toxicity/summary/<uuid:molecule_id>')
api.add_resource(BulkToxicityResource, '/toxicity/bulk')
api.add_resource(ToxicityScoreResource, '/toxicity/score/<uuid:molecule_id>')
api.add_resource(SimilarToxicityResource, '/toxicity/similar/<uuid:molecule_id>')
```

## 3. Implementation Guide

### 3.1 Step-by-Step Approach

1. **Schema Optimization**
   - Create and apply the new toxicity schema
   - Migrate existing data to the new schema
   - Enable RLS on all new tables

2. **Performance Indexes**
   - Create indexes for all tables
   - Test query performance
   - Adjust indexes based on usage patterns

3. **Materialized Views**
   - Create materialized views for common data access patterns
   - Set up refresh schedules
   - Create indexes on materialized views

4. **Query Optimization**
   - Create stored functions for complex calculations
   - Optimize API queries to use materialized views
   - Implement efficient filtering and joining

5. **API Optimization**
   - Implement caching and ETags
   - Create bulk endpoints
   - Update API documentation

### 3.2 Implementation Script

```sql
-- Begin transaction
BEGIN;

-- Step 1: Create schema
-- (toxicity_schema.sql contents)

-- Step 2: Create indexes
-- (toxicity_indexes.sql contents)

-- Step 3: Create materialized views
-- (toxicity_materialized_views.sql contents)

-- Migrate existing data if needed
-- This assumes there's an existing toxicity_data_old table
INSERT INTO toxicity_data (
    molecule_id, source, source_id, toxicity_type, species,
    route_of_administration, value, unit, uncertainty,
    confidence_score, is_predicted, prediction_method,
    reference_doi, created_at, created_by, updated_at
)
SELECT 
    molecule_id, source, source_id, toxicity_type, species,
    route_of_administration, value, unit, uncertainty,
    confidence_score, is_predicted, prediction_method,
    reference_doi, created_at, created_by, updated_at
FROM toxicity_data_old;

-- Commit transaction
COMMIT;

-- Refresh materialized views
SELECT refresh_toxicity_materialized_views();
```

### 3.3 Testing Recommendations

1. **Performance Testing**
   - Benchmark query performance before and after optimization
   - Test with varying data volumes
   - Test with different query patterns
   - Measure response times for API endpoints

2. **Data Integrity Testing**
   - Verify data migration is complete and accurate
   - Test RLS policy effectiveness
   - Verify materialized view data is accurate
   - Test stored functions return correct results

3. **API Testing**
   - Test API endpoint response times
   - Verify caching works correctly
   - Test ETag functionality
   - Test bulk endpoints with varying payloads

## 4. Expected Performance Improvements

### Before Optimization
- Average query time for toxicity data: 300-500ms
- Bulk retrieval (100 molecules): 3-5 seconds
- Database CPU utilization during toxicity queries: 50-70%
- Caching effectiveness: Minimal

### After Optimization
- Average query time for toxicity data: 30-50ms (85-90% improvement)
- Bulk retrieval (100 molecules): 300-500ms (90% improvement)
- Database CPU utilization during toxicity queries: 10-20%
- Caching effectiveness: High (cache hit rate >80%)

## 5. Performance Monitoring

### Key Metrics to Monitor
- Query execution time for toxicity-related queries
- Materialized view refresh time
- Cache hit rate for toxicity endpoints
- Database CPU utilization during peak load
- API response times for toxicity endpoints

### Monitoring Implementation

```python
# monitoring_toxicity.py
import time
import logging
import functools
from flask import g, request, current_app

logger = logging.getLogger(__name__)

# Query performance tracking decorator
def track_query_performance(name):
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            result = f(*args, **kwargs)
            execution_time = time.time() - start_time
            
            # Log query performance
            logger.info(f"Query {name} executed in {execution_time:.4f}s")
            
            # Record in Prometheus metrics (if available)
            if hasattr(current_app, 'metrics'):
                current_app.metrics.query_execution_time.labels(
                    query_name=name
                ).observe(execution_time)
                
            return result
        return wrapper
    return decorator

# API endpoint performance tracking middleware
def track_api_performance():
    start_time = time.time()
    
    # Register teardown callback
    def teardown(exception):
        execution_time = time.time() - start_time
        endpoint = request.endpoint
        method = request.method
        
        # Only track toxicity endpoints
        if endpoint and 'toxicity' in endpoint:
            # Log performance
            logger.info(f"API {method} {endpoint} responded in {execution_time:.4f}s")
            
            # Record in Prometheus metrics (if available)
            if hasattr(current_app, 'metrics'):
                current_app.metrics.api_response_time.labels(
                    endpoint=endpoint,
                    method=method
                ).observe(execution_time)
                
                # Track cache hits
                if request.headers.get('If-None-Match') and getattr(g, 'cache_hit', False):
                    current_app.metrics.cache_hits.labels(
                        endpoint=endpoint
                    ).inc()
    
    return teardown
```

## 6. Implementation Risks and Mitigations

### Risks
- **Data migration issues**: Data loss or corruption during schema migration
- **Query regression**: Some queries might perform worse with new schema
- **RLS policy errors**: Access control issues with new tables
- **Materialized view staleness**: Views may become outdated between refreshes

### Mitigations
- Create thorough backups before migration
- Test all query patterns before and after changes
- Verify RLS policies thoroughly
- Set appropriate refresh schedules for views
- Implement monitoring for early detection of issues

## 7. Conclusion

Optimizing the toxicity data schema, creating efficient indexes, implementing materialized views, and enhancing API endpoints will significantly improve the performance and scalability of toxicity data handling in the CryoProtect system. These optimizations will enable faster toxicity screening, more efficient data retrieval, and better user experience when working with toxicity data.