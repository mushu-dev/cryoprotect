-- Migration: 030_scientific_models_schema
-- Description: Creates tables and indexes for Phase 4 scientific models

-- Create mixture_optimizations table for storing optimization results
CREATE TABLE IF NOT EXISTS mixture_optimizations (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    components JSONB NOT NULL, -- molecules and their concentrations
    optimization_parameters JSONB NOT NULL, -- parameters used for optimization
    effectiveness_score DOUBLE PRECISION,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by TEXT,
    notes TEXT,
    
    -- Add constraint to ensure components is a valid JSON array
    CONSTRAINT valid_components_json CHECK (jsonb_typeof(components) = 'array')
);

-- Create index on components for faster lookup
CREATE INDEX idx_mixture_optimizations_components ON mixture_optimizations USING GIN (components);

-- Create concentration_models table for storing concentration-dependent model data
CREATE TABLE IF NOT EXISTS concentration_models (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES molecules(id),
    model_type TEXT NOT NULL, -- type of concentration model
    parameters JSONB NOT NULL, -- model parameters
    valid_range NUMRANGE, -- valid concentration range
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by TEXT,
    
    -- Ensure only one model of each type per molecule
    CONSTRAINT unique_molecule_model_type UNIQUE (molecule_id, model_type)
);

-- Create index on molecule_id for faster lookup
CREATE INDEX idx_concentration_models_molecule_id ON concentration_models(molecule_id);

-- Create temperature_models table for storing temperature-dependent model data
CREATE TABLE IF NOT EXISTS temperature_models (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES molecules(id),
    model_type TEXT NOT NULL, -- type of temperature model
    parameters JSONB NOT NULL, -- model parameters
    temperature_range NUMRANGE, -- valid temperature range
    glass_transition_temp DOUBLE PRECISION, -- Tg value in Kelvin
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by TEXT,
    
    -- Ensure only one model of each type per molecule
    CONSTRAINT unique_molecule_temp_model_type UNIQUE (molecule_id, model_type)
);

-- Create index on molecule_id for faster lookup
CREATE INDEX idx_temperature_models_molecule_id ON temperature_models(molecule_id);

-- Create tissue_compatibility table for storing tissue-specific compatibility data
CREATE TABLE IF NOT EXISTS tissue_compatibility (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES molecules(id),
    tissue_type TEXT NOT NULL, -- type of tissue
    compatibility_score DOUBLE PRECISION,
    penetration_rate DOUBLE PRECISION, -- rate of cellular penetration
    toxicity_score DOUBLE PRECISION, -- tissue-specific toxicity
    data_source TEXT, -- source of the data (experimental, predicted, etc.)
    notes TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by TEXT,
    
    -- Ensure only one entry per molecule-tissue combination
    CONSTRAINT unique_molecule_tissue UNIQUE (molecule_id, tissue_type)
);

-- Create index on molecule_id for faster lookup
CREATE INDEX idx_tissue_compatibility_molecule_id ON tissue_compatibility(molecule_id);
CREATE INDEX idx_tissue_compatibility_tissue_type ON tissue_compatibility(tissue_type);

-- Create view that combines scientific model data for a molecule
CREATE OR REPLACE VIEW molecule_scientific_data AS
SELECT 
    m.id as molecule_id,
    m.name as molecule_name,
    m.smiles,
    m.pubchem_cid,
    m.chembl_id,
    jsonb_agg(DISTINCT cm.model_type) FILTER (WHERE cm.id IS NOT NULL) as concentration_models,
    jsonb_agg(DISTINCT tm.model_type) FILTER (WHERE tm.id IS NOT NULL) as temperature_models,
    jsonb_agg(DISTINCT tc.tissue_type) FILTER (WHERE tc.id IS NOT NULL) as compatible_tissues,
    CASE 
        WHEN count(DISTINCT tc.id) > 0 THEN 
            avg(tc.compatibility_score) FILTER (WHERE tc.id IS NOT NULL)
        ELSE NULL
    END as avg_tissue_compatibility,
    CASE
        WHEN count(DISTINCT tm.id) > 0 AND tm.glass_transition_temp IS NOT NULL THEN
            min(tm.glass_transition_temp) FILTER (WHERE tm.id IS NOT NULL)
        ELSE NULL
    END as min_glass_transition_temp
FROM 
    molecules m
LEFT JOIN 
    concentration_models cm ON m.id = cm.molecule_id
LEFT JOIN 
    temperature_models tm ON m.id = tm.molecule_id
LEFT JOIN 
    tissue_compatibility tc ON m.id = tc.molecule_id
GROUP BY 
    m.id, m.name, m.smiles, m.pubchem_cid, m.chembl_id;

-- Create materialized view for efficient scientific data retrieval
CREATE MATERIALIZED VIEW mv_molecule_scientific_properties AS
SELECT 
    m.id as molecule_id,
    m.name as molecule_name,
    m.smiles,
    m.pubchem_cid,
    m.chembl_id,
    m.formula,
    m.molecular_weight,
    (SELECT jsonb_agg(jsonb_build_object(
        'model_type', cm.model_type,
        'valid_range', cm.valid_range::text,
        'parameters', cm.parameters
    )) FROM concentration_models cm WHERE cm.molecule_id = m.id) as concentration_models,
    
    (SELECT jsonb_agg(jsonb_build_object(
        'model_type', tm.model_type,
        'temperature_range', tm.temperature_range::text,
        'glass_transition_temp', tm.glass_transition_temp,
        'parameters', tm.parameters
    )) FROM temperature_models tm WHERE tm.molecule_id = m.id) as temperature_models,
    
    (SELECT jsonb_agg(jsonb_build_object(
        'tissue_type', tc.tissue_type,
        'compatibility_score', tc.compatibility_score,
        'penetration_rate', tc.penetration_rate,
        'toxicity_score', tc.toxicity_score
    )) FROM tissue_compatibility tc WHERE tc.molecule_id = m.id) as tissue_compatibility,
    
    (SELECT count(*) > 0 FROM mixture_optimizations mo WHERE mo.components @> jsonb_build_array(jsonb_build_object('molecule_id', m.id::text))) as used_in_optimizations
FROM 
    molecules m;

-- Create index on the materialized view for faster lookups
CREATE UNIQUE INDEX idx_mv_molecule_scientific_properties_id ON mv_molecule_scientific_properties(molecule_id);

-- Create view that expands mixture components with their properties
CREATE OR REPLACE VIEW mv_mixture_components_expanded AS
SELECT 
    mix.id as mixture_id,
    mix.name as mixture_name,
    comp.molecule_id,
    m.name as molecule_name,
    m.smiles,
    comp.concentration,
    comp.concentration_units,
    msp.concentration_models,
    msp.temperature_models,
    msp.tissue_compatibility
FROM 
    mixtures mix
JOIN 
    mixture_components comp ON mix.id = comp.mixture_id
JOIN 
    molecules m ON comp.molecule_id = m.id
LEFT JOIN 
    mv_molecule_scientific_properties msp ON m.id = msp.molecule_id;

-- Create function to refresh materialized view automatically
CREATE OR REPLACE FUNCTION refresh_molecule_scientific_properties()
RETURNS TRIGGER AS $$
BEGIN
    REFRESH MATERIALIZED VIEW CONCURRENTLY mv_molecule_scientific_properties;
    RETURN NULL;
END;
$$ LANGUAGE plpgsql;

-- Create triggers to refresh the materialized view when data changes
DROP TRIGGER IF EXISTS refresh_mv_on_concentration_models ON concentration_models;
CREATE TRIGGER refresh_mv_on_concentration_models
AFTER INSERT OR UPDATE OR DELETE ON concentration_models
FOR EACH STATEMENT EXECUTE FUNCTION refresh_molecule_scientific_properties();

DROP TRIGGER IF EXISTS refresh_mv_on_temperature_models ON temperature_models;
CREATE TRIGGER refresh_mv_on_temperature_models
AFTER INSERT OR UPDATE OR DELETE ON temperature_models
FOR EACH STATEMENT EXECUTE FUNCTION refresh_molecule_scientific_properties();

DROP TRIGGER IF EXISTS refresh_mv_on_tissue_compatibility ON tissue_compatibility;
CREATE TRIGGER refresh_mv_on_tissue_compatibility
AFTER INSERT OR UPDATE OR DELETE ON tissue_compatibility
FOR EACH STATEMENT EXECUTE FUNCTION refresh_molecule_scientific_properties();

-- Create custom types for scientific models
CREATE TYPE concentration_model_type AS ENUM ('linear', 'exponential', 'sigmoid', 'custom');
CREATE TYPE temperature_model_type AS ENUM ('linear', 'arrhenius', 'piecewise', 'custom');
CREATE TYPE tissue_type AS ENUM (
    'blood', 'liver', 'kidney', 'heart', 'lung', 
    'brain', 'skin', 'muscle', 'bone', 'pancreas',
    'cell_culture', 'sperm', 'oocyte', 'embryo', 'other'
);

-- Add comment to explain the tables and their relationships
COMMENT ON TABLE mixture_optimizations IS 
'Stores the results of mixture optimization algorithms, including the components and their concentrations';

COMMENT ON TABLE concentration_models IS 
'Contains concentration-dependent models for molecules, describing how properties change with concentration';

COMMENT ON TABLE temperature_models IS 
'Contains temperature-dependent models for molecules, including glass transition temperatures';

COMMENT ON TABLE tissue_compatibility IS 
'Stores tissue-specific compatibility data for molecules, including penetration rates and toxicity';

COMMENT ON MATERIALIZED VIEW mv_molecule_scientific_properties IS 
'Provides efficient access to all scientific properties of molecules, including concentration, temperature, and tissue models';

-- Post-migration verification
DO $$
BEGIN
    RAISE NOTICE 'Scientific models schema migration complete';
    RAISE NOTICE 'Tables created: mixture_optimizations, concentration_models, temperature_models, tissue_compatibility';
    RAISE NOTICE 'Views created: molecule_scientific_data, mv_mixture_components_expanded';
    RAISE NOTICE 'Materialized view created: mv_molecule_scientific_properties';
    RAISE NOTICE 'Function created: refresh_molecule_scientific_properties';
    RAISE NOTICE 'Custom types created: concentration_model_type, temperature_model_type, tissue_type';
    
    -- Check if all tables were created
    IF EXISTS (SELECT 1 FROM pg_tables WHERE tablename = 'mixture_optimizations') AND
       EXISTS (SELECT 1 FROM pg_tables WHERE tablename = 'concentration_models') AND
       EXISTS (SELECT 1 FROM pg_tables WHERE tablename = 'temperature_models') AND
       EXISTS (SELECT 1 FROM pg_tables WHERE tablename = 'tissue_compatibility') THEN
        RAISE NOTICE 'All tables created successfully';
    ELSE
        RAISE WARNING 'Some tables may not have been created correctly';
    END IF;
END $$;