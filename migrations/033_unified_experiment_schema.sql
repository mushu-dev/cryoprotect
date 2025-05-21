-- Migration 033: Unified Experiment Schema
-- Implements the enhanced experimental data schema

-- Enable Row Level Security (RLS)
ALTER DATABASE postgres SET row_security TO ON;

--------------------------------------------------------------------------
-- EXPERIMENTAL SCHEMA ENHANCEMENT
--------------------------------------------------------------------------

-- 1. Create protocols table
CREATE TABLE IF NOT EXISTS protocols (
    id UUID PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    steps JSONB NOT NULL,
    parameters JSONB,
    version INTEGER NOT NULL DEFAULT 1,
    parent_id UUID REFERENCES protocols(id),
    created_by UUID REFERENCES user_profile(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

COMMENT ON TABLE protocols IS 'Template-based protocol definitions with versioning support';

-- 2. Enhance experiment_types table (if it exists, add missing columns)
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'experiment_types') THEN
        -- Add protocol_details column if it doesn't exist
        IF NOT EXISTS (SELECT FROM information_schema.columns 
                      WHERE table_name = 'experiment_types' AND column_name = 'protocol_details') THEN
            ALTER TABLE experiment_types ADD COLUMN protocol_details JSONB;
        END IF;
    ELSE
        -- Create experiment_types table
        CREATE TABLE experiment_types (
            id UUID PRIMARY KEY,
            name VARCHAR(255) NOT NULL,
            description TEXT,
            protocol_details JSONB,
            created_by UUID REFERENCES user_profile(id),
            created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
        );
    END IF;
END
$$;

COMMENT ON TABLE experiment_types IS 'Categorizes different types of cryopreservation experiments';

-- 3. Enhance tissue_types table (if it exists, add missing columns)
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'tissue_types') THEN
        -- Add properties column if it doesn't exist
        IF NOT EXISTS (SELECT FROM information_schema.columns 
                      WHERE table_name = 'tissue_types' AND column_name = 'properties') THEN
            ALTER TABLE tissue_types ADD COLUMN properties JSONB;
        END IF;
    ELSE
        -- Create tissue_types table
        CREATE TABLE tissue_types (
            id UUID PRIMARY KEY,
            name VARCHAR(255) NOT NULL,
            description TEXT,
            species VARCHAR(255),
            taxonomy_id INTEGER,
            properties JSONB,
            created_by UUID REFERENCES user_profile(id),
            created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
        );
    END IF;
END
$$;

COMMENT ON TABLE tissue_types IS 'Biological samples used in experiments';

-- 4. Create or enhance experiments table
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'experiments') THEN
        -- Add new columns if they don't exist
        IF NOT EXISTS (SELECT FROM information_schema.columns 
                      WHERE table_name = 'experiments' AND column_name = 'protocol_id') THEN
            ALTER TABLE experiments ADD COLUMN protocol_id UUID REFERENCES protocols(id);
        END IF;
        
        IF NOT EXISTS (SELECT FROM information_schema.columns 
                      WHERE table_name = 'experiments' AND column_name = 'version') THEN
            ALTER TABLE experiments ADD COLUMN version INTEGER NOT NULL DEFAULT 1;
        END IF;
        
        IF NOT EXISTS (SELECT FROM information_schema.columns 
                      WHERE table_name = 'experiments' AND column_name = 'provenance') THEN
            ALTER TABLE experiments ADD COLUMN provenance JSONB;
        END IF;
        
        -- Rename parameters column if it exists with a different name
        IF NOT EXISTS (SELECT FROM information_schema.columns 
                      WHERE table_name = 'experiments' AND column_name = 'parameters') 
           AND EXISTS (SELECT FROM information_schema.columns 
                      WHERE table_name = 'experiments' AND column_name = 'additional_parameters') THEN
            ALTER TABLE experiments RENAME COLUMN additional_parameters TO parameters;
        END IF;
    ELSE
        -- Create experiments table
        CREATE TABLE experiments (
            id UUID PRIMARY KEY,
            experiment_type_id UUID NOT NULL REFERENCES experiment_types(id),
            title VARCHAR(255) NOT NULL,
            description TEXT,
            protocol_id UUID REFERENCES protocols(id),
            date_performed DATE,
            temperature NUMERIC,
            cooling_rate NUMERIC,
            thawing_rate NUMERIC,
            parameters JSONB,
            version INTEGER NOT NULL DEFAULT 1,
            provenance JSONB,
            created_by UUID REFERENCES user_profile(id),
            created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
        );
    END IF;
END
$$;

COMMENT ON TABLE experiments IS 'Core experiment records';

-- 5. Create or enhance experiment_results table
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'experiment_results') THEN
        -- Add new columns if they don't exist
        IF NOT EXISTS (SELECT FROM information_schema.columns 
                      WHERE table_name = 'experiment_results' AND column_name = 'uncertainty') THEN
            ALTER TABLE experiment_results ADD COLUMN uncertainty JSONB;
        END IF;
        
        IF NOT EXISTS (SELECT FROM information_schema.columns 
                      WHERE table_name = 'experiment_results' AND column_name = 'provenance') THEN
            ALTER TABLE experiment_results ADD COLUMN provenance JSONB;
        END IF;
    ELSE
        -- Create experiment_results table
        CREATE TABLE experiment_results (
            id UUID PRIMARY KEY,
            experiment_id UUID NOT NULL REFERENCES experiments(id),
            molecule_id UUID REFERENCES molecules(id),
            mixture_id UUID REFERENCES mixtures(id),
            tissue_type_id UUID NOT NULL REFERENCES tissue_types(id),
            concentration NUMERIC,
            concentration_unit VARCHAR(50),
            viability_percentage NUMERIC,
            recovery_rate NUMERIC,
            functionality_score NUMERIC,
            result_details JSONB,
            uncertainty JSONB,
            provenance JSONB,
            notes TEXT,
            created_by UUID REFERENCES user_profile(id),
            created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            CONSTRAINT check_molecule_or_mixture CHECK (
                (molecule_id IS NOT NULL AND mixture_id IS NULL) OR
                (molecule_id IS NULL AND mixture_id IS NOT NULL)
            )
        );
    END IF;
END
$$;

COMMENT ON TABLE experiment_results IS 'Results from experiments on specific molecules or mixtures';

-- 6. Create equipment-related tables
CREATE TABLE IF NOT EXISTS equipment_types (
    id UUID PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    manufacturer VARCHAR(255),
    created_by UUID REFERENCES user_profile(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

COMMENT ON TABLE equipment_types IS 'Categories of laboratory equipment';

CREATE TABLE IF NOT EXISTS equipment (
    id UUID PRIMARY KEY,
    equipment_type_id UUID NOT NULL REFERENCES equipment_types(id),
    name VARCHAR(255) NOT NULL,
    model VARCHAR(255),
    description TEXT,
    parameters JSONB,
    created_by UUID REFERENCES user_profile(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

COMMENT ON TABLE equipment IS 'Specific equipment instances';

CREATE TABLE IF NOT EXISTS experiment_equipment (
    id UUID PRIMARY KEY,
    experiment_id UUID NOT NULL REFERENCES experiments(id),
    equipment_id UUID NOT NULL REFERENCES equipment(id),
    parameters JSONB,
    created_by UUID REFERENCES user_profile(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

COMMENT ON TABLE experiment_equipment IS 'Links experiments to equipment used';

-- 7. Create time-series tables
CREATE TABLE IF NOT EXISTS time_series (
    id UUID PRIMARY KEY,
    experiment_id UUID NOT NULL REFERENCES experiments(id),
    name VARCHAR(255) NOT NULL,
    description TEXT,
    property_type VARCHAR(255) NOT NULL,
    unit VARCHAR(50),
    metadata JSONB,
    created_by UUID REFERENCES user_profile(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

COMMENT ON TABLE time_series IS 'Time-series data definitions';

CREATE TABLE IF NOT EXISTS time_series_data (
    id UUID PRIMARY KEY,
    time_series_id UUID NOT NULL REFERENCES time_series(id),
    timestamp TIMESTAMPTZ NOT NULL,
    value NUMERIC NOT NULL,
    uncertainty NUMERIC,
    provenance JSONB,
    created_by UUID REFERENCES user_profile(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

COMMENT ON TABLE time_series_data IS 'Individual time-series data points';

-- 8. Create validation_rules table
CREATE TABLE IF NOT EXISTS validation_rules (
    id UUID PRIMARY KEY,
    property_type VARCHAR(255) NOT NULL,
    rule_type VARCHAR(50) NOT NULL,
    parameters JSONB NOT NULL,
    description TEXT,
    created_by UUID REFERENCES user_profile(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

COMMENT ON TABLE validation_rules IS 'Rules for validating experimental data';

--------------------------------------------------------------------------
-- CREATE INDEXES
--------------------------------------------------------------------------

-- Indexes on experiment_types
CREATE INDEX IF NOT EXISTS idx_experiment_types_name ON experiment_types (name);

-- Indexes on protocols
CREATE INDEX IF NOT EXISTS idx_protocols_name ON protocols (name);
CREATE INDEX IF NOT EXISTS idx_protocols_parent_id ON protocols (parent_id);
CREATE INDEX IF NOT EXISTS idx_protocols_version ON protocols (version);

-- Indexes on tissue_types
CREATE INDEX IF NOT EXISTS idx_tissue_types_name ON tissue_types (name);
CREATE INDEX IF NOT EXISTS idx_tissue_types_species ON tissue_types (species);

-- Indexes on experiments
CREATE INDEX IF NOT EXISTS idx_experiments_experiment_type_id ON experiments (experiment_type_id);
CREATE INDEX IF NOT EXISTS idx_experiments_protocol_id ON experiments (protocol_id);
CREATE INDEX IF NOT EXISTS idx_experiments_date_performed ON experiments (date_performed);
CREATE INDEX IF NOT EXISTS idx_experiments_title ON experiments (title);

-- Indexes on experiment_results
CREATE INDEX IF NOT EXISTS idx_experiment_results_experiment_id ON experiment_results (experiment_id);
CREATE INDEX IF NOT EXISTS idx_experiment_results_molecule_id ON experiment_results (molecule_id);
CREATE INDEX IF NOT EXISTS idx_experiment_results_mixture_id ON experiment_results (mixture_id);
CREATE INDEX IF NOT EXISTS idx_experiment_results_tissue_type_id ON experiment_results (tissue_type_id);
CREATE INDEX IF NOT EXISTS idx_experiment_results_viability ON experiment_results (viability_percentage);

-- Indexes on equipment-related tables
CREATE INDEX IF NOT EXISTS idx_equipment_types_name ON equipment_types (name);
CREATE INDEX IF NOT EXISTS idx_equipment_equipment_type_id ON equipment (equipment_type_id);
CREATE INDEX IF NOT EXISTS idx_equipment_name ON equipment (name);
CREATE INDEX IF NOT EXISTS idx_experiment_equipment_experiment_id ON experiment_equipment (experiment_id);
CREATE INDEX IF NOT EXISTS idx_experiment_equipment_equipment_id ON experiment_equipment (equipment_id);

-- Indexes on time-series tables
CREATE INDEX IF NOT EXISTS idx_time_series_experiment_id ON time_series (experiment_id);
CREATE INDEX IF NOT EXISTS idx_time_series_property_type ON time_series (property_type);
CREATE INDEX IF NOT EXISTS idx_time_series_data_time_series_id ON time_series_data (time_series_id);
CREATE INDEX IF NOT EXISTS idx_time_series_data_timestamp ON time_series_data (timestamp);

-- Indexes on validation_rules
CREATE INDEX IF NOT EXISTS idx_validation_rules_property_type ON validation_rules (property_type);
CREATE INDEX IF NOT EXISTS idx_validation_rules_rule_type ON validation_rules (rule_type);

-- JSON indexes for faster querying of JSON fields
CREATE INDEX IF NOT EXISTS idx_protocols_steps ON protocols USING GIN (steps);
CREATE INDEX IF NOT EXISTS idx_experiments_parameters ON experiments USING GIN (parameters);
CREATE INDEX IF NOT EXISTS idx_experiments_provenance ON experiments USING GIN (provenance);
CREATE INDEX IF NOT EXISTS idx_experiment_results_result_details ON experiment_results USING GIN (result_details);
CREATE INDEX IF NOT EXISTS idx_experiment_results_uncertainty ON experiment_results USING GIN (uncertainty);

--------------------------------------------------------------------------
-- CREATE VIEWS FOR COMPATIBILITY AND SIMPLIFIED QUERIES
--------------------------------------------------------------------------

-- View for experiment summary
CREATE OR REPLACE VIEW experiment_summary AS
SELECT 
    e.id AS experiment_id,
    e.title AS experiment_title,
    et.name AS experiment_type,
    e.date_performed,
    e.temperature,
    e.cooling_rate,
    e.thawing_rate,
    COUNT(er.id) AS result_count,
    AVG(er.viability_percentage) AS avg_viability,
    MIN(er.viability_percentage) AS min_viability,
    MAX(er.viability_percentage) AS max_viability
FROM 
    experiments e
JOIN 
    experiment_types et ON e.experiment_type_id = et.id
LEFT JOIN 
    experiment_results er ON e.id = er.experiment_id
GROUP BY 
    e.id, e.title, et.name, e.date_performed, e.temperature, e.cooling_rate, e.thawing_rate;

-- View for molecule experiment results
CREATE OR REPLACE VIEW molecule_experiment_results AS
SELECT 
    er.id,
    e.title AS experiment_title,
    et.name AS experiment_type,
    m.name AS molecule_name,
    m.smiles AS molecule_smiles,
    tt.name AS tissue_type,
    tt.species,
    er.concentration,
    er.concentration_unit,
    er.viability_percentage,
    er.recovery_rate,
    er.functionality_score,
    er.result_details,
    e.temperature,
    e.cooling_rate,
    e.thawing_rate,
    e.date_performed
FROM 
    experiment_results er
JOIN 
    experiments e ON er.experiment_id = e.id
JOIN 
    experiment_types et ON e.experiment_type_id = et.id
JOIN 
    tissue_types tt ON er.tissue_type_id = tt.id
JOIN 
    molecules m ON er.molecule_id = m.id
WHERE 
    er.molecule_id IS NOT NULL;

-- View for mixture experiment results
CREATE OR REPLACE VIEW mixture_experiment_results AS
SELECT 
    er.id,
    e.title AS experiment_title,
    et.name AS experiment_type,
    mx.name AS mixture_name,
    tt.name AS tissue_type,
    tt.species,
    er.concentration,
    er.concentration_unit,
    er.viability_percentage,
    er.recovery_rate,
    er.functionality_score,
    er.result_details,
    e.temperature,
    e.cooling_rate,
    e.thawing_rate,
    e.date_performed
FROM 
    experiment_results er
JOIN 
    experiments e ON er.experiment_id = e.id
JOIN 
    experiment_types et ON e.experiment_type_id = et.id
JOIN 
    tissue_types tt ON er.tissue_type_id = tt.id
JOIN 
    mixtures mx ON er.mixture_id = mx.id
WHERE 
    er.mixture_id IS NOT NULL;

-- View for protocol templates
CREATE OR REPLACE VIEW protocol_templates AS
SELECT 
    p.id,
    p.name,
    p.description,
    p.steps,
    p.parameters,
    p.version,
    COUNT(e.id) AS usage_count
FROM 
    protocols p
LEFT JOIN 
    experiments e ON p.id = e.protocol_id
GROUP BY 
    p.id, p.name, p.description, p.steps, p.parameters, p.version;

--------------------------------------------------------------------------
-- ROW LEVEL SECURITY POLICIES
--------------------------------------------------------------------------

-- Enable RLS on all new tables
ALTER TABLE protocols ENABLE ROW LEVEL SECURITY;
ALTER TABLE equipment_types ENABLE ROW LEVEL SECURITY;
ALTER TABLE equipment ENABLE ROW LEVEL SECURITY;
ALTER TABLE experiment_equipment ENABLE ROW LEVEL SECURITY;
ALTER TABLE time_series ENABLE ROW LEVEL SECURITY;
ALTER TABLE time_series_data ENABLE ROW LEVEL SECURITY;
ALTER TABLE validation_rules ENABLE ROW LEVEL SECURITY;

-- Make sure RLS is enabled on existing tables if they exist
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'experiment_types') THEN
        ALTER TABLE experiment_types ENABLE ROW LEVEL SECURITY;
    END IF;
    
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'tissue_types') THEN
        ALTER TABLE tissue_types ENABLE ROW LEVEL SECURITY;
    END IF;
    
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'experiments') THEN
        ALTER TABLE experiments ENABLE ROW LEVEL SECURITY;
    END IF;
    
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'experiment_results') THEN
        ALTER TABLE experiment_results ENABLE ROW LEVEL SECURITY;
    END IF;
END
$$;

-- Create RLS policies for public access (read-only)
CREATE POLICY protocols_public_read ON protocols
    FOR SELECT USING (true);

CREATE POLICY experiment_types_public_read ON experiment_types
    FOR SELECT USING (true);

CREATE POLICY tissue_types_public_read ON tissue_types
    FOR SELECT USING (true);

CREATE POLICY experiments_public_read ON experiments
    FOR SELECT USING (true);

CREATE POLICY experiment_results_public_read ON experiment_results
    FOR SELECT USING (true);

CREATE POLICY equipment_types_public_read ON equipment_types
    FOR SELECT USING (true);

CREATE POLICY equipment_public_read ON equipment
    FOR SELECT USING (true);

CREATE POLICY experiment_equipment_public_read ON experiment_equipment
    FOR SELECT USING (true);

CREATE POLICY time_series_public_read ON time_series
    FOR SELECT USING (true);

CREATE POLICY time_series_data_public_read ON time_series_data
    FOR SELECT USING (true);

CREATE POLICY validation_rules_public_read ON validation_rules
    FOR SELECT USING (true);

-- Create RLS policies for authenticated users (CRUD on own records)
-- Note: These assume a user_profile table with auth_user_id column linking to auth.users
CREATE POLICY protocols_auth_all ON protocols
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY experiment_types_auth_all ON experiment_types
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY tissue_types_auth_all ON tissue_types
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY experiments_auth_all ON experiments
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY experiment_results_auth_all ON experiment_results
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY equipment_types_auth_all ON equipment_types
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY equipment_auth_all ON equipment
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY experiment_equipment_auth_all ON experiment_equipment
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY time_series_auth_all ON time_series
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY time_series_data_auth_all ON time_series_data
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

CREATE POLICY validation_rules_auth_all ON validation_rules
    USING (created_by IN (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));

-- Create triggers for updated_at timestamps
CREATE OR REPLACE FUNCTION update_timestamp()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Apply triggers to all new tables
CREATE TRIGGER update_protocols_timestamp
BEFORE UPDATE ON protocols
FOR EACH ROW EXECUTE PROCEDURE update_timestamp();

CREATE TRIGGER update_equipment_types_timestamp
BEFORE UPDATE ON equipment_types
FOR EACH ROW EXECUTE PROCEDURE update_timestamp();

CREATE TRIGGER update_equipment_timestamp
BEFORE UPDATE ON equipment
FOR EACH ROW EXECUTE PROCEDURE update_timestamp();

CREATE TRIGGER update_experiment_equipment_timestamp
BEFORE UPDATE ON experiment_equipment
FOR EACH ROW EXECUTE PROCEDURE update_timestamp();

CREATE TRIGGER update_time_series_timestamp
BEFORE UPDATE ON time_series
FOR EACH ROW EXECUTE PROCEDURE update_timestamp();

CREATE TRIGGER update_time_series_data_timestamp
BEFORE UPDATE ON time_series_data
FOR EACH ROW EXECUTE PROCEDURE update_timestamp();

CREATE TRIGGER update_validation_rules_timestamp
BEFORE UPDATE ON validation_rules
FOR EACH ROW EXECUTE PROCEDURE update_timestamp();

-- Apply triggers to existing tables if they exist
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'experiment_types') 
       AND NOT EXISTS (SELECT FROM information_schema.triggers WHERE trigger_name = 'update_experiment_types_timestamp') THEN
        CREATE TRIGGER update_experiment_types_timestamp
        BEFORE UPDATE ON experiment_types
        FOR EACH ROW EXECUTE PROCEDURE update_timestamp();
    END IF;
    
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'tissue_types')
       AND NOT EXISTS (SELECT FROM information_schema.triggers WHERE trigger_name = 'update_tissue_types_timestamp') THEN
        CREATE TRIGGER update_tissue_types_timestamp
        BEFORE UPDATE ON tissue_types
        FOR EACH ROW EXECUTE PROCEDURE update_timestamp();
    END IF;
    
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'experiments')
       AND NOT EXISTS (SELECT FROM information_schema.triggers WHERE trigger_name = 'update_experiments_timestamp') THEN
        CREATE TRIGGER update_experiments_timestamp
        BEFORE UPDATE ON experiments
        FOR EACH ROW EXECUTE PROCEDURE update_timestamp();
    END IF;
    
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'experiment_results')
       AND NOT EXISTS (SELECT FROM information_schema.triggers WHERE trigger_name = 'update_experiment_results_timestamp') THEN
        CREATE TRIGGER update_experiment_results_timestamp
        BEFORE UPDATE ON experiment_results
        FOR EACH ROW EXECUTE PROCEDURE update_timestamp();
    END IF;
END
$$;