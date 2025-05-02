-- CryoProtect Analyzer: Toxicity Database Schema Migration
-- Migration: 012_toxicity_schema.sql
-- This migration implements the schema for toxicity data integration from Tox21.
-- It is idempotent and safe to run multiple times.

-- =========================
-- 1. Toxicity Data Source Table
-- =========================
CREATE TABLE IF NOT EXISTS toxicity_data_source (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT NOT NULL UNIQUE,
    description TEXT,
    version TEXT,
    url TEXT,
    last_updated TIMESTAMPTZ DEFAULT now(),
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

-- Insert standard toxicity data sources
INSERT INTO toxicity_data_source (name, description, url)
VALUES
    ('Tox21', 'Toxicology in the 21st Century federal collaboration', 'https://ntp.niehs.nih.gov/whatwestudy/tox21/index.html')
ON CONFLICT (name) DO NOTHING;

-- =========================
-- 2. Toxicity Assay Table
-- =========================
CREATE TABLE IF NOT EXISTS toxicity_assay (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    source_id uuid REFERENCES toxicity_data_source(id),
    assay_name TEXT NOT NULL,
    assay_id TEXT NOT NULL,
    description TEXT,
    assay_type TEXT,
    organism TEXT,
    tissue TEXT,
    cell_line TEXT,
    assay_target TEXT,
    assay_function_type TEXT,
    detection_technology TEXT,
    key_assay_reagent TEXT,
    assay_footprint TEXT,
    assay_format TEXT,
    assay_design_type TEXT,
    biological_process_target TEXT,
    detection_technology_type TEXT,
    signal_direction TEXT,
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now(),
    UNIQUE (source_id, assay_id)
);

CREATE INDEX IF NOT EXISTS idx_toxicity_assay_source_id ON toxicity_assay (source_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_assay_assay_id ON toxicity_assay (assay_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_assay_assay_name ON toxicity_assay (assay_name);

-- =========================
-- 3. Molecule Identifier Mapping Table
-- =========================
CREATE TABLE IF NOT EXISTS molecule_identifier_mapping (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    molecule_id uuid REFERENCES molecule(id),
    source_id uuid REFERENCES toxicity_data_source(id),
    external_id TEXT NOT NULL,
    external_id_type TEXT NOT NULL, -- e.g., 'CASRN', 'DTXSID', 'ChEMBL_ID'
    confidence_score DECIMAL(5,4), -- 0-1 score for mapping confidence
    mapping_method TEXT, -- e.g., 'exact_match', 'inchikey_match', 'name_match'
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now(),
    UNIQUE (molecule_id, source_id, external_id_type, external_id)
);

CREATE INDEX IF NOT EXISTS idx_molecule_identifier_mapping_molecule_id ON molecule_identifier_mapping (molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecule_identifier_mapping_external_id ON molecule_identifier_mapping (external_id);
CREATE INDEX IF NOT EXISTS idx_molecule_identifier_mapping_source_id ON molecule_identifier_mapping (source_id);

-- =========================
-- 4. Toxicity Data Table
-- =========================
CREATE TABLE IF NOT EXISTS toxicity_data (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    molecule_id uuid REFERENCES molecule(id),
    assay_id uuid REFERENCES toxicity_assay(id),
    mapping_id uuid REFERENCES molecule_identifier_mapping(id),
    activity_value DECIMAL(18,8),
    activity_unit TEXT,
    activity_type TEXT, -- e.g., 'AC50', 'IC50', 'EC50'
    hit_call BOOLEAN,
    significance DECIMAL(5,4), -- p-value or similar
    reliability_score DECIMAL(5,4), -- 0-1 score for data reliability
    data_quality_comment TEXT,
    citation TEXT,
    created_by uuid REFERENCES user_profile(id),
    data_source TEXT,
    version INTEGER DEFAULT 1,
    modification_history JSONB DEFAULT '[]',
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS idx_toxicity_data_molecule_id ON toxicity_data (molecule_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_assay_id ON toxicity_data (assay_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_data_mapping_id ON toxicity_data (mapping_id);

-- =========================
-- 5. Toxicity Score Table (for aggregated/calculated toxicity scores)
-- =========================
CREATE TABLE IF NOT EXISTS toxicity_score (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    molecule_id uuid REFERENCES molecule(id),
    score_type TEXT NOT NULL, -- e.g., 'overall_toxicity', 'developmental_toxicity', 'carcinogenicity'
    score_value DECIMAL(18,8) NOT NULL,
    score_unit TEXT,
    confidence DECIMAL(5,4),
    method_id uuid REFERENCES calculation_method(id),
    assays_used JSONB, -- List of assay IDs used in calculation
    created_by uuid REFERENCES user_profile(id),
    data_source TEXT,
    version INTEGER DEFAULT 1,
    modification_history JSONB DEFAULT '[]',
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now(),
    UNIQUE (molecule_id, score_type, method_id)
);

CREATE INDEX IF NOT EXISTS idx_toxicity_score_molecule_id ON toxicity_score (molecule_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_score_score_type ON toxicity_score (score_type);

-- =========================
-- 6. Toxicity Endpoint Table (for categorized toxicity endpoints)
-- =========================
CREATE TABLE IF NOT EXISTS toxicity_endpoint (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT NOT NULL UNIQUE,
    description TEXT,
    category TEXT, -- e.g., 'developmental', 'carcinogenicity', 'endocrine'
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

-- Insert common toxicity endpoints
INSERT INTO toxicity_endpoint (name, description, category)
VALUES 
    ('Acute Toxicity', 'Adverse effects occurring within a short time of exposure', 'systemic'),
    ('Carcinogenicity', 'Ability to cause cancer', 'carcinogenicity'),
    ('Developmental Toxicity', 'Adverse effects on developing organism', 'developmental'),
    ('Endocrine Disruption', 'Interference with hormone systems', 'endocrine'),
    ('Genotoxicity', 'Damage to genetic material', 'genetic'),
    ('Hepatotoxicity', 'Liver toxicity', 'organ'),
    ('Neurotoxicity', 'Nervous system toxicity', 'organ'),
    ('Reproductive Toxicity', 'Adverse effects on reproductive function', 'reproductive'),
    ('Skin Sensitization', 'Allergic response following skin contact', 'dermal'),
    ('Cytotoxicity', 'Toxicity to cells', 'cellular')
ON CONFLICT (name) DO NOTHING;

-- =========================
-- 7. Assay-Endpoint Mapping Table
-- =========================
CREATE TABLE IF NOT EXISTS assay_endpoint_mapping (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    assay_id uuid REFERENCES toxicity_assay(id),
    endpoint_id uuid REFERENCES toxicity_endpoint(id),
    relevance_score DECIMAL(5,4), -- 0-1 score for how relevant this assay is to the endpoint
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now(),
    UNIQUE (assay_id, endpoint_id)
);

CREATE INDEX IF NOT EXISTS idx_assay_endpoint_mapping_assay_id ON assay_endpoint_mapping (assay_id);
CREATE INDEX IF NOT EXISTS idx_assay_endpoint_mapping_endpoint_id ON assay_endpoint_mapping (endpoint_id);

-- =========================
-- 8. Toxicity Import Job Table (for tracking data imports)
-- =========================
CREATE TABLE IF NOT EXISTS toxicity_import_job (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    source_id uuid REFERENCES toxicity_data_source(id),
    status TEXT NOT NULL, -- 'pending', 'in_progress', 'completed', 'failed'
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ,
    records_processed INTEGER DEFAULT 0,
    records_imported INTEGER DEFAULT 0,
    records_failed INTEGER DEFAULT 0,
    error_log JSONB DEFAULT '[]',
    created_by uuid REFERENCES user_profile(id),
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS idx_toxicity_import_job_source_id ON toxicity_import_job (source_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_import_job_status ON toxicity_import_job (status);

-- =========================
-- 9. Add toxicity_score column to molecule table
-- =========================
ALTER TABLE molecule
    ADD COLUMN IF NOT EXISTS toxicity_score DECIMAL(5,2);

CREATE INDEX IF NOT EXISTS idx_molecule_toxicity_score ON molecule (toxicity_score);

-- =========================
-- 10. Add toxicity_data_available flag to molecule table
-- =========================
ALTER TABLE molecule
    ADD COLUMN IF NOT EXISTS toxicity_data_available BOOLEAN DEFAULT FALSE;

CREATE INDEX IF NOT EXISTS idx_molecule_toxicity_data_available ON molecule (toxicity_data_available);

-- =========================
-- 11. Add calculation method for toxicity scoring
-- =========================
INSERT INTO calculation_method (name, description, method_type, reference)
VALUES (
    'Tox21 Score',
    'Toxicity score based on Tox21 assay data',
    'computational',
    'Tox21 program'
)
ON CONFLICT (name) DO NOTHING;

-- =========================
-- 12. RLS Policies
-- =========================
-- Enable RLS on all toxicity tables
ALTER TABLE toxicity_data_source ENABLE ROW LEVEL SECURITY;
ALTER TABLE toxicity_assay ENABLE ROW LEVEL SECURITY;
ALTER TABLE molecule_identifier_mapping ENABLE ROW LEVEL SECURITY;
ALTER TABLE toxicity_data ENABLE ROW LEVEL SECURITY;
ALTER TABLE toxicity_score ENABLE ROW LEVEL SECURITY;
ALTER TABLE toxicity_endpoint ENABLE ROW LEVEL SECURITY;
ALTER TABLE assay_endpoint_mapping ENABLE ROW LEVEL SECURITY;
ALTER TABLE toxicity_import_job ENABLE ROW LEVEL SECURITY;

-- Create policies for toxicity_data_source (reference data, viewable by all authenticated users)
CREATE POLICY toxicity_data_source_select_policy ON toxicity_data_source
    FOR SELECT USING (auth.role() = 'authenticated');

-- Create policies for toxicity_assay (reference data, viewable by all authenticated users)
CREATE POLICY toxicity_assay_select_policy ON toxicity_assay
    FOR SELECT USING (auth.role() = 'authenticated');

-- Create policies for molecule_identifier_mapping (follows molecule access)
CREATE POLICY molecule_identifier_mapping_select_policy ON molecule_identifier_mapping
    FOR SELECT USING (
        auth.role() = 'authenticated' AND
        EXISTS (
            SELECT 1 FROM molecule m
            WHERE m.id = molecule_identifier_mapping.molecule_id
            AND (
                m.created_by = auth.uid() OR
                m.team_id IN (SELECT team_id FROM team_member WHERE user_id = auth.uid()) OR
                m.project_id IN (SELECT project_id FROM project_member WHERE user_id = auth.uid())
            )
        )
    );

-- Create policies for toxicity_data (follows molecule access)
CREATE POLICY toxicity_data_select_policy ON toxicity_data
    FOR SELECT USING (
        auth.role() = 'authenticated' AND
        EXISTS (
            SELECT 1 FROM molecule m
            WHERE m.id = toxicity_data.molecule_id
            AND (
                m.created_by = auth.uid() OR
                m.team_id IN (SELECT team_id FROM team_member WHERE user_id = auth.uid()) OR
                m.project_id IN (SELECT project_id FROM project_member WHERE user_id = auth.uid())
            )
        )
    );

-- Create policies for toxicity_score (follows molecule access)
CREATE POLICY toxicity_score_select_policy ON toxicity_score
    FOR SELECT USING (
        auth.role() = 'authenticated' AND
        EXISTS (
            SELECT 1 FROM molecule m
            WHERE m.id = toxicity_score.molecule_id
            AND (
                m.created_by = auth.uid() OR
                m.team_id IN (SELECT team_id FROM team_member WHERE user_id = auth.uid()) OR
                m.project_id IN (SELECT project_id FROM project_member WHERE user_id = auth.uid())
            )
        )
    );

-- Create policies for toxicity_endpoint (reference data, viewable by all authenticated users)
CREATE POLICY toxicity_endpoint_select_policy ON toxicity_endpoint
    FOR SELECT USING (auth.role() = 'authenticated');

-- Create policies for assay_endpoint_mapping (reference data, viewable by all authenticated users)
CREATE POLICY assay_endpoint_mapping_select_policy ON assay_endpoint_mapping
    FOR SELECT USING (auth.role() = 'authenticated');

-- Create policies for toxicity_import_job (admin only)
CREATE POLICY toxicity_import_job_select_policy ON toxicity_import_job
    FOR SELECT USING (auth.role() = 'authenticated' AND auth.uid() IN (SELECT auth_user_id FROM user_profile WHERE is_admin = TRUE));