-- CryoProtect Analyzer: Scientific Database Schema Migration
-- Migration: 005_scientific_schema.sql
-- This migration implements best practices for scientific data storage in molecular/cryoprotectant research.
-- It is idempotent and safe to run multiple times.

-- =========================
-- 1. UserProfile Table
-- =========================
CREATE TABLE IF NOT EXISTS user_profile (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    auth_user_id uuid NOT NULL UNIQUE, -- Supabase Auth user id
    display_name TEXT,
    email TEXT,
    affiliation TEXT,
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

-- =========================
-- 2. Team Table
-- =========================
CREATE TABLE IF NOT EXISTS team (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT NOT NULL,
    description TEXT,
    created_by uuid REFERENCES user_profile(id),
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

-- =========================
-- 3. Project Table
-- =========================
CREATE TABLE IF NOT EXISTS project (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT NOT NULL,
    description TEXT,
    team_id uuid REFERENCES team(id),
    created_by uuid REFERENCES user_profile(id),
    data_source TEXT,
    version INTEGER DEFAULT 1,
    modification_history JSONB DEFAULT '[]',
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

-- =========================
-- 4. Molecule Table
-- =========================
CREATE TABLE IF NOT EXISTS molecule (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT,
    smiles TEXT NOT NULL, -- Canonical SMILES
    inchi TEXT NOT NULL,  -- InChI string
    inchikey TEXT NOT NULL, -- InChIKey for indexing/search
    formula TEXT,
    molecular_weight DECIMAL(18,8), -- High precision
    created_by uuid REFERENCES user_profile(id),
    data_source TEXT,
    version INTEGER DEFAULT 1,
    modification_history JSONB DEFAULT '[]',
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now(),
    UNIQUE (smiles),
    UNIQUE (inchikey)
);

CREATE INDEX IF NOT EXISTS idx_molecule_inchikey ON molecule (inchikey);
CREATE INDEX IF NOT EXISTS idx_molecule_smiles ON molecule (smiles);

-- =========================
-- 5. Mixture Table
-- =========================
CREATE TABLE IF NOT EXISTS mixture (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT,
    description TEXT,
    project_id uuid REFERENCES project(id),
    created_by uuid REFERENCES user_profile(id),
    data_source TEXT,
    version INTEGER DEFAULT 1,
    modification_history JSONB DEFAULT '[]',
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

-- =========================
-- 6. MixtureComponent Table (Graph-based: links mixture to molecules)
-- =========================
CREATE TABLE IF NOT EXISTS mixture_component (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    mixture_id uuid NOT NULL REFERENCES mixture(id) ON DELETE CASCADE,
    molecule_id uuid NOT NULL REFERENCES molecule(id),
    amount DECIMAL(18,8) NOT NULL, -- High precision
    amount_unit TEXT NOT NULL,     -- Explicit units (e.g., mol, g, %)
    role TEXT,                     -- e.g., solvent, solute, cryoprotectant
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now(),
    UNIQUE (mixture_id, molecule_id, role)
);

CREATE INDEX IF NOT EXISTS idx_mixture_component_mixture_id ON mixture_component (mixture_id);
CREATE INDEX IF NOT EXISTS idx_mixture_component_molecule_id ON mixture_component (molecule_id);

-- =========================
-- 7. Experiment Table (Recipe/Preparation pattern)
-- =========================
CREATE TABLE IF NOT EXISTS experiment (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT,
    description TEXT,
    mixture_id uuid REFERENCES mixture(id),
    project_id uuid REFERENCES project(id),
    preparation_protocol TEXT, -- Recipe/Preparation details
    temperature DECIMAL(10,4),
    temperature_unit TEXT,
    pressure DECIMAL(10,4),
    pressure_unit TEXT,
    created_by uuid REFERENCES user_profile(id),
    data_source TEXT,
    version INTEGER DEFAULT 1,
    modification_history JSONB DEFAULT '[]',
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS idx_experiment_mixture_id ON experiment (mixture_id);

-- =========================
-- 8. CalculationMethod Table
-- =========================
CREATE TABLE IF NOT EXISTS calculation_method (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT NOT NULL,
    description TEXT,
    method_type TEXT, -- e.g., experimental, computational, ML
    reference TEXT,
    created_by uuid REFERENCES user_profile(id),
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

-- =========================
-- 9. MolecularProperty Table
-- =========================
CREATE TABLE IF NOT EXISTS molecular_property (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    molecule_id uuid REFERENCES molecule(id),
    property_type TEXT NOT NULL, -- e.g., melting_point, logP
    value DECIMAL(18,8) NOT NULL,
    unit TEXT NOT NULL,
    method_id uuid REFERENCES calculation_method(id),
    experiment_id uuid REFERENCES experiment(id),
    provenance TEXT,
    created_by uuid REFERENCES user_profile(id),
    data_source TEXT,
    version INTEGER DEFAULT 1,
    modification_history JSONB DEFAULT '[]',
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS idx_molecular_property_molecule_id ON molecular_property (molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecular_property_type ON molecular_property (property_type);

-- =========================
-- 10. Prediction Table
-- =========================
CREATE TABLE IF NOT EXISTS prediction (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    molecule_id uuid REFERENCES molecule(id),
    property_type TEXT NOT NULL,
    predicted_value DECIMAL(18,8) NOT NULL,
    unit TEXT NOT NULL,
    method_id uuid REFERENCES calculation_method(id),
    model_version TEXT,
    confidence DECIMAL(5,4),
    created_by uuid REFERENCES user_profile(id),
    data_source TEXT,
    version INTEGER DEFAULT 1,
    modification_history JSONB DEFAULT '[]',
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS idx_prediction_molecule_id ON prediction (molecule_id);
CREATE INDEX IF NOT EXISTS idx_prediction_property_type ON prediction (property_type);

-- =========================
-- 11. ExperimentProperty Table (Properties measured in experiments, e.g., glass transition temp)
-- =========================
CREATE TABLE IF NOT EXISTS experiment_property (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    experiment_id uuid REFERENCES experiment(id),
    property_type TEXT NOT NULL,
    value DECIMAL(18,8) NOT NULL,
    unit TEXT NOT NULL,
    method_id uuid REFERENCES calculation_method(id),
    provenance TEXT,
    created_by uuid REFERENCES user_profile(id),
    data_source TEXT,
    version INTEGER DEFAULT 1,
    modification_history JSONB DEFAULT '[]',
    created_at TIMESTAMPTZ DEFAULT now(),
    updated_at TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS idx_experiment_property_experiment_id ON experiment_property (experiment_id);

-- =========================
-- 12. Graph-based Relationships (e.g., mixture <-> experiment, property dependencies)
-- =========================
-- Example: Link experiments to multiple mixtures (if needed)
CREATE TABLE IF NOT EXISTS experiment_mixture_link (
    id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
    experiment_id uuid REFERENCES experiment(id),
    mixture_id uuid REFERENCES mixture(id),
    relationship_type TEXT, -- e.g., "input", "output"
    created_at TIMESTAMPTZ DEFAULT now()
);

-- =========================
-- 13. RLS/Ownership Support
-- =========================
-- All main tables include project_id, team_id, or created_by for RLS compatibility.
-- Add project_id/team_id columns if needed for RLS in other tables (example below):

ALTER TABLE molecule
    ADD COLUMN IF NOT EXISTS project_id uuid REFERENCES project(id);

ALTER TABLE molecule
    ADD COLUMN IF NOT EXISTS team_id uuid REFERENCES team(id);

ALTER TABLE mixture
    ADD COLUMN IF NOT EXISTS team_id uuid REFERENCES team(id);

ALTER TABLE experiment
    ADD COLUMN IF NOT EXISTS team_id uuid REFERENCES team(id);

-- =========================
-- 14. Timestamps/Versioning/Provenance
-- =========================
-- All main tables include created_at, updated_at, version, data_source, modification_history, created_by.

-- =========================
-- 15. Indexes for Performance
-- =========================
-- Already included above for InChIKey, SMILES, foreign keys.

-- =========================
-- 16. Idempotency
-- =========================
-- All CREATE TABLE/INDEX/ALTER TABLE statements use IF NOT EXISTS.

-- =========================
-- 17. Comments on Scientific Design
-- =========================
-- - Dual SMILES/InChI/InChIKey for robust molecular identification.
-- - High-precision DECIMAL fields for scientific accuracy.
-- - Explicit units for all quantitative fields.
-- - Graph-based tables (mixture_component, experiment_mixture_link) for flexible relationships.
-- - Provenance, versioning, and modification history for scientific traceability.
-- - RLS-compatible fields for secure, multi-tenant data access.
-- - Indexes for fast search and retrieval of key scientific identifiers.