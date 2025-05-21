-- Migration: 033_add_performance_indexes.sql
-- Purpose: Add performance-optimizing indexes to frequently queried columns
-- This migration adds indexes to improve query performance

-- Molecules table indexes
-- For molecule lookup by ID (already has primary key index)
-- For molecule lookup by name (case-insensitive search)
CREATE INDEX idx_molecules_name_lower ON molecules USING btree (lower(name));

-- For pattern matching on molecule name
CREATE INDEX idx_molecules_name_trgm ON molecules USING gin (name gin_trgm_ops);

-- For molecule lookup by formula
CREATE INDEX idx_molecules_formula ON molecules (formula);

-- For molecule lookup by type
CREATE INDEX idx_molecules_type ON molecules (type);

-- For searching by pubchem_cid
CREATE INDEX idx_molecules_pubchem_cid ON molecules (pubchem_cid);

-- For filtering by public/private status
CREATE INDEX idx_molecules_is_public ON molecules (is_public);

-- For ordering by creation date
CREATE INDEX idx_molecules_created_at ON molecules (created_at);

-- Molecular properties table indexes
-- For property lookup by molecule
CREATE INDEX idx_molecular_properties_molecule_id ON molecular_properties (molecule_id);

-- For property lookup by name
CREATE INDEX idx_molecular_properties_property_name ON molecular_properties (property_name);

-- For property lookup by numeric value ranges
CREATE INDEX idx_molecular_properties_numeric_value ON molecular_properties (numeric_value);

-- For combined molecule+property queries
CREATE INDEX idx_molecular_properties_molecule_property ON molecular_properties (molecule_id, property_name);

-- For filtering properties by type
CREATE INDEX idx_molecular_properties_property_type_id ON molecular_properties (property_type_id);

-- For unit conversion queries
CREATE INDEX idx_molecular_properties_unit_id ON molecular_properties (unit_id);

-- Consolidated molecules table indexes
-- For lookup by primary ID
CREATE INDEX idx_consolidated_molecules_primary_id ON consolidated_molecules (primary_id);

-- For finding duplicates
CREATE INDEX idx_consolidated_molecules_duplicate_of ON consolidated_molecules (duplicate_of);

-- Property types table indexes
-- For lookups by property name
CREATE INDEX idx_property_types_name ON property_types (name);

-- For lookups by value type
CREATE INDEX idx_property_types_value_type ON property_types (value_type);

-- Units table indexes
-- For lookups by unit name
CREATE INDEX idx_units_name ON units (name);

-- For lookups by category
CREATE INDEX idx_units_category ON units (category);

-- For lookups by symbol
CREATE INDEX idx_units_symbol ON units (symbol);

-- Mixtures table indexes (if exists)
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables 
               WHERE table_schema = 'public' 
               AND table_name = 'mixtures') THEN
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_mixtures_name ON mixtures (name)';
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_mixtures_created_at ON mixtures (created_at)';
    END IF;
END $$;

-- Mixture components table indexes (if exists)
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables 
               WHERE table_schema = 'public' 
               AND table_name = 'mixture_components') THEN
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_mixture_components_mixture_id ON mixture_components (mixture_id)';
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_mixture_components_molecule_id ON mixture_components (molecule_id)';
    END IF;
END $$;

-- Experiments table indexes (if exists)
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables 
               WHERE table_schema = 'public' 
               AND table_name = 'experiments') THEN
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_experiments_name ON experiments (name)';
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_experiments_created_by ON experiments (created_by)';
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_experiments_created_at ON experiments (created_at)';
    END IF;
END $$;

-- User related indexes (if user-related tables exist)
DO $$
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables 
               WHERE table_schema = 'public' 
               AND table_name = 'users') THEN
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_users_email ON users (email)';
    END IF;
    
    IF EXISTS (SELECT FROM information_schema.tables 
               WHERE table_schema = 'public' 
               AND table_name = 'teams') THEN
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_teams_name ON teams (name)';
    END IF;
    
    IF EXISTS (SELECT FROM information_schema.tables 
               WHERE table_schema = 'public' 
               AND table_name = 'team_members') THEN
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_team_members_team_id ON team_members (team_id)';
        EXECUTE 'CREATE INDEX IF NOT EXISTS idx_team_members_user_id ON team_members (user_id)';
    END IF;
END $$;