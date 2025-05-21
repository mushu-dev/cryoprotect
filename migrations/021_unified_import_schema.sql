-- CryoProtect Analyzer: Unified Import Schema Migration
-- Migration: 021_unified_import_schema.sql
-- Author: Claude
-- Date: 2025-05-11
--
-- This migration implements schema changes for the unified chemical data importer
-- which can handle multiple data sources including PubChem, ChEMBL, etc.
-- It is idempotent and safe to run multiple times.

-- Begin transaction to ensure all changes are applied atomically
BEGIN;

-- =========================
-- 1. Create data source tracking table
-- =========================
CREATE TABLE IF NOT EXISTS public.data_sources (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT UNIQUE NOT NULL,
    description TEXT,
    url TEXT,
    api_base_url TEXT,
    version TEXT,
    status TEXT DEFAULT 'active', -- 'active', 'deprecated', 'offline'
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- Insert default data sources
INSERT INTO public.data_sources (name, description, url, api_base_url, version)
VALUES 
    ('PubChem', 'National Center for Biotechnology Information PubChem database', 'https://pubchem.ncbi.nlm.nih.gov', 'https://pubchem.ncbi.nlm.nih.gov/rest/pug', 'current'),
    ('ChEMBL', 'European Molecular Biology Laboratory ChEMBL database', 'https://www.ebi.ac.uk/chembl', 'https://www.ebi.ac.uk/chembl/api', 'current'),
    ('Unified', 'CryoProtect Unified Import System', 'internal', 'internal', '1.0')
ON CONFLICT (name) DO NOTHING;

-- =========================
-- 2. Create molecule synonym/identifier table for better tracking
-- =========================
CREATE TABLE IF NOT EXISTS public.molecule_synonyms (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
    data_source_id UUID REFERENCES public.data_sources(id),
    synonym_type TEXT NOT NULL, -- 'name', 'iupac_name', 'common_name', 'trade_name', etc.
    synonym_value TEXT NOT NULL,
    is_preferred BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    UNIQUE(molecule_id, data_source_id, synonym_type, synonym_value)
);

-- Create indexes for synonym searches
CREATE INDEX IF NOT EXISTS idx_molecule_synonyms_molecule_id ON public.molecule_synonyms(molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecule_synonyms_value ON public.molecule_synonyms(synonym_value);
CREATE INDEX IF NOT EXISTS idx_molecule_synonyms_type_value ON public.molecule_synonyms(synonym_type, synonym_value);

-- Apply RLS to molecule_synonyms table
ALTER TABLE public.molecule_synonyms ENABLE ROW LEVEL SECURITY;

-- Create policies for molecule_synonyms (viewable by all, managed by authenticated users)
DO $$
BEGIN
    IF NOT EXISTS (SELECT FROM pg_policies WHERE tablename = 'molecule_synonyms' AND policyname = 'Molecule synonyms are viewable by everyone') THEN
        CREATE POLICY "Molecule synonyms are viewable by everyone" 
        ON public.molecule_synonyms FOR SELECT USING (true);
    END IF;
    
    IF NOT EXISTS (SELECT FROM pg_policies WHERE tablename = 'molecule_synonyms' AND policyname = 'Molecule synonyms can be inserted by authenticated users') THEN
        CREATE POLICY "Molecule synonyms can be inserted by authenticated users" 
        ON public.molecule_synonyms FOR INSERT WITH CHECK (auth.role() = 'authenticated');
    END IF;
    
    IF NOT EXISTS (SELECT FROM pg_policies WHERE tablename = 'molecule_synonyms' AND policyname = 'Molecule synonyms can be updated by authenticated users') THEN
        CREATE POLICY "Molecule synonyms can be updated by authenticated users" 
        ON public.molecule_synonyms FOR UPDATE USING (auth.role() = 'authenticated');
    END IF;
END $$;

-- Create trigger for updated_at timestamp on molecule_synonyms
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_molecule_synonyms') THEN
        CREATE TRIGGER set_timestamp_molecule_synonyms
        BEFORE UPDATE ON public.molecule_synonyms
        FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    END IF;
END $$;

-- =========================
-- 3. Add and modify columns on molecules table to support unified import
-- =========================

-- Add default source column if it doesn't exist
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                   WHERE table_schema = 'public' 
                   AND table_name = 'molecules' 
                   AND column_name = 'data_source_id') THEN
        ALTER TABLE public.molecules 
        ADD COLUMN data_source_id UUID REFERENCES public.data_sources(id);
        
        -- Create index on data_source_id
        CREATE INDEX idx_molecules_data_source_id ON public.molecules(data_source_id);
    END IF;
END $$;

-- Add unified import fields
DO $$
BEGIN
    -- Add import_source column (can store JSON with source info)
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                   WHERE table_schema = 'public' 
                   AND table_name = 'molecules' 
                   AND column_name = 'import_source') THEN
        ALTER TABLE public.molecules 
        ADD COLUMN import_source JSONB DEFAULT '{}';
    END IF;
    
    -- Add import_timestamp for tracking when data was imported
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                   WHERE table_schema = 'public' 
                   AND table_name = 'molecules' 
                   AND column_name = 'import_timestamp') THEN
        ALTER TABLE public.molecules 
        ADD COLUMN import_timestamp TIMESTAMPTZ DEFAULT NOW();
    END IF;
    
    -- Add last_updated_source to track which source last updated this molecule
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                   WHERE table_schema = 'public' 
                   AND table_name = 'molecules' 
                   AND column_name = 'last_updated_source') THEN
        ALTER TABLE public.molecules 
        ADD COLUMN last_updated_source TEXT;
    END IF;
    
    -- Add InChI and InChIKey if they don't exist (for new unified schema)
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                  WHERE table_schema = 'public' 
                  AND table_name = 'molecules' 
                  AND column_name = 'inchi') THEN
        ALTER TABLE public.molecules 
        ADD COLUMN inchi TEXT;
    END IF;
    
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                  WHERE table_schema = 'public' 
                  AND table_name = 'molecules' 
                  AND column_name = 'inchikey') THEN
        ALTER TABLE public.molecules 
        ADD COLUMN inchikey TEXT;
        
        -- Create unique index on inchikey (allowing nulls for backward compatibility)
        CREATE UNIQUE INDEX IF NOT EXISTS idx_molecules_inchikey_unique 
        ON public.molecules(inchikey) 
        WHERE inchikey IS NOT NULL;
    END IF;
END $$;

-- =========================
-- 4. Add source tracking to molecular_properties table
-- =========================
DO $$
BEGIN
    -- Add source column if it doesn't exist
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                   WHERE table_schema = 'public' 
                   AND table_name = 'molecular_properties' 
                   AND column_name = 'source') THEN
        ALTER TABLE public.molecular_properties 
        ADD COLUMN source TEXT;
    END IF;
    
    -- Add data_source_id to link to the data_sources table
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                   WHERE table_schema = 'public' 
                   AND table_name = 'molecular_properties' 
                   AND column_name = 'data_source_id') THEN
        ALTER TABLE public.molecular_properties 
        ADD COLUMN data_source_id UUID REFERENCES public.data_sources(id);
    END IF;
    
    -- Add confidence score for property value
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                   WHERE table_schema = 'public' 
                   AND table_name = 'molecular_properties' 
                   AND column_name = 'confidence') THEN
        ALTER TABLE public.molecular_properties 
        ADD COLUMN confidence NUMERIC CHECK (confidence >= 0 AND confidence <= 1);
    END IF;
END $$;

-- =========================
-- 5. Create import job history table
-- =========================
CREATE TABLE IF NOT EXISTS public.import_jobs (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    data_source_id UUID REFERENCES public.data_sources(id),
    job_type TEXT NOT NULL, -- 'initial', 'update', 'batch', etc.
    status TEXT NOT NULL, -- 'queued', 'running', 'completed', 'failed'
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ,
    records_processed INTEGER DEFAULT 0,
    records_inserted INTEGER DEFAULT 0,
    records_updated INTEGER DEFAULT 0,
    records_skipped INTEGER DEFAULT 0,
    records_failed INTEGER DEFAULT 0,
    error_log JSONB DEFAULT '[]',
    options JSONB DEFAULT '{}',
    created_by UUID REFERENCES auth.users(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- Index on status and data_source_id for filtering jobs
CREATE INDEX IF NOT EXISTS idx_import_jobs_status ON public.import_jobs(status);
CREATE INDEX IF NOT EXISTS idx_import_jobs_data_source_id ON public.import_jobs(data_source_id);

-- Apply RLS to import_jobs table
ALTER TABLE public.import_jobs ENABLE ROW LEVEL SECURITY;

-- Create policies for import_jobs
DO $$
BEGIN
    IF NOT EXISTS (SELECT FROM pg_policies WHERE tablename = 'import_jobs' AND policyname = 'Import jobs are viewable by authenticated users') THEN
        CREATE POLICY "Import jobs are viewable by authenticated users" 
        ON public.import_jobs FOR SELECT USING (auth.role() = 'authenticated');
    END IF;
    
    IF NOT EXISTS (SELECT FROM pg_policies WHERE tablename = 'import_jobs' AND policyname = 'Import jobs can be managed by authenticated users') THEN
        CREATE POLICY "Import jobs can be managed by authenticated users" 
        ON public.import_jobs FOR ALL USING (auth.role() = 'authenticated');
    END IF;
END $$;

-- Create trigger for updated_at timestamp on import_jobs
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_import_jobs') THEN
        CREATE TRIGGER set_timestamp_import_jobs
        BEFORE UPDATE ON public.import_jobs
        FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    END IF;
END $$;

-- =========================
-- 6. Create molecule cross-reference table
-- =========================
CREATE TABLE IF NOT EXISTS public.molecule_cross_references (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
    source_data_source_id UUID NOT NULL REFERENCES public.data_sources(id),
    target_data_source_id UUID NOT NULL REFERENCES public.data_sources(id),
    source_identifier TEXT NOT NULL,
    target_identifier TEXT NOT NULL,
    confidence NUMERIC CHECK (confidence >= 0 AND confidence <= 1),
    mapping_method TEXT, -- 'exact_match', 'inchikey_match', 'name_match', etc.
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    UNIQUE(source_data_source_id, source_identifier, target_data_source_id, target_identifier)
);

-- Create indexes for cross-reference lookups
CREATE INDEX IF NOT EXISTS idx_molecule_cross_references_molecule_id 
ON public.molecule_cross_references(molecule_id);

CREATE INDEX IF NOT EXISTS idx_molecule_cross_references_source 
ON public.molecule_cross_references(source_data_source_id, source_identifier);

CREATE INDEX IF NOT EXISTS idx_molecule_cross_references_target 
ON public.molecule_cross_references(target_data_source_id, target_identifier);

-- Apply RLS to molecule_cross_references table
ALTER TABLE public.molecule_cross_references ENABLE ROW LEVEL SECURITY;

-- Create policies for molecule_cross_references
DO $$
BEGIN
    IF NOT EXISTS (SELECT FROM pg_policies WHERE tablename = 'molecule_cross_references' AND policyname = 'Molecule cross references are viewable by authenticated users') THEN
        CREATE POLICY "Molecule cross references are viewable by authenticated users" 
        ON public.molecule_cross_references FOR SELECT USING (auth.role() = 'authenticated');
    END IF;
    
    IF NOT EXISTS (SELECT FROM pg_policies WHERE tablename = 'molecule_cross_references' AND policyname = 'Molecule cross references can be managed by authenticated users') THEN
        CREATE POLICY "Molecule cross references can be managed by authenticated users" 
        ON public.molecule_cross_references FOR ALL USING (auth.role() = 'authenticated');
    END IF;
END $$;

-- Create trigger for updated_at timestamp on molecule_cross_references
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_molecule_cross_references') THEN
        CREATE TRIGGER set_timestamp_molecule_cross_references
        BEFORE UPDATE ON public.molecule_cross_references
        FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    END IF;
END $$;

-- =========================
-- 7. Create functions for unified import system
-- =========================

-- Function to find a molecule by any identifier across datasources
CREATE OR REPLACE FUNCTION public.find_molecule_by_identifier(
    identifier_value TEXT,
    identifier_type TEXT DEFAULT NULL
) RETURNS TABLE (
    molecule_id UUID,
    match_type TEXT,
    confidence NUMERIC
) AS $$
BEGIN
    -- First try direct molecule table lookups for common identifiers
    RETURN QUERY 
    SELECT 
        m.id AS molecule_id,
        'direct_match' AS match_type,
        1.0 AS confidence
    FROM 
        public.molecules m
    WHERE 
        m.cid::TEXT = identifier_value OR
        m.chembl_id = identifier_value OR
        m.inchikey = identifier_value
    LIMIT 1;
    
    -- If no results, try molecule_synonyms table
    IF NOT FOUND THEN
        RETURN QUERY
        SELECT 
            ms.molecule_id,
            'synonym_match' AS match_type,
            CASE WHEN ms.is_preferred THEN 0.9 ELSE 0.8 END AS confidence
        FROM 
            public.molecule_synonyms ms
        WHERE 
            ms.synonym_value = identifier_value
            AND (identifier_type IS NULL OR ms.synonym_type = identifier_type)
        ORDER BY 
            ms.is_preferred DESC
        LIMIT 1;
    END IF;
    
    -- If still no results, try cross-references
    IF NOT FOUND THEN
        RETURN QUERY
        SELECT 
            mcr.molecule_id,
            'cross_reference_match' AS match_type,
            mcr.confidence
        FROM 
            public.molecule_cross_references mcr
        WHERE 
            mcr.source_identifier = identifier_value OR
            mcr.target_identifier = identifier_value
        ORDER BY 
            mcr.confidence DESC
        LIMIT 1;
    END IF;
    
    -- Still no results, try fuzzy name matching with molecule_identifier_mapping as fallback
    IF NOT FOUND AND identifier_type IS NULL THEN
        RETURN QUERY
        SELECT 
            mim.molecule_id,
            'mapping_match' AS match_type,
            mim.confidence_score AS confidence
        FROM 
            public.molecule_identifier_mapping mim
        WHERE 
            mim.external_id = identifier_value
        ORDER BY 
            mim.confidence_score DESC
        LIMIT 1;
    END IF;
    
    RETURN;
END;
$$ LANGUAGE plpgsql;

-- Function to merge molecule data from multiple sources
CREATE OR REPLACE FUNCTION public.merge_molecule_data(
    primary_molecule_id UUID, 
    secondary_molecule_id UUID, 
    merge_properties BOOLEAN DEFAULT TRUE,
    merge_synonyms BOOLEAN DEFAULT TRUE,
    remove_secondary BOOLEAN DEFAULT FALSE
) RETURNS BOOLEAN AS $$
DECLARE
    success BOOLEAN := TRUE;
BEGIN
    -- Merge properties if requested
    IF merge_properties THEN
        -- Copy properties from secondary to primary if they don't exist
        INSERT INTO public.molecular_properties (
            molecule_id, property_type_id, numeric_value, text_value, 
            boolean_value, data_source_id, source, confidence
        )
        SELECT 
            primary_molecule_id, mp.property_type_id, mp.numeric_value, mp.text_value,
            mp.boolean_value, mp.data_source_id, mp.source, mp.confidence
        FROM 
            public.molecular_properties mp
        WHERE 
            mp.molecule_id = secondary_molecule_id
            -- Only import if primary doesn't have this property type yet
            AND NOT EXISTS (
                SELECT 1 FROM public.molecular_properties 
                WHERE molecule_id = primary_molecule_id 
                AND property_type_id = mp.property_type_id
            );
    END IF;
    
    -- Merge synonyms if requested
    IF merge_synonyms THEN
        -- Copy synonyms from secondary to primary if they don't exist
        INSERT INTO public.molecule_synonyms (
            molecule_id, data_source_id, synonym_type, synonym_value, is_preferred
        )
        SELECT 
            primary_molecule_id, ms.data_source_id, ms.synonym_type, 
            ms.synonym_value, ms.is_preferred
        FROM 
            public.molecule_synonyms ms
        WHERE 
            ms.molecule_id = secondary_molecule_id
            -- Only import if primary doesn't have this synonym yet
            AND NOT EXISTS (
                SELECT 1 FROM public.molecule_synonyms 
                WHERE molecule_id = primary_molecule_id 
                AND synonym_type = ms.synonym_type
                AND synonym_value = ms.synonym_value
            );
    END IF;
    
    -- Add cross-reference between the two molecule IDs
    INSERT INTO public.molecule_cross_references (
        molecule_id, source_data_source_id, target_data_source_id,
        source_identifier, target_identifier, confidence, mapping_method
    )
    SELECT 
        primary_molecule_id, 
        m1.data_source_id, 
        m2.data_source_id,
        COALESCE(m1.cid::TEXT, m1.chembl_id, m1.id::TEXT), 
        COALESCE(m2.cid::TEXT, m2.chembl_id, m2.id::TEXT),
        1.0, 
        'direct_merge'
    FROM 
        public.molecules m1, 
        public.molecules m2
    WHERE 
        m1.id = primary_molecule_id AND m2.id = secondary_molecule_id
    ON CONFLICT (source_data_source_id, source_identifier, target_data_source_id, target_identifier) 
    DO NOTHING;
    
    -- If requested, remove the secondary molecule
    IF remove_secondary THEN
        BEGIN
            DELETE FROM public.molecules WHERE id = secondary_molecule_id;
        EXCEPTION WHEN OTHERS THEN
            success := FALSE;
        END;
    END IF;
    
    RETURN success;
END;
$$ LANGUAGE plpgsql;

-- Function to update InChI representation data when SMILES is updated
CREATE OR REPLACE FUNCTION public.update_inchi_from_smiles()
RETURNS TRIGGER AS $$
BEGIN
    -- This is just a placeholder for a function that would be implemented in application code.
    -- Since we can't actually generate InChI from SMILES directly in SQL,
    -- we'll just ensure the trigger is set up so the application can handle it.
    -- In practice, this would be done by having RDKit or other tools generate InChI from SMILES.
    
    -- Dummy implementation that doesn't actually change anything:
    NEW.inchi := NEW.inchi;
    NEW.inchikey := NEW.inchikey;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Create or replace the trigger for InChI updates if needed
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_trigger WHERE tgname = 'update_inchi_data') THEN
        CREATE TRIGGER update_inchi_data
        BEFORE UPDATE OF smiles ON public.molecules
        FOR EACH ROW
        WHEN (OLD.smiles IS DISTINCT FROM NEW.smiles)
        EXECUTE PROCEDURE public.update_inchi_from_smiles();
    END IF;
END $$;

-- =========================
-- 8. Create the needed property types if they don't exist
-- =========================

-- Insert common property types that will be used by the unified importer
INSERT INTO public.property_types (name, description, units, data_type)
VALUES
    ('molecular_weight', 'Molecular Weight', 'g/mol', 'numeric'),
    ('logp', 'Octanol-Water Partition Coefficient (LogP)', NULL, 'numeric'),
    ('tpsa', 'Topological Polar Surface Area', 'Å²', 'numeric'),
    ('h_bond_donors', 'Number of Hydrogen Bond Donors', NULL, 'numeric'),
    ('h_bond_acceptors', 'Number of Hydrogen Bond Acceptors', NULL, 'numeric'),
    ('ro5_violations', 'Number of Lipinski Rule of 5 Violations', NULL, 'numeric'),
    ('med_chem_friendly', 'Medchem Friendliness Flag', NULL, 'text'),
    ('rotatable_bonds', 'Number of Rotatable Bonds', NULL, 'numeric'),
    ('heavy_atoms', 'Number of Heavy (Non-hydrogen) Atoms', NULL, 'numeric'),
    ('iupac_name', 'IUPAC Name', NULL, 'text'),
    ('canonical_smiles', 'Canonical SMILES Representation', NULL, 'text'),
    ('inchi_string', 'International Chemical Identifier (InChI)', NULL, 'text'),
    ('inchikey', 'InChIKey Identifier', NULL, 'text')
ON CONFLICT (name) DO NOTHING;

-- =========================
-- 9. Add RLS Policies to new data_sources table
-- =========================
ALTER TABLE public.data_sources ENABLE ROW LEVEL SECURITY;

DO $$
BEGIN
    -- Create policies for data_sources
    IF NOT EXISTS (SELECT FROM pg_policies WHERE tablename = 'data_sources' AND policyname = 'Data sources are viewable by everyone') THEN
        CREATE POLICY "Data sources are viewable by everyone" 
        ON public.data_sources FOR SELECT USING (true);
    END IF;
    
    IF NOT EXISTS (SELECT FROM pg_policies WHERE tablename = 'data_sources' AND policyname = 'Data sources can be managed by authenticated users') THEN
        CREATE POLICY "Data sources can be managed by authenticated users" 
        ON public.data_sources FOR ALL USING (auth.role() = 'authenticated');
    END IF;
END $$;

-- Create trigger for updated_at timestamp on data_sources
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_data_sources') THEN
        CREATE TRIGGER set_timestamp_data_sources
        BEFORE UPDATE ON public.data_sources
        FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    END IF;
END $$;

-- =========================
-- 10. Link existing molecules to their data sources
-- =========================
DO $$
DECLARE
    pubchem_id UUID;
    chembl_id UUID;
BEGIN
    -- Get data source IDs
    SELECT id INTO pubchem_id FROM public.data_sources WHERE name = 'PubChem';
    SELECT id INTO chembl_id FROM public.data_sources WHERE name = 'ChEMBL';
    
    -- Update molecules with PubChem as data source if they have a CID
    IF pubchem_id IS NOT NULL THEN
        UPDATE public.molecules
        SET data_source_id = pubchem_id
        WHERE cid IS NOT NULL 
        AND (data_source_id IS NULL OR data_source = 'PubChem')
        AND data_source_id IS DISTINCT FROM pubchem_id;
    END IF;
    
    -- Update molecules with ChEMBL as data source if they have a chembl_id
    IF chembl_id IS NOT NULL THEN
        UPDATE public.molecules
        SET data_source_id = chembl_id
        WHERE chembl_id IS NOT NULL 
        AND (data_source_id IS NULL OR data_source = 'ChEMBL')
        AND data_source_id IS DISTINCT FROM chembl_id;
    END IF;
END $$;

-- =========================
-- 11. Create indexes to improve import performance
-- =========================

-- Create indexes on molecule identifiers if they don't exist
CREATE INDEX IF NOT EXISTS idx_molecules_smiles ON public.molecules(smiles);
CREATE INDEX IF NOT EXISTS idx_molecules_chemical_formula ON public.molecules(molecular_formula);
CREATE INDEX IF NOT EXISTS idx_molecules_name_trgm ON public.molecules USING gin (name gin_trgm_ops);

-- Enable pg_trgm extension if not enabled
DO $$
BEGIN
    CREATE EXTENSION IF NOT EXISTS pg_trgm;
EXCEPTION
    WHEN insufficient_privilege THEN
        RAISE NOTICE 'Skipping pg_trgm extension creation due to insufficient privileges';
END $$;

-- =========================
-- 12. Add any missing columns to ensure backward compatibility
-- =========================
DO $$
BEGIN
    -- Add data_source text column if it doesn't exist (for backward compatibility)
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                   WHERE table_schema = 'public' 
                   AND table_name = 'molecules' 
                   AND column_name = 'data_source') THEN
        ALTER TABLE public.molecules 
        ADD COLUMN data_source TEXT;
    END IF;
END $$;

-- =========================
-- 13. Migration comment
-- =========================
COMMENT ON TABLE public.data_sources IS 'Tracks different chemical databases and data sources used in the unified importer';
COMMENT ON TABLE public.molecule_synonyms IS 'Stores alternative names and identifiers for molecules';
COMMENT ON TABLE public.import_jobs IS 'Tracks chemical data import job history and statistics';
COMMENT ON TABLE public.molecule_cross_references IS 'Maps molecule identifiers between different data sources';
COMMENT ON FUNCTION public.find_molecule_by_identifier IS 'Looks up a molecule by any of its identifiers across all data sources';
COMMENT ON FUNCTION public.merge_molecule_data IS 'Merges data from two molecule records, optionally removing the second record';

-- Commit transaction
COMMIT;