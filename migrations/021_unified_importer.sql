-- Migration: 021_unified_importer.sql
-- Description: Schema changes to support the unified chemical data importer
--
-- This migration adds tables and fields needed for the unified chemical data importer,
-- including data source tracking, molecule synonyms, import job history, and
-- performance optimizations.

-- Only run statements if the migration hasn't been applied
DO $$
BEGIN

-- Check if migration has already been applied
IF NOT EXISTS (
    SELECT 1 FROM pg_tables 
    WHERE tablename = 'data_sources'
) THEN

    ------------------------------------------------------
    -- Create data sources table
    ------------------------------------------------------
    
    CREATE TABLE IF NOT EXISTS data_sources (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        name TEXT NOT NULL UNIQUE,
        description TEXT,
        url TEXT,
        created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
        updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
    );
    
    -- Add comment to table
    COMMENT ON TABLE data_sources IS 'Sources of molecular data such as PubChem, ChEMBL, etc.';
    
    -- Insert standard data sources
    INSERT INTO data_sources (name, description, url) VALUES
        ('PubChem', 'PubChem Database', 'https://pubchem.ncbi.nlm.nih.gov'),
        ('ChEMBL', 'ChEMBL Database', 'https://www.ebi.ac.uk/chembl'),
        ('CryoProtect', 'CryoProtect Database (internal)', NULL);
    
    ------------------------------------------------------
    -- Molecule table enhancements
    ------------------------------------------------------
    
    -- Add data source fields to molecules table
    ALTER TABLE molecules
        ADD COLUMN IF NOT EXISTS data_source_id UUID REFERENCES data_sources(id),
        ADD COLUMN IF NOT EXISTS external_id TEXT,
        ADD COLUMN IF NOT EXISTS inchi TEXT,
        ADD COLUMN IF NOT EXISTS inchikey TEXT,
        ADD COLUMN IF NOT EXISTS import_timestamp TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
        ADD COLUMN IF NOT EXISTS pubchem_cid TEXT,
        ADD COLUMN IF NOT EXISTS chembl_id TEXT;
    
    -- Set existing molecules to CryoProtect source
    UPDATE molecules 
    SET data_source_id = (SELECT id FROM data_sources WHERE name = 'CryoProtect') 
    WHERE data_source_id IS NULL;
    
    -- Add indexes for new fields
    CREATE INDEX IF NOT EXISTS idx_molecules_data_source_id ON molecules (data_source_id);
    CREATE INDEX IF NOT EXISTS idx_molecules_external_id ON molecules (external_id);
    CREATE INDEX IF NOT EXISTS idx_molecules_pubchem_cid ON molecules (pubchem_cid);
    CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON molecules (chembl_id);
    CREATE INDEX IF NOT EXISTS idx_molecules_inchikey ON molecules (inchikey);
    
    -- Add trigram extension for name searching if it doesn't exist
    CREATE EXTENSION IF NOT EXISTS pg_trgm;
    
    -- Add trigram index for name searching
    CREATE INDEX IF NOT EXISTS idx_molecules_name_trigram ON molecules USING GIN (name gin_trgm_ops);
    
    ------------------------------------------------------
    -- Create molecule synonyms table
    ------------------------------------------------------
    
    CREATE TABLE IF NOT EXISTS molecule_synonyms (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        molecule_id UUID NOT NULL REFERENCES molecules(id) ON DELETE CASCADE,
        name TEXT NOT NULL,
        synonym_type TEXT,
        is_preferred BOOLEAN DEFAULT FALSE,
        data_source_id UUID REFERENCES data_sources(id),
        created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
        updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
    );
    
    -- Add comment to table
    COMMENT ON TABLE molecule_synonyms IS 'Alternative names, identifiers, and synonyms for molecules';
    
    -- Add indexes
    CREATE INDEX IF NOT EXISTS idx_molecule_synonyms_molecule_id ON molecule_synonyms(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_molecule_synonyms_name ON molecule_synonyms(name);
    CREATE INDEX IF NOT EXISTS idx_molecule_synonyms_name_trigram ON molecule_synonyms USING GIN (name gin_trgm_ops);
    
    ------------------------------------------------------
    -- Enhance molecular properties table
    ------------------------------------------------------
    
    -- Add data source fields to properties table
    ALTER TABLE molecular_properties
        ADD COLUMN IF NOT EXISTS data_source_id UUID REFERENCES data_sources(id),
        ADD COLUMN IF NOT EXISTS confidence_score NUMERIC,
        ADD COLUMN IF NOT EXISTS data_quality TEXT,
        ADD COLUMN IF NOT EXISTS external_reference TEXT,
        ADD COLUMN IF NOT EXISTS import_timestamp TIMESTAMP WITH TIME ZONE DEFAULT NOW();
    
    -- Set existing properties to CryoProtect source
    UPDATE molecular_properties 
    SET data_source_id = (SELECT id FROM data_sources WHERE name = 'CryoProtect') 
    WHERE data_source_id IS NULL;
    
    -- Add indexes
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_data_source_id ON molecular_properties(data_source_id);
    
    ------------------------------------------------------
    -- Create molecule cross-references table
    ------------------------------------------------------
    
    CREATE TABLE IF NOT EXISTS molecule_cross_references (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        molecule_id UUID NOT NULL REFERENCES molecules(id) ON DELETE CASCADE,
        reference_source_id UUID NOT NULL REFERENCES data_sources(id),
        reference_id TEXT NOT NULL,
        confidence_score NUMERIC,
        relationship_type TEXT,
        created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
        updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
        UNIQUE(molecule_id, reference_source_id, reference_id)
    );
    
    -- Add comment to table
    COMMENT ON TABLE molecule_cross_references IS 'Links between molecules across different data sources';
    
    -- Add indexes
    CREATE INDEX IF NOT EXISTS idx_molecule_cross_references_molecule_id ON molecule_cross_references(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_molecule_cross_references_reference_source_id_reference_id ON molecule_cross_references(reference_source_id, reference_id);
    
    ------------------------------------------------------
    -- Create import job history table
    ------------------------------------------------------
    
    CREATE TABLE IF NOT EXISTS import_jobs (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        data_source_id UUID NOT NULL REFERENCES data_sources(id),
        job_type TEXT NOT NULL, -- e.g., 'import', 'update', 'verify'
        start_time TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
        end_time TIMESTAMP WITH TIME ZONE,
        total_records INTEGER,
        success_count INTEGER,
        error_count INTEGER,
        skipped_count INTEGER,
        query_parameters JSONB,
        checkpoint_file TEXT,
        log_file TEXT,
        status TEXT NOT NULL -- e.g., 'running', 'completed', 'failed'
    );
    
    -- Add comment to table
    COMMENT ON TABLE import_jobs IS 'History of chemical data import jobs';
    
    -- Add indexes
    CREATE INDEX IF NOT EXISTS idx_import_jobs_data_source_id ON import_jobs(data_source_id);
    CREATE INDEX IF NOT EXISTS idx_import_jobs_status ON import_jobs(status);
    
    ------------------------------------------------------
    -- Create import job errors table
    ------------------------------------------------------
    
    CREATE TABLE IF NOT EXISTS import_job_errors (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        import_job_id UUID NOT NULL REFERENCES import_jobs(id) ON DELETE CASCADE,
        record_identifier TEXT,
        error_type TEXT,
        error_message TEXT,
        error_details JSONB,
        created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
    );
    
    -- Add comment to table
    COMMENT ON TABLE import_job_errors IS 'Error details for failed import job records';
    
    -- Add indexes
    CREATE INDEX IF NOT EXISTS idx_import_job_errors_import_job_id ON import_job_errors(import_job_id);
    
    ------------------------------------------------------
    -- Standardize property types for unified importer
    ------------------------------------------------------
    
    -- Create a type lookup table if it doesn't exist
    CREATE TABLE IF NOT EXISTS property_types (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        name TEXT UNIQUE NOT NULL,
        description TEXT,
        data_type TEXT, -- 'numeric', 'text', 'boolean'
        category TEXT, -- 'physicochemical', 'structural', 'biological', etc.
        standard_unit TEXT,
        created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
        updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
    );
    
    -- Add comment to table
    COMMENT ON TABLE property_types IS 'Standardized property types for molecular data';
    
    -- Insert standard property types if table is empty
    INSERT INTO property_types (name, description, data_type, category, standard_unit)
    SELECT * FROM (
        VALUES
            ('LogP', 'Octanol-water partition coefficient', 'numeric', 'physicochemical', 'log units'),
            ('MolecularWeight', 'Molecular weight', 'numeric', 'physicochemical', 'g/mol'),
            ('TPSA', 'Topological polar surface area', 'numeric', 'physicochemical', 'Å²'),
            ('HBondDonorCount', 'Number of hydrogen bond donors', 'numeric', 'structural', 'count'),
            ('HBondAcceptorCount', 'Number of hydrogen bond acceptors', 'numeric', 'structural', 'count'),
            ('RotatableBondCount', 'Number of rotatable bonds', 'numeric', 'structural', 'count'),
            ('HeavyAtomCount', 'Number of non-hydrogen atoms', 'numeric', 'structural', 'count'),
            ('pKa', 'Acid dissociation constant', 'numeric', 'physicochemical', 'pKa units'),
            ('SolubilityAqueous', 'Aqueous solubility', 'numeric', 'physicochemical', 'mg/mL'),
            ('MeltingPoint', 'Melting point', 'numeric', 'physicochemical', '°C'),
            ('BoilingPoint', 'Boiling point', 'numeric', 'physicochemical', '°C'),
            ('Complexity', 'Molecular complexity score', 'numeric', 'structural', 'score')
    ) AS v(name, description, data_type, category, standard_unit)
    WHERE NOT EXISTS (SELECT 1 FROM property_types);
    
    ------------------------------------------------------
    -- Create helpful functions for the importer
    ------------------------------------------------------
    
    -- Function to find a molecule by any identifier
    CREATE OR REPLACE FUNCTION find_molecule_by_identifier(
        p_identifier TEXT
    ) RETURNS UUID AS $$
    DECLARE
        v_molecule_id UUID;
    BEGIN
        -- Try exact match on id
        BEGIN
            v_molecule_id := p_identifier::UUID;
            IF EXISTS (SELECT 1 FROM molecules WHERE id = v_molecule_id) THEN
                RETURN v_molecule_id;
            END IF;
        EXCEPTION WHEN OTHERS THEN
            -- Not a UUID, continue with other checks
        END;
        
        -- Check molecule name
        SELECT id INTO v_molecule_id FROM molecules WHERE name = p_identifier LIMIT 1;
        IF v_molecule_id IS NOT NULL THEN
            RETURN v_molecule_id;
        END IF;
        
        -- Check inchikey
        SELECT id INTO v_molecule_id FROM molecules WHERE inchikey = p_identifier LIMIT 1;
        IF v_molecule_id IS NOT NULL THEN
            RETURN v_molecule_id;
        END IF;
        
        -- Check PubChem CID
        SELECT id INTO v_molecule_id FROM molecules WHERE pubchem_cid = p_identifier LIMIT 1;
        IF v_molecule_id IS NOT NULL THEN
            RETURN v_molecule_id;
        END IF;
        
        -- Check ChEMBL ID
        SELECT id INTO v_molecule_id FROM molecules WHERE chembl_id = p_identifier LIMIT 1;
        IF v_molecule_id IS NOT NULL THEN
            RETURN v_molecule_id;
        END IF;
        
        -- Check external_id
        SELECT id INTO v_molecule_id FROM molecules WHERE external_id = p_identifier LIMIT 1;
        IF v_molecule_id IS NOT NULL THEN
            RETURN v_molecule_id;
        END IF;
        
        -- Check synonyms
        SELECT molecule_id INTO v_molecule_id FROM molecule_synonyms WHERE name = p_identifier LIMIT 1;
        IF v_molecule_id IS NOT NULL THEN
            RETURN v_molecule_id;
        END IF;
        
        -- Check cross references
        SELECT molecule_id INTO v_molecule_id FROM molecule_cross_references WHERE reference_id = p_identifier LIMIT 1;
        IF v_molecule_id IS NOT NULL THEN
            RETURN v_molecule_id;
        END IF;
        
        -- No match found
        RETURN NULL;
    END;
    $$ LANGUAGE plpgsql;
    
    -- Function to merge molecule data
    CREATE OR REPLACE FUNCTION merge_molecule_data(
        p_source_id UUID,
        p_target_id UUID
    ) RETURNS BOOLEAN AS $$
    DECLARE
        v_result BOOLEAN := FALSE;
    BEGIN
        -- Don't merge if IDs are the same
        IF p_source_id = p_target_id THEN
            RETURN FALSE;
        END IF;
        
        -- Move properties
        UPDATE molecular_properties
        SET molecule_id = p_target_id
        WHERE molecule_id = p_source_id;
        
        -- Move synonyms
        UPDATE molecule_synonyms
        SET molecule_id = p_target_id
        WHERE molecule_id = p_source_id;
        
        -- Move cross references
        UPDATE molecule_cross_references
        SET molecule_id = p_target_id
        WHERE molecule_id = p_source_id;
        
        -- Create a cross reference between the two molecules
        INSERT INTO molecule_cross_references (molecule_id, reference_source_id, reference_id, confidence_score)
        VALUES (
            p_target_id,
            (SELECT data_source_id FROM molecules WHERE id = p_source_id),
            (SELECT external_id FROM molecules WHERE id = p_source_id),
            1.0
        )
        ON CONFLICT DO NOTHING;
        
        -- Mark source as deleted or archive it
        -- (Not actually deleting to maintain references)
        UPDATE molecules
        SET name = 'MERGED: ' || name
        WHERE id = p_source_id;
        
        RETURN TRUE;
    END;
    $$ LANGUAGE plpgsql;
    
    -- Function to update InChI from SMILES
    CREATE OR REPLACE FUNCTION update_inchi_data() RETURNS TRIGGER AS $$
    BEGIN
        -- In a real implementation, this would call an RDKit function
        -- or other chemistry toolkit to convert SMILES to InChI and InChIKey
        -- For now, we just set placeholders
        
        -- This is a placeholder - in production, use RDKit or similar
        IF NEW.smiles IS NOT NULL AND (NEW.inchi IS NULL OR NEW.inchikey IS NULL) THEN
            -- In real implementation: 
            -- NEW.inchi := rdkit.mol_to_inchi(NEW.smiles);
            -- NEW.inchikey := rdkit.inchi_to_inchikey(NEW.inchi);
            
            -- Placeholder for now
            NEW.inchi := 'InChI derived from ' || NEW.smiles;
            NEW.inchikey := 'InChIKey derived from ' || NEW.smiles;
        END IF;
        
        RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    
    -- Create trigger for InChI updates
    DROP TRIGGER IF EXISTS trigger_update_inchi_data ON molecules;
    CREATE TRIGGER trigger_update_inchi_data
    BEFORE INSERT OR UPDATE OF smiles ON molecules
    FOR EACH ROW
    EXECUTE FUNCTION update_inchi_data();
    
    ------------------------------------------------------
    -- Add RLS policies for new tables
    ------------------------------------------------------
    
    -- Only add these if RLS is enabled
    DO $$
    BEGIN
        -- Check if RLS is enabled
        IF EXISTS (SELECT 1 FROM pg_tables WHERE tablename = 'molecules' AND rowsecurity = true) THEN
            -- RLS for data_sources table
            ALTER TABLE data_sources ENABLE ROW LEVEL SECURITY;
            
            -- Everyone can view data sources
            CREATE POLICY data_sources_select_policy 
            ON data_sources FOR SELECT 
            USING (true);
            
            -- Only authenticated users can modify data sources
            CREATE POLICY data_sources_insert_policy 
            ON data_sources FOR INSERT 
            WITH CHECK (auth.role() = 'authenticated');
            
            CREATE POLICY data_sources_update_policy 
            ON data_sources FOR UPDATE 
            USING (auth.role() = 'authenticated');
            
            -- RLS for molecule_synonyms table
            ALTER TABLE molecule_synonyms ENABLE ROW LEVEL SECURITY;
            
            -- Everyone can view synonyms
            CREATE POLICY molecule_synonyms_select_policy 
            ON molecule_synonyms FOR SELECT 
            USING (true);
            
            -- Only authenticated users can modify synonyms
            CREATE POLICY molecule_synonyms_insert_policy 
            ON molecule_synonyms FOR INSERT 
            WITH CHECK (auth.role() = 'authenticated');
            
            CREATE POLICY molecule_synonyms_update_policy 
            ON molecule_synonyms FOR UPDATE 
            USING (auth.role() = 'authenticated');
            
            CREATE POLICY molecule_synonyms_delete_policy 
            ON molecule_synonyms FOR DELETE 
            USING (auth.role() = 'authenticated');
            
            -- RLS for molecule_cross_references table
            ALTER TABLE molecule_cross_references ENABLE ROW LEVEL SECURITY;
            
            -- Everyone can view cross references
            CREATE POLICY molecule_cross_references_select_policy 
            ON molecule_cross_references FOR SELECT 
            USING (true);
            
            -- Only authenticated users can modify cross references
            CREATE POLICY molecule_cross_references_insert_policy 
            ON molecule_cross_references FOR INSERT 
            WITH CHECK (auth.role() = 'authenticated');
            
            CREATE POLICY molecule_cross_references_update_policy 
            ON molecule_cross_references FOR UPDATE 
            USING (auth.role() = 'authenticated');
            
            CREATE POLICY molecule_cross_references_delete_policy 
            ON molecule_cross_references FOR DELETE 
            USING (auth.role() = 'authenticated');
            
            -- RLS for import_jobs table
            ALTER TABLE import_jobs ENABLE ROW LEVEL SECURITY;
            
            -- Everyone can view import jobs
            CREATE POLICY import_jobs_select_policy 
            ON import_jobs FOR SELECT 
            USING (true);
            
            -- Only authenticated users can modify import jobs
            CREATE POLICY import_jobs_insert_policy 
            ON import_jobs FOR INSERT 
            WITH CHECK (auth.role() = 'authenticated');
            
            CREATE POLICY import_jobs_update_policy 
            ON import_jobs FOR UPDATE 
            USING (auth.role() = 'authenticated');
            
            -- RLS for import_job_errors table
            ALTER TABLE import_job_errors ENABLE ROW LEVEL SECURITY;
            
            -- Everyone can view import job errors
            CREATE POLICY import_job_errors_select_policy 
            ON import_job_errors FOR SELECT 
            USING (true);
            
            -- Only authenticated users can modify import job errors
            CREATE POLICY import_job_errors_insert_policy 
            ON import_job_errors FOR INSERT 
            WITH CHECK (auth.role() = 'authenticated');
            
            -- RLS for property_types table
            ALTER TABLE property_types ENABLE ROW LEVEL SECURITY;
            
            -- Everyone can view property types
            CREATE POLICY property_types_select_policy 
            ON property_types FOR SELECT 
            USING (true);
            
            -- Only authenticated users can modify property types
            CREATE POLICY property_types_insert_policy 
            ON property_types FOR INSERT 
            WITH CHECK (auth.role() = 'authenticated');
            
            CREATE POLICY property_types_update_policy 
            ON property_types FOR UPDATE 
            USING (auth.role() = 'authenticated');
        END IF;
    END $$;

END IF;

END $$;