-- CryoProtect Analyzer Database Schema Migration
-- Compatible with Supabase and PostgreSQL

-- Enable UUID extension
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Enable RLS
ALTER DATABASE postgres SET "app.jwt_secret" TO 'your-jwt-secret';
ALTER DATABASE postgres SET "app.jwt_exp" TO 3600;

-- Create schemas
CREATE SCHEMA IF NOT EXISTS public;
CREATE SCHEMA IF NOT EXISTS private;

-- Set default privileges
ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT ALL ON TABLES TO postgres, anon, authenticated, service_role;
ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT ALL ON FUNCTIONS TO postgres, anon, authenticated, service_role;
ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT ALL ON SEQUENCES TO postgres, anon, authenticated, service_role;

-- Create tables

-- 1. Molecules table (stores basic molecule information)
CREATE TABLE public.molecules (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    cid INTEGER UNIQUE NOT NULL,  -- PubChem Compound ID
    name TEXT,
    molecular_formula TEXT,
    smiles TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id),
    pubchem_link TEXT GENERATED ALWAYS AS ('https://pubchem.ncbi.nlm.nih.gov/compound/' || cid) STORED
);

-- 2. Property types table (defines different types of properties that can be measured)
CREATE TABLE public.property_types (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT UNIQUE NOT NULL,
    description TEXT,
    units TEXT,
    data_type TEXT NOT NULL,  -- 'numeric', 'text', 'boolean', etc.
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- 3. Molecular properties table (stores properties of molecules)
CREATE TABLE public.molecular_properties (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
    property_type_id UUID NOT NULL REFERENCES public.property_types(id),
    numeric_value NUMERIC,
    text_value TEXT,
    boolean_value BOOLEAN,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id),
    UNIQUE(molecule_id, property_type_id)
);

-- 4. Mixtures table (stores information about mixtures of molecules)
CREATE TABLE public.mixtures (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT,
    description TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id)
);

-- 5. Mixture components table (stores which molecules are in which mixtures)
CREATE TABLE public.mixture_components (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    mixture_id UUID NOT NULL REFERENCES public.mixtures(id) ON DELETE CASCADE,
    molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
    concentration NUMERIC NOT NULL,  -- Percentage or molarity
    concentration_unit TEXT NOT NULL,  -- '%', 'M', etc.
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id),
    UNIQUE(mixture_id, molecule_id)
);

-- 6. Calculation methods table (stores different methods used for predictions)
CREATE TABLE public.calculation_methods (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT UNIQUE NOT NULL,
    description TEXT,
    version TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- 7. Predictions table (stores computational predictions of mixture properties)
CREATE TABLE public.predictions (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    mixture_id UUID NOT NULL REFERENCES public.mixtures(id) ON DELETE CASCADE,
    property_type_id UUID NOT NULL REFERENCES public.property_types(id),
    calculation_method_id UUID NOT NULL REFERENCES public.calculation_methods(id),
    numeric_value NUMERIC,
    text_value TEXT,
    boolean_value BOOLEAN,
    confidence NUMERIC,  -- Confidence level (0-1)
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id),
    UNIQUE(mixture_id, property_type_id, calculation_method_id)
);

-- 8. Experiments table (stores experimental data)
CREATE TABLE public.experiments (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    mixture_id UUID NOT NULL REFERENCES public.mixtures(id) ON DELETE CASCADE,
    property_type_id UUID NOT NULL REFERENCES public.property_types(id),
    numeric_value NUMERIC,
    text_value TEXT,
    boolean_value BOOLEAN,
    experimental_conditions TEXT,
    date_performed DATE,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id)
);

-- Create indexes for performance
CREATE INDEX idx_molecules_cid ON public.molecules(cid);
CREATE INDEX idx_molecular_properties_molecule_id ON public.molecular_properties(molecule_id);
CREATE INDEX idx_molecular_properties_property_type_id ON public.molecular_properties(property_type_id);
CREATE INDEX idx_mixture_components_mixture_id ON public.mixture_components(mixture_id);
CREATE INDEX idx_mixture_components_molecule_id ON public.mixture_components(molecule_id);
CREATE INDEX idx_predictions_mixture_id ON public.predictions(mixture_id);
CREATE INDEX idx_predictions_property_type_id ON public.predictions(property_type_id);
CREATE INDEX idx_experiments_mixture_id ON public.experiments(mixture_id);
CREATE INDEX idx_experiments_property_type_id ON public.experiments(property_type_id);

-- Create triggers for updated_at timestamps
CREATE OR REPLACE FUNCTION public.trigger_set_timestamp()
RETURNS TRIGGER AS $$
BEGIN
  NEW.updated_at = NOW();
  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Apply the trigger to all tables with updated_at
CREATE TRIGGER set_timestamp_molecules
BEFORE UPDATE ON public.molecules
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

CREATE TRIGGER set_timestamp_property_types
BEFORE UPDATE ON public.property_types
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

CREATE TRIGGER set_timestamp_molecular_properties
BEFORE UPDATE ON public.molecular_properties
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

CREATE TRIGGER set_timestamp_mixtures
BEFORE UPDATE ON public.mixtures
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

CREATE TRIGGER set_timestamp_mixture_components
BEFORE UPDATE ON public.mixture_components
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

CREATE TRIGGER set_timestamp_calculation_methods
BEFORE UPDATE ON public.calculation_methods
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

CREATE TRIGGER set_timestamp_predictions
BEFORE UPDATE ON public.predictions
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

CREATE TRIGGER set_timestamp_experiments
BEFORE UPDATE ON public.experiments
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

-- Enable Row Level Security (RLS)
ALTER TABLE public.molecules ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.property_types ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.molecular_properties ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.mixtures ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.mixture_components ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.calculation_methods ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.predictions ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.experiments ENABLE ROW LEVEL SECURITY;

-- Create RLS Policies

-- Molecules: everyone can view, only authenticated users can insert, only creators can update/delete
CREATE POLICY "Molecules are viewable by everyone" 
ON public.molecules FOR SELECT USING (true);

CREATE POLICY "Molecules can be inserted by authenticated users" 
ON public.molecules FOR INSERT WITH CHECK (auth.role() = 'authenticated');

CREATE POLICY "Molecules can be updated by the creator" 
ON public.molecules FOR UPDATE USING (auth.uid() = created_by);

CREATE POLICY "Molecules can be deleted by the creator" 
ON public.molecules FOR DELETE USING (auth.uid() = created_by);

-- Property Types: everyone can view, only authenticated users can insert/update/delete
CREATE POLICY "Property types are viewable by everyone" 
ON public.property_types FOR SELECT USING (true);

CREATE POLICY "Property types can be managed by authenticated users" 
ON public.property_types FOR ALL USING (auth.role() = 'authenticated');

-- Molecular Properties: everyone can view, only authenticated users can insert, only creators can update/delete
CREATE POLICY "Molecular properties are viewable by everyone" 
ON public.molecular_properties FOR SELECT USING (true);

CREATE POLICY "Molecular properties can be inserted by authenticated users" 
ON public.molecular_properties FOR INSERT WITH CHECK (auth.role() = 'authenticated');

CREATE POLICY "Molecular properties can be updated by the creator" 
ON public.molecular_properties FOR UPDATE USING (auth.uid() = created_by);

CREATE POLICY "Molecular properties can be deleted by the creator" 
ON public.molecular_properties FOR DELETE USING (auth.uid() = created_by);

-- Mixtures: everyone can view, only authenticated users can insert, only creators can update/delete
CREATE POLICY "Mixtures are viewable by everyone" 
ON public.mixtures FOR SELECT USING (true);

CREATE POLICY "Mixtures can be inserted by authenticated users" 
ON public.mixtures FOR INSERT WITH CHECK (auth.role() = 'authenticated');

CREATE POLICY "Mixtures can be updated by the creator" 
ON public.mixtures FOR UPDATE USING (auth.uid() = created_by);

CREATE POLICY "Mixtures can be deleted by the creator" 
ON public.mixtures FOR DELETE USING (auth.uid() = created_by);

-- Mixture Components: everyone can view, only authenticated users can insert, only creators can update/delete
CREATE POLICY "Mixture components are viewable by everyone" 
ON public.mixture_components FOR SELECT USING (true);

CREATE POLICY "Mixture components can be inserted by authenticated users" 
ON public.mixture_components FOR INSERT WITH CHECK (auth.role() = 'authenticated');

CREATE POLICY "Mixture components can be updated by the creator" 
ON public.mixture_components FOR UPDATE USING (auth.uid() = created_by);

CREATE POLICY "Mixture components can be deleted by the creator" 
ON public.mixture_components FOR DELETE USING (auth.uid() = created_by);

-- Calculation Methods: everyone can view, only authenticated users can insert/update/delete
CREATE POLICY "Calculation methods are viewable by everyone" 
ON public.calculation_methods FOR SELECT USING (true);

CREATE POLICY "Calculation methods can be managed by authenticated users" 
ON public.calculation_methods FOR ALL USING (auth.role() = 'authenticated');

-- Predictions: everyone can view, only authenticated users can insert, only creators can update/delete
CREATE POLICY "Predictions are viewable by everyone" 
ON public.predictions FOR SELECT USING (true);

CREATE POLICY "Predictions can be inserted by authenticated users" 
ON public.predictions FOR INSERT WITH CHECK (auth.role() = 'authenticated');

CREATE POLICY "Predictions can be updated by the creator" 
ON public.predictions FOR UPDATE USING (auth.uid() = created_by);

CREATE POLICY "Predictions can be deleted by the creator" 
ON public.predictions FOR DELETE USING (auth.uid() = created_by);

-- Experiments: everyone can view, only authenticated users can insert, only creators can update/delete
CREATE POLICY "Experiments are viewable by everyone" 
ON public.experiments FOR SELECT USING (true);

CREATE POLICY "Experiments can be inserted by authenticated users" 
ON public.experiments FOR INSERT WITH CHECK (auth.role() = 'authenticated');

CREATE POLICY "Experiments can be updated by the creator" 
ON public.experiments FOR UPDATE USING (auth.uid() = created_by);

CREATE POLICY "Experiments can be deleted by the creator" 
ON public.experiments FOR DELETE USING (auth.uid() = created_by);

-- Insert initial property types based on the Python script
INSERT INTO public.property_types (name, description, units, data_type) VALUES
('Molecular Weight', 'The mass of a molecule', 'g/mol', 'numeric'),
('LogP', 'Octanol-water partition coefficient', NULL, 'numeric'),
('TPSA', 'Topological Polar Surface Area', 'Å²', 'numeric'),
('H-Bond Donors', 'Number of hydrogen bond donor groups', NULL, 'numeric'),
('H-Bond Acceptors', 'Number of hydrogen bond acceptor groups', NULL, 'numeric'),
('Toxicity', 'Toxicity information', NULL, 'text'),
('Stability', 'Chemical stability information', NULL, 'text'),
('Environmental Safety', 'Environmental impact information', NULL, 'text'),
('Hydrogen Bonding Score', 'Score for hydrogen bonding capability', NULL, 'numeric'),
('Solubility Score', 'Score for solubility and polarity', NULL, 'numeric'),
('Membrane Permeability Score', 'Score for membrane permeability', NULL, 'numeric'),
('Toxicity Score', 'Score for toxicity and biocompatibility', NULL, 'numeric'),
('Stability Score', 'Score for stability and reactivity', NULL, 'numeric'),
('Environmental Score', 'Score for environmental safety', NULL, 'numeric'),
('Total Score', 'Overall cryoprotectant score', NULL, 'numeric');

-- Insert initial calculation methods
INSERT INTO public.calculation_methods (name, description, version) VALUES
('PubChem Properties', 'Properties fetched directly from PubChem', '1.0'),
('CryoProtect Scoring', 'Scoring algorithm based on weighted properties', '1.0'),
('Experimental Validation', 'Results from laboratory experiments', '1.0');

-- Create a view for easy querying of molecules with their properties
CREATE OR REPLACE VIEW public.molecule_with_properties AS
SELECT 
    m.id,
    m.cid,
    m.name,
    m.molecular_formula,
    m.smiles,
    m.pubchem_link,
    m.created_at,
    m.updated_at,
    m.created_by,
    jsonb_object_agg(pt.name, 
        CASE 
            WHEN pt.data_type = 'numeric' THEN to_jsonb(mp.numeric_value)
            WHEN pt.data_type = 'text' THEN to_jsonb(mp.text_value)
            WHEN pt.data_type = 'boolean' THEN to_jsonb(mp.boolean_value)
            ELSE NULL
        END
    ) AS properties
FROM 
    public.molecules m
LEFT JOIN 
    public.molecular_properties mp ON m.id = mp.molecule_id
LEFT JOIN 
    public.property_types pt ON mp.property_type_id = pt.id
GROUP BY 
    m.id;

-- Create a view for mixtures with their components
CREATE OR REPLACE VIEW public.mixture_with_components AS
SELECT 
    mix.id,
    mix.name,
    mix.description,
    mix.created_at,
    mix.updated_at,
    mix.created_by,
    jsonb_agg(
        jsonb_build_object(
            'molecule_id', mol.id,
            'cid', mol.cid,
            'name', mol.name,
            'concentration', mc.concentration,
            'concentration_unit', mc.concentration_unit
        )
    ) AS components
FROM 
    public.mixtures mix
LEFT JOIN 
    public.mixture_components mc ON mix.id = mc.mixture_id
LEFT JOIN 
    public.molecules mol ON mc.molecule_id = mol.id
GROUP BY 
    mix.id;

-- Create a function to import molecules from PubChem
CREATE OR REPLACE FUNCTION public.import_molecule_from_pubchem(
    p_cid INTEGER,
    p_user_id UUID
) RETURNS UUID AS $$
DECLARE
    v_molecule_id UUID;
    v_property_type_id UUID;
    v_pubchem_data JSONB;
    v_property_value JSONB;
BEGIN
    -- Check if molecule already exists
    SELECT id INTO v_molecule_id FROM public.molecules WHERE cid = p_cid;
    
    IF v_molecule_id IS NULL THEN
        -- Insert new molecule
        INSERT INTO public.molecules (cid, created_by)
        VALUES (p_cid, p_user_id)
        RETURNING id INTO v_molecule_id;
        
        -- TODO: Fetch data from PubChem API and insert properties
        -- This would be implemented in a separate function or trigger
    END IF;
    
    RETURN v_molecule_id;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Create a function to calculate mixture scores
CREATE OR REPLACE FUNCTION public.calculate_mixture_score(
    p_mixture_id UUID
) RETURNS NUMERIC AS $$
DECLARE
    v_score NUMERIC := 0;
    v_component RECORD;
    v_molecule_score NUMERIC;
    v_concentration NUMERIC;
    v_total_concentration NUMERIC := 0;
BEGIN
    -- Calculate weighted average of component scores
    FOR v_component IN 
        SELECT mc.molecule_id, mc.concentration
        FROM public.mixture_components mc
        WHERE mc.mixture_id = p_mixture_id
    LOOP
        -- Get molecule score
        SELECT mp.numeric_value INTO v_molecule_score
        FROM public.molecular_properties mp
        JOIN public.property_types pt ON mp.property_type_id = pt.id
        WHERE mp.molecule_id = v_component.molecule_id
        AND pt.name = 'Total Score';
        
        IF v_molecule_score IS NOT NULL THEN
            v_score := v_score + (v_molecule_score * v_component.concentration);
            v_total_concentration := v_total_concentration + v_component.concentration;
        END IF;
    END LOOP;
    
    -- Calculate final score
    IF v_total_concentration > 0 THEN
        v_score := v_score / v_total_concentration;
    ELSE
        v_score := 0;
    END IF;
    
    RETURN v_score;
END;
$$ LANGUAGE plpgsql;

-- Create a function to compare prediction with experiment
CREATE OR REPLACE FUNCTION public.compare_prediction_with_experiment(
    p_mixture_id UUID,
    p_property_type_id UUID
) RETURNS JSONB AS $$
DECLARE
    v_prediction RECORD;
    v_experiment RECORD;
    v_result JSONB;
BEGIN
    -- Get prediction
    SELECT 
        p.numeric_value, p.text_value, p.boolean_value, p.confidence,
        cm.name AS method_name
    INTO v_prediction
    FROM public.predictions p
    JOIN public.calculation_methods cm ON p.calculation_method_id = cm.id
    WHERE p.mixture_id = p_mixture_id
    AND p.property_type_id = p_property_type_id
    ORDER BY p.created_at DESC
    LIMIT 1;
    
    -- Get experiment
    SELECT 
        e.numeric_value, e.text_value, e.boolean_value,
        e.experimental_conditions, e.date_performed
    INTO v_experiment
    FROM public.experiments e
    WHERE e.mixture_id = p_mixture_id
    AND e.property_type_id = p_property_type_id
    ORDER BY e.date_performed DESC, e.created_at DESC
    LIMIT 1;
    
    -- Build result
    v_result := jsonb_build_object(
        'prediction', jsonb_build_object(
            'numeric_value', v_prediction.numeric_value,
            'text_value', v_prediction.text_value,
            'boolean_value', v_prediction.boolean_value,
            'confidence', v_prediction.confidence,
            'method', v_prediction.method_name
        ),
        'experiment', jsonb_build_object(
            'numeric_value', v_experiment.numeric_value,
            'text_value', v_experiment.text_value,
            'boolean_value', v_experiment.boolean_value,
            'conditions', v_experiment.experimental_conditions,
            'date', v_experiment.date_performed
        )
    );
    
    -- Add comparison metrics for numeric values
    IF v_prediction.numeric_value IS NOT NULL AND v_experiment.numeric_value IS NOT NULL THEN
        v_result := v_result || jsonb_build_object(
            'difference', v_experiment.numeric_value - v_prediction.numeric_value,
            'percent_error', 
                CASE 
                    WHEN v_experiment.numeric_value = 0 THEN NULL
                    ELSE ABS((v_experiment.numeric_value - v_prediction.numeric_value) / v_experiment.numeric_value) * 100
                END
        );
    END IF;
    
    RETURN v_result;
END;
$$ LANGUAGE plpgsql;

COMMENT ON TABLE public.molecules IS 'Stores basic information about molecules from PubChem';
COMMENT ON TABLE public.property_types IS 'Defines different types of properties that can be measured or predicted';
COMMENT ON TABLE public.molecular_properties IS 'Stores properties of individual molecules';
COMMENT ON TABLE public.mixtures IS 'Stores information about mixtures of molecules';
COMMENT ON TABLE public.mixture_components IS 'Stores which molecules are in which mixtures and their concentrations';
COMMENT ON TABLE public.calculation_methods IS 'Stores different methods used for predictions';
COMMENT ON TABLE public.predictions IS 'Stores computational predictions of mixture properties';
COMMENT ON TABLE public.experiments IS 'Stores experimental data for mixtures';