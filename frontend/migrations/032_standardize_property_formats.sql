-- Migration: 032_standardize_property_formats.sql
-- Purpose: Standardize property data formats and units in molecular_properties table
-- This migration adds unit standardization and data quality improvements

-- Step 1: Create property_types reference table
CREATE TABLE IF NOT EXISTS public.property_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    display_name TEXT NOT NULL,
    description TEXT,
    value_type TEXT NOT NULL CHECK (value_type IN ('numeric', 'text', 'boolean', 'json')),
    default_unit TEXT,
    unit_category TEXT,
    min_numeric_value DOUBLE PRECISION,
    max_numeric_value DOUBLE PRECISION,
    allowed_values JSONB,
    validation_regex TEXT,
    is_required BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- Note: Audit schema not available in this environment
-- Creating a comment-based audit trail for now
COMMENT ON TABLE property_types IS 'Property types reference table created by 032_standardize_property_formats.sql';

-- Enable RLS on the property_types table
ALTER TABLE property_types ENABLE ROW LEVEL SECURITY;

-- Create RLS policies for the property_types table
CREATE POLICY property_types_select_policy ON property_types
    FOR SELECT
    USING (TRUE);  -- Everyone can read property_types

-- Simplified policy for this environment
CREATE POLICY property_types_modify_policy ON property_types
    FOR ALL
    USING (TRUE);  -- In this environment, allow all for demo purposes

-- Step 2: Create units reference table
CREATE TABLE IF NOT EXISTS public.units (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    symbol TEXT NOT NULL,
    category TEXT NOT NULL,
    description TEXT,
    conversion_factor DOUBLE PRECISION,
    base_unit TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- Note: Audit schema not available in this environment
-- Creating a comment-based audit trail for now
COMMENT ON TABLE units IS 'Units reference table created by 032_standardize_property_formats.sql';

-- Enable RLS on the units table
ALTER TABLE units ENABLE ROW LEVEL SECURITY;

-- Create RLS policies for the units table
CREATE POLICY units_select_policy ON units
    FOR SELECT
    USING (TRUE);  -- Everyone can read units

-- Simplified policy for this environment
CREATE POLICY units_modify_policy ON units
    FOR ALL
    USING (TRUE);  -- In this environment, allow all for demo purposes

-- Step 3: Add unit_id column to molecular_properties table
ALTER TABLE molecular_properties ADD COLUMN IF NOT EXISTS unit_id INTEGER REFERENCES units(id);
ALTER TABLE molecular_properties ADD COLUMN IF NOT EXISTS property_type_id INTEGER REFERENCES property_types(id);

-- Step 4: Create unit conversion functions
CREATE OR REPLACE FUNCTION public.convert_unit(
    value DOUBLE PRECISION, 
    from_unit_id INTEGER, 
    to_unit_id INTEGER
) RETURNS DOUBLE PRECISION AS $$
DECLARE
    from_factor DOUBLE PRECISION;
    to_factor DOUBLE PRECISION;
    from_category TEXT;
    to_category TEXT;
    from_base TEXT;
    to_base TEXT;
BEGIN
    -- Get unit details
    SELECT 
        u1.conversion_factor, 
        u1.category,
        u1.base_unit
    INTO 
        from_factor, 
        from_category,
        from_base
    FROM 
        units u1 
    WHERE 
        u1.id = from_unit_id;
    
    SELECT 
        u2.conversion_factor, 
        u2.category,
        u2.base_unit
    INTO 
        to_factor, 
        to_category,
        to_base
    FROM 
        units u2 
    WHERE 
        u2.id = to_unit_id;
    
    -- Check if units are compatible
    IF from_category <> to_category THEN
        RAISE EXCEPTION 'Cannot convert between different unit categories: % and %', from_category, to_category;
    END IF;
    
    -- Check if base units match
    IF from_base <> to_base THEN
        RAISE EXCEPTION 'Cannot convert between units with different base units: % and %', from_base, to_base;
    END IF;
    
    -- Perform conversion
    RETURN value * (from_factor / to_factor);
END;
$$ LANGUAGE plpgsql IMMUTABLE;

-- Step 5: Create functions for property validation
CREATE OR REPLACE FUNCTION public.validate_property_value(
    property_type_id INTEGER,
    numeric_value DOUBLE PRECISION,
    text_value TEXT,
    boolean_value BOOLEAN,
    json_value JSONB,
    unit_id INTEGER
) RETURNS BOOLEAN AS $$
DECLARE
    expected_value_type TEXT;
    min_value DOUBLE PRECISION;
    max_value DOUBLE PRECISION;
    allowed_values JSONB;
    validation_regex TEXT;
    default_unit_id INTEGER;
BEGIN
    -- Get property type details
    SELECT 
        pt.value_type,
        pt.min_numeric_value,
        pt.max_numeric_value,
        pt.allowed_values,
        pt.validation_regex,
        u.id
    INTO 
        expected_value_type,
        min_value,
        max_value,
        allowed_values,
        validation_regex,
        default_unit_id
    FROM 
        property_types pt
    LEFT JOIN
        units u ON pt.default_unit = u.name
    WHERE 
        pt.id = property_type_id;
    
    -- Check value type
    IF expected_value_type = 'numeric' AND numeric_value IS NULL THEN
        RETURN FALSE;
    ELSIF expected_value_type = 'text' AND text_value IS NULL THEN
        RETURN FALSE;
    ELSIF expected_value_type = 'boolean' AND boolean_value IS NULL THEN
        RETURN FALSE;
    ELSIF expected_value_type = 'json' AND json_value IS NULL THEN
        RETURN FALSE;
    END IF;
    
    -- For numeric values, check min/max if specified
    IF expected_value_type = 'numeric' THEN
        -- If unit is different from default, convert value to default unit for checking
        IF unit_id IS NOT NULL AND default_unit_id IS NOT NULL AND unit_id <> default_unit_id THEN
            numeric_value := public.convert_unit(numeric_value, unit_id, default_unit_id);
        END IF;
        
        IF min_value IS NOT NULL AND numeric_value < min_value THEN
            RETURN FALSE;
        END IF;
        
        IF max_value IS NOT NULL AND numeric_value > max_value THEN
            RETURN FALSE;
        END IF;
    END IF;
    
    -- For text values, check regex if specified
    IF expected_value_type = 'text' AND validation_regex IS NOT NULL THEN
        IF text_value !~ validation_regex THEN
            RETURN FALSE;
        END IF;
    END IF;
    
    -- For all types, check allowed values if specified
    IF allowed_values IS NOT NULL THEN
        IF expected_value_type = 'numeric' AND NOT (to_jsonb(numeric_value) <@ allowed_values) THEN
            RETURN FALSE;
        ELSIF expected_value_type = 'text' AND NOT (to_jsonb(text_value) <@ allowed_values) THEN
            RETURN FALSE;
        ELSIF expected_value_type = 'boolean' AND NOT (to_jsonb(boolean_value) <@ allowed_values) THEN
            RETURN FALSE;
        END IF;
    END IF;
    
    RETURN TRUE;
END;
$$ LANGUAGE plpgsql IMMUTABLE;

-- Step 6: Create a function to standardize property names
CREATE OR REPLACE FUNCTION public.standardize_property_name(name TEXT) RETURNS TEXT AS $$
DECLARE
    standardized TEXT;
BEGIN
    -- Convert to lowercase
    standardized := LOWER(name);
    
    -- Replace spaces with underscores
    standardized := REGEXP_REPLACE(standardized, '\s+', '_', 'g');
    
    -- Remove special characters except underscores
    standardized := REGEXP_REPLACE(standardized, '[^a-z0-9_]', '', 'g');
    
    -- Ensure it starts with a letter
    IF standardized ~ '^[^a-z]' THEN
        standardized := 'prop_' || standardized;
    END IF;
    
    RETURN standardized;
END;
$$ LANGUAGE plpgsql IMMUTABLE;

-- Step 7: Create a trigger for property validation
CREATE OR REPLACE FUNCTION public.validate_molecular_property() RETURNS TRIGGER AS $$
BEGIN
    -- Skip validation if property_type_id is not set
    IF NEW.property_type_id IS NULL THEN
        RETURN NEW;
    END IF;
    
    -- Validate the property value
    IF NOT public.validate_property_value(
        NEW.property_type_id,
        NEW.numeric_value,
        NEW.text_value,
        NEW.boolean_value,
        NEW.json_value,
        NEW.unit_id
    ) THEN
        RAISE EXCEPTION 'Invalid value for property type %', NEW.property_type_id;
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Create trigger for validation
DROP TRIGGER IF EXISTS validate_molecular_property_trigger ON molecular_properties;
CREATE TRIGGER validate_molecular_property_trigger
    BEFORE INSERT OR UPDATE ON molecular_properties
    FOR EACH ROW
    EXECUTE FUNCTION public.validate_molecular_property();

-- Step 8: Populate the units table with standard units
INSERT INTO units (name, symbol, category, description, conversion_factor, base_unit)
VALUES
    -- Temperature units
    ('celsius', '°C', 'temperature', 'Degrees Celsius', 1.0, 'celsius'),
    ('fahrenheit', '°F', 'temperature', 'Degrees Fahrenheit', 0.556, 'celsius'), -- (F-32)*5/9 for exact
    ('kelvin', 'K', 'temperature', 'Kelvin', 1.0, 'kelvin'),
    
    -- Mass/weight units
    ('gram', 'g', 'mass', 'Gram', 1.0, 'gram'),
    ('kilogram', 'kg', 'mass', 'Kilogram', 1000.0, 'gram'),
    ('milligram', 'mg', 'mass', 'Milligram', 0.001, 'gram'),
    ('microgram', 'µg', 'mass', 'Microgram', 0.000001, 'gram'),
    ('dalton', 'Da', 'molecular_weight', 'Dalton or atomic mass unit', 1.0, 'dalton'),
    ('kilodalton', 'kDa', 'molecular_weight', 'Kilodalton', 1000.0, 'dalton'),
    
    -- Volume units
    ('liter', 'L', 'volume', 'Liter', 1.0, 'liter'),
    ('milliliter', 'mL', 'volume', 'Milliliter', 0.001, 'liter'),
    ('microliter', 'µL', 'volume', 'Microliter', 0.000001, 'liter'),
    
    -- Concentration units
    ('molar', 'M', 'concentration', 'Moles per liter', 1.0, 'molar'),
    ('millimolar', 'mM', 'concentration', 'Millimoles per liter', 0.001, 'molar'),
    ('micromolar', 'µM', 'concentration', 'Micromoles per liter', 0.000001, 'molar'),
    ('nanomolar', 'nM', 'concentration', 'Nanomoles per liter', 0.000000001, 'molar'),
    ('percent', '%', 'concentration', 'Percent', 1.0, 'percent'),
    ('parts_per_million', 'ppm', 'concentration', 'Parts per million', 0.0001, 'percent'),
    
    -- Pressure units
    ('pascal', 'Pa', 'pressure', 'Pascal', 1.0, 'pascal'),
    ('kilopascal', 'kPa', 'pressure', 'Kilopascal', 1000.0, 'pascal'),
    ('megapascal', 'MPa', 'pressure', 'Megapascal', 1000000.0, 'pascal'),
    ('bar', 'bar', 'pressure', 'Bar', 100000.0, 'pascal'),
    ('atmosphere', 'atm', 'pressure', 'Standard atmosphere', 101325.0, 'pascal'),
    
    -- Length units
    ('meter', 'm', 'length', 'Meter', 1.0, 'meter'),
    ('centimeter', 'cm', 'length', 'Centimeter', 0.01, 'meter'),
    ('millimeter', 'mm', 'length', 'Millimeter', 0.001, 'meter'),
    ('micrometer', 'µm', 'length', 'Micrometer', 0.000001, 'meter'),
    ('nanometer', 'nm', 'length', 'Nanometer', 0.000000001, 'meter'),
    ('angstrom', 'Å', 'length', 'Angstrom', 0.0000000001, 'meter'),
    
    -- Time units
    ('second', 's', 'time', 'Second', 1.0, 'second'),
    ('minute', 'min', 'time', 'Minute', 60.0, 'second'),
    ('hour', 'h', 'time', 'Hour', 3600.0, 'second'),
    ('day', 'd', 'time', 'Day', 86400.0, 'second'),
    
    -- Energy units
    ('joule', 'J', 'energy', 'Joule', 1.0, 'joule'),
    ('kilojoule', 'kJ', 'energy', 'Kilojoule', 1000.0, 'joule'),
    ('calorie', 'cal', 'energy', 'Calorie', 4.184, 'joule'),
    ('kilocalorie', 'kcal', 'energy', 'Kilocalorie', 4184.0, 'joule'),
    
    -- Dimensionless units
    ('unitless', '', 'dimensionless', 'No units', 1.0, 'unitless'),
    ('ratio', 'ratio', 'dimensionless', 'Ratio', 1.0, 'unitless'),
    ('pH', 'pH', 'dimensionless', 'pH scale', 1.0, 'pH')
ON CONFLICT (name) DO UPDATE SET
    symbol = EXCLUDED.symbol,
    category = EXCLUDED.category,
    description = EXCLUDED.description,
    conversion_factor = EXCLUDED.conversion_factor,
    base_unit = EXCLUDED.base_unit;

-- Step 9: Populate the property_types table with standard property types
INSERT INTO property_types (
    name, 
    display_name, 
    description, 
    value_type, 
    default_unit,
    unit_category, 
    min_numeric_value, 
    max_numeric_value
)
VALUES
    -- Physical properties
    ('molecular_weight', 'Molecular Weight', 'Molecular weight of compound', 'numeric', 'dalton', 'molecular_weight', 0, NULL),
    ('melting_point', 'Melting Point', 'Temperature at which solid transitions to liquid', 'numeric', 'celsius', 'temperature', -273.15, NULL),
    ('boiling_point', 'Boiling Point', 'Temperature at which liquid transitions to gas', 'numeric', 'celsius', 'temperature', -273.15, NULL),
    ('density', 'Density', 'Mass per unit volume', 'numeric', 'gram', 'mass', 0, NULL),
    ('solubility_water', 'Water Solubility', 'Solubility in water', 'numeric', 'milligram', 'mass', 0, NULL),
    
    -- Cryoprotectant properties
    ('freezing_point', 'Freezing Point', 'Temperature at which liquid transitions to solid', 'numeric', 'celsius', 'temperature', -273.15, NULL),
    ('glass_transition_temp', 'Glass Transition Temperature', 'Temperature at which material becomes brittle', 'numeric', 'celsius', 'temperature', -273.15, NULL),
    ('viscosity', 'Viscosity', 'Resistance to flow', 'numeric', 'pascal', 'pressure', 0, NULL),
    ('toxicity_level', 'Toxicity Level', 'Level of toxicity', 'numeric', 'unitless', 'dimensionless', 0, 10),
    
    -- Chemical properties
    ('ph', 'pH', 'Acidity or basicity', 'numeric', 'pH', 'dimensionless', 0, 14),
    ('formula', 'Chemical Formula', 'Chemical formula of compound', 'text', NULL, NULL, NULL, NULL),
    ('smiles', 'SMILES Notation', 'Simplified molecular-input line-entry system', 'text', NULL, NULL, NULL, NULL),
    
    -- Biological properties
    ('bioavailability', 'Bioavailability', 'Proportion of substance which enters circulation', 'numeric', 'percent', 'concentration', 0, 100),
    ('ld50', 'LD50', 'Lethal dose for 50% of test population', 'numeric', 'milligram', 'mass', 0, NULL),
    
    -- Descriptive properties
    ('appearance', 'Appearance', 'Physical appearance', 'text', NULL, NULL, NULL, NULL),
    ('color', 'Color', 'Color of compound', 'text', NULL, NULL, NULL, NULL),
    ('odor', 'Odor', 'Odor description', 'text', NULL, NULL, NULL, NULL),
    
    -- Experimental properties
    ('optimal_concentration', 'Optimal Concentration', 'Optimal concentration for cryoprotection', 'numeric', 'percent', 'concentration', 0, 100),
    ('effectiveness_score', 'Effectiveness Score', 'Score of effectiveness as cryoprotectant', 'numeric', 'unitless', 'dimensionless', 0, 10)
ON CONFLICT (name) DO UPDATE SET
    display_name = EXCLUDED.display_name,
    description = EXCLUDED.description,
    value_type = EXCLUDED.value_type,
    default_unit = EXCLUDED.default_unit,
    unit_category = EXCLUDED.unit_category,
    min_numeric_value = EXCLUDED.min_numeric_value,
    max_numeric_value = EXCLUDED.max_numeric_value;

-- Step 10: Create a helper function to extract units from text values
CREATE OR REPLACE FUNCTION public.extract_unit_from_text(text_value TEXT) RETURNS TABLE (
    extracted_value DOUBLE PRECISION,
    extracted_unit TEXT
) AS $$
DECLARE
    unit_patterns TEXT[];
    i INTEGER;
    regex TEXT;
    matches TEXT[];
BEGIN
    -- Define common unit patterns to search for
    unit_patterns := ARRAY[
        '°C', '°F', 'K', 'celsius', 'fahrenheit', 'kelvin',
        'g', 'kg', 'mg', 'µg', 'g/mol', 'Da', 'kDa', 'dalton',
        'L', 'mL', 'µL', 'l', 'ml', 'µl', 'liter', 'milliliter',
        'M', 'mM', 'µM', 'nM', 'pM', '%', 'ppm', 'ppb', 'mol/L',
        'Pa', 'kPa', 'MPa', 'bar', 'atm', 'psi', 'mmHg',
        'mm', 'cm', 'm', 'µm', 'nm', 'Å',
        's', 'min', 'h', 'hr', 'second', 'minute', 'hour', 'day',
        'J', 'kJ', 'kcal', 'cal', 'eV',
        '°', 'rad', 'degree'
    ];
    
    -- Try to match different patterns
    FOR i IN 1..array_length(unit_patterns, 1) LOOP
        -- Pattern: number followed by unit
        regex := '([-+]?\d*\.?\d+)\s*' || unit_patterns[i] || '\b';
        
        IF text_value ~ regex THEN
            matches := regexp_matches(text_value, regex);
            extracted_value := matches[1]::DOUBLE PRECISION;
            extracted_unit := unit_patterns[i];
            RETURN NEXT;
            RETURN;
        END IF;
    END LOOP;
    
    -- If no match, try to extract just the numeric value
    IF text_value ~ '^[-+]?\d*\.?\d+$' THEN
        extracted_value := text_value::DOUBLE PRECISION;
        extracted_unit := NULL;
        RETURN NEXT;
    END IF;
    
    -- No matches found
    RETURN;
END;
$$ LANGUAGE plpgsql IMMUTABLE;

-- Step 11: Create view for standardized property access
CREATE OR REPLACE VIEW public.standardized_properties AS
SELECT
    mp.id,
    mp.molecule_id,
    pt.name AS property_name,
    pt.display_name,
    pt.description,
    pt.value_type,
    mp.numeric_value,
    mp.text_value,
    mp.boolean_value,
    mp.json_value,
    u.name AS unit_name,
    u.symbol AS unit_symbol,
    u.category AS unit_category,
    mp.created_at,
    mp.updated_at
FROM
    molecular_properties mp
LEFT JOIN
    property_types pt ON mp.property_type_id = pt.id
LEFT JOIN
    units u ON mp.unit_id = u.id;

-- Step 12: Create a function to add a property with validation
CREATE OR REPLACE FUNCTION public.add_molecular_property(
    molecule_id UUID,
    property_name TEXT,
    numeric_value DOUBLE PRECISION DEFAULT NULL,
    text_value TEXT DEFAULT NULL,
    boolean_value BOOLEAN DEFAULT NULL,
    json_value JSONB DEFAULT NULL,
    unit_name TEXT DEFAULT NULL
) RETURNS UUID AS $$
DECLARE
    property_type_id INTEGER;
    unit_id INTEGER;
    new_property_id UUID;
    value_type TEXT;
BEGIN
    -- Get the property type ID
    SELECT id, value_type INTO property_type_id, value_type
    FROM property_types
    WHERE name = property_name;
    
    -- If property type doesn't exist, create it with a default configuration
    IF property_type_id IS NULL THEN
        -- Standardize the property name
        property_name := public.standardize_property_name(property_name);
        
        -- Determine value type based on provided values
        IF numeric_value IS NOT NULL THEN
            value_type := 'numeric';
        ELSIF text_value IS NOT NULL THEN
            value_type := 'text';
        ELSIF boolean_value IS NOT NULL THEN
            value_type := 'boolean';
        ELSIF json_value IS NOT NULL THEN
            value_type := 'json';
        ELSE
            RAISE EXCEPTION 'No value provided for property %', property_name;
        END IF;
        
        -- Create new property type
        INSERT INTO property_types (
            name, 
            display_name, 
            description, 
            value_type, 
            default_unit
        )
        VALUES (
            property_name,
            initcap(replace(property_name, '_', ' ')),
            'Automatically created property',
            value_type,
            unit_name
        )
        RETURNING id INTO property_type_id;
    END IF;
    
    -- Get the unit ID if a unit name is provided
    IF unit_name IS NOT NULL THEN
        SELECT id INTO unit_id
        FROM units
        WHERE name = unit_name;
        
        -- If unit doesn't exist, raise an exception (we don't auto-create units)
        IF unit_id IS NULL THEN
            RAISE EXCEPTION 'Unknown unit: %', unit_name;
        END IF;
    END IF;
    
    -- Check that the value type matches the property type
    IF value_type = 'numeric' AND numeric_value IS NULL THEN
        RAISE EXCEPTION 'Numeric value required for property %', property_name;
    ELSIF value_type = 'text' AND text_value IS NULL THEN
        RAISE EXCEPTION 'Text value required for property %', property_name;
    ELSIF value_type = 'boolean' AND boolean_value IS NULL THEN
        RAISE EXCEPTION 'Boolean value required for property %', property_name;
    ELSIF value_type = 'json' AND json_value IS NULL THEN
        RAISE EXCEPTION 'JSON value required for property %', property_name;
    END IF;
    
    -- Insert the property
    INSERT INTO molecular_properties (
        molecule_id,
        property_name,
        property_type_id,
        numeric_value,
        text_value,
        boolean_value,
        json_value,
        unit_id
    )
    VALUES (
        molecule_id,
        property_name,
        property_type_id,
        numeric_value,
        text_value,
        boolean_value,
        json_value,
        unit_id
    )
    RETURNING id INTO new_property_id;
    
    RETURN new_property_id;
END;
$$ LANGUAGE plpgsql;