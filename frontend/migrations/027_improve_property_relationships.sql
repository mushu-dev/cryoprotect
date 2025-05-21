-- Migration: 027_improve_property_relationships.sql
-- Purpose: Improve the relationships between molecules and their properties

-- Step 1: Create an index on molecule properties by property name
-- This will help with queries that filter by property name
CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_name 
ON molecular_properties(property_name);

-- Step 2: Create an index on molecule properties for faster text search
CREATE INDEX IF NOT EXISTS idx_molecular_properties_text_search
ON molecular_properties USING gin(to_tsvector('english', COALESCE(property_name, '') || ' ' || COALESCE(property_value::text, '')));

-- Step 3: Create a view for easier querying of molecule properties
CREATE OR REPLACE VIEW molecule_with_properties AS
SELECT
    m.id,
    m.name,
    m.smiles,
    m.formula,
    m.pubchem_cid,
    m.molecular_weight,
    m.is_public,
    m.created_at,
    m.updated_at,
    m.created_by,
    (
        SELECT jsonb_object_agg(mp.property_name, 
            CASE 
                WHEN mp.numeric_value IS NOT NULL THEN to_jsonb(mp.numeric_value)
                WHEN mp.text_value IS NOT NULL THEN to_jsonb(mp.text_value)
                WHEN mp.boolean_value IS NOT NULL THEN to_jsonb(mp.boolean_value)
                ELSE NULL
            END)
        FROM molecular_properties mp
        WHERE mp.molecule_id = m.id
    ) AS properties
FROM
    molecules m;

-- Add a comment to the view
COMMENT ON VIEW molecule_with_properties IS 'Consolidated view of molecules with all their properties as a JSONB object';

-- Step 4: Create functions to simplify working with molecule properties

-- Function to get a property value for a molecule
CREATE OR REPLACE FUNCTION get_molecule_property(
    molecule_id UUID,
    property_name TEXT
) RETURNS TEXT AS $$
DECLARE
    result TEXT;
BEGIN
    SELECT
        CASE 
            WHEN mp.numeric_value IS NOT NULL THEN mp.numeric_value::TEXT
            WHEN mp.text_value IS NOT NULL THEN mp.text_value
            WHEN mp.boolean_value IS NOT NULL THEN mp.boolean_value::TEXT
            ELSE NULL
        END INTO result
    FROM
        molecular_properties mp
    WHERE
        mp.molecule_id = get_molecule_property.molecule_id
        AND mp.property_name = get_molecule_property.property_name
    LIMIT 1;
    
    RETURN result;
END;
$$ LANGUAGE plpgsql;

-- Function to set a property value for a molecule
CREATE OR REPLACE FUNCTION set_molecule_property(
    molecule_id UUID,
    property_name TEXT,
    property_value TEXT,
    property_type TEXT DEFAULT NULL
) RETURNS VOID AS $$
DECLARE
    existing_id UUID;
    property_type_id UUID;
    data_type TEXT;
BEGIN
    -- Check if property already exists
    SELECT id INTO existing_id
    FROM molecular_properties
    WHERE molecule_id = set_molecule_property.molecule_id
      AND property_name = set_molecule_property.property_name;
      
    -- Determine property type if not provided
    IF property_type IS NULL THEN
        -- Try to get it from existing record
        IF existing_id IS NOT NULL THEN
            SELECT mp.property_type INTO property_type
            FROM molecular_properties mp
            WHERE mp.id = existing_id;
        END IF;
        
        -- If still null, try to guess from value format
        IF property_type IS NULL THEN
            -- Is it numeric?
            IF property_value ~ '^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$' THEN
                property_type := 'numeric';
            -- Is it boolean?
            ELSIF property_value IN ('true', 'false', 't', 'f', 'yes', 'no', 'y', 'n', '1', '0') THEN
                property_type := 'boolean';
            -- Default to text
            ELSE
                property_type := 'text';
            END IF;
        END IF;
    END IF;
    
    -- Get property type ID (create if it doesn't exist)
    SELECT id INTO property_type_id
    FROM property_types
    WHERE name = property_name;
    
    IF property_type_id IS NULL THEN
        INSERT INTO property_types (name, data_type, created_at, updated_at)
        VALUES (property_name, property_type, NOW(), NOW())
        RETURNING id INTO property_type_id;
    END IF;
    
    -- Update or insert the property
    IF existing_id IS NOT NULL THEN
        UPDATE molecular_properties
        SET 
            numeric_value = CASE WHEN property_type = 'numeric' THEN property_value::numeric ELSE NULL END,
            text_value = CASE WHEN property_type = 'text' THEN property_value ELSE NULL END,
            boolean_value = CASE WHEN property_type = 'boolean' THEN (
                property_value IN ('true', 't', 'yes', 'y', '1')
            ) ELSE NULL END,
            property_value = property_value,
            property_type = property_type,
            updated_at = NOW()
        WHERE id = existing_id;
    ELSE
        INSERT INTO molecular_properties (
            molecule_id, 
            property_type_id, 
            numeric_value,
            text_value,
            boolean_value,
            property_name,
            property_value,
            property_type,
            created_at,
            updated_at
        )
        VALUES (
            molecule_id,
            property_type_id,
            CASE WHEN property_type = 'numeric' THEN property_value::numeric ELSE NULL END,
            CASE WHEN property_type = 'text' THEN property_value ELSE NULL END,
            CASE WHEN property_type = 'boolean' THEN (
                property_value IN ('true', 't', 'yes', 'y', '1')
            ) ELSE NULL END,
            property_name,
            property_value,
            property_type,
            NOW(),
            NOW()
        );
    END IF;
END;
$$ LANGUAGE plpgsql;

-- Step 5: Create a function to search molecules by property values
CREATE OR REPLACE FUNCTION search_molecules_by_property(
    property_name TEXT,
    property_value TEXT
) RETURNS TABLE (
    molecule_id UUID,
    molecule_name TEXT,
    molecule_smiles TEXT,
    molecule_formula TEXT,
    value_found TEXT
) AS $$
BEGIN
    RETURN QUERY
    SELECT
        m.id AS molecule_id,
        m.name AS molecule_name,
        m.smiles AS molecule_smiles,
        m.formula AS molecule_formula,
        CASE 
            WHEN mp.numeric_value IS NOT NULL THEN mp.numeric_value::TEXT
            WHEN mp.text_value IS NOT NULL THEN mp.text_value
            WHEN mp.boolean_value IS NOT NULL THEN mp.boolean_value::TEXT
            ELSE mp.property_value
        END AS value_found
    FROM
        molecules m
    JOIN
        molecular_properties mp ON m.id = mp.molecule_id
    WHERE
        mp.property_name = search_molecules_by_property.property_name
        AND (
            (mp.numeric_value IS NOT NULL AND mp.numeric_value::TEXT = property_value)
            OR (mp.text_value IS NOT NULL AND mp.text_value ILIKE '%' || property_value || '%')
            OR (mp.boolean_value IS NOT NULL AND mp.boolean_value::TEXT = property_value)
            OR (mp.property_value ILIKE '%' || property_value || '%')
        );
END;
$$ LANGUAGE plpgsql;

-- Add comments to the functions
COMMENT ON FUNCTION get_molecule_property IS 'Get a property value for a molecule by ID and property name';
COMMENT ON FUNCTION set_molecule_property IS 'Set a property value for a molecule, creating it if it doesn''t exist';
COMMENT ON FUNCTION search_molecules_by_property IS 'Search for molecules by property name and value';