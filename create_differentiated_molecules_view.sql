-- View for accessing differentiated molecule groups
-- This view makes it easier to query molecules that have structural differences
-- but similar names or identifiers

-- Create helper function to get all molecules in a differentiation group
CREATE OR REPLACE FUNCTION get_differentiation_group_members(group_id TEXT)
RETURNS TABLE (
    id UUID,
    name TEXT,
    molecular_formula TEXT,
    smiles TEXT,
    pubchem_cid INTEGER,
    differentiation_info JSONB
) AS $$
BEGIN
    RETURN QUERY
    SELECT m.id, m.name, m.molecular_formula, m.smiles, m.pubchem_cid,
           m.properties->'differentiation' AS differentiation_info
    FROM molecules m
    WHERE m.properties->'differentiation'->>'differentiation_group' = group_id;
END;
$$ LANGUAGE plpgsql;

-- Create view for differentiated molecules
CREATE OR REPLACE VIEW differentiated_molecules AS
SELECT 
    m.id,
    m.name,
    m.molecular_formula,
    m.smiles,
    m.pubchem_cid,
    m.properties->'differentiation'->>'differentiation_group' AS differentiation_group,
    (SELECT COUNT(*) 
     FROM molecules 
     WHERE properties->'differentiation'->>'differentiation_group' = 
           m.properties->'differentiation'->>'differentiation_group'
    ) AS group_size
FROM 
    molecules m
WHERE 
    m.properties->'differentiation' IS NOT NULL;

-- Create view that groups differentiated molecules
CREATE OR REPLACE VIEW differentiation_groups AS
WITH groups AS (
    SELECT 
        properties->'differentiation'->>'differentiation_group' AS group_id,
        COUNT(*) AS molecule_count
    FROM 
        molecules
    WHERE 
        properties->'differentiation' IS NOT NULL
    GROUP BY 
        properties->'differentiation'->>'differentiation_group'
)
SELECT 
    g.group_id,
    g.molecule_count,
    json_agg(
        json_build_object(
            'id', m.id,
            'name', m.name,
            'molecular_formula', m.molecular_formula,
            'pubchem_cid', m.pubchem_cid
        )
    ) AS molecules
FROM 
    groups g
JOIN 
    molecules m ON g.group_id = m.properties->'differentiation'->>'differentiation_group'
GROUP BY 
    g.group_id, g.molecule_count;

-- Example queries:

-- 1. Find all differentiated molecules
-- SELECT * FROM differentiated_molecules;

-- 2. Find all molecules in a specific differentiation group
-- SELECT * FROM get_differentiation_group_members('1d968532-913b-4744-b67a-65ee497767a8');

-- 3. Get summary of all differentiation groups
-- SELECT * FROM differentiation_groups;

-- 4. Find molecules with the same name but different structures
-- SELECT 
--     name, 
--     COUNT(*) AS variant_count,
--     json_agg(json_build_object('id', id, 'formula', molecular_formula, 'smiles', smiles)) AS variants
-- FROM 
--     differentiated_molecules
-- GROUP BY 
--     name
-- HAVING 
--     COUNT(*) > 1;