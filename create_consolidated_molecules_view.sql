-- Create a view for consolidated molecules
-- This view provides a unified way to query molecules, automatically handling consolidation

CREATE OR REPLACE VIEW consolidated_molecules AS
WITH consolidated_mapping AS (
  SELECT 
    m.id,
    CASE 
      WHEN m.properties->>'consolidated_to' IS NOT NULL THEN 
        m.properties->>'consolidated_to'
      ELSE 
        m.id::text
    END AS primary_id,
    m.properties->>'consolidated_to' IS NOT NULL AS is_consolidated
  FROM 
    molecules m
)
SELECT 
  m.id,
  m.name,
  m.smiles,
  m.inchi,
  m.inchikey,
  m.formula,
  m.molecular_weight,
  m.molecular_formula,
  m.pubchem_cid,
  m.chembl_id,
  m.pubchem_link,
  m.is_public,
  m.data_source,
  m.created_by,
  m.created_at,
  m.updated_at,
  cm.is_consolidated,
  cm.primary_id::uuid AS primary_molecule_id,
  pm.name AS primary_molecule_name,
  pm.pubchem_cid AS primary_pubchem_cid,
  pm.molecular_formula AS primary_molecular_formula,
  CASE 
    WHEN cm.is_consolidated THEN 'Consolidated'
    WHEN EXISTS (
      SELECT 1 FROM consolidated_mapping 
      WHERE primary_id = m.id::text AND id != m.id
    ) THEN 'Primary'
    ELSE 'Unique'
  END AS molecule_status
FROM 
  molecules m
JOIN 
  consolidated_mapping cm ON m.id = cm.id
LEFT JOIN 
  molecules pm ON pm.id = cm.primary_id::uuid
ORDER BY 
  molecule_status, name;

COMMENT ON VIEW consolidated_molecules IS 
'View that provides consolidated molecule information, automatically mapping secondary molecules to their primary molecules.
- primary_molecule_id: The ID of the primary molecule (same as id for primaries)
- is_consolidated: Whether this molecule has been consolidated to another
- molecule_status: "Primary", "Consolidated", or "Unique"
';

-- Create a view for primary molecules only
CREATE OR REPLACE VIEW primary_molecules AS
SELECT *
FROM consolidated_molecules
WHERE molecule_status IN ('Primary', 'Unique');

COMMENT ON VIEW primary_molecules IS 
'View that shows only primary and unique molecules, filtering out consolidated molecules.
Use this view when you want to avoid duplicate molecules in your queries.
';

-- Create a function to get the primary molecule ID
CREATE OR REPLACE FUNCTION get_primary_molecule_id(molecule_uuid uuid)
RETURNS uuid AS $$
DECLARE
  primary_id uuid;
BEGIN
  SELECT 
    CASE 
      WHEN properties->>'consolidated_to' IS NOT NULL THEN 
        (properties->>'consolidated_to')::uuid
      ELSE 
        id
    END INTO primary_id
  FROM 
    molecules
  WHERE 
    id = molecule_uuid;
    
  RETURN primary_id;
END;
$$ LANGUAGE plpgsql;

COMMENT ON FUNCTION get_primary_molecule_id(uuid) IS 
'Function that returns the primary molecule ID for any molecule.
If the input molecule has been consolidated, returns its primary molecule ID.
Otherwise returns the input ID unchanged.
';

-- Create a function to check if a molecule is a primary
CREATE OR REPLACE FUNCTION is_primary_molecule(molecule_uuid uuid)
RETURNS boolean AS $$
DECLARE
  result boolean;
BEGIN
  SELECT 
    CASE 
      WHEN EXISTS (
        SELECT 1 FROM molecules
        WHERE properties->>'consolidated_to' = molecule_uuid::text
      ) THEN true
      WHEN EXISTS (
        SELECT 1 FROM molecules
        WHERE id = molecule_uuid AND properties->>'consolidated_to' IS NULL
      ) THEN true
      ELSE false
    END INTO result;
    
  RETURN result;
END;
$$ LANGUAGE plpgsql;

COMMENT ON FUNCTION is_primary_molecule(uuid) IS 
'Function that checks if a molecule is a primary molecule.
Returns true if:
1. The molecule is referenced as a primary by other molecules, or
2. The molecule exists and has not been consolidated
Returns false otherwise.
';