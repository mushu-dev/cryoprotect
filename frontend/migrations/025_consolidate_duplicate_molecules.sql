-- Migration: 025_consolidate_duplicate_molecules.sql
-- Purpose: Clean up duplicate molecules in the consolidated_molecules table
-- and ensure data consistency

-- Step 1: Create a temporary table to identify the primary records for duplicate molecules
CREATE TEMP TABLE molecule_consolidation AS
WITH ranked_duplicates AS (
  SELECT 
    id,
    name,
    smiles,
    formula,
    molecular_weight,
    pubchem_cid,
    primary_molecule_id,
    ROW_NUMBER() OVER (
      PARTITION BY smiles 
      ORDER BY 
        -- Prefer records with pubchem_cid
        (pubchem_cid IS NOT NULL) DESC,
        -- Prefer records with more complete data
        (formula IS NOT NULL) DESC,
        (molecular_weight IS NOT NULL) DESC,
        -- Prefer records with canonical naming
        (CASE WHEN name LIKE 'TEST_%' THEN 0 ELSE 1 END) DESC,
        -- Prefer records that are already marked primary
        (primary_molecule_id IS NULL) DESC
    ) as duplicate_rank
  FROM 
    consolidated_molecules
  WHERE 
    smiles IS NOT NULL
)
SELECT 
  id,
  name,
  smiles,
  formula,
  molecular_weight,
  pubchem_cid,
  duplicate_rank,
  CASE WHEN duplicate_rank = 1 THEN true ELSE false END as is_primary
FROM 
  ranked_duplicates;

-- Step 2: Create a map of molecules to be merged
CREATE TEMP TABLE consolidation_map AS
SELECT 
  d.id as duplicate_id,
  p.id as primary_id
FROM 
  molecule_consolidation d
JOIN 
  molecule_consolidation p ON d.smiles = p.smiles AND p.is_primary = true
WHERE 
  d.is_primary = false;

-- Step 3: Update the molecules table to point duplicates to their primary molecules
UPDATE consolidated_molecules
SET 
  primary_molecule_id = cm.primary_id,
  primary_molecule_name = (SELECT name FROM consolidated_molecules WHERE id = cm.primary_id),
  primary_pubchem_cid = (SELECT pubchem_cid FROM consolidated_molecules WHERE id = cm.primary_id),
  primary_molecular_formula = (SELECT formula FROM consolidated_molecules WHERE id = cm.primary_id),
  molecule_status = 'duplicate',
  updated_at = NOW()
FROM 
  consolidation_map cm
WHERE 
  consolidated_molecules.id = cm.duplicate_id;

-- Step 4: Fill in missing data in primary molecules where possible
WITH molecule_data_aggregation AS (
  SELECT
    cm.primary_id,
    COALESCE(
      MAX(CASE WHEN m.formula IS NOT NULL THEN m.formula ELSE NULL END),
      NULL
    ) as formula,
    COALESCE(
      MAX(CASE WHEN m.molecular_weight IS NOT NULL THEN m.molecular_weight ELSE NULL END),
      NULL
    ) as molecular_weight,
    COALESCE(
      MAX(CASE WHEN m.pubchem_cid IS NOT NULL THEN m.pubchem_cid ELSE NULL END),
      NULL
    ) as pubchem_cid
  FROM
    consolidation_map cm
  JOIN
    consolidated_molecules m ON cm.duplicate_id = m.id
  GROUP BY
    cm.primary_id
)
UPDATE consolidated_molecules
SET
  formula = COALESCE(consolidated_molecules.formula, agg.formula),
  molecular_weight = COALESCE(consolidated_molecules.molecular_weight, agg.molecular_weight),
  pubchem_cid = COALESCE(consolidated_molecules.pubchem_cid, agg.pubchem_cid),
  molecule_status = 'primary',
  is_consolidated = true,
  updated_at = NOW()
FROM
  molecule_data_aggregation agg
WHERE
  consolidated_molecules.id = agg.primary_id;

-- Step 5: Update molecular_properties table to ensure all properties 
-- are connected to the primary molecules
UPDATE molecular_properties
SET
  molecule_id = cm.primary_id,
  updated_at = NOW()
FROM
  consolidation_map cm
WHERE
  molecular_properties.molecule_id = cm.duplicate_id;

-- Step 6: Create index on consolidated_molecules for better performance 
-- when querying by SMILES
CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_smiles ON consolidated_molecules(smiles);

-- Step 7: Fill in any remaining missing formulas by computing them from SMILES
-- This could be done in the future with RDKit or other chemical libraries
-- For now, we'll just note that this step exists

-- Step 8: Add a unique constraint on smiles to prevent future duplicates
-- Note: We'll do this in a separate migration after verifying the consolidation worked
-- ALTER TABLE consolidated_molecules ADD CONSTRAINT uq_consolidated_molecules_smiles UNIQUE (smiles);

-- Drop the temporary tables
DROP TABLE IF EXISTS molecule_consolidation;
DROP TABLE IF EXISTS consolidation_map;