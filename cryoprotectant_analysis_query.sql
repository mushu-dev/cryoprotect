-- Cryoprotectant Analysis Query
-- This query extracts molecules with their calculated properties for visualization and analysis

WITH property_types AS (
    -- Get all property type ids we're interested in
    SELECT id, name 
    FROM property_types 
    WHERE name IN (
        'molecular_weight', 'logp', 'tpsa', 'h_donors', 'h_acceptors',
        'rotatable_bonds', 'ring_count', 'cryoprotectant_type',
        'h_bond_donor_acceptor_ratio', 'membrane_interaction_score',
        'ice_interaction_potential', 'vitrification_potential', 
        'estimated_toxicity', 'cryoprotectant_score'
    )
),
molecule_properties AS (
    -- Extract all relevant properties for molecules
    SELECT 
        m.id AS molecule_id,
        m.name AS molecule_name,
        m.cid,
        m.smiles,
        pt.name AS property_name,
        COALESCE(mp.numeric_value, 0) AS numeric_value,
        mp.text_value
    FROM molecules m
    CROSS JOIN property_types pt
    LEFT JOIN molecular_properties mp ON m.id = mp.molecule_id AND pt.id = mp.property_type_id
),
molecule_property_pivot AS (
    -- Pivot the properties to get one row per molecule with all properties as columns
    SELECT
        molecule_id,
        molecule_name,
        cid,
        smiles,
        MAX(CASE WHEN property_name = 'molecular_weight' THEN numeric_value END) AS molecular_weight,
        MAX(CASE WHEN property_name = 'logp' THEN numeric_value END) AS logp,
        MAX(CASE WHEN property_name = 'tpsa' THEN numeric_value END) AS tpsa,
        MAX(CASE WHEN property_name = 'h_donors' THEN numeric_value END) AS h_donors,
        MAX(CASE WHEN property_name = 'h_acceptors' THEN numeric_value END) AS h_acceptors,
        MAX(CASE WHEN property_name = 'rotatable_bonds' THEN numeric_value END) AS rotatable_bonds,
        MAX(CASE WHEN property_name = 'ring_count' THEN numeric_value END) AS ring_count,
        MAX(CASE WHEN property_name = 'h_bond_donor_acceptor_ratio' THEN numeric_value END) AS h_bond_ratio,
        MAX(CASE WHEN property_name = 'membrane_interaction_score' THEN numeric_value END) AS membrane_score,
        MAX(CASE WHEN property_name = 'ice_interaction_potential' THEN numeric_value END) AS ice_interaction,
        MAX(CASE WHEN property_name = 'vitrification_potential' THEN numeric_value END) AS vitrification,
        MAX(CASE WHEN property_name = 'estimated_toxicity' THEN numeric_value END) AS toxicity,
        MAX(CASE WHEN property_name = 'cryoprotectant_score' THEN numeric_value END) AS cryo_score,
        MAX(CASE WHEN property_name = 'cryoprotectant_type' THEN text_value END) AS cryo_type
    FROM molecule_properties
    GROUP BY molecule_id, molecule_name, cid, smiles
),
classification AS (
    -- Classify molecules based on their cryoprotectant score
    SELECT
        *,
        CASE 
            WHEN cryo_score >= 8.0 THEN 'Excellent'
            WHEN cryo_score >= 6.5 THEN 'Very Good'
            WHEN cryo_score >= 5.0 THEN 'Good'
            WHEN cryo_score >= 3.5 THEN 'Moderate'
            WHEN cryo_score >= 2.0 THEN 'Poor'
            ELSE 'Very Poor'
        END AS effectiveness_class,
        -- If we have a known cryoprotectant type, use it; otherwise make a prediction
        COALESCE(
            cryo_type,
            CASE 
                WHEN molecular_weight < 150 AND logp > -2 AND logp < 1 THEN 'PENETRATING'
                WHEN molecular_weight >= 150 AND h_donors + h_acceptors > 8 THEN 'NON_PENETRATING'
                ELSE NULL
            END
        ) AS cryoprotectant_class
    FROM molecule_property_pivot
)

-- Final result
SELECT
    c.molecule_id,
    c.molecule_name,
    c.cid,
    c.smiles,
    c.molecular_weight,
    c.logp,
    c.tpsa,
    c.h_donors,
    c.h_acceptors,
    c.rotatable_bonds,
    c.ring_count,
    c.h_bond_ratio,
    c.membrane_score,
    c.ice_interaction,
    c.vitrification,
    c.toxicity,
    c.cryo_score,
    c.cryoprotectant_class,
    c.effectiveness_class,
    -- Add a curated visualization URL for the molecule, either via PubChem or RDKit
    CASE 
        WHEN c.cid IS NOT NULL THEN 
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' || c.cid || '/PNG'
        ELSE NULL
    END AS structure_image_url
FROM classification c
WHERE cryo_score > 0  -- Only include molecules that have been scored
ORDER BY 
    CASE 
        WHEN c.cryo_type IS NOT NULL THEN 0  -- Known cryoprotectants first
        ELSE 1
    END,
    c.cryo_score DESC  -- Then by descending score
LIMIT 1000;  -- Limit to prevent overwhelming result sets