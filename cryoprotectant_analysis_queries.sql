-- SQL Queries for Cryoprotectant Analysis

-- 1. Find top-scoring cryoprotectant candidates
SELECT 
    m.name,
    m.smiles,
    m.molecular_formula,
    m.pubchem_cid,
    cp.numeric_value AS cryoprotectant_score
FROM 
    molecules m
JOIN 
    molecular_properties cp ON m.id = cp.molecule_id
JOIN 
    property_types pt ON cp.property_type_id = pt.id
WHERE 
    pt.name = 'cryoprotectant_score'
ORDER BY 
    cp.numeric_value DESC
LIMIT 50;

-- 2. Get complete property profile for a molecule by name
WITH target_molecule AS (
    SELECT id FROM molecules WHERE name ILIKE $1 LIMIT 1
)
SELECT 
    pt.name AS property_name,
    mp.numeric_value,
    pt.units
FROM 
    target_molecule tm
JOIN 
    molecular_properties mp ON tm.id = mp.molecule_id
JOIN 
    property_types pt ON mp.property_type_id = pt.id
WHERE 
    mp.numeric_value IS NOT NULL
ORDER BY 
    pt.name;

-- 3. Compare known cryoprotectants
WITH known_cryoprotectants AS (
    SELECT id, name 
    FROM molecules 
    WHERE 
        name ILIKE '%glycerol%' OR 
        name ILIKE '%dmso%' OR 
        name ILIKE '%dimethyl sulfoxide%' OR
        name ILIKE '%ethylene glycol%' OR
        name ILIKE '%propylene glycol%' OR
        name ILIKE '%trehalose%' OR
        name ILIKE '%sucrose%'
    LIMIT 20
),
property_values AS (
    SELECT 
        k.name AS molecule_name,
        pt.name AS property_name,
        mp.numeric_value
    FROM 
        known_cryoprotectants k
    JOIN 
        molecular_properties mp ON k.id = mp.molecule_id
    JOIN 
        property_types pt ON mp.property_type_id = pt.id
    WHERE 
        pt.name IN (
            'molecular_weight', 'logp', 'tpsa', 'h_donors', 'h_acceptors',
            'membrane_interaction_score', 'ice_interaction_potential',
            'vitrification_potential', 'estimated_toxicity', 'cryoprotectant_score'
        )
)
SELECT 
    molecule_name,
    MAX(CASE WHEN property_name = 'molecular_weight' THEN numeric_value END) AS molecular_weight,
    MAX(CASE WHEN property_name = 'logp' THEN numeric_value END) AS logp,
    MAX(CASE WHEN property_name = 'tpsa' THEN numeric_value END) AS tpsa,
    MAX(CASE WHEN property_name = 'h_donors' THEN numeric_value END) AS h_donors,
    MAX(CASE WHEN property_name = 'h_acceptors' THEN numeric_value END) AS h_acceptors,
    MAX(CASE WHEN property_name = 'membrane_interaction_score' THEN numeric_value END) AS membrane_score,
    MAX(CASE WHEN property_name = 'ice_interaction_potential' THEN numeric_value END) AS ice_interaction,
    MAX(CASE WHEN property_name = 'vitrification_potential' THEN numeric_value END) AS vitrification,
    MAX(CASE WHEN property_name = 'estimated_toxicity' THEN numeric_value END) AS toxicity,
    MAX(CASE WHEN property_name = 'cryoprotectant_score' THEN numeric_value END) AS cryo_score
FROM 
    property_values
GROUP BY 
    molecule_name
ORDER BY 
    MAX(CASE WHEN property_name = 'cryoprotectant_score' THEN numeric_value END) DESC;

-- 4. Find molecules with similar properties to known cryoprotectants
WITH known_cryo AS (
    -- Get average properties of top known cryoprotectants
    SELECT 
        AVG(CASE WHEN pt.name = 'molecular_weight' THEN mp.numeric_value END) AS avg_mw,
        AVG(CASE WHEN pt.name = 'logp' THEN mp.numeric_value END) AS avg_logp,
        AVG(CASE WHEN pt.name = 'tpsa' THEN mp.numeric_value END) AS avg_tpsa,
        AVG(CASE WHEN pt.name = 'h_donors' THEN mp.numeric_value END) AS avg_hd,
        AVG(CASE WHEN pt.name = 'h_acceptors' THEN mp.numeric_value END) AS avg_ha
    FROM 
        molecules m
    JOIN 
        molecular_properties mp ON m.id = mp.molecule_id
    JOIN 
        property_types pt ON mp.property_type_id = pt.id
    WHERE 
        (m.name ILIKE '%glycerol%' OR 
         m.name ILIKE '%dmso%' OR 
         m.name ILIKE '%dimethyl sulfoxide%' OR
         m.name ILIKE '%ethylene glycol%')
    AND pt.name IN ('molecular_weight', 'logp', 'tpsa', 'h_donors', 'h_acceptors')
),
molecule_properties AS (
    -- Get properties for all molecules
    SELECT 
        m.id,
        m.name,
        m.smiles,
        MAX(CASE WHEN pt.name = 'molecular_weight' THEN mp.numeric_value END) AS mw,
        MAX(CASE WHEN pt.name = 'logp' THEN mp.numeric_value END) AS logp,
        MAX(CASE WHEN pt.name = 'tpsa' THEN mp.numeric_value END) AS tpsa,
        MAX(CASE WHEN pt.name = 'h_donors' THEN mp.numeric_value END) AS hd,
        MAX(CASE WHEN pt.name = 'h_acceptors' THEN mp.numeric_value END) AS ha,
        MAX(CASE WHEN pt.name = 'cryoprotectant_score' THEN mp.numeric_value END) AS score
    FROM 
        molecules m
    JOIN 
        molecular_properties mp ON m.id = mp.molecule_id
    JOIN 
        property_types pt ON mp.property_type_id = pt.id
    WHERE 
        pt.name IN ('molecular_weight', 'logp', 'tpsa', 'h_donors', 'h_acceptors', 'cryoprotectant_score')
    GROUP BY 
        m.id, m.name, m.smiles
)
SELECT 
    mp.name,
    mp.smiles,
    mp.mw AS molecular_weight,
    mp.logp,
    mp.tpsa,
    mp.hd AS h_donors,
    mp.ha AS h_acceptors,
    mp.score AS cryoprotectant_score,
    -- Calculate similarity score (lower is more similar)
    (
        ABS(mp.mw - kc.avg_mw) / NULLIF(kc.avg_mw, 0) +
        ABS(mp.logp - kc.avg_logp) / NULLIF(ABS(kc.avg_logp), 0) +
        ABS(mp.tpsa - kc.avg_tpsa) / NULLIF(kc.avg_tpsa, 0) +
        ABS(mp.hd - kc.avg_hd) / NULLIF(kc.avg_hd, 0) +
        ABS(mp.ha - kc.avg_ha) / NULLIF(kc.avg_ha, 0)
    ) AS property_distance
FROM 
    molecule_properties mp
CROSS JOIN 
    known_cryo kc
WHERE 
    -- Exclude known cryoprotectants from the results
    mp.name NOT ILIKE '%glycerol%' AND
    mp.name NOT ILIKE '%dmso%' AND
    mp.name NOT ILIKE '%dimethyl sulfoxide%' AND
    mp.name NOT ILIKE '%ethylene glycol%' AND
    mp.name NOT ILIKE '%propylene glycol%' AND
    mp.name NOT ILIKE '%trehalose%' AND
    mp.name NOT ILIKE '%sucrose%' AND
    mp.name NOT ILIKE '%mannitol%'
    -- Ensure we have all properties
    AND mp.mw IS NOT NULL 
    AND mp.logp IS NOT NULL 
    AND mp.tpsa IS NOT NULL 
    AND mp.hd IS NOT NULL 
    AND mp.ha IS NOT NULL
ORDER BY 
    property_distance ASC
LIMIT 50;

-- 5. Find molecules with specific desirable cryoprotectant properties
SELECT 
    m.name,
    m.smiles,
    m.molecular_formula,
    MAX(CASE WHEN pt.name = 'molecular_weight' THEN mp.numeric_value END) AS molecular_weight,
    MAX(CASE WHEN pt.name = 'logp' THEN mp.numeric_value END) AS logp,
    MAX(CASE WHEN pt.name = 'h_donors' THEN mp.numeric_value END) AS h_donors,
    MAX(CASE WHEN pt.name = 'h_acceptors' THEN mp.numeric_value END) AS h_acceptors,
    MAX(CASE WHEN pt.name = 'ice_interaction_potential' THEN mp.numeric_value END) AS ice_interaction,
    MAX(CASE WHEN pt.name = 'estimated_toxicity' THEN mp.numeric_value END) AS toxicity,
    MAX(CASE WHEN pt.name = 'cryoprotectant_score' THEN mp.numeric_value END) AS cryo_score
FROM 
    molecules m
JOIN 
    molecular_properties mp ON m.id = mp.molecule_id
JOIN 
    property_types pt ON mp.property_type_id = pt.id
WHERE 
    pt.name IN (
        'molecular_weight', 'logp', 'h_donors', 'h_acceptors', 
        'ice_interaction_potential', 'estimated_toxicity', 'cryoprotectant_score'
    )
GROUP BY 
    m.id, m.name, m.smiles, m.molecular_formula
HAVING 
    -- Low molecular weight (for penetration)
    MAX(CASE WHEN pt.name = 'molecular_weight' THEN mp.numeric_value END) < 150
    -- Balanced LogP (not too hydrophobic or hydrophilic)
    AND MAX(CASE WHEN pt.name = 'logp' THEN mp.numeric_value END) BETWEEN -2 AND 1
    -- Good hydrogen bonding capability
    AND MAX(CASE WHEN pt.name = 'h_donors' THEN mp.numeric_value END) >= 2
    AND MAX(CASE WHEN pt.name = 'h_acceptors' THEN mp.numeric_value END) >= 2
    -- Low toxicity
    AND MAX(CASE WHEN pt.name = 'estimated_toxicity' THEN mp.numeric_value END) < 0.6
    -- Good ice interaction
    AND MAX(CASE WHEN pt.name = 'ice_interaction_potential' THEN mp.numeric_value END) > 0.7
ORDER BY 
    MAX(CASE WHEN pt.name = 'cryoprotectant_score' THEN mp.numeric_value END) DESC
LIMIT 50;