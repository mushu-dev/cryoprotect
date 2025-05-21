-- Migration: 022_toxicity_data_migration.sql
-- Purpose: Migrate data from the old toxicity schema to the new optimized schema.

BEGIN;

-- =======================================
-- Migrate toxicity data to new structure
-- =======================================

-- Step 1: Migrate basic toxicity data
INSERT INTO toxicity_data_new (
    molecule_id,
    source,
    source_id,
    toxicity_type,
    species,
    route_of_administration,
    value,
    unit,
    uncertainty,
    confidence_score,
    is_predicted,
    prediction_method,
    reference_doi,
    created_at,
    created_by,
    updated_at
)
SELECT 
    td.molecule_id,
    tas.name AS source,
    ta.assay_id AS source_id,
    CASE
        WHEN aem.endpoint_id IS NOT NULL THEN te.name
        ELSE 'Assay Result'
    END AS toxicity_type,
    ta.organism AS species,
    NULL AS route_of_administration,
    COALESCE(td.activity_value, 0) AS value,
    COALESCE(td.activity_unit, 'score') AS unit,
    NULL AS uncertainty,
    COALESCE(td.reliability_score, 0) AS confidence_score,
    false AS is_predicted,
    NULL AS prediction_method,
    td.citation AS reference_doi,
    td.created_at,
    td.created_by,
    td.updated_at
FROM 
    toxicity_data td
JOIN 
    toxicity_assay ta ON td.assay_id = ta.id
JOIN 
    toxicity_data_source tas ON ta.source_id = tas.id
LEFT JOIN 
    assay_endpoint_mapping aem ON ta.id = aem.assay_id
LEFT JOIN 
    toxicity_endpoint te ON aem.endpoint_id = te.id
ON CONFLICT DO NOTHING;

-- Step 2: Migrate LD50 data
-- We assume LD50-specific data is stored in activity_value with appropriate units
WITH ld50_data AS (
    INSERT INTO toxicity_data_new (
        molecule_id,
        source,
        source_id,
        toxicity_type,
        species,
        route_of_administration,
        value,
        unit,
        uncertainty,
        confidence_score,
        is_predicted,
        prediction_method,
        reference_doi,
        created_at,
        created_by,
        updated_at
    )
    SELECT 
        ts.molecule_id,
        ts.data_source AS source,
        NULL AS source_id,
        'LD50' AS toxicity_type,
        'Rat' AS species, -- Assuming rat as default species for now
        'Oral' AS route_of_administration, -- Assuming oral as default route for now
        ts.score_value AS value,
        COALESCE(ts.score_unit, 'mg/kg') AS unit,
        NULL AS uncertainty,
        ts.confidence AS confidence_score,
        false AS is_predicted,
        NULL AS prediction_method,
        NULL AS reference_doi,
        ts.created_at,
        ts.created_by,
        ts.updated_at
    FROM 
        toxicity_score ts
    WHERE 
        ts.score_type = 'LD50' OR ts.score_type = 'Acute Toxicity'
    RETURNING id, molecule_id
)
INSERT INTO toxicity_ld50 (
    toxicity_data_id,
    organism_strain,
    age_group,
    gender,
    duration,
    duration_unit,
    observation_period,
    observation_period_unit,
    dosing_regime,
    number_of_animals,
    death_count,
    additional_effects
)
SELECT 
    ld.id AS toxicity_data_id,
    NULL AS organism_strain,
    NULL AS age_group,
    NULL AS gender,
    14 AS duration, -- Assuming standard 14-day observation period
    'days' AS duration_unit,
    14 AS observation_period,
    'days' AS observation_period_unit,
    'Single dose' AS dosing_regime,
    NULL AS number_of_animals,
    NULL AS death_count,
    NULL AS additional_effects
FROM 
    ld50_data ld;

-- Step 3: Migrate Tox21 data
WITH tox21_data AS (
    SELECT
        td.id,
        td.molecule_id,
        ta.assay_id,
        ta.assay_target,
        CASE
            WHEN td.hit_call = true THEN 'Active'
            ELSE 'Inactive'
        END AS activity_outcome,
        td.activity_value AS activity_score,
        NULL AS curve_description,
        ta.biological_process_target AS intended_target_family,
        ta.assay_function_type AS intended_target_type,
        ta.assay_format,
        ta.assay_design_type AS assay_format_type
    FROM
        toxicity_data td
    JOIN
        toxicity_assay ta ON td.assay_id = ta.id
    JOIN
        toxicity_data_source tas ON ta.source_id = tas.id
    WHERE
        tas.name = 'Tox21'
)
INSERT INTO toxicity_tox21 (
    toxicity_data_id,
    assay_id,
    assay_target,
    activity_outcome,
    activity_score,
    curve_description,
    intended_target_family,
    intended_target_type,
    assay_format,
    assay_format_type
)
SELECT
    tdn.id AS toxicity_data_id,
    t21.assay_id,
    t21.assay_target,
    t21.activity_outcome,
    t21.activity_score,
    t21.curve_description,
    t21.intended_target_family,
    t21.intended_target_type,
    t21.assay_format,
    t21.assay_format_type
FROM
    tox21_data t21
JOIN
    toxicity_data_new tdn ON tdn.molecule_id = t21.molecule_id
                         AND tdn.source = 'Tox21'
                         AND tdn.source_id = t21.assay_id;

-- Step 4: Migrate classification data
-- Assuming this is derived from toxicity scores with classification scores
INSERT INTO toxicity_classification (
    molecule_id,
    classification_system,
    hazard_class,
    hazard_category,
    hazard_statement,
    hazard_code,
    pictogram,
    signal_word,
    source,
    is_predicted,
    confidence_score,
    created_at,
    created_by,
    updated_at
)
SELECT
    ts.molecule_id,
    'GHS' AS classification_system,
    te.name AS hazard_class,
    te.category AS hazard_category,
    'Toxic' AS hazard_statement, -- Placeholder
    NULL AS hazard_code,
    NULL AS pictogram,
    CASE
        WHEN ts.score_value > 7 THEN 'Danger'
        WHEN ts.score_value > 4 THEN 'Warning'
        ELSE 'Caution'
    END AS signal_word,
    COALESCE(ts.data_source, 'Calculated') AS source,
    CASE
        WHEN cm.method_type = 'computational' THEN true
        ELSE false
    END AS is_predicted,
    ts.confidence AS confidence_score,
    ts.created_at,
    ts.created_by,
    ts.updated_at
FROM
    toxicity_score ts
JOIN
    toxicity_endpoint te ON ts.score_type = te.name
LEFT JOIN
    calculation_method cm ON ts.method_id = cm.id
WHERE
    ts.score_type <> 'LD50' AND ts.score_type <> 'Acute Toxicity';
    
-- =======================================
-- Update molecule toxicity flags
-- =======================================

-- Update toxicity_data_available flag
UPDATE molecule m
SET toxicity_data_available = EXISTS (
    SELECT 1 FROM toxicity_data_new td WHERE td.molecule_id = m.id
);

-- Update toxicity_score based on our new calculation function
UPDATE molecule m
SET toxicity_score = (
    SELECT calculate_toxicity_score(m.id)
)
WHERE toxicity_data_available = true;

-- =======================================
-- Perform initial refresh of materialized views
-- =======================================

SELECT refresh_toxicity_materialized_views();

-- Insert a migration completion record
INSERT INTO migration_history (name, applied_at, description) 
VALUES (
    '022_toxicity_data_migration', 
    NOW(), 
    'Migrated data from old toxicity schema to the new optimized schema'
);

COMMIT;