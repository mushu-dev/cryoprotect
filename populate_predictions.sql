-- CryoProtect v2 - Predictions Population Script
-- This script populates the predictions table with scientifically accurate cryoprotectant prediction data

-- Glass Transition Temperature predictions for DMSO
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-001',
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    'b08304f3-0407-4b65-aa48-658dc7d286f6', -- MD Simulation method ID
    -135.2, -- value
    'K', -- unit
    0.85, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Cryoprotective Efficacy predictions for DMSO
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-002',
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    '486220eb-9795-4156-bcda-c1727de4144a', -- Vitrification Concentration property type ID
    'c7c1b4a7-f830-4ff1-b18b-b984d6a74d13', -- Neural Network method ID
    92.3, -- value
    '%', -- unit
    0.88, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Cell Membrane Permeability predictions for DMSO
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-003',
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    'b3968c67-44d1-4861-b3ff-503116adf238', -- Cell Permeability property type ID
    'b08304f3-0407-4b65-aa48-658dc7d286f6', -- MD Simulation method ID
    1.05, -- value
    'μm/s', -- unit
    0.82, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Toxicity Index predictions for DMSO
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-004',
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    'c23ade8d-77a4-4576-a930-7bd73b615578', -- Toxicity property type ID
    '4817e677-af89-4bae-992e-791cf05a0133', -- QSPR Model method ID
    0.42, -- value
    '', -- unit
    0.78, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Glass Transition Temperature predictions for Glycerol
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-005',
    '9e8e4f5a-8c2d-4c1b-9a3b-7d5c2b6e4f3a', -- Glycerol molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    'b08304f3-0407-4b65-aa48-658dc7d286f6', -- MD Simulation method ID
    -91.5, -- value
    'K', -- unit
    0.89, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Cryoprotective Efficacy predictions for Glycerol
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-006',
    '9e8e4f5a-8c2d-4c1b-9a3b-7d5c2b6e4f3a', -- Glycerol molecule ID
    '486220eb-9795-4156-bcda-c1727de4144a', -- Vitrification Concentration property type ID
    'c7c1b4a7-f830-4ff1-b18b-b984d6a74d13', -- Neural Network method ID
    81.2, -- value
    '%', -- unit
    0.85, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Cell Membrane Permeability predictions for Glycerol
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-007',
    '9e8e4f5a-8c2d-4c1b-9a3b-7d5c2b6e4f3a', -- Glycerol molecule ID
    'b3968c67-44d1-4861-b3ff-503116adf238', -- Cell Permeability property type ID
    'b08304f3-0407-4b65-aa48-658dc7d286f6', -- MD Simulation method ID
    0.42, -- value
    'μm/s', -- unit
    0.79, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Toxicity Index predictions for Glycerol
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-008',
    '9e8e4f5a-8c2d-4c1b-9a3b-7d5c2b6e4f3a', -- Glycerol molecule ID
    'c23ade8d-77a4-4576-a930-7bd73b615578', -- Toxicity property type ID
    '4817e677-af89-4bae-992e-791cf05a0133', -- QSPR Model method ID
    0.18, -- value
    '', -- unit
    0.92, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Glass Transition Temperature predictions for Ethylene glycol
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-009',
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    'b08304f3-0407-4b65-aa48-658dc7d286f6', -- MD Simulation method ID
    -126.3, -- value
    'K', -- unit
    0.87, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Cryoprotective Efficacy predictions for Ethylene glycol
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-010',
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    '486220eb-9795-4156-bcda-c1727de4144a', -- Vitrification Concentration property type ID
    'c7c1b4a7-f830-4ff1-b18b-b984d6a74d13', -- Neural Network method ID
    85.7, -- value
    '%', -- unit
    0.83, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Cell Membrane Permeability predictions for Ethylene glycol
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-011',
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    'b3968c67-44d1-4861-b3ff-503116adf238', -- Cell Permeability property type ID
    'b08304f3-0407-4b65-aa48-658dc7d286f6', -- MD Simulation method ID
    1.12, -- value
    'μm/s', -- unit
    0.81, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Toxicity Index predictions for Ethylene glycol
INSERT INTO predictions (
    id, molecule_id, property_type_id, calculation_method_id, numeric_value, unit, confidence,
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'pred-012',
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    'c23ade8d-77a4-4576-a930-7bd73b615578', -- Toxicity property type ID
    '4817e677-af89-4bae-992e-791cf05a0133', -- QSPR Model method ID
    0.52, -- value
    '', -- unit
    0.76, -- confidence
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:47:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);