-- CryoProtect v2 - Experiment Properties Population Script
-- This script populates the experiment_properties table with scientifically accurate data

-- Properties for Slow Freezing Protocol - DMSO/EG
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-001-01',
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    'c23ade8d-77a4-4576-a930-7bd73b615578', -- Toxicity property type ID
    78.5, -- Cell Viability (%)
    '%',
    'Experimental: Slow Freezing Protocol - DMSO/EG',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-001-02',
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    '3d5de8ff-3790-4417-a853-3e7e2bf63919', -- Ice Crystal Growth Inhibition property type ID
    12.3, -- Ice Crystal Formation (%)
    '%',
    'Experimental: Slow Freezing Protocol - DMSO/EG',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-001-03',
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    'b3968c67-44d1-4861-b3ff-503116adf238', -- Cell Permeability property type ID
    65.2, -- Recovery Rate (%)
    '%',
    'Experimental: Slow Freezing Protocol - DMSO/EG',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Vitrification Protocol - Trehalose/Glycerol
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-002-01',
    'c3d4e5f6-a7b8-4c5d-9e7f-8a9b0c1d2e3f', -- Trehalose molecule ID
    'c23ade8d-77a4-4576-a930-7bd73b615578', -- Toxicity property type ID
    92.1, -- Cell Viability (%)
    '%',
    'Experimental: Vitrification Protocol - Trehalose/Glycerol',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-002-02',
    'c3d4e5f6-a7b8-4c5d-9e7f-8a9b0c1d2e3f', -- Trehalose molecule ID
    '3d5de8ff-3790-4417-a853-3e7e2bf63919', -- Ice Crystal Growth Inhibition property type ID
    0.5, -- Ice Crystal Formation (%)
    '%',
    'Experimental: Vitrification Protocol - Trehalose/Glycerol',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-002-03',
    'c3d4e5f6-a7b8-4c5d-9e7f-8a9b0c1d2e3f', -- Trehalose molecule ID
    'b3968c67-44d1-4861-b3ff-503116adf238', -- Cell Permeability property type ID
    88.7, -- Recovery Rate (%)
    '%',
    'Experimental: Vitrification Protocol - Trehalose/Glycerol',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Controlled Rate Freezing - EAFS
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-003-01',
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    'c23ade8d-77a4-4576-a930-7bd73b615578', -- Toxicity property type ID
    85.7, -- Cell Viability (%)
    '%',
    'Experimental: Controlled Rate Freezing - EAFS',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-003-02',
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    '3d5de8ff-3790-4417-a853-3e7e2bf63919', -- Ice Crystal Growth Inhibition property type ID
    5.2, -- Ice Crystal Formation (%)
    '%',
    'Experimental: Controlled Rate Freezing - EAFS',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-003-03',
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    'b3968c67-44d1-4861-b3ff-503116adf238', -- Cell Permeability property type ID
    79.3, -- Recovery Rate (%)
    '%',
    'Experimental: Controlled Rate Freezing - EAFS',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Oocyte Cryopreservation - PROH/Sucrose
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-004-01',
    'b8c9d0e1-f2a3-4b8c-9d0e-1f2a3b4c5d6e', -- 1,2-Propanediol molecule ID
    'c23ade8d-77a4-4576-a930-7bd73b615578', -- Toxicity property type ID
    81.3, -- Cell Viability (%)
    '%',
    'Experimental: Oocyte Cryopreservation - PROH/Sucrose',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-004-02',
    'b8c9d0e1-f2a3-4b8c-9d0e-1f2a3b4c5d6e', -- 1,2-Propanediol molecule ID
    '3d5de8ff-3790-4417-a853-3e7e2bf63919', -- Ice Crystal Growth Inhibition property type ID
    3.7, -- Ice Crystal Formation (%)
    '%',
    'Experimental: Oocyte Cryopreservation - PROH/Sucrose',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-004-03',
    'b8c9d0e1-f2a3-4b8c-9d0e-1f2a3b4c5d6e', -- 1,2-Propanediol molecule ID
    'b3968c67-44d1-4861-b3ff-503116adf238', -- Cell Permeability property type ID
    72.8, -- Recovery Rate (%)
    '%',
    'Experimental: Oocyte Cryopreservation - PROH/Sucrose',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Fish Sperm Cryopreservation - Methanol/Glycerol
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-005-01',
    'e5f6a7b8-c9d0-4e5f-9a7b-8c9d0e1f2a3b', -- Methanol molecule ID
    'b3968c67-44d1-4861-b3ff-503116adf238', -- Cell Permeability property type ID
    73.2, -- Motility (%)
    '%',
    'Experimental: Fish Sperm Cryopreservation - Methanol/Glycerol',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-005-02',
    'e5f6a7b8-c9d0-4e5f-9a7b-8c9d0e1f2a3b', -- Methanol molecule ID
    'c23ade8d-77a4-4576-a930-7bd73b615578', -- Toxicity property type ID
    68.9, -- Membrane Integrity (%)
    '%',
    'Experimental: Fish Sperm Cryopreservation - Methanol/Glycerol',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'ep-005-03',
    'e5f6a7b8-c9d0-4e5f-9a7b-8c9d0e1f2a3b', -- Methanol molecule ID
    'b3968c67-44d1-4861-b3ff-503116adf238', -- Cell Permeability property type ID
    65.4, -- Fertilization Capacity (%)
    '%',
    'Experimental: Fish Sperm Cryopreservation - Methanol/Glycerol',
    1,
    '[{"timestamp": "2025-04-18T19:46:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);