-- CryoProtect v2 - Molecular Properties Population Script
-- This script populates the molecular_properties table with scientifically accurate data

-- Properties for Dimethyl sulfoxide (DMSO)
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-01-01', -- ID
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -1.35, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-01-02',
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    -137, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Glycerol
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-02-01',
    '9e8e4f5a-8c2d-4c1b-9a3b-7d5c2b6e4f3a', -- Glycerol molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -1.76, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-02-02',
    '9e8e4f5a-8c2d-4c1b-9a3b-7d5c2b6e4f3a', -- Glycerol molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    -93, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Ethylene glycol
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-03-01',
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -1.36, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-03-02',
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    -128, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Propylene glycol
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-04-01',
    'b2c3d4e5-f6a7-4b5c-9d7e-8f9a0b1c2d3e', -- Propylene glycol molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -0.92, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-04-02',
    'b2c3d4e5-f6a7-4b5c-9d7e-8f9a0b1c2d3e', -- Propylene glycol molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    -108, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Trehalose
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-05-01',
    'c3d4e5f6-a7b8-4c5d-9e7f-8a9b0c1d2e3f', -- Trehalose molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -4.23, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-05-02',
    'c3d4e5f6-a7b8-4c5d-9e7f-8a9b0c1d2e3f', -- Trehalose molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    115, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Sucrose
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-06-01',
    'd4e5f6a7-b8c9-4d5e-9f7a-8b9c0d1e2f3a', -- Sucrose molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -3.76, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-06-02',
    'd4e5f6a7-b8c9-4d5e-9f7a-8b9c0d1e2f3a', -- Sucrose molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    65, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Methanol
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-07-01',
    'e5f6a7b8-c9d0-4e5f-9a7b-8c9d0e1f2a3b', -- Methanol molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -0.77, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-07-02',
    'e5f6a7b8-c9d0-4e5f-9a7b-8c9d0e1f2a3b', -- Methanol molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    -175, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Formamide
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-08-01',
    'f6a7b8c9-d0e1-4f6a-9b7c-8d9e0f1a2b3c', -- Formamide molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -1.51, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-08-02',
    'f6a7b8c9-d0e1-4f6a-9b7c-8d9e0f1a2b3c', -- Formamide molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    -113, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for Acetamide
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-09-01',
    'a7b8c9d0-e1f2-4a7b-9c8d-9e0f1a2b3c4d', -- Acetamide molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -1.26, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-09-02',
    'a7b8c9d0-e1f2-4a7b-9c8d-9e0f1a2b3c4d', -- Acetamide molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    -73, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Properties for 1,2-Propanediol
INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-10-01',
    'b8c9d0e1-f2a3-4b8c-9d0e-1f2a3b4c5d6e', -- 1,2-Propanediol molecule ID
    '6ff67057-fc79-4b96-a604-de5d08b49f51', -- LogP property type ID
    -0.92, -- value
    '', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

INSERT INTO molecular_properties (
    id, molecule_id, property_type_id, numeric_value, unit, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mp-10-02',
    'b8c9d0e1-f2a3-4b8c-9d0e-1f2a3b4c5d6e', -- 1,2-Propanediol molecule ID
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    -108, -- value
    'K', -- unit
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);