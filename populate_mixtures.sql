-- CryoProtect v2 - Mixtures Population Script
-- This script populates the mixtures table with scientifically accurate cryoprotectant mixture data

-- DMSO/EG Vitrification Solution
INSERT INTO mixtures (
    id, name, description, is_public, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mix-001',
    'DMSO/EG Vitrification Solution',
    '1:1 molar mixture of DMSO and ethylene glycol, widely used for oocyte and embryo vitrification. [Source: Fahy et al., Cryobiology 21.4 (1984): 407-426; Rall & Fahy, Nature 313 (1985): 573-575.]',
    true,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:45:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Trehalose/Glycerol Mixture
INSERT INTO mixtures (
    id, name, description, is_public, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mix-002',
    'Trehalose/Glycerol Mixture',
    '2:1 mass ratio of trehalose to glycerol, used for slow freezing of cells. [Source: Best, B. P., Rejuvenation Research 18.5 (2015): 422-436.]',
    true,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:45:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- EAFS Solution
INSERT INTO mixtures (
    id, name, description, is_public, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mix-003',
    'EAFS Solution',
    'EAFS (ethylene glycol, acetamide, Ficoll, sucrose) solution for embryo vitrification. [Source: Kasai et al., Biology of Reproduction 28.3 (1983): 687-693.]',
    true,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:45:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- PROH/Sucrose Solution
INSERT INTO mixtures (
    id, name, description, is_public, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mix-004',
    'PROH/Sucrose Solution',
    '1.5 mol/L 1,2-propanediol (PROH) with 0.2 mol/L sucrose, used for oocyte cryopreservation. [Source: Fabbri et al., Human Reproduction 16.7 (2001): 1469-1476.]',
    true,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:45:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Methanol/Glycerol Fish Cryoprotectant
INSERT INTO mixtures (
    id, name, description, is_public, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'mix-005',
    'Methanol/Glycerol Fish Cryoprotectant',
    '1.3 mol/L methanol and 0.5 mol/L glycerol, used for fish sperm cryopreservation. [Source: Cabrita et al., Aquaculture 249.1-4 (2005): 533-546.]',
    true,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:45:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);