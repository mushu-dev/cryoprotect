-- CryoProtect v2 - Mixture Components Population Script
-- This script populates the mixture_components table with scientifically accurate data

-- Components for DMSO/EG Vitrification Solution
INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-001-01',
    'mix-001', -- DMSO/EG Vitrification Solution
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c', -- DMSO molecule ID
    7.5,
    'mol/L',
    'cryoprotectant',
    NOW(),
    NOW()
);

INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-001-02',
    'mix-001', -- DMSO/EG Vitrification Solution
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    7.5,
    'mol/L',
    'cryoprotectant',
    NOW(),
    NOW()
);

-- Components for Trehalose/Glycerol Mixture
INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-002-01',
    'mix-002', -- Trehalose/Glycerol Mixture
    'c3d4e5f6-a7b8-4c5d-9e7f-8a9b0c1d2e3f', -- Trehalose molecule ID
    0.2,
    'g/mL',
    'cryoprotectant',
    NOW(),
    NOW()
);

INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-002-02',
    'mix-002', -- Trehalose/Glycerol Mixture
    '9e8e4f5a-8c2d-4c1b-9a3b-7d5c2b6e4f3a', -- Glycerol molecule ID
    0.1,
    'g/mL',
    'cryoprotectant',
    NOW(),
    NOW()
);

-- Components for EAFS Solution
INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-003-01',
    'mix-003', -- EAFS Solution
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d', -- Ethylene glycol molecule ID
    3.4,
    'mol/L',
    'cryoprotectant',
    NOW(),
    NOW()
);

INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-003-02',
    'mix-003', -- EAFS Solution
    'a7b8c9d0-e1f2-4a7b-9c8d-9e0f1a2b3c4d', -- Acetamide molecule ID
    2.0,
    'mol/L',
    'cryoprotectant',
    NOW(),
    NOW()
);

INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-003-03',
    'mix-003', -- EAFS Solution
    'd4e5f6a7-b8c9-4d5e-9f7a-8b9c0d1e2f3a', -- Sucrose molecule ID
    0.3,
    'mol/L',
    'osmoprotectant',
    NOW(),
    NOW()
);

-- Components for PROH/Sucrose Solution
INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-004-01',
    'mix-004', -- PROH/Sucrose Solution
    'b8c9d0e1-f2a3-4b8c-9d0e-1f2a3b4c5d6e', -- 1,2-Propanediol molecule ID
    1.5,
    'mol/L',
    'cryoprotectant',
    NOW(),
    NOW()
);

INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-004-02',
    'mix-004', -- PROH/Sucrose Solution
    'd4e5f6a7-b8c9-4d5e-9f7a-8b9c0d1e2f3a', -- Sucrose molecule ID
    0.2,
    'mol/L',
    'osmoprotectant',
    NOW(),
    NOW()
);

-- Components for Methanol/Glycerol Fish Cryoprotectant
INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-005-01',
    'mix-005', -- Methanol/Glycerol Fish Cryoprotectant
    'e5f6a7b8-c9d0-4e5f-9a7b-8c9d0e1f2a3b', -- Methanol molecule ID
    1.3,
    'mol/L',
    'cryoprotectant',
    NOW(),
    NOW()
);

INSERT INTO mixture_components (
    id, mixture_id, molecule_id, concentration, concentration_unit, role,
    created_at, updated_at
) VALUES (
    'mc-005-02',
    'mix-005', -- Methanol/Glycerol Fish Cryoprotectant
    '9e8e4f5a-8c2d-4c1b-9a3b-7d5c2b6e4f3a', -- Glycerol molecule ID
    0.5,
    'mol/L',
    'cryoprotectant',
    NOW(),
    NOW()
);