-- CryoProtect v2 - Experiments Population Script
-- This script populates the experiments table with scientifically accurate cryoprotectant experiment data

-- Slow Freezing Protocol - DMSO/EG
INSERT INTO experiments (
    id, mixture_id, property_type_id, experimental_conditions, 
    created_at, updated_at
) VALUES (
    'exp-001',
    'mix-001', -- DMSO/EG Vitrification Solution
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    '1. Equilibrate cells in 10% DMSO/EG solution for 15 min at 4°C
2. Cool at 1°C/min to -80°C
3. Plunge into liquid nitrogen
Temperature: -196°C
Pressure: 101.3 kPa',
    NOW(),
    NOW()
);

-- Vitrification Protocol - Trehalose/Glycerol
INSERT INTO experiments (
    id, mixture_id, property_type_id, experimental_conditions, 
    created_at, updated_at
) VALUES (
    'exp-002',
    'mix-002', -- Trehalose/Glycerol Mixture
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    '1. Equilibrate cells in 20% Trehalose/Glycerol solution for 5 min at room temperature
2. Plunge directly into liquid nitrogen
Temperature: -196°C
Pressure: 101.3 kPa',
    NOW(),
    NOW()
);

-- Controlled Rate Freezing - EAFS
INSERT INTO experiments (
    id, mixture_id, property_type_id, experimental_conditions, 
    created_at, updated_at
) VALUES (
    'exp-003',
    'mix-003', -- EAFS Solution
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    '1. Equilibrate cells in EAFS solution for 10 min at 4°C
2. Cool at 0.5°C/min to -40°C
3. Cool at 10°C/min to -80°C
4. Plunge into liquid nitrogen
Temperature: -196°C
Pressure: 101.3 kPa',
    NOW(),
    NOW()
);

-- Oocyte Cryopreservation - PROH/Sucrose
INSERT INTO experiments (
    id, mixture_id, property_type_id, experimental_conditions, 
    created_at, updated_at
) VALUES (
    'exp-004',
    'mix-004', -- PROH/Sucrose Solution
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    '1. Expose oocytes to 1.5M PROH for 10 min
2. Transfer to 1.5M PROH + 0.2M sucrose for 5 min
3. Plunge into liquid nitrogen
Temperature: -196°C
Pressure: 101.3 kPa',
    NOW(),
    NOW()
);

-- Fish Sperm Cryopreservation - Methanol/Glycerol
INSERT INTO experiments (
    id, mixture_id, property_type_id, experimental_conditions, 
    created_at, updated_at
) VALUES (
    'exp-005',
    'mix-005', -- Methanol/Glycerol Fish Cryoprotectant
    '4c3e067b-58ca-4224-80c8-540ee21587f3', -- Glass Transition Temperature property type ID
    '1. Mix sperm with Methanol/Glycerol solution at 1:3 ratio
2. Equilibrate for 2 min at 4°C
3. Freeze in 0.5mL straws at 30°C/min to -80°C
4. Store in liquid nitrogen
Temperature: -196°C
Pressure: 101.3 kPa',
    NOW(),
    NOW()
);