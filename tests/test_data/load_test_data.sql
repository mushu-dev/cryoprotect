-- CryoProtect v2 Test Data Loading Script
-- This script provides an alternative method to load the test datasets directly via SQL
-- Note: This is meant to be run on a development/test database only

-- Set variables
SET @default_user_id = '00000000-0000-0000-0000-000000000001';
SET @timestamp = NOW();

-- Function to generate UUIDs (for databases that don't support uuid_generate_v4())
DELIMITER //
CREATE OR REPLACE FUNCTION generate_uuid() RETURNS CHAR(36)
BEGIN
    RETURN UUID();
END //
DELIMITER ;

-- Clear existing test data (optional - uncomment if needed)
/*
DELETE FROM mixture_components WHERE created_by = @default_user_id;
DELETE FROM mixtures WHERE created_by = @default_user_id;
DELETE FROM molecular_properties WHERE created_by = @default_user_id;
DELETE FROM molecules WHERE created_by = @default_user_id;
*/

-- Insert core cryoprotectants
-- DMSO
SET @dmso_id = generate_uuid();
INSERT INTO molecules (id, name, smiles, molecular_formula, created_at, updated_at, created_by)
VALUES (@dmso_id, 'Dimethyl sulfoxide', 'CS(=O)C', 'C2H6OS', @timestamp, @timestamp, @default_user_id);

-- Glycerol
SET @glycerol_id = generate_uuid();
INSERT INTO molecules (id, name, smiles, molecular_formula, created_at, updated_at, created_by)
VALUES (@glycerol_id, 'Glycerol', 'C(C(CO)O)O', 'C3H8O3', @timestamp, @timestamp, @default_user_id);

-- Ethylene glycol
SET @eg_id = generate_uuid();
INSERT INTO molecules (id, name, smiles, molecular_formula, created_at, updated_at, created_by)
VALUES (@eg_id, 'Ethylene glycol', 'C(CO)O', 'C2H6O2', @timestamp, @timestamp, @default_user_id);

-- Propylene glycol
SET @pg_id = generate_uuid();
INSERT INTO molecules (id, name, smiles, molecular_formula, created_at, updated_at, created_by)
VALUES (@pg_id, 'Propylene glycol', 'CC(O)CO', 'C3H8O2', @timestamp, @timestamp, @default_user_id);

-- Trehalose
SET @trehalose_id = generate_uuid();
INSERT INTO molecules (id, name, smiles, molecular_formula, created_at, updated_at, created_by)
VALUES (@trehalose_id, 'Trehalose', 'C(C1C(C(C(C(O1)O)O)O)O)O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O', 'C12H22O11', @timestamp, @timestamp, @default_user_id);

-- Insert molecular properties (example for DMSO)
-- Get property type IDs
SET @mw_property_id = (SELECT id FROM property_types WHERE name = 'Molecular Weight' LIMIT 1);
SET @logp_property_id = (SELECT id FROM property_types WHERE name = 'LogP' LIMIT 1);
SET @hbd_property_id = (SELECT id FROM property_types WHERE name = 'H-Bond Donors' LIMIT 1);
SET @hba_property_id = (SELECT id FROM property_types WHERE name = 'H-Bond Acceptors' LIMIT 1);

-- Insert DMSO properties
INSERT INTO molecular_properties (id, molecule_id, property_type_id, numeric_value, created_at, updated_at, created_by)
VALUES 
  (generate_uuid(), @dmso_id, @mw_property_id, 78.133, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @dmso_id, @logp_property_id, -1.35, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @dmso_id, @hbd_property_id, 0, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @dmso_id, @hba_property_id, 1, @timestamp, @timestamp, @default_user_id);

-- Insert Glycerol properties
INSERT INTO molecular_properties (id, molecule_id, property_type_id, numeric_value, created_at, updated_at, created_by)
VALUES 
  (generate_uuid(), @glycerol_id, @mw_property_id, 92.094, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @glycerol_id, @logp_property_id, -1.76, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @glycerol_id, @hbd_property_id, 3, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @glycerol_id, @hba_property_id, 3, @timestamp, @timestamp, @default_user_id);

-- Insert Ethylene glycol properties
INSERT INTO molecular_properties (id, molecule_id, property_type_id, numeric_value, created_at, updated_at, created_by)
VALUES 
  (generate_uuid(), @eg_id, @mw_property_id, 62.068, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @eg_id, @logp_property_id, -1.36, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @eg_id, @hbd_property_id, 2, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @eg_id, @hba_property_id, 2, @timestamp, @timestamp, @default_user_id);

-- Insert example mixture (DMSO/EG)
SET @mixture_id = generate_uuid();
INSERT INTO mixtures (id, name, description, created_at, updated_at, created_by)
VALUES (@mixture_id, 'DMSO/EG Mixture', 'Standard 1:1 DMSO/Ethylene glycol mixture for vitrification', @timestamp, @timestamp, @default_user_id);

-- Insert mixture components
INSERT INTO mixture_components (id, mixture_id, molecule_id, concentration, concentration_unit, created_at, updated_at, created_by)
VALUES 
  (generate_uuid(), @mixture_id, @dmso_id, 50.0, '%', @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @mixture_id, @eg_id, 50.0, '%', @timestamp, @timestamp, @default_user_id);

-- Insert example edge case (Water)
SET @water_id = generate_uuid();
INSERT INTO molecules (id, name, smiles, molecular_formula, created_at, updated_at, created_by)
VALUES (@water_id, 'Water', 'O', 'H2O', @timestamp, @timestamp, @default_user_id);

-- Insert Water properties
INSERT INTO molecular_properties (id, molecule_id, property_type_id, numeric_value, created_at, updated_at, created_by)
VALUES 
  (generate_uuid(), @water_id, @mw_property_id, 18.015, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @water_id, @logp_property_id, -0.77, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @water_id, @hbd_property_id, 1, @timestamp, @timestamp, @default_user_id),
  (generate_uuid(), @water_id, @hba_property_id, 1, @timestamp, @timestamp, @default_user_id);

-- Note: This SQL script provides a simplified version of the data loading
-- For complete data loading with all properties and calculations, use the Python script load_test_data.py