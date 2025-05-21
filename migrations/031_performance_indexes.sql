-- Migration: Performance Indexes
-- Created: 2025-05-14T00:31:43.675229
-- Description: Adds performance indexes for common query patterns

BEGIN;

-- Performance indexes for frequently queried columns
CREATE INDEX idx_molecular_properties_unit ON molecular_properties (unit);
CREATE INDEX idx_molecular_properties_data_source ON molecular_properties (data_source);
CREATE INDEX idx_molecular_properties_property_name ON molecular_properties (property_name);
CREATE INDEX idx_molecular_properties_property_value ON molecular_properties (property_value);
CREATE INDEX idx_molecular_properties_property_type ON molecular_properties (property_type);
CREATE INDEX idx_molecular_properties_source ON molecular_properties (source);
CREATE INDEX idx_molecular_properties_version ON molecular_properties (version);
CREATE INDEX idx_molecules_molecule_id ON molecules (molecule_id);
CREATE INDEX idx_molecules_smiles ON molecules (smiles);
CREATE INDEX idx_molecules_inchi ON molecules (inchi);
CREATE INDEX idx_molecules_data_source ON molecules (data_source);
CREATE INDEX idx_molecules_pubchem_link ON molecules (pubchem_link);
CREATE INDEX idx_molecules_molecular_formula ON molecules (molecular_formula);
CREATE INDEX idx_molecules_molecular_weight ON molecules (molecular_weight);
CREATE INDEX idx_molecules_version ON molecules (version);
CREATE INDEX idx_property_calculation_queue_molecule_id ON property_calculation_queue (molecule_id);
CREATE INDEX idx_property_calculation_queue_smiles ON property_calculation_queue (smiles);
CREATE INDEX idx_property_calculation_queue_status ON property_calculation_queue (status);
CREATE INDEX idx_property_calculation_queue_error ON property_calculation_queue (error);
CREATE INDEX idx_property_calculation_queue_properties_added ON property_calculation_queue (properties_added);
CREATE INDEX idx_scientific_data_audit_operation ON scientific_data_audit (operation);
CREATE INDEX idx_scientific_data_audit_ip_address ON scientific_data_audit (ip_address);
CREATE INDEX idx_scientific_data_audit_application_context ON scientific_data_audit (application_context);
CREATE INDEX idx_property_types_property_type_id ON property_types (property_type_id);

-- End of migration
COMMIT;
