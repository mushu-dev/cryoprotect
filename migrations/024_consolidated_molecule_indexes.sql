-- Migration: 024_consolidated_molecule_indexes
-- Description: Adds indexes for optimizing queries involving consolidated molecules

-- Add index for consolidated_to lookups
CREATE INDEX IF NOT EXISTS idx_molecules_consolidated_to 
ON molecules(consolidated_to)
WHERE consolidated_to IS NOT NULL;

-- Add index for differentiation group lookups
CREATE INDEX IF NOT EXISTS idx_molecular_properties_differentiation_group
ON molecular_properties(property_value, molecule_id)
WHERE property_type_id IN (
    SELECT id FROM property_types WHERE name = 'differentiationGroup'
);

-- Create index on the property type to speed up lookups
CREATE INDEX IF NOT EXISTS idx_property_types_name
ON property_types(name)
WHERE name IN ('differentiationGroup', 'differentiationDescription');

-- Comment on indexes to document their purpose
COMMENT ON INDEX idx_molecules_consolidated_to IS 
'Index for efficiently querying molecules by their primary molecule (consolidated_to)';

COMMENT ON INDEX idx_molecular_properties_differentiation_group IS
'Index for efficiently querying molecules in the same differentiation group';

COMMENT ON INDEX idx_property_types_name IS
'Index for efficiently looking up differentiation-related property types';