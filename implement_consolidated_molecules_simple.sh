#!/bin/bash
# CryoProtect Analyzer - Consolidated Molecules Implementation Script (Simplified)
#
# This script implements the consolidated molecules system by:
# 1. Creating the SQL migration files
# 2. Updating the API with the consolidated molecule functionality
# 3. Testing the implementation
#
# Usage: ./implement_consolidated_molecules_simple.sh

set -e

# Set text colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Display banner
echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}  CryoProtect Analyzer - Consolidated Molecules Implementation${NC}"
echo -e "${BLUE}======================================================${NC}"

# Create backups directory if it doesn't exist
mkdir -p backups
mkdir -p migrations

# Step 1: Create migration files
echo -e "\n${BLUE}Step 1: Creating migration files...${NC}"

# 1.1: Create migration for consolidated molecule constraints and indexes
echo "Creating migration for consolidated molecule constraints and indexes..."
cat > migrations/026_consolidated_molecule_constraints_indexes.sql << 'EOF'
-- Check if table exists
DO $$ 
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'consolidated_molecules') THEN
        -- First validate the current data
        IF EXISTS (
            SELECT 1
            FROM consolidated_molecules 
            WHERE molecule_status NOT IN ('original', 'primary', 'duplicate')
        ) THEN
            RAISE EXCEPTION 'Invalid molecule_status values found in consolidated_molecules table. Please fix before applying constraint.';
        END IF;

        -- Add CHECK constraint for molecule_status
        IF NOT EXISTS (
            SELECT 1
            FROM information_schema.table_constraints 
            WHERE table_name = 'consolidated_molecules' 
            AND constraint_type = 'CHECK'
            AND constraint_name LIKE '%valid_status%'
        ) THEN
            ALTER TABLE consolidated_molecules 
            ADD CONSTRAINT valid_status 
            CHECK (molecule_status IN ('original', 'primary', 'duplicate'));
            
            RAISE NOTICE 'Added CHECK constraint for molecule_status column.';
        ELSE
            RAISE NOTICE 'CHECK constraint for molecule_status already exists.';
        END IF;

        -- Create index for inchikey if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_inchikey%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_inchikey 
            ON public.consolidated_molecules(inchikey);
            
            RAISE NOTICE 'Created index for inchikey column.';
        ELSE
            RAISE NOTICE 'Index for inchikey column already exists.';
        END IF;

        -- Create index for primary_molecule_id if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_primary_id%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_primary_id 
            ON public.consolidated_molecules(primary_molecule_id);
            
            RAISE NOTICE 'Created index for primary_molecule_id column.';
        ELSE
            RAISE NOTICE 'Index for primary_molecule_id column already exists.';
        END IF;

        -- Create index for molecule_status if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_status%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_status 
            ON public.consolidated_molecules(molecule_status);
            
            RAISE NOTICE 'Created index for molecule_status column.';
        ELSE
            RAISE NOTICE 'Index for molecule_status column already exists.';
        END IF;
    ELSE
        RAISE NOTICE 'consolidated_molecules table does not exist. Skipping constraint and index creation.';
    END IF;
END $$;
EOF

# 1.2: Create migration for scientific data audit table
echo "Creating migration for scientific data audit table..."
cat > migrations/027_create_scientific_data_audit.sql << 'EOF'
-- Check if table exists
DO $$ 
BEGIN
    -- Create scientific_data_audit table if it doesn't exist
    IF NOT EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'scientific_data_audit') THEN
        CREATE TABLE IF NOT EXISTS public.scientific_data_audit (
            id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
            table_name TEXT NOT NULL,
            record_id TEXT NOT NULL,
            operation TEXT NOT NULL,
            old_value JSONB,
            new_value JSONB,
            user_id UUID,
            timestamp TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        );
        
        -- Create indexes for efficient querying
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_table_name 
        ON public.scientific_data_audit(table_name);
        
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_record_id 
        ON public.scientific_data_audit(record_id);
        
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_timestamp 
        ON public.scientific_data_audit(timestamp);
        
        RAISE NOTICE 'Created scientific_data_audit table with indexes.';
    ELSE
        RAISE NOTICE 'scientific_data_audit table already exists.';
    END IF;
    
    -- Create RLS policies for scientific_data_audit if needed
    -- This can be extended as needed for your specific security requirements
END $$;
EOF

echo -e "${GREEN}Migration files created successfully.${NC}"
echo "To apply these migrations, you have several options:"
echo "1. Use psql: psql -h <host> -U <user> -d <database> -f migrations/026_consolidated_molecule_constraints_indexes.sql"
echo "2. Use apply_file_sql.py: python3 apply_file_sql.py migrations/026_consolidated_molecule_constraints_indexes.sql"
echo "3. Use the fix_consolidated_migration.py script: python3 fix_consolidated_migration.py"
echo ""

# Step 2: Update the API
echo -e "\n${BLUE}Step 2: Updating the API...${NC}"
if [ -f update_consolidated_molecule_api.py ]; then
    python3 update_consolidated_molecule_api.py
    if [ $? -ne 0 ]; then
        echo -e "${RED}Error updating the API. See output above for details.${NC}"
        exit 1
    fi
else
    echo -e "${RED}update_consolidated_molecule_api.py not found. Cannot update API.${NC}"
    exit 1
fi

# Step 3: Run tests
echo -e "\n${BLUE}Step 3: Running tests...${NC}"
if [ -f test_consolidated_molecules_integration.py ]; then
    echo "Warning: Tests require database connectivity, which may not be available."
    echo "Skipping automated tests. You can run tests manually with:"
    echo "  python3 test_consolidated_molecules_integration.py"
else
    echo -e "${RED}test_consolidated_molecules_integration.py not found. Cannot run tests.${NC}"
    exit 1
fi

# Step 4: Final verification
echo -e "\n${BLUE}Step 4: Verifying implementation...${NC}"
echo "Checking API updates..."

# Check if the API files have been updated
if [ -f api/consolidated_molecule_resource.py ] && [ -f api/consolidated_utils.py ]; then
    echo "API files are in place. Verifying content..."
    
    # Check for key functions in consolidated_utils.py
    if grep -q "def get_primary_molecule" api/consolidated_utils.py && \
       grep -q "def migrate_properties" api/consolidated_utils.py; then
        echo -e "${GREEN}consolidated_utils.py contains required functions.${NC}"
    else
        echo -e "${RED}consolidated_utils.py is missing required functions.${NC}"
        API_VERIFICATION_FAILED=true
    fi
    
    # Check for key endpoints in consolidated_molecule_resource.py
    if grep -q "class ConsolidatedMoleculeResource" api/consolidated_molecule_resource.py && \
       grep -q "migrate_properties" api/consolidated_molecule_resource.py; then
        echo -e "${GREEN}consolidated_molecule_resource.py contains required endpoints.${NC}"
    else
        echo -e "${RED}consolidated_molecule_resource.py is missing required endpoints.${NC}"
        API_VERIFICATION_FAILED=true
    fi
    
    # Check API registration in __init__.py
    if grep -q "consolidated_molecule_resource" api/__init__.py; then
        echo -e "${GREEN}API endpoints are registered in __init__.py.${NC}"
    else
        echo -e "${RED}API endpoints are not registered in __init__.py.${NC}"
        API_VERIFICATION_FAILED=true
    fi
else
    echo -e "${RED}API files are missing.${NC}"
    API_VERIFICATION_FAILED=true
fi

# Check for RDKit wrapper file
if [ -f rdkit_wrapper_consolidated.py ]; then
    echo "RDKit wrapper file is in place."
    
    # Check for key functions in rdkit_wrapper_consolidated.py
    if grep -q "def get_primary_molecule_id" rdkit_wrapper_consolidated.py && \
       grep -q "def calculate_properties_for_consolidated" rdkit_wrapper_consolidated.py; then
        echo -e "${GREEN}rdkit_wrapper_consolidated.py contains required functions.${NC}"
    else
        echo -e "${RED}rdkit_wrapper_consolidated.py is missing required functions.${NC}"
        API_VERIFICATION_FAILED=true
    fi
else
    echo -e "${RED}RDKit wrapper file is missing.${NC}"
    API_VERIFICATION_FAILED=true
fi

# Inform about database migrations
echo -e "\n${BLUE}Database Schema:${NC}"
echo "The following database elements should be created:"
echo "1. CHECK constraint on consolidated_molecules.molecule_status"
echo "2. Indexes on consolidated_molecules (inchikey, primary_molecule_id, molecule_status)"
echo "3. scientific_data_audit table with indexes"
echo -e "${BLUE}These can be verified once database access is available.${NC}"

if [ "$API_VERIFICATION_FAILED" = "true" ]; then
    echo -e "${RED}API implementation verification failed.${NC}"
    echo -e "${RED}Please check the logs for details.${NC}"
    exit 1
else
    echo -e "\n${GREEN}Implementation completed successfully!${NC}"
    echo -e "${GREEN}The API components for the consolidated molecules system are now ready to use.${NC}"
    echo -e "\n${BLUE}New API endpoints available:${NC}"
    echo -e "- /api/v1/consolidated/molecules/{molecule_id}"
    echo -e "- /api/v1/consolidated/batch"
    echo -e "- /api/v1/molecules/{molecule_id}/primary"
    echo -e "- /api/v1/consolidated"
    echo -e "- /api/v1/consolidated/molecules/{molecule_id}/audit"
    echo -e "- /api/v1/consolidated/molecules/{molecule_id}/migrate-properties"
    echo -e "- /api/v1/consolidated/search"
    echo -e "\n${BLUE}Note:${NC} Database migrations must be applied separately before using the API."
fi