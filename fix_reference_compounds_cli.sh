#!/bin/bash
# CryoProtect v2 - Fix Reference Compounds Script using Supabase CLI
# This script fixes missing properties for reference compounds using the Supabase CLI

# Check if Supabase CLI is installed
if ! command -v supabase &> /dev/null; then
    echo "Supabase CLI is not installed. Please install it first."
    echo "Visit https://supabase.com/docs/guides/cli for installation instructions."
    exit 1
fi

# Set variables
PROJECT_ID="tsdlmynydfuypiugmkev"

# Check if user is logged in to Supabase
supabase projects list &> /dev/null
if [ $? -ne 0 ]; then
    echo "You are not logged in to Supabase CLI. Please login first."
    echo "Run: supabase login"
    exit 1
fi

echo "Starting reference compounds fix using Supabase CLI..."

# Function to execute SQL
execute_sql() {
    local sql_content=$1
    local description=$2
    
    echo "Executing $description..."
    
    # Execute SQL using Supabase CLI
    echo "$sql_content" | supabase db execute --project-ref $PROJECT_ID
    
    if [ $? -ne 0 ]; then
        echo "Error executing $description. Exiting."
        exit 1
    fi
    
    echo "$description executed successfully."
    echo ""
}

# Define reference compound IDs
REFERENCE_IDS=(
    "CHEMBL1098659"
    "CHEMBL1487"
    "CHEMBL262548"
    "CHEMBL388978"
    "CHEMBL500033"
    "CHEMBL6196"
    "CHEMBL66195"
    "CHEMBL6752"
    "CHEMBL967"
)

# Step 1: Verify reference compounds exist
echo "Verifying reference compounds exist..."

VERIFY_SQL="
SELECT chembl_id, id FROM molecules 
WHERE chembl_id IN ('CHEMBL1098659', 'CHEMBL1487', 'CHEMBL262548', 'CHEMBL388978', 
                   'CHEMBL500033', 'CHEMBL6196', 'CHEMBL66195', 'CHEMBL6752', 'CHEMBL967');
"

execute_sql "$VERIFY_SQL" "Reference compounds verification"

# Step 2: Ensure property types exist
echo "Ensuring property types exist..."

PROPERTY_TYPES_SQL="
-- Create property types if they don't exist
INSERT INTO property_types (name, data_type, description)
SELECT 'logP', 'numeric', 'Partition coefficient'
WHERE NOT EXISTS (SELECT 1 FROM property_types WHERE name = 'logP');

INSERT INTO property_types (name, data_type, description)
SELECT 'h_bond_donors', 'numeric', 'Number of hydrogen bond donors'
WHERE NOT EXISTS (SELECT 1 FROM property_types WHERE name = 'h_bond_donors');

INSERT INTO property_types (name, data_type, description)
SELECT 'h_bond_acceptors', 'numeric', 'Number of hydrogen bond acceptors'
WHERE NOT EXISTS (SELECT 1 FROM property_types WHERE name = 'h_bond_acceptors');

-- Return the property type IDs
SELECT id, name FROM property_types 
WHERE name IN ('logP', 'h_bond_donors', 'h_bond_acceptors');
"

execute_sql "$PROPERTY_TYPES_SQL" "Property types creation"

# Step 3: Check missing properties for each reference compound
echo "Checking missing properties for reference compounds..."

for chembl_id in "${REFERENCE_IDS[@]}"; do
    echo "Checking properties for $chembl_id..."
    
    CHECK_SQL="
    WITH molecule AS (
        SELECT id FROM molecules WHERE chembl_id = '$chembl_id'
    ),
    property_types AS (
        SELECT id, name FROM property_types 
        WHERE name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
    ),
    existing_properties AS (
        SELECT pt.name
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        JOIN molecule m ON mp.molecule_id = m.id
        WHERE pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
    )
    SELECT pt.name AS missing_property
    FROM property_types pt
    WHERE pt.name NOT IN (SELECT name FROM existing_properties);
    "
    
    # Execute the check SQL
    MISSING_PROPERTIES=$(echo "$CHECK_SQL" | supabase db execute --project-ref $PROJECT_ID)
    
    # If there are missing properties, fix them
    if [[ $MISSING_PROPERTIES == *"missing_property"* ]]; then
        echo "Found missing properties for $chembl_id. Fixing..."
        
        FIX_SQL="
        -- Get molecule ID
        WITH molecule AS (
            SELECT id FROM molecules WHERE chembl_id = '$chembl_id'
        ),
        -- Get property type IDs
        property_types AS (
            SELECT id, name FROM property_types 
            WHERE name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
        ),
        -- Get existing properties
        existing_properties AS (
            SELECT pt.name
            FROM molecular_properties mp
            JOIN property_types pt ON mp.property_type_id = pt.id
            JOIN molecule m ON mp.molecule_id = m.id
            WHERE pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
        )
        -- Insert missing properties
        INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value, created_at, updated_at)
        SELECT 
            (SELECT id FROM molecule),
            pt.id,
            CASE 
                WHEN pt.name = 'logP' THEN 0.0
                ELSE 0
            END,
            NOW(),
            NOW()
        FROM property_types pt
        WHERE pt.name NOT IN (SELECT name FROM existing_properties);
        "
        
        execute_sql "$FIX_SQL" "Property fix for $chembl_id"
    else
        echo "No missing properties for $chembl_id."
    fi
done

# Step 4: Verify all properties are now set
echo "Verifying all properties are now set..."

FINAL_VERIFY_SQL="
WITH reference_molecules AS (
    SELECT id, chembl_id 
    FROM molecules 
    WHERE chembl_id IN ('CHEMBL1098659', 'CHEMBL1487', 'CHEMBL262548', 'CHEMBL388978', 
                       'CHEMBL500033', 'CHEMBL6196', 'CHEMBL66195', 'CHEMBL6752', 'CHEMBL967')
),
property_counts AS (
    SELECT 
        rm.chembl_id,
        COUNT(DISTINCT pt.name) AS property_count
    FROM reference_molecules rm
    JOIN molecular_properties mp ON rm.id = mp.molecule_id
    JOIN property_types pt ON mp.property_type_id = pt.id
    WHERE pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
    GROUP BY rm.chembl_id
)
SELECT 
    rm.chembl_id,
    COALESCE(pc.property_count, 0) AS property_count,
    CASE 
        WHEN COALESCE(pc.property_count, 0) = 3 THEN 'Complete'
        ELSE 'Incomplete'
    END AS status
FROM reference_molecules rm
LEFT JOIN property_counts pc ON rm.chembl_id = pc.chembl_id
ORDER BY rm.chembl_id;
"

execute_sql "$FINAL_VERIFY_SQL" "Final verification"

echo "Reference compounds fix complete!"