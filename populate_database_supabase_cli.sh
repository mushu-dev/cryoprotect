#!/bin/bash
# CryoProtect v2 - Database Population Script using Supabase CLI
# This script populates the CryoProtect database with initial data using the Supabase CLI

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

echo "Starting database population using Supabase CLI..."

# Function to execute SQL script
execute_sql_script() {
    local script_name=$1
    echo "Executing $script_name..."
    
    # Read SQL file
    SQL_CONTENT=$(cat $script_name)
    
    # Remove transaction commands and psql-specific commands
    SQL_CONTENT=$(echo "$SQL_CONTENT" | grep -v "^BEGIN;" | grep -v "^COMMIT;" | grep -v "^\\\\echo" | grep -v "^\\\\i")
    
    # Execute SQL using Supabase CLI
    echo "$SQL_CONTENT" | supabase db execute --project-ref $PROJECT_ID
    
    if [ $? -ne 0 ]; then
        echo "Error executing $script_name. Exiting."
        exit 1
    fi
    
    echo "$script_name executed successfully."
    echo ""
}

# Execute scripts in order
execute_sql_script "populate_molecules.sql"
execute_sql_script "populate_molecular_properties.sql"
execute_sql_script "populate_mixtures.sql"
execute_sql_script "populate_mixture_components.sql"
execute_sql_script "populate_experiments.sql"
execute_sql_script "populate_experiment_properties.sql"
execute_sql_script "populate_predictions.sql"

echo "Database population completed successfully!"

# Verify data was inserted correctly
echo "Verifying data was inserted correctly..."

VERIFICATION_SQL="
SELECT 'molecules' AS table_name, COUNT(*) AS record_count FROM molecules
UNION ALL
SELECT 'molecular_properties' AS table_name, COUNT(*) AS record_count FROM molecular_properties
UNION ALL
SELECT 'mixtures' AS table_name, COUNT(*) AS record_count FROM mixtures
UNION ALL
SELECT 'mixture_components' AS table_name, COUNT(*) AS record_count FROM mixture_components
UNION ALL
SELECT 'experiments' AS table_name, COUNT(*) AS record_count FROM experiments
UNION ALL
SELECT 'predictions' AS table_name, COUNT(*) AS record_count FROM predictions
ORDER BY table_name;
"

echo "$VERIFICATION_SQL" | supabase db execute --project-ref $PROJECT_ID

echo "Database population verification complete!"