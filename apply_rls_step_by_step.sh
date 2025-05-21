#!/bin/bash
# apply_rls_step_by_step.sh
# Script to apply RLS policy improvements step by step

# Ensure script execution stops if any command fails
set -e

# Ensure we're in the project root directory
cd "$(dirname "$0")"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "=== RLS Policy Improvements Step-by-Step Applier ==="
echo "This script will apply improved RLS policies based on verification findings."
echo ""

# Source .env file if it exists
if [ -f .env ]; then
    echo "Loading environment variables from .env file"
    export $(grep -v '^#' .env | xargs)
fi

# Check for required environment variables
if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_USER" ] || [ -z "$SUPABASE_DB_PASSWORD" ]; then
    echo "ERROR: Missing required environment variables"
    echo "Please set the following environment variables or create a .env file:"
    echo "  SUPABASE_DB_HOST"
    echo "  SUPABASE_DB_PORT (defaults to 5432)"
    echo "  SUPABASE_DB_NAME (defaults to postgres)"
    echo "  SUPABASE_DB_USER"
    echo "  SUPABASE_DB_PASSWORD"
    echo "  SUPABASE_DB_SSLMODE (defaults to require)"
    exit 1
fi

echo "Environment variables set correctly"
echo "Applying RLS policy improvements step by step..."

# Ensure the migrations directory exists
mkdir -p "migrations/rls_helpers"

# Apply each SQL file in the rls_helpers directory
SQL_FILES_DIR="migrations/rls_helpers"
for sql_file in "$SQL_FILES_DIR"/*.sql; do
    echo "Applying $sql_file"
    python apply_file_sql.py "$sql_file"
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to apply $sql_file"
        exit 1
    fi
    echo "Successfully applied $sql_file"
    echo ""
done

# Record the migration
echo "Recording the migration in the migrations table"
python - << EOF
import psycopg2
import os
from datetime import datetime

# Get database connection parameters from environment variables
db_params = {
    "host": os.environ.get("SUPABASE_DB_HOST"),
    "port": os.environ.get("SUPABASE_DB_PORT", "5432"),
    "database": os.environ.get("SUPABASE_DB_NAME", "postgres"),
    "user": os.environ.get("SUPABASE_DB_USER"),
    "password": os.environ.get("SUPABASE_DB_PASSWORD"),
    "sslmode": os.environ.get("SUPABASE_DB_SSLMODE", "require")
}

try:
    # Create a connection to the database
    conn = psycopg2.connect(**db_params)
    conn.autocommit = True
    cursor = conn.cursor()
    
    # Create migrations table if it doesn't exist
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS migrations (
            id SERIAL PRIMARY KEY,
            name VARCHAR(255) NOT NULL,
            applied_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        );
    """)
    
    # Record the migration
    cursor.execute("""
        INSERT INTO migrations (name, applied_at)
        VALUES (%s, NOW());
    """, ('020_rls_improvements.sql',))
    
    print("Migration recorded successfully")
    
    cursor.close()
    conn.close()
except Exception as e:
    print(f"ERROR: Failed to record migration: {e}")
    exit(1)
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to record migration"
    exit 1
fi

echo ""
echo "RLS policy improvements applied successfully!"
echo "You can now run the verification tests to confirm the changes using:"
echo "  python run_rls_verification.py"