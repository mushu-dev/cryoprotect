#!/bin/bash
# apply_rls_policy_improvements.sh
# Script to apply RLS policy improvements to the database

# Ensure script execution stops if any command fails
set -e

# Ensure we're in the project root directory
cd "$(dirname "$0")"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "=== RLS Policy Improvements Application Tool ==="
echo "This script will apply improved RLS policies based on verification findings."
echo "Applying the improved team-based RLS policies will enhance security and performance."
echo ""

# Check for command line arguments
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --dry-run     Check what would be applied without making changes"
    echo "  --mcp         Execute using the Supabase MCP tool"
    echo "  --help, -h    Show this help message"
    exit 0
fi

# Check if this is a dry run
if [ "$1" == "--dry-run" ]; then
    echo "DRY RUN MODE - No changes will be made to the database"
    python apply_rls_policy_improvements.py --dry-run
    exit 0
fi

# Check if we should use the MCP tool
if [ "$1" == "--mcp" ]; then
    echo "Using Supabase MCP tool to apply migration..."
    
    if ! command -v python3 &> /dev/null; then
        echo "Error: python3 is required but could not be found"
        exit 1
    fi
    
    # Create a temporary Python script to use the MCP tool
    TEMP_SCRIPT=$(mktemp)
    cat << 'EOF' > "$TEMP_SCRIPT"
#!/usr/bin/env python3
"""
Script to apply RLS policy improvements using the Supabase MCP tool.
"""
import os
import sys
import json
from datetime import datetime

# Add necessary paths
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

try:
    # Import MCP-specific modules if available
    from supabase_mcp_tools import execute_sql_with_mcp
    
    # Read the migration file
    migration_path = "migrations/020_rls_policy_improvements.sql"
    with open(migration_path, 'r') as file:
        sql_content = file.read()
    
    # Get the MCP project ID
    project_id = os.environ.get("SUPABASE_PROJECT_ID")
    if not project_id:
        print("ERROR: SUPABASE_PROJECT_ID environment variable is required")
        sys.exit(1)
    
    # Execute the SQL using MCP
    print(f"Executing SQL with MCP for project {project_id}...")
    response = execute_sql_with_mcp(project_id, sql_content)
    
    # Process the response
    if response.get("status") == "success":
        print("Migration applied successfully using MCP")
        
        # Record the migration
        record_sql = """
            CREATE TABLE IF NOT EXISTS migrations (
                id SERIAL PRIMARY KEY,
                name VARCHAR(255) NOT NULL,
                applied_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
            );
            
            INSERT INTO migrations (name, applied_at)
            VALUES ('020_rls_policy_improvements.sql', NOW());
        """
        
        record_response = execute_sql_with_mcp(project_id, record_sql)
        
        if record_response.get("status") == "success":
            print("Migration recorded successfully")
        else:
            print(f"WARNING: Migration applied but failed to record: {record_response.get('error')}")
            
    else:
        print(f"ERROR: Failed to apply migration: {response.get('error')}")
        sys.exit(1)
        
except ImportError:
    print("ERROR: Could not import MCP tools. Make sure the MCP tooling is properly set up.")
    sys.exit(1)
except Exception as e:
    print(f"ERROR: An unexpected error occurred: {str(e)}")
    sys.exit(1)
EOF
    
    # Make the script executable and run it
    chmod +x "$TEMP_SCRIPT"
    python3 "$TEMP_SCRIPT"
    
    # Clean up
    rm "$TEMP_SCRIPT"
    
    exit $?
fi

# Default case: Apply directly using the Python script
echo "Checking environment variables..."

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
echo "Applying RLS policy improvements..."

# Ensure the Python script is executable
chmod +x apply_rls_policy_improvements.py

# Run the Python script
python apply_rls_policy_improvements.py

# Check for success
if [ $? -eq 0 ]; then
    echo "RLS policy improvements applied successfully!"
    echo "You can now run the verification tests to confirm the changes using:"
    echo "  python run_rls_verification.py"
else
    echo "ERROR: Failed to apply RLS policy improvements"
    exit 1
fi