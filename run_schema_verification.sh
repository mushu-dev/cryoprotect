#!/bin/bash

# Run Database Schema Verification Tests
# This script runs the database schema verification tests

# Exit on error
set -e

# Print header
echo "==================================================="
echo "Database Schema Verification Tests"
echo "==================================================="

# Set up environment
echo -e "\nSetting up environment variables..."
if [ -f .env ]; then
    echo "Loading .env file"
    export $(grep -v '^#' .env | xargs)
else
    echo "No .env file found, please make sure database connection parameters are set in the environment"
    read -p "Do you want to continue? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]
    then
        exit 1
    fi
fi

# Check if required variables exist
if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_NAME" ] || [ -z "$SUPABASE_DB_USER" ] || [ -z "$SUPABASE_DB_PASSWORD" ]; then
    echo "Error: Required database connection variables are not set"
    echo "Please set SUPABASE_DB_HOST, SUPABASE_DB_NAME, SUPABASE_DB_USER, and SUPABASE_DB_PASSWORD"
    exit 1
fi

# Create Python virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo -e "\nCreating virtual environment..."
    python -m venv venv
fi

# Activate virtual environment
echo -e "\nActivating virtual environment..."
source venv/bin/activate

# Install required packages
echo -e "\nInstalling required packages..."
pip install -r requirements.txt
pip install python-dotenv psycopg2-binary

# Create reports directory if it doesn't exist
mkdir -p reports

# Timestamp for reports
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Step 1: Check database connection
echo -e "\n==================================================="
echo "Step 1: Checking database connection"
echo "==================================================="
echo "Testing connection to Supabase database..."
python test_db_connection.py
if [ $? -ne 0 ]; then
    echo "Database connection test failed. Please fix the connection issues before continuing."
    exit 1
fi
echo "Database connection test passed!"

# Step 2: Generate schema snapshot
echo -e "\n==================================================="
echo "Step 2: Generating database schema snapshot"
echo "==================================================="
echo "Creating snapshot of current database schema..."

# Create schema snapshot script
cat > generate_schema_snapshot.py << EOF
"""
Generate Database Schema Snapshot
"""
import os
import sys
import json
import datetime
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Get database connection parameters
DB_CONNECTION_PARAMS = {
    'host': os.environ.get('SUPABASE_DB_HOST'),
    'port': os.environ.get('SUPABASE_DB_PORT', '5432'),
    'dbname': os.environ.get('SUPABASE_DB_NAME'),
    'user': os.environ.get('SUPABASE_DB_USER'),
    'password': os.environ.get('SUPABASE_DB_PASSWORD'),
    'sslmode': os.environ.get('SUPABASE_DB_SSLMODE', 'require')
}

def generate_schema_snapshot():
    """Generate a snapshot of the current database schema"""
    conn = None
    try:
        # Connect to database
        conn = psycopg2.connect(**DB_CONNECTION_PARAMS)
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        
        # Get all tables
        cursor.execute("""
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = 'public'
            AND table_type = 'BASE TABLE'
            ORDER BY table_name
        """)
        tables = [table['table_name'] for table in cursor.fetchall()]
        
        # Create schema snapshot
        schema = {
            'database': DB_CONNECTION_PARAMS['dbname'],
            'timestamp': datetime.datetime.now().isoformat(),
            'tables': {}
        }
        
        # Get columns for each table
        for table in tables:
            cursor.execute(f"""
                SELECT column_name, data_type, is_nullable, column_default
                FROM information_schema.columns
                WHERE table_name = '{table}'
                AND table_schema = 'public'
                ORDER BY ordinal_position
            """)
            schema['tables'][table] = {'columns': cursor.fetchall()}
            
            # Get constraints for table
            cursor.execute(f"""
                SELECT
                    tc.constraint_name,
                    tc.constraint_type,
                    kcu.column_name,
                    ccu.table_name AS foreign_table_name,
                    ccu.column_name AS foreign_column_name
                FROM information_schema.table_constraints tc
                LEFT JOIN information_schema.key_column_usage kcu
                    ON tc.constraint_name = kcu.constraint_name
                    AND tc.table_schema = kcu.table_schema
                LEFT JOIN information_schema.constraint_column_usage ccu
                    ON ccu.constraint_name = tc.constraint_name
                    AND ccu.table_schema = tc.table_schema
                WHERE tc.table_name = '{table}'
                AND tc.table_schema = 'public'
                ORDER BY tc.constraint_name, kcu.column_name
            """)
            schema['tables'][table]['constraints'] = cursor.fetchall()
            
            # Get indexes for table
            cursor.execute(f"""
                SELECT
                    i.relname AS index_name,
                    array_agg(a.attname) AS column_names,
                    ix.indisunique AS is_unique,
                    ix.indisprimary AS is_primary
                FROM pg_class t
                JOIN pg_index ix ON t.oid = ix.indrelid
                JOIN pg_class i ON i.oid = ix.indexrelid
                JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(ix.indkey)
                JOIN pg_namespace n ON n.oid = t.relnamespace
                WHERE t.relkind = 'r'
                AND n.nspname = 'public'
                AND t.relname = '{table}'
                GROUP BY i.relname, ix.indisunique, ix.indisprimary
                ORDER BY i.relname
            """)
            schema['tables'][table]['indexes'] = cursor.fetchall()
        
        # Add RLS policies
        cursor.execute("""
            SELECT
                schemaname,
                tablename,
                policyname,
                roles,
                cmd,
                qual,
                with_check
            FROM pg_policies
            WHERE schemaname = 'public'
            ORDER BY tablename, policyname
        """)
        schema['rls_policies'] = cursor.fetchall()
        
        # Add security definer functions
        cursor.execute("""
            SELECT
                p.proname AS function_name,
                l.lanname AS language,
                p.prosecdef AS security_definer
            FROM pg_proc p
            JOIN pg_namespace n ON p.pronamespace = n.oid
            JOIN pg_language l ON p.prolang = l.oid
            WHERE n.nspname = 'public'
            AND p.prosecdef = true
            ORDER BY p.proname
        """)
        schema['security_definer_functions'] = cursor.fetchall()
        
        # Save schema snapshot to file
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        os.makedirs('reports', exist_ok=True)
        
        snapshot_file = f'reports/database_schema_snapshot_{timestamp}.json'
        with open(snapshot_file, 'w') as f:
            json.dump(schema, f, indent=2, default=str)
        
        print(f"Schema snapshot saved to {snapshot_file}")
        return snapshot_file
    
    except Exception as e:
        print(f"Error generating schema snapshot: {e}")
        return None
    
    finally:
        if conn:
            conn.close()

if __name__ == '__main__':
    generate_schema_snapshot()
EOF

# Run schema snapshot script
python generate_schema_snapshot.py
SNAPSHOT_RESULT=$?
if [ $SNAPSHOT_RESULT -ne 0 ]; then
    echo "Error generating schema snapshot. Continuing with tests anyway."
fi

# Step 3: Run schema verification tests
echo -e "\n==================================================="
echo "Step 3: Running database schema verification tests"
echo "==================================================="
echo "Executing schema verification tests..."
python -m tests.test_database_schema_verification
TEST_RESULT=$?
echo "Schema verification tests completed with exit code: $TEST_RESULT"

# Generate comprehensive report
echo -e "\n==================================================="
echo "Generating comprehensive schema verification report"
echo "==================================================="

# Create a markdown report
REPORT_FILE="reports/schema_verification_report_${TIMESTAMP}.md"
cat > "$REPORT_FILE" << EOF
# Database Schema Verification Report

Date: $(date)
Database: $SUPABASE_DB_NAME
Host: $SUPABASE_DB_HOST

## Summary

This report summarizes the results of the database schema verification tests.

| Test Category | Status | Notes |
|---------------|--------|-------|
| Table Structure | $([ $TEST_RESULT -eq 0 ] && echo "✅ Passed" || echo "❌ Failed") | See test logs for details |
| Foreign Key Relationships | $([ $TEST_RESULT -eq 0 ] && echo "✅ Passed" || echo "❌ Failed") | See test logs for details |
| Index Coverage | $([ $TEST_RESULT -eq 0 ] && echo "✅ Passed" || echo "❌ Failed") | See test logs for details |

## Detailed Test Logs

### Database Schema Tests

\`\`\`
$(cat database_schema_tests.log 2>/dev/null || echo "Log file not found")
\`\`\`

## Recommendations

Based on the test results, the following recommendations are made:

1. $([ $TEST_RESULT -eq 0 ] && echo "The database schema is correctly implemented according to requirements." || echo "Fix issues with the database schema as indicated in the test logs.")
2. Maintain proper foreign key relationships for all tables to ensure data integrity.
3. Ensure indexes exist for all columns used in RLS policies to optimize performance.
4. Continue to verify database schema after any migrations to ensure integrity.

## Next Steps

1. Run performance tests against the database with the current schema.
2. Monitor query performance in production to identify any needed additional indexes.
3. Document the database schema for future reference.
EOF

echo "Comprehensive report generated: $REPORT_FILE"

# Final status
if [ $TEST_RESULT -eq 0 ]; then
    echo -e "\nDatabase schema verification tests passed successfully! 🎉"
    exit 0
else
    echo -e "\nDatabase schema verification tests failed. Please check the logs and fix any issues."
    exit 1
fi