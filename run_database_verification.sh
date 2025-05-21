#!/bin/bash

# Database Verification Suite
# This script orchestrates the execution of all database verification tests

# Exit on error
set -e

# Print header
echo "==================================================="
echo "CryoProtect Database Verification Suite"
echo "==================================================="
echo "Running a comprehensive series of tests to verify:"
echo "  1. RLS Policy Implementation"
echo "  2. Connection Pool Performance"
echo "  3. Database Schema Integrity"
echo "  4. Migration Framework"
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

# Create logs directory if it doesn't exist
mkdir -p logs

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

# Step 2: Analyze current RLS implementation
echo -e "\n==================================================="
echo "Step 2: Analyzing RLS implementation"
echo "==================================================="
echo "Generating RLS verification report..."
python verify_optimized_rls.py
if [ $? -ne 0 ]; then
    echo "Warning: RLS analysis encountered issues. Proceeding with tests anyway."
fi
echo "RLS analysis completed."

# Step 3: Run RLS Policy Verification Tests
echo -e "\n==================================================="
echo "Step 3: Running RLS Policy Verification Tests"
echo "==================================================="
echo "Executing RLS policy tests..."
python -m tests.test_rls_policies_verification
TEST_RESULT_RLS=$?
echo "RLS policy tests completed with exit code: $TEST_RESULT_RLS"

# Step 4: Run Connection Pool Verification Tests
echo -e "\n==================================================="
echo "Step 4: Running Connection Pool Verification Tests"
echo "==================================================="
echo "Executing connection pool tests..."
python -m tests.test_connection_pool_verification
TEST_RESULT_POOL=$?
echo "Connection pool tests completed with exit code: $TEST_RESULT_POOL"

# Step 5: Run Database Schema Verification (if available)
echo -e "\n==================================================="
echo "Step 5: Running Database Schema Verification"
echo "==================================================="
if [ -f "tests/test_database_schema_verification.py" ]; then
    echo "Executing database schema tests..."
    python -m tests.test_database_schema_verification
    TEST_RESULT_SCHEMA=$?
    echo "Database schema tests completed with exit code: $TEST_RESULT_SCHEMA"
else
    echo "Database schema verification tests not found. Creating validation script..."
    
    cat > tests/test_database_schema_verification.py << EOF
"""
Database Schema Verification Tests
Based on the test cases defined in DATABASE_VERIFICATION_PLAN.md
"""
import os
import sys
import unittest
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Configure logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename='database_schema_tests.log'
)
logger = logging.getLogger(__name__)

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

class DatabaseSchemaTests(unittest.TestCase):
    """Tests for Database Schema (TC-DB-001 to TC-DB-005)"""
    
    @classmethod
    def setUpClass(cls):
        """Set up database connection"""
        cls.conn = psycopg2.connect(**DB_CONNECTION_PARAMS)
    
    @classmethod
    def tearDownClass(cls):
        """Close database connection"""
        cls.conn.close()
    
    def test_molecule_table_schema(self):
        """TC-DB-001: Verify molecule table schema"""
        cursor = self.conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute("""
            SELECT column_name, data_type, is_nullable
            FROM information_schema.columns
            WHERE table_name = 'molecule'
            ORDER BY ordinal_position
        """)
        columns = cursor.fetchall()
        
        # Verify essential columns
        essential_columns = ['id', 'project_id', 'name', 'smiles']
        for column in essential_columns:
            self.assertTrue(
                any(col['column_name'] == column for col in columns),
                f"Column {column} not found in molecule table"
            )
        
        logger.info("Molecule table schema verified")
    
    def test_foreign_key_relationships(self):
        """TC-DB-010 to TC-DB-015: Verify foreign key relationships"""
        cursor = self.conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute("""
            SELECT
                tc.table_name,
                kcu.column_name,
                ccu.table_name AS foreign_table_name,
                ccu.column_name AS foreign_column_name
            FROM information_schema.table_constraints tc
            JOIN information_schema.key_column_usage kcu
                ON tc.constraint_name = kcu.constraint_name
            JOIN information_schema.constraint_column_usage ccu
                ON ccu.constraint_name = tc.constraint_name
            WHERE tc.constraint_type = 'FOREIGN KEY'
            AND tc.table_name IN ('molecule', 'mixture', 'mixture_component',
                                 'molecular_property', 'experiment', 'experiment_property')
        """)
        foreign_keys = cursor.fetchall()
        
        # Expected relationships
        expected_relationships = [
            {'table': 'molecule', 'column': 'project_id', 'foreign_table': 'project', 'foreign_column': 'id'},
            {'table': 'mixture', 'column': 'project_id', 'foreign_table': 'project', 'foreign_column': 'id'},
            {'table': 'mixture_component', 'column': 'mixture_id', 'foreign_table': 'mixture', 'foreign_column': 'id'},
            {'table': 'molecular_property', 'column': 'molecule_id', 'foreign_table': 'molecule', 'foreign_column': 'id'},
            {'table': 'experiment', 'column': 'project_id', 'foreign_table': 'project', 'foreign_column': 'id'},
            {'table': 'experiment_property', 'column': 'experiment_id', 'foreign_table': 'experiment', 'foreign_column': 'id'}
        ]
        
        for expected in expected_relationships:
            self.assertTrue(
                any(
                    fk['table_name'] == expected['table'] and
                    fk['column_name'] == expected['column'] and
                    fk['foreign_table_name'] == expected['foreign_table'] and
                    fk['foreign_column_name'] == expected['foreign_column']
                    for fk in foreign_keys
                ),
                f"Foreign key relationship not found: {expected['table']}.{expected['column']} -> {expected['foreign_table']}.{expected['foreign_column']}"
            )
        
        logger.info("Foreign key relationships verified")
    
    def test_rls_policy_indexes(self):
        """TC-DB-020: Verify indexes for RLS policy columns"""
        cursor = self.conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute("""
            SELECT
                t.relname AS table_name,
                i.relname AS index_name,
                a.attname AS column_name
            FROM pg_class t
            JOIN pg_index ix ON t.oid = ix.indrelid
            JOIN pg_class i ON i.oid = ix.indexrelid
            JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(ix.indkey)
            JOIN pg_namespace n ON n.oid = t.relnamespace
            WHERE t.relkind = 'r'
            AND n.nspname = 'public'
            AND a.attname IN ('project_id', 'user_id')
            ORDER BY t.relname, i.relname
        """)
        indexes = cursor.fetchall()
        
        # Expected indexes for RLS policy columns
        expected_rls_indexes = [
            {'table': 'molecule', 'column': 'project_id'},
            {'table': 'mixture', 'column': 'project_id'},
            {'table': 'experiment', 'column': 'project_id'},
            {'table': 'user_profile', 'column': 'user_id'},
            {'table': 'user_profile', 'column': 'project_id'}
        ]
        
        for expected in expected_rls_indexes:
            self.assertTrue(
                any(
                    idx['table_name'] == expected['table'] and
                    idx['column_name'] == expected['column']
                    for idx in indexes
                ),
                f"RLS index not found: {expected['table']}.{expected['column']}"
            )
        
        logger.info("RLS policy indexes verified")

if __name__ == '__main__':
    unittest.main()
EOF
    
    echo "Running newly created database schema tests..."
    python -m tests.test_database_schema_verification
    TEST_RESULT_SCHEMA=$?
    echo "Database schema tests completed with exit code: $TEST_RESULT_SCHEMA"
fi

# Step 6: Generate comprehensive report
echo -e "\n==================================================="
echo "Step 6: Generating Comprehensive Verification Report"
echo "==================================================="

# Create a markdown report
REPORT_FILE="reports/database_verification_${TIMESTAMP}.md"
cat > "$REPORT_FILE" << EOF
# Database Verification Report

Date: $(date)
Database: $SUPABASE_DB_NAME
Host: $SUPABASE_DB_HOST

## Summary

This report summarizes the results of the database verification tests.

| Test Type | Status | Notes |
|-----------|--------|-------|
| Database Connection | ✅ Passed | Connection to Supabase database successful |
| RLS Implementation | $([ $TEST_RESULT_RLS -eq 0 ] && echo "✅ Passed" || echo "❌ Failed") | See test logs for details |
| Connection Pool | $([ $TEST_RESULT_POOL -eq 0 ] && echo "✅ Passed" || echo "❌ Failed") | See test logs for details |
| Database Schema | $([ $TEST_RESULT_SCHEMA -eq 0 ] && echo "✅ Passed" || echo "❌ Failed") | See test logs for details |

## Detailed Test Results

### RLS Policy Verification

$(cat rls_verification_tests.log 2>/dev/null || echo "Log file not found")

### Connection Pool Verification

$(cat connection_pool_tests.log 2>/dev/null || echo "Log file not found")

### Database Schema Verification

$(cat database_schema_tests.log 2>/dev/null || echo "Log file not found")

## Recommendations

Based on the test results, the following recommendations are made:

1. $([ $TEST_RESULT_RLS -eq 0 ] && echo "No action needed for RLS policies" || echo "Fix issues with RLS policy implementation")
2. $([ $TEST_RESULT_POOL -eq 0 ] && echo "No action needed for connection pool" || echo "Optimize connection pool implementation")
3. $([ $TEST_RESULT_SCHEMA -eq 0 ] && echo "No action needed for database schema" || echo "Fix issues with database schema")

## Next Steps

1. Address any failed tests and rerun the verification suite
2. Implement performance optimizations based on the test results
3. Document the verified database architecture
EOF

echo "Comprehensive report generated: $REPORT_FILE"

# Final summary
echo -e "\n==================================================="
echo "Database Verification Suite Completed"
echo "==================================================="
echo "Results Summary:"
echo "  - RLS Policy Tests: $([ $TEST_RESULT_RLS -eq 0 ] && echo "PASSED" || echo "FAILED")"
echo "  - Connection Pool Tests: $([ $TEST_RESULT_POOL -eq 0 ] && echo "PASSED" || echo "FAILED")"
echo "  - Database Schema Tests: $([ $TEST_RESULT_SCHEMA -eq 0 ] && echo "PASSED" || echo "FAILED")"
echo -e "\nDetailed report: $REPORT_FILE"
echo "==================================================="

# Calculate overall status
if [ $TEST_RESULT_RLS -eq 0 ] && [ $TEST_RESULT_POOL -eq 0 ] && [ $TEST_RESULT_SCHEMA -eq 0 ]; then
    echo "All tests passed successfully! 🎉"
    exit 0
else
    echo "Some tests failed. Please address the issues before proceeding."
    exit 1
fi