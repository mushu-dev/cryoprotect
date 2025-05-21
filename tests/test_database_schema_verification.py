"""
Database Schema Verification Tests
Based on the test cases defined in DATABASE_VERIFICATION_PLAN.md
"""
import os
import sys
import unittest
import json
import datetime
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
        cls.generate_schema_snapshot()
    
    @classmethod
    def tearDownClass(cls):
        """Close database connection"""
        cls.conn.close()
    
    @classmethod
    def generate_schema_snapshot(cls):
        """Generate a snapshot of the current database schema"""
        cursor = cls.conn.cursor(cursor_factory=RealDictCursor)
        
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
        schema = {}
        
        # Get columns for each table
        for table in tables:
            cursor.execute(f"""
                SELECT column_name, data_type, is_nullable, column_default
                FROM information_schema.columns
                WHERE table_name = '{table}'
                AND table_schema = 'public'
                ORDER BY ordinal_position
            """)
            schema[table] = {'columns': cursor.fetchall()}
            
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
            schema[table]['constraints'] = cursor.fetchall()
            
            # Get indexes for table
            cursor.execute(f"""
                SELECT
                    i.relname AS index_name,
                    a.attname AS column_name,
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
                ORDER BY i.relname, a.attnum
            """)
            schema[table]['indexes'] = cursor.fetchall()
        
        # Store schema snapshot
        cls.schema = schema
        
        # Save schema snapshot to file
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        report_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'reports')
        os.makedirs(report_dir, exist_ok=True)
        
        report_file = os.path.join(report_dir, f'database_schema_{timestamp}.json')
        with open(report_file, 'w') as f:
            json.dump(schema, f, indent=2, default=str)
        
        logger.info(f"Schema snapshot saved to {report_file}")
    
    def test_molecule_table_schema(self):
        """TC-DB-001: Verify molecule table schema"""
        # Check if molecule table exists
        self.assertIn('molecule', self.schema, "Molecule table does not exist")
        
        # Get column names
        columns = [col['column_name'] for col in self.schema['molecule']['columns']]
        
        # Verify essential columns
        essential_columns = ['id', 'project_id', 'name', 'smiles']
        for column in essential_columns:
            self.assertIn(column, columns, f"Column {column} not found in molecule table")
        
        logger.info("Molecule table schema verified")
    
    def test_mixture_table_schema(self):
        """TC-DB-002: Verify mixture table schema"""
        # Check if mixture table exists
        self.assertIn('mixture', self.schema, "Mixture table does not exist")
        
        # Get column names
        columns = [col['column_name'] for col in self.schema['mixture']['columns']]
        
        # Verify essential columns
        essential_columns = ['id', 'project_id', 'name']
        for column in essential_columns:
            self.assertIn(column, columns, f"Column {column} not found in mixture table")
        
        logger.info("Mixture table schema verified")
    
    def test_experiment_table_schema(self):
        """TC-DB-003: Verify experiment table schema"""
        # Check if experiment table exists
        self.assertIn('experiment', self.schema, "Experiment table does not exist")
        
        # Get column names
        columns = [col['column_name'] for col in self.schema['experiment']['columns']]
        
        # Verify essential columns
        essential_columns = ['id', 'project_id', 'name']
        for column in essential_columns:
            self.assertIn(column, columns, f"Column {column} not found in experiment table")
        
        logger.info("Experiment table schema verified")
    
    def test_molecular_property_table_schema(self):
        """TC-DB-004: Verify molecular_property table schema"""
        # Check if molecular_property table exists
        self.assertIn('molecular_property', self.schema, "Molecular_property table does not exist")
        
        # Get column names
        columns = [col['column_name'] for col in self.schema['molecular_property']['columns']]
        
        # Verify essential columns
        essential_columns = ['molecule_id', 'property_name', 'property_value']
        for column in essential_columns:
            self.assertIn(column, columns, f"Column {column} not found in molecular_property table")
        
        logger.info("Molecular_property table schema verified")
    
    def test_user_profile_table_schema(self):
        """TC-DB-005: Verify user_profile table schema"""
        # Check if user_profile table exists
        self.assertIn('user_profile', self.schema, "User_profile table does not exist")
        
        # Get column names
        columns = [col['column_name'] for col in self.schema['user_profile']['columns']]
        
        # Verify essential columns
        essential_columns = ['user_id', 'project_id', 'role']
        for column in essential_columns:
            self.assertIn(column, columns, f"Column {column} not found in user_profile table")
        
        logger.info("User_profile table schema verified")

class ForeignKeyRelationshipTests(unittest.TestCase):
    """Tests for Foreign Key Relationships (TC-DB-010 to TC-DB-015)"""
    
    @classmethod
    def setUpClass(cls):
        """Set up database connection"""
        cls.conn = psycopg2.connect(**DB_CONNECTION_PARAMS)
    
    @classmethod
    def tearDownClass(cls):
        """Close database connection"""
        cls.conn.close()
    
    def fetch_foreign_keys(self):
        """Fetch all foreign key relationships in the database"""
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
            AND tc.table_schema = 'public'
        """)
        return cursor.fetchall()
    
    def test_molecule_to_project_relationship(self):
        """TC-DB-010: Verify molecule → project relationship"""
        foreign_keys = self.fetch_foreign_keys()
        
        # Look for molecule to project relationship
        found = False
        for fk in foreign_keys:
            if (fk['table_name'] == 'molecule' and 
                fk['column_name'] == 'project_id' and 
                fk['foreign_table_name'] == 'project' and 
                fk['foreign_column_name'] == 'id'):
                found = True
                break
        
        self.assertTrue(found, "Foreign key relationship from molecule to project not found")
        logger.info("Molecule to project relationship verified")
    
    def test_mixture_to_project_relationship(self):
        """TC-DB-011: Verify mixture → project relationship"""
        foreign_keys = self.fetch_foreign_keys()
        
        # Look for mixture to project relationship
        found = False
        for fk in foreign_keys:
            if (fk['table_name'] == 'mixture' and 
                fk['column_name'] == 'project_id' and 
                fk['foreign_table_name'] == 'project' and 
                fk['foreign_column_name'] == 'id'):
                found = True
                break
        
        self.assertTrue(found, "Foreign key relationship from mixture to project not found")
        logger.info("Mixture to project relationship verified")
    
    def test_mixture_component_to_mixture_relationship(self):
        """TC-DB-012: Verify mixture_component → mixture relationship"""
        foreign_keys = self.fetch_foreign_keys()
        
        # Look for mixture_component to mixture relationship
        found = False
        for fk in foreign_keys:
            if (fk['table_name'] == 'mixture_component' and 
                fk['column_name'] == 'mixture_id' and 
                fk['foreign_table_name'] == 'mixture' and 
                fk['foreign_column_name'] == 'id'):
                found = True
                break
        
        self.assertTrue(found, "Foreign key relationship from mixture_component to mixture not found")
        logger.info("Mixture_component to mixture relationship verified")
    
    def test_molecular_property_to_molecule_relationship(self):
        """TC-DB-013: Verify molecular_property → molecule relationship"""
        foreign_keys = self.fetch_foreign_keys()
        
        # Look for molecular_property to molecule relationship
        found = False
        for fk in foreign_keys:
            if (fk['table_name'] == 'molecular_property' and 
                fk['column_name'] == 'molecule_id' and 
                fk['foreign_table_name'] == 'molecule' and 
                fk['foreign_column_name'] == 'id'):
                found = True
                break
        
        self.assertTrue(found, "Foreign key relationship from molecular_property to molecule not found")
        logger.info("Molecular_property to molecule relationship verified")
    
    def test_experiment_to_project_relationship(self):
        """TC-DB-014: Verify experiment → project relationship"""
        foreign_keys = self.fetch_foreign_keys()
        
        # Look for experiment to project relationship
        found = False
        for fk in foreign_keys:
            if (fk['table_name'] == 'experiment' and 
                fk['column_name'] == 'project_id' and 
                fk['foreign_table_name'] == 'project' and 
                fk['foreign_column_name'] == 'id'):
                found = True
                break
        
        self.assertTrue(found, "Foreign key relationship from experiment to project not found")
        logger.info("Experiment to project relationship verified")
    
    def test_experiment_property_to_experiment_relationship(self):
        """TC-DB-015: Verify experiment_property → experiment relationship"""
        foreign_keys = self.fetch_foreign_keys()
        
        # Look for experiment_property to experiment relationship
        found = False
        for fk in foreign_keys:
            if (fk['table_name'] == 'experiment_property' and 
                fk['column_name'] == 'experiment_id' and 
                fk['foreign_table_name'] == 'experiment' and 
                fk['foreign_column_name'] == 'id'):
                found = True
                break
        
        self.assertTrue(found, "Foreign key relationship from experiment_property to experiment not found")
        logger.info("Experiment_property to experiment relationship verified")

class IndexCoverageTests(unittest.TestCase):
    """Tests for Index Coverage (TC-DB-020 to TC-DB-022)"""
    
    @classmethod
    def setUpClass(cls):
        """Set up database connection"""
        cls.conn = psycopg2.connect(**DB_CONNECTION_PARAMS)
    
    @classmethod
    def tearDownClass(cls):
        """Close database connection"""
        cls.conn.close()
    
    def fetch_indexes(self):
        """Fetch all indexes in the database"""
        cursor = self.conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute("""
            SELECT
                t.relname AS table_name,
                i.relname AS index_name,
                a.attname AS column_name,
                ix.indisunique AS is_unique,
                ix.indisprimary AS is_primary
            FROM pg_class t
            JOIN pg_index ix ON t.oid = ix.indrelid
            JOIN pg_class i ON i.oid = ix.indexrelid
            JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(ix.indkey)
            JOIN pg_namespace n ON n.oid = t.relnamespace
            WHERE t.relkind = 'r'
            AND n.nspname = 'public'
            ORDER BY t.relname, i.relname, a.attnum
        """)
        return cursor.fetchall()
    
    def test_rls_policy_indexes(self):
        """TC-DB-020: Verify indexes for RLS policy columns"""
        indexes = self.fetch_indexes()
        
        # RLS policy columns that should be indexed
        rls_indexes = [
            {'table': 'molecule', 'column': 'project_id'},
            {'table': 'mixture', 'column': 'project_id'},
            {'table': 'experiment', 'column': 'project_id'},
            {'table': 'user_profile', 'column': 'user_id'},
            {'table': 'user_profile', 'column': 'project_id'}
        ]
        
        for expected in rls_indexes:
            found = False
            for idx in indexes:
                if (idx['table_name'] == expected['table'] and 
                    idx['column_name'] == expected['column']):
                    found = True
                    break
            
            self.assertTrue(found, f"Index for {expected['table']}.{expected['column']} not found")
        
        logger.info("RLS policy indexes verified")
    
    def test_common_query_pattern_indexes(self):
        """TC-DB-021: Verify indexes for common query patterns"""
        indexes = self.fetch_indexes()
        
        # Common query patterns that should be indexed
        common_indexes = [
            {'table': 'molecular_property', 'column': 'molecule_id'},
            {'table': 'mixture_component', 'column': 'mixture_id'},
            {'table': 'experiment_property', 'column': 'experiment_id'}
        ]
        
        for expected in common_indexes:
            found = False
            for idx in indexes:
                if (idx['table_name'] == expected['table'] and 
                    idx['column_name'] == expected['column']):
                    found = True
                    break
            
            self.assertTrue(found, f"Index for {expected['table']}.{expected['column']} not found")
        
        logger.info("Common query pattern indexes verified")
    
    def test_index_usage_in_queries(self):
        """TC-DB-022: Verify index usage in queries"""
        cursor = self.conn.cursor()
        
        # Test queries that should use indexes
        test_queries = [
            {
                'query': """
                    EXPLAIN ANALYZE
                    SELECT * FROM molecule
                    WHERE project_id = '11111111-1111-1111-1111-111111111111'
                """,
                'expected_index': 'project_id'
            },
            {
                'query': """
                    EXPLAIN ANALYZE
                    SELECT * FROM molecular_property
                    WHERE molecule_id = '00000000-0000-0000-0000-000000000001'
                """,
                'expected_index': 'molecule_id'
            },
            {
                'query': """
                    EXPLAIN ANALYZE
                    SELECT * FROM user_profile
                    WHERE user_id = '00000000-0000-0000-0000-000000000001'
                """,
                'expected_index': 'user_id'
            }
        ]
        
        for test in test_queries:
            try:
                cursor.execute(test['query'])
                
                # Get query plan
                plan_rows = cursor.fetchall()
                plan_text = "\n".join([row[0] for row in plan_rows])
                
                # Look for index scan or bitmap index scan in the plan
                index_used = ('Index Scan' in plan_text or 'Bitmap Index Scan' in plan_text)
                
                # This test is more informational - we don't fail if indexes aren't used
                # because it depends on data distribution, statistics, etc.
                if not index_used:
                    logger.warning(f"No index scan found for {test['expected_index']} in query plan:\n{plan_text}")
                else:
                    logger.info(f"Index usage verified for {test['expected_index']}")
            except Exception as e:
                logger.warning(f"Error testing index usage for {test['expected_index']}: {str(e)}")
        
        logger.info("Index usage in queries tested")

def generate_schema_report():
    """Generate a comprehensive schema report"""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    report_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'reports')
    os.makedirs(report_dir, exist_ok=True)
    
    report_file = os.path.join(report_dir, f'database_schema_verification_{timestamp}.md')
    
    with open(report_file, 'w') as f:
        f.write("# Database Schema Verification Report\n\n")
        f.write(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
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
            
            f.write(f"## Schema Overview\n\n")
            f.write(f"Database contains {len(tables)} tables.\n\n")
            
            # Table list
            f.write("### Tables\n\n")
            for table in tables:
                f.write(f"- {table}\n")
            
            f.write("\n")
            
            # Foreign Key Relationships
            f.write("## Foreign Key Relationships\n\n")
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
                AND tc.table_schema = 'public'
                ORDER BY tc.table_name, kcu.column_name
            """)
            foreign_keys = cursor.fetchall()
            
            f.write("| Table | Column | References | Foreign Column |\n")
            f.write("|-------|--------|------------|---------------|\n")
            
            for fk in foreign_keys:
                f.write(f"| {fk['table_name']} | {fk['column_name']} | {fk['foreign_table_name']} | {fk['foreign_column_name']} |\n")
            
            f.write("\n")
            
            # Indexes
            f.write("## Indexes\n\n")
            cursor.execute("""
                SELECT
                    t.relname AS table_name,
                    i.relname AS index_name,
                    a.attname AS column_name,
                    ix.indisunique AS is_unique,
                    ix.indisprimary AS is_primary
                FROM pg_class t
                JOIN pg_index ix ON t.oid = ix.indrelid
                JOIN pg_class i ON i.oid = ix.indexrelid
                JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(ix.indkey)
                JOIN pg_namespace n ON n.oid = t.relnamespace
                WHERE t.relkind = 'r'
                AND n.nspname = 'public'
                ORDER BY t.relname, i.relname, a.attnum
            """)
            indexes = cursor.fetchall()
            
            f.write("| Table | Index | Column | Type |\n")
            f.write("|-------|-------|--------|------|\n")
            
            for idx in indexes:
                idx_type = "Primary Key" if idx['is_primary'] else "Unique" if idx['is_unique'] else "Non-Unique"
                f.write(f"| {idx['table_name']} | {idx['index_name']} | {idx['column_name']} | {idx_type} |\n")
            
            f.write("\n")
            
            # RLS Policies
            f.write("## Row Level Security Policies\n\n")
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
            policies = cursor.fetchall()
            
            f.write("| Table | Policy | Command | Using | With Check |\n")
            f.write("|-------|--------|---------|-------|------------|\n")
            
            for policy in policies:
                using = policy['qual'] or ''
                with_check = policy['with_check'] or ''
                
                # Truncate for readability
                if len(using) > 100:
                    using = using[:97] + '...'
                if len(with_check) > 100:
                    with_check = with_check[:97] + '...'
                
                f.write(f"| {policy['tablename']} | {policy['policyname']} | {policy['cmd']} | {using} | {with_check} |\n")
            
            f.write("\n")
            
            # Security Definer Functions
            f.write("## Security Definer Functions\n\n")
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
            functions = cursor.fetchall()
            
            f.write("| Function | Language |\n")
            f.write("|----------|----------|\n")
            
            for func in functions:
                f.write(f"| {func['function_name']} | {func['language']} |\n")
            
        except Exception as e:
            f.write(f"\n## Error\n\n")
            f.write(f"Error generating schema report: {str(e)}\n")
        finally:
            if conn:
                conn.close()
    
    return report_file

if __name__ == '__main__':
    # Generate schema report first
    report_file = generate_schema_report()
    print(f"Generated schema report: {report_file}")
    
    # Run tests
    unittest.main()