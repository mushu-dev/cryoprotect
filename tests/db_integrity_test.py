#!/usr/bin/env python3
"""
Database integrity test suite for CryoProtect.

This script validates data integrity before and after database schema changes.
It ensures that critical relationships and data are preserved during migrations.
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
import unittest
from datetime import datetime

class DatabaseIntegrityTest(unittest.TestCase):
    """Test case for database integrity validation."""
    
    @classmethod
    def setUpClass(cls):
        """Set up database connection and snapshot the current state."""
        # Get database connection parameters from environment
        cls.db_params = {
            'host': os.environ.get('DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
            'port': os.environ.get('DB_PORT', '5432'),
            'dbname': os.environ.get('DB_NAME', 'postgres'),
            'user': os.environ.get('DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
            'password': os.environ.get('DB_PASSWORD'),
            'sslmode': 'require'
        }
        
        # Ensure password is provided
        if not cls.db_params['password']:
            raise ValueError("DB_PASSWORD environment variable must be set")
        
        # Connect to the database
        cls.conn = psycopg2.connect(**cls.db_params)
        cls.conn.autocommit = False
        
        # Take an initial snapshot
        cls.initial_snapshot = cls._get_database_snapshot()
        
        print(f"Initial database snapshot taken at {datetime.now().isoformat()}")
    
    @classmethod
    def tearDownClass(cls):
        """Close database connection."""
        if cls.conn:
            cls.conn.close()
    
    @classmethod
    def _get_database_snapshot(cls):
        """Get a snapshot of critical database tables and counts."""
        snapshot = {}
        
        with cls.conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Get table counts
            cursor.execute("""
                SELECT
                    table_name,
                    (SELECT count(*) FROM information_schema.columns 
                     WHERE table_schema = 'public' AND columns.table_name = tables.table_name) as column_count,
                    pg_total_relation_size('public.' || table_name) as table_size
                FROM
                    information_schema.tables tables
                WHERE
                    table_schema = 'public'
                    AND table_type = 'BASE TABLE'
                ORDER BY
                    table_name;
            """)
            
            snapshot['tables'] = {}
            for row in cursor.fetchall():
                table_name = row['table_name']
                snapshot['tables'][table_name] = {
                    'column_count': row['column_count'],
                    'size': row['table_size']
                }
                
                # Get row count for each table
                cursor.execute(f"SELECT count(*) as row_count FROM {table_name}")
                count_result = cursor.fetchone()
                snapshot['tables'][table_name]['row_count'] = count_result['row_count']
            
            # Get specific critical counts
            critical_counts = {
                'molecules': "SELECT COUNT(*) FROM molecules",
                'molecular_properties': "SELECT COUNT(*) FROM molecular_properties",
                'property_types': "SELECT COUNT(*) FROM property_types",
                'predictions': "SELECT COUNT(*) FROM predictions",
                'mixtures': "SELECT COUNT(*) FROM mixtures",
                'mixture_components': "SELECT COUNT(*) FROM mixture_components"
            }
            
            snapshot['critical_counts'] = {}
            for name, query in critical_counts.items():
                cursor.execute(query)
                result = cursor.fetchone()
                snapshot['critical_counts'][name] = result['count']
                
            # Check specific relationships
            if 'calculation_method' in snapshot['tables'] and 'calculation_methods' in snapshot['tables']:
                # Check calculation method data
                cursor.execute("SELECT id, name FROM calculation_method")
                snapshot['calculation_method_data'] = cursor.fetchall()
                
                cursor.execute("SELECT id, name FROM calculation_methods")
                snapshot['calculation_methods_data'] = cursor.fetchall()
                
                # Check foreign key usage
                cursor.execute("""
                    SELECT 
                        table_name, 
                        column_name 
                    FROM 
                        information_schema.constraint_column_usage 
                    WHERE 
                        table_schema = 'public' 
                        AND table_name = 'calculation_method'
                """)
                snapshot['calculation_method_references'] = cursor.fetchall()
                
                cursor.execute("""
                    SELECT 
                        table_name, 
                        column_name 
                    FROM 
                        information_schema.constraint_column_usage 
                    WHERE 
                        table_schema = 'public' 
                        AND table_name = 'calculation_methods'
                """)
                snapshot['calculation_methods_references'] = cursor.fetchall()
                
                # Check if predictions table references calculation_method or calculation_methods
                cursor.execute("""
                    SELECT 
                        ccu.table_name AS referenced_table,
                        ccu.column_name AS referenced_column,
                        tc.constraint_name
                    FROM 
                        information_schema.table_constraints AS tc 
                    JOIN 
                        information_schema.key_column_usage AS kcu
                        ON tc.constraint_name = kcu.constraint_name
                        AND tc.table_schema = kcu.table_schema
                    JOIN 
                        information_schema.constraint_column_usage AS ccu
                        ON ccu.constraint_name = tc.constraint_name
                        AND ccu.table_schema = tc.table_schema
                    WHERE 
                        tc.constraint_type = 'FOREIGN KEY' 
                        AND tc.table_schema = 'public'
                        AND tc.table_name = 'predictions'
                        AND (ccu.table_name = 'calculation_method' OR ccu.table_name = 'calculation_methods');
                """)
                snapshot['predictions_references_calculation'] = cursor.fetchall()
            
        return snapshot
    
    def test_table_existence(self):
        """Test that all tables exist with expected row counts."""
        current_snapshot = self._get_database_snapshot()
        
        # Check that all tables from initial snapshot still exist
        for table_name in self.initial_snapshot['tables']:
            self.assertIn(table_name, current_snapshot['tables'], 
                        f"Table {table_name} no longer exists")
    
    def test_critical_row_counts(self):
        """Test that critical tables have the same number of rows."""
        current_snapshot = self._get_database_snapshot()
        
        # Check critical counts
        for table_name, count in self.initial_snapshot['critical_counts'].items():
            if table_name in current_snapshot['critical_counts']:
                self.assertEqual(count, current_snapshot['critical_counts'][table_name],
                               f"Row count for {table_name} has changed. Expected {count}, got {current_snapshot['critical_counts'][table_name]}")
    
    def test_calculation_methods_consolidation(self):
        """Test that calculation_methods table contains all data from calculation_method."""
        # Skip if we don't have both tables in the initial snapshot
        if ('calculation_method' not in self.initial_snapshot['tables'] or 
            'calculation_methods' not in self.initial_snapshot['tables']):
            self.skipTest("One or both calculation tables don't exist")
        
        current_snapshot = self._get_database_snapshot()
        
        # If calculation_method still exists, check data integrity
        if 'calculation_method' in current_snapshot['tables']:
            # Check that all calculation_method entries exist in calculation_methods
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                for method in self.initial_snapshot['calculation_method_data']:
                    method_id = method['id']
                    cursor.execute("""
                        SELECT COUNT(*) as count FROM calculation_methods WHERE id = %s
                    """, (method_id,))
                    result = cursor.fetchone()
                    self.assertGreater(result['count'], 0, 
                                     f"Calculation method with ID {method_id} not found in calculation_methods")
        
        # If calculation_method has been dropped, check that all its data was migrated
        if 'calculation_method' not in current_snapshot['tables']:
            for method in self.initial_snapshot['calculation_method_data']:
                method_id = method['id']
                method_name = method['name']
                
                # Find this method in the current calculation_methods data
                found = False
                for current_method in current_snapshot.get('calculation_methods_data', []):
                    if current_method['id'] == method_id or current_method['name'] == method_name:
                        found = True
                        break
                
                self.assertTrue(found, f"Method {method_name} (ID: {method_id}) was not migrated to calculation_methods")
    
    def test_foreign_key_updates(self):
        """Test that foreign keys were properly updated during consolidation."""
        # Skip if there were no foreign keys in the initial snapshot
        if not self.initial_snapshot.get('predictions_references_calculation'):
            self.skipTest("No foreign key references to check")
        
        current_snapshot = self._get_database_snapshot()
        
        # If calculation_method has been dropped, check that predictions now references calculation_methods
        if ('calculation_method' not in current_snapshot['tables'] and 
            self.initial_snapshot.get('predictions_references_calculation')):
            # Check that predictions table now references calculation_methods
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute("""
                    SELECT 
                        ccu.table_name AS referenced_table,
                        ccu.column_name AS referenced_column
                    FROM 
                        information_schema.table_constraints AS tc 
                    JOIN 
                        information_schema.key_column_usage AS kcu
                        ON tc.constraint_name = kcu.constraint_name
                        AND tc.table_schema = kcu.table_schema
                    JOIN 
                        information_schema.constraint_column_usage AS ccu
                        ON ccu.constraint_name = tc.constraint_name
                        AND ccu.table_schema = tc.table_schema
                    WHERE 
                        tc.constraint_type = 'FOREIGN KEY' 
                        AND tc.table_schema = 'public'
                        AND tc.table_name = 'predictions'
                        AND ccu.table_name = 'calculation_methods';
                """)
                references = cursor.fetchall()
                
                self.assertGreater(len(references), 0, 
                                 "predictions table should reference calculation_methods after consolidation")

def main():
    unittest.main()

if __name__ == "__main__":
    main()