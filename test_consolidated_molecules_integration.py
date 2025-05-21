#!/usr/bin/env python3
"""
CryoProtect Analyzer - Consolidated Molecules Integration Test

This script tests the consolidated molecules integration, including:
1. Database schema verification
2. Consolidated molecule API functionality
3. RDKit integration with consolidated molecules
4. Property migration between consolidated molecules

Run this script after applying the consolidated molecule updates.
"""

import os
import sys
import logging
import json
import unittest
import requests
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple
import uuid

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(f'consolidated_test_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)

# Add the parent directory to the path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

# Import database adapter
try:
    from database.adapter import get_connection
except ImportError:
    logger.error("Could not import database adapter. Make sure you're running from the project root.")
    sys.exit(1)

# Try to import the consolidated RDKit wrapper
try:
    from rdkit_wrapper_consolidated import (
        get_primary_molecule_id,
        calculate_properties_for_consolidated,
        migrate_properties
    )
    RDKIT_WRAPPER_AVAILABLE = True
except ImportError:
    logger.warning("Could not import consolidated RDKit wrapper. Some tests will be skipped.")
    RDKIT_WRAPPER_AVAILABLE = False

class ConsolidatedMoleculesTest(unittest.TestCase):
    """Test case for consolidated molecules functionality."""
    
    def setUp(self):
        """Set up test database connection."""
        self.conn = get_connection()
        self.cursor = self.conn.cursor()
        
        # Get some test molecule IDs
        self.cursor.execute("""
        SELECT id FROM molecules LIMIT 10
        """)
        self.test_molecules = [row[0] for row in self.cursor.fetchall()]
        
        # Check for any existing consolidated molecules
        self.cursor.execute("""
        SELECT 
            cm.id, 
            cm.molecule_status,
            cm.primary_molecule_id
        FROM consolidated_molecules cm
        WHERE cm.molecule_status = 'duplicate'
        LIMIT 5
        """)
        
        self.consolidated_molecules = []
        results = self.cursor.fetchall()
        
        for row in results:
            self.consolidated_molecules.append({
                'id': row[0],
                'status': row[1],
                'primary_id': row[2]
            })
        
        # Log test setup
        logger.info(f"Test setup complete. Found {len(self.test_molecules)} test molecules and {len(self.consolidated_molecules)} consolidated molecules.")
    
    def tearDown(self):
        """Clean up after tests."""
        if self.conn:
            self.conn.close()
    
    def test_consolidated_molecules_schema(self):
        """Test that the consolidated_molecules table has the correct structure."""
        self.cursor.execute("""
        SELECT column_name, data_type
        FROM information_schema.columns
        WHERE table_name = 'consolidated_molecules'
        """)
        
        columns = {row[0]: row[1] for row in self.cursor.fetchall()}
        
        # Check for required columns
        self.assertIn('id', columns, "id column is missing from consolidated_molecules table")
        self.assertIn('molecule_status', columns, "molecule_status column is missing from consolidated_molecules table")
        self.assertIn('primary_molecule_id', columns, "primary_molecule_id column is missing from consolidated_molecules table")
        
        # Check column types
        self.assertEqual(columns.get('molecule_status'), 'character varying', "molecule_status should be of type character varying")
        
        # Check for constraints
        self.cursor.execute("""
        SELECT constraint_name, constraint_type
        FROM information_schema.table_constraints
        WHERE table_name = 'consolidated_molecules' AND constraint_type = 'CHECK'
        """)
        
        constraints = {row[0]: row[1] for row in self.cursor.fetchall()}
        self.assertTrue(any('status' in name.lower() for name in constraints.keys()), 
                      "No CHECK constraint found for molecule_status column")
        
        # Check for indexes
        self.cursor.execute("""
        SELECT indexname
        FROM pg_indexes
        WHERE tablename = 'consolidated_molecules'
        """)
        
        indexes = [row[0] for row in self.cursor.fetchall()]
        self.assertTrue(any('inchikey' in name.lower() for name in indexes), 
                      "No index found for inchikey column")
        self.assertTrue(any('primary' in name.lower() for name in indexes), 
                      "No index found for primary_molecule_id column")
        
        logger.info("Database schema verification passed.")
    
    def test_scientific_data_audit_schema(self):
        """Test that the scientific_data_audit table has the correct structure."""
        self.cursor.execute("""
        SELECT column_name, data_type
        FROM information_schema.columns
        WHERE table_name = 'scientific_data_audit'
        """)
        
        columns = {row[0]: row[1] for row in self.cursor.fetchall()}
        
        # Check for required columns
        self.assertIn('id', columns, "id column is missing from scientific_data_audit table")
        self.assertIn('table_name', columns, "table_name column is missing from scientific_data_audit table")
        self.assertIn('record_id', columns, "record_id column is missing from scientific_data_audit table")
        self.assertIn('operation', columns, "operation column is missing from scientific_data_audit table")
        self.assertIn('old_value', columns, "old_value column is missing from scientific_data_audit table")
        self.assertIn('new_value', columns, "new_value column is missing from scientific_data_audit table")
        self.assertIn('timestamp', columns, "timestamp column is missing from scientific_data_audit table")
        
        # Check column types
        self.assertEqual(columns.get('old_value'), 'jsonb', "old_value should be of type jsonb")
        self.assertEqual(columns.get('new_value'), 'jsonb', "new_value should be of type jsonb")
        
        logger.info("Scientific data audit schema verification passed.")
    
    @unittest.skipIf(not RDKIT_WRAPPER_AVAILABLE, "Consolidated RDKit wrapper not available")
    def test_get_primary_molecule_id(self):
        """Test the get_primary_molecule_id function."""
        # Skip if no consolidated molecules found
        if not self.consolidated_molecules:
            logger.warning("No consolidated molecules found, skipping test_get_primary_molecule_id")
            return
        
        for molecule in self.consolidated_molecules:
            duplicate_id = molecule['id']
            expected_primary_id = molecule['primary_id']
            
            # Test the function
            actual_primary_id = get_primary_molecule_id(duplicate_id)
            self.assertEqual(expected_primary_id, actual_primary_id,
                           f"get_primary_molecule_id returned {actual_primary_id} but expected {expected_primary_id}")
            
            # Test with a non-consolidated molecule ID
            original_id = self.test_molecules[0]
            self.assertEqual(original_id, get_primary_molecule_id(original_id),
                           f"get_primary_molecule_id should return the original ID for non-consolidated molecules")
        
        logger.info("get_primary_molecule_id test passed.")
    
    @unittest.skipIf(not RDKIT_WRAPPER_AVAILABLE, "Consolidated RDKit wrapper not available")
    def test_property_migration(self):
        """Test property migration between consolidated molecules."""
        # Skip if no consolidated molecules found
        if not self.consolidated_molecules:
            logger.warning("No consolidated molecules found, skipping test_property_migration")
            return
        
        for molecule in self.consolidated_molecules:
            duplicate_id = molecule['id']
            primary_id = molecule['primary_id']
            
            # Create a test property
            test_property_type = "test_property_" + datetime.now().strftime("%Y%m%d%H%M%S")
            test_property_value = "test_value_" + str(uuid.uuid4())
            
            # First, create the property type
            self.cursor.execute("""
            INSERT INTO property_types (name, data_type, description, created_at, updated_at)
            VALUES (%s, 'string', 'Test property for migration', NOW(), NOW())
            RETURNING id
            """, (test_property_type,))
            
            property_type_id = self.cursor.fetchone()[0]
            
            # Add the property to the duplicate molecule
            self.cursor.execute("""
            INSERT INTO molecular_properties (molecule_id, property_type_id, property_value, calculation_method, created_at, updated_at)
            VALUES (%s, %s, %s, 'test', NOW(), NOW())
            """, (duplicate_id, property_type_id, test_property_value))
            
            self.conn.commit()
            
            # Now migrate the property
            result = migrate_properties(duplicate_id, primary_id, [test_property_type])
            
            # Check that migration was successful
            self.assertEqual(result['source_id'], duplicate_id, "Source ID should match")
            self.assertEqual(result['target_id'], primary_id, "Target ID should match")
            self.assertEqual(len(result['migrated_properties']), 1, "Should have migrated exactly one property")
            self.assertEqual(result['migrated_properties'][0]['property_type'], test_property_type, 
                           "Migrated property type should match")
            self.assertEqual(result['migrated_properties'][0]['value'], test_property_value, 
                           "Migrated property value should match")
            
            # Verify in the database
            self.cursor.execute("""
            SELECT property_value
            FROM molecular_properties
            WHERE molecule_id = %s AND property_type_id = %s
            """, (primary_id, property_type_id))
            
            db_value = self.cursor.fetchone()
            self.assertIsNotNone(db_value, "Property should exist in the primary molecule")
            self.assertEqual(db_value[0], test_property_value, "Property value should match")
            
            # Clean up
            self.cursor.execute("""
            DELETE FROM molecular_properties WHERE property_type_id = %s
            """, (property_type_id,))
            
            self.cursor.execute("""
            DELETE FROM property_types WHERE id = %s
            """, (property_type_id,))
            
            self.conn.commit()
        
        logger.info("Property migration test passed.")
    
    def test_api_endpoints(self):
        """Test the API endpoints for consolidated molecules."""
        # This requires the API to be running, so skip if not available
        api_url = os.environ.get('API_URL', 'http://localhost:5000')
        
        try:
            # Check if API is available
            response = requests.get(f"{api_url}/api/v1/health")
            if response.status_code != 200:
                logger.warning(f"API not available at {api_url}, skipping API tests")
                return
            
            # Test consolidated molecules endpoint
            if self.consolidated_molecules:
                duplicate_id = self.consolidated_molecules[0]['id']
                primary_id = self.consolidated_molecules[0]['primary_id']
                
                # Test getting a duplicate molecule
                response = requests.get(f"{api_url}/api/v1/consolidated/molecules/{duplicate_id}")
                self.assertEqual(response.status_code, 200, "Should get successful response for duplicate molecule")
                data = response.json()
                self.assertIn('data', data, "Response should contain data field")
                self.assertIn('is_consolidated', data['data'], "Response should include is_consolidated field")
                self.assertTrue(data['data']['is_consolidated'], "is_consolidated should be True for duplicate molecule")
                self.assertEqual(data['data']['primary_molecule_id'], primary_id, 
                               "primary_molecule_id should match expected value")
                
                # Test getting a primary molecule
                response = requests.get(f"{api_url}/api/v1/consolidated/molecules/{primary_id}")
                self.assertEqual(response.status_code, 200, "Should get successful response for primary molecule")
                data = response.json()
                self.assertIn('data', data, "Response should contain data field")
                self.assertIn('molecule_status', data['data'], "Response should include molecule_status field")
                self.assertEqual(data['data']['molecule_status'], 'primary', 
                               "molecule_status should be 'primary' for primary molecule")
                
                # Test getting audit history
                response = requests.get(f"{api_url}/api/v1/consolidated/molecules/{duplicate_id}/audit")
                self.assertEqual(response.status_code, 200, "Should get successful response for audit history")
                
                # Test search endpoint
                response = requests.get(f"{api_url}/api/v1/consolidated/search?status=duplicate&limit=5")
                self.assertEqual(response.status_code, 200, "Should get successful response for search")
                data = response.json()
                self.assertIn('data', data, "Response should contain data field")
                self.assertIn('molecules', data['data'], "Response should include molecules field")
                self.assertGreater(len(data['data']['molecules']), 0, "Should find at least one duplicate molecule")
            
            logger.info("API endpoint tests passed.")
        except requests.exceptions.ConnectionError:
            logger.warning(f"Could not connect to API at {api_url}, skipping API tests")

def run_full_test():
    """Run all tests."""
    runner = unittest.TextTestRunner(verbosity=2)
    suite = unittest.TestLoader().loadTestsFromTestCase(ConsolidatedMoleculesTest)
    result = runner.run(suite)
    
    # Generate report
    success = result.wasSuccessful()
    test_count = result.testsRun
    failure_count = len(result.failures)
    error_count = len(result.errors)
    
    logger.info("=" * 60)
    logger.info(f"Test Results: {'SUCCESS' if success else 'FAILURE'}")
    logger.info(f"Tests Run: {test_count}")
    logger.info(f"Failures: {failure_count}")
    logger.info(f"Errors: {error_count}")
    logger.info("=" * 60)
    
    if not success:
        logger.error("Tests failed. See details above.")
        return False
    
    return True

if __name__ == "__main__":
    run_full_test()