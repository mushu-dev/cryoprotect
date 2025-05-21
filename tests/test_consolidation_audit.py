"""
Tests for the consolidated molecule audit trail functionality.

This module tests the database triggers, API endpoints, and overall
functionality of the consolidation audit system.
"""

import unittest
import uuid
from datetime import datetime, timedelta
import json

from flask import url_for
from database.adapter import get_database_connection
from api.utils import get_supabase_client
import app as flask_app

class ConsolidationAuditTestCase(unittest.TestCase):
    """Test case for consolidated molecule audit functionality."""
    
    def setUp(self):
        """Set up test environment before each test."""
        self.app = flask_app.create_app('test')
        self.app_context = self.app.app_context()
        self.app_context.push()
        self.client = self.app.test_client()
        
        # Get database connection
        self.db = get_database_connection()
        
        # Create test molecules for consolidation testing
        self.primary_id = str(uuid.uuid4())
        self.secondary_id = str(uuid.uuid4())
        
        # Insert test molecules
        self.db.execute(
            """
            INSERT INTO molecules (id, name, smiles)
            VALUES (%s, %s, %s)
            """,
            (self.primary_id, "Test Primary Molecule", "C")
        )
        
        self.db.execute(
            """
            INSERT INTO molecules (id, name, smiles)
            VALUES (%s, %s, %s)
            """,
            (self.secondary_id, "Test Secondary Molecule", "CC")
        )
        
        # Create a test admin user
        self.admin_user_id = str(uuid.uuid4())
        self.admin_token = "test_admin_token"
        
        # Create mock user data in the app context
        self.app.config['TEST_USER'] = {
            'id': self.admin_user_id,
            'email': 'admin@test.com',
            'app_metadata': {
                'roles': ['admin']
            }
        }
        self.app.config['TEST_TOKEN'] = self.admin_token
    
    def tearDown(self):
        """Clean up after each test."""
        # Remove test molecules
        self.db.execute(
            "DELETE FROM molecules WHERE id IN (%s, %s)",
            (self.primary_id, self.secondary_id)
        )
        
        # Remove any audit records created by tests
        self.db.execute(
            "DELETE FROM molecule_consolidation_audit WHERE primary_molecule_id = %s OR secondary_molecule_id = %s",
            (self.primary_id, self.secondary_id)
        )
        
        self.app_context.pop()
    
    def test_audit_trigger_on_consolidation(self):
        """Test that the audit trigger fires when a molecule is consolidated."""
        # Consolidate the secondary molecule to the primary
        result = self.db.execute(
            """
            UPDATE molecules 
            SET consolidated_to = %s
            WHERE id = %s
            RETURNING id
            """,
            (self.primary_id, self.secondary_id)
        )
        
        # Check if update was successful
        self.assertIsNotNone(result.fetchone())
        
        # Verify audit record was created
        result = self.db.execute(
            """
            SELECT * FROM molecule_consolidation_audit
            WHERE primary_molecule_id = %s AND secondary_molecule_id = %s
            """,
            (self.primary_id, self.secondary_id)
        )
        
        audit_record = result.fetchone()
        self.assertIsNotNone(audit_record)
        self.assertEqual(audit_record['operation_type'], 'CONSOLIDATE')
    
    def test_audit_trigger_on_deconsolidation(self):
        """Test that the audit trigger fires when a molecule is deconsolidated."""
        # First consolidate the molecule
        self.db.execute(
            """
            UPDATE molecules 
            SET consolidated_to = %s
            WHERE id = %s
            """,
            (self.primary_id, self.secondary_id)
        )
        
        # Then deconsolidate it
        result = self.db.execute(
            """
            UPDATE molecules 
            SET consolidated_to = NULL
            WHERE id = %s
            RETURNING id
            """,
            (self.secondary_id,)
        )
        
        # Check if update was successful
        self.assertIsNotNone(result.fetchone())
        
        # Verify audit records were created (should be two - one for consolidate, one for deconsolidate)
        result = self.db.execute(
            """
            SELECT * FROM molecule_consolidation_audit
            WHERE secondary_molecule_id = %s
            ORDER BY performed_at DESC
            """,
            (self.secondary_id,)
        )
        
        audit_records = result.fetchall()
        self.assertEqual(len(audit_records), 2)
        self.assertEqual(audit_records[0]['operation_type'], 'DECONSOLIDATE')
    
    def test_audit_api_endpoint(self):
        """Test the consolidation audit API endpoint."""
        # Perform a consolidation to create an audit record
        self.db.execute(
            """
            UPDATE molecules 
            SET consolidated_to = %s
            WHERE id = %s
            """,
            (self.primary_id, self.secondary_id)
        )
        
        # Call the API endpoint with admin token
        response = self.client.get(
            '/api/v1/admin/consolidation-audit',
            headers={
                'Authorization': f'Bearer {self.admin_token}'
            }
        )
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        
        # Basic validation of response structure
        self.assertIn('data', data)
        self.assertIn('audit_records', data['data'])
        self.assertIn('pagination', data)
        
        # Verify our test record is in the results
        found = False
        for record in data['data']['audit_records']:
            if (record['primary_molecule_id'] == self.primary_id and 
                record['secondary_molecule_id'] == self.secondary_id):
                found = True
                self.assertEqual(record['operation_type'], 'CONSOLIDATE')
                self.assertEqual(record['primary_molecule_name'], 'Test Primary Molecule')
                self.assertEqual(record['secondary_molecule_name'], 'Test Secondary Molecule')
                break
        
        self.assertTrue(found, "Test audit record not found in API response")
    
    def test_audit_api_permissions(self):
        """Test that non-admin users cannot access the audit endpoint."""
        # Create a non-admin user context
        self.app.config['TEST_USER'] = {
            'id': str(uuid.uuid4()),
            'email': 'user@test.com',
            'app_metadata': {
                'roles': ['user']
            }
        }
        self.app.config['TEST_TOKEN'] = "test_user_token"
        
        # Call the API endpoint with non-admin token
        response = self.client.get(
            '/api/v1/admin/consolidation-audit',
            headers={
                'Authorization': f'Bearer test_user_token'
            }
        )
        
        # Check response - should be forbidden
        self.assertEqual(response.status_code, 403)
    
    def test_audit_api_filtering(self):
        """Test filtering of audit records in the API."""
        # Create multiple audit records by consolidating and deconsolidating
        self.db.execute(
            """
            UPDATE molecules 
            SET consolidated_to = %s
            WHERE id = %s
            """,
            (self.primary_id, self.secondary_id)
        )
        
        # Wait a moment to ensure timestamps are different
        import time
        time.sleep(1)
        
        self.db.execute(
            """
            UPDATE molecules 
            SET consolidated_to = NULL
            WHERE id = %s
            """,
            (self.secondary_id,)
        )
        
        # Test filtering by operation_type
        response = self.client.get(
            '/api/v1/admin/consolidation-audit?operation_type=DECONSOLIDATE',
            headers={
                'Authorization': f'Bearer {self.admin_token}'
            }
        )
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        
        # Verify only DECONSOLIDATE operations are returned
        for record in data['data']['audit_records']:
            self.assertEqual(record['operation_type'], 'DECONSOLIDATE')
        
        # Test filtering by molecule_id
        response = self.client.get(
            f'/api/v1/admin/consolidation-audit?molecule_id={self.secondary_id}',
            headers={
                'Authorization': f'Bearer {self.admin_token}'
            }
        )
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        
        # Verify only records for secondary_id are returned
        for record in data['data']['audit_records']:
            self.assertTrue(
                record['primary_molecule_id'] == self.secondary_id or 
                record['secondary_molecule_id'] == self.secondary_id
            )

if __name__ == '__main__':
    unittest.main()