"""
API Batch Endpoints Integration Tests for CryoProtect v2

Covers:
- /api/v1/scoring/batch
- /api/v1/batch/<operation_id> (GET: status, DELETE: cancel)
- Error handling, authentication, and edge cases
"""

import os
import sys
import json
import uuid
import unittest
from unittest.mock import patch, MagicMock

# Add parent directory to path to import from api
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import MockSupabaseBaseTestCase
from app import create_app

class TestAPIBatchEndpoints(MockSupabaseBaseTestCase):
    """Integration tests for batch API endpoints."""

    def setUp(self):
        super().setUp()
        self.client = self.app.test_client()
        self.auth_token = str(uuid.uuid4())
        self.headers = {'Authorization': f'Bearer {self.auth_token}'}

    @patch('api.scoring.batch_score_molecules')
    def test_scoring_batch_endpoint_success(self, mock_batch_score):
        """Test POST /api/v1/scoring/batch with valid data."""
        mock_batch_score.return_value = [
            {"id": "mol-1", "score": 85.0},
            {"id": "mol-2", "score": 72.5}
        ]
        payload = {"molecule_ids": ["mol-1", "mol-2"]}
        response = self.client.post('/api/v1/scoring/batch', json=payload, headers=self.headers)
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 2)
        self.assertIn("score", data[0])

    def test_scoring_batch_endpoint_missing_auth(self):
        """Test POST /api/v1/scoring/batch without authentication."""
        payload = {"molecule_ids": ["mol-1"]}
        response = self.client.post('/api/v1/scoring/batch', json=payload)
        self.assertIn(response.status_code, [401, 403])

    def test_scoring_batch_endpoint_invalid_payload(self):
        """Test POST /api/v1/scoring/batch with invalid payload."""
        response = self.client.post('/api/v1/scoring/batch', json={}, headers=self.headers)
        self.assertIn(response.status_code, [400, 422])

    @patch('api.batch_resources.get_supabase_client')
    @patch('api.batch_resources.get_user_id')
    def test_batch_operation_status_success(self, mock_get_user_id, mock_get_supabase_client):
        """Test GET /api/v1/batch/<operation_id> for a valid operation."""
        operation_id = str(uuid.uuid4())
        mock_get_user_id.return_value = "user-1"
        mock_supabase = MagicMock()
        mock_supabase.table().select().eq().execute.return_value.data = [{
            "id": operation_id,
            "status": "SUCCESS",
            "operation": "property_calculation",
            "entity_type": "molecule",
            "created_by": "user-1",
            "created_at": "2025-04-21T00:00:00",
            "completed_at": "2025-04-21T00:01:00",
            "result_count": 2,
            "error_count": 0,
            "is_cancelled": False
        }]
        mock_get_supabase_client.return_value = mock_supabase

        response = self.client.get(f'/api/v1/batch/{operation_id}', headers=self.headers)
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["operation_id"], operation_id)
        self.assertEqual(data["status"], "SUCCESS")

    @patch('api.batch_resources.get_supabase_client')
    @patch('api.batch_resources.get_user_id')
    def test_batch_operation_status_not_found(self, mock_get_user_id, mock_get_supabase_client):
        """Test GET /api/v1/batch/<operation_id> for a nonexistent operation."""
        operation_id = str(uuid.uuid4())
        mock_get_user_id.return_value = "user-1"
        mock_supabase = MagicMock()
        mock_supabase.table().select().eq().execute.return_value.data = []
        mock_get_supabase_client.return_value = mock_supabase

        response = self.client.get(f'/api/v1/batch/{operation_id}', headers=self.headers)
        self.assertEqual(response.status_code, 404)
        data = json.loads(response.data)
        self.assertEqual(data["status"], "ERROR")

    @patch('api.batch_resources.get_supabase_client')
    @patch('api.batch_resources.get_user_id')
    def test_batch_operation_cancel_success(self, mock_get_user_id, mock_get_supabase_client):
        """Test DELETE /api/v1/batch/<operation_id> for a valid operation."""
        operation_id = str(uuid.uuid4())
        mock_get_user_id.return_value = "user-1"
        mock_supabase = MagicMock()
        # Simulate operation found and not completed
        mock_supabase.table().select().eq().execute.return_value.data = [{
            "id": operation_id,
            "status": "RUNNING",
            "created_by": "user-1"
        }]
        mock_supabase.table().update().eq().execute.return_value.error = None
        mock_get_supabase_client.return_value = mock_supabase

        response = self.client.delete(f'/api/v1/batch/{operation_id}', headers=self.headers)
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["status"], "SUCCESS")

    @patch('api.batch_resources.get_supabase_client')
    @patch('api.batch_resources.get_user_id')
    def test_batch_operation_cancel_already_completed(self, mock_get_user_id, mock_get_supabase_client):
        """Test DELETE /api/v1/batch/<operation_id> for an already completed operation."""
        operation_id = str(uuid.uuid4())
        mock_get_user_id.return_value = "user-1"
        mock_supabase = MagicMock()
        # Simulate operation found and already completed
        mock_supabase.table().select().eq().execute.return_value.data = [{
            "id": operation_id,
            "status": "SUCCESS",
            "created_by": "user-1"
        }]
        mock_get_supabase_client.return_value = mock_supabase

        response = self.client.delete(f'/api/v1/batch/{operation_id}', headers=self.headers)
        self.assertEqual(response.status_code, 400)
        data = json.loads(response.data)
        self.assertEqual(data["status"], "ERROR")

    def test_batch_operation_status_missing_auth(self):
        """Test GET /api/v1/batch/<operation_id> without authentication."""
        operation_id = str(uuid.uuid4())
        response = self.client.get(f'/api/v1/batch/{operation_id}')
        self.assertIn(response.status_code, [401, 403])

    def test_batch_operation_cancel_missing_auth(self):
        """Test DELETE /api/v1/batch/<operation_id> without authentication."""
        operation_id = str(uuid.uuid4())
        response = self.client.delete(f'/api/v1/batch/{operation_id}')
        self.assertIn(response.status_code, [401, 403])

if __name__ == '__main__':
    unittest.main()