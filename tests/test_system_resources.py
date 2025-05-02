"""
API Tests for System Resources Endpoints in CryoProtect v2

Covers:
- /api/v1/system/status
- /api/v1/system/logs
- /api/v1/system/metrics
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

class TestSystemResourcesEndpoints(MockSupabaseBaseTestCase):
    """Integration tests for system resources API endpoints."""

    def setUp(self):
        super().setUp()
        self.client = self.app.test_client()
        self.auth_token = str(uuid.uuid4())
        self.headers = {'Authorization': f'Bearer {self.auth_token}'}
        
        # Mock admin user
        self.admin_user_id = str(uuid.uuid4())
        self.admin_headers = {'Authorization': f'Bearer admin-{self.auth_token}'}
        
        # Mock regular user
        self.regular_user_id = str(uuid.uuid4())
        self.regular_headers = {'Authorization': f'Bearer regular-{self.auth_token}'}

    @patch('api.system_resources.get_system_status')
    @patch('api.rbac.UserRoleManager.get_user_roles')
    @patch('api.system_resources.get_user_id')
    def test_system_status_endpoint_admin(self, mock_get_user_id, mock_get_user_roles, mock_get_system_status):
        """Test GET /api/v1/system/status with admin user."""
        # Mock admin user
        mock_get_user_id.return_value = self.admin_user_id
        mock_get_user_roles.return_value = [{"name": "admin"}]
        
        # Mock system status
        mock_get_system_status.return_value = {
            "status": "healthy",
            "version": "2.0.0",
            "uptime": "3d 12h 45m",
            "database_connection": "connected",
            "api_server": "running",
            "rdkit_integration": "available",
            "memory_usage": "45%",
            "cpu_usage": "32%",
            "disk_usage": "68%"
        }
        
        # Make request
        response = self.client.get('/api/v1/system/status', headers=self.admin_headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["status"], "healthy")
        self.assertIn("version", data)
        self.assertIn("uptime", data)
        self.assertIn("database_connection", data)
        self.assertIn("memory_usage", data)

    @patch('api.rbac.UserRoleManager.get_user_roles')
    @patch('api.system_resources.get_user_id')
    def test_system_status_endpoint_regular_user(self, mock_get_user_id, mock_get_user_roles):
        """Test GET /api/v1/system/status with regular user (should be restricted)."""
        # Mock regular user
        mock_get_user_id.return_value = self.regular_user_id
        mock_get_user_roles.return_value = [{"name": "user"}]
        
        # Make request
        response = self.client.get('/api/v1/system/status', headers=self.regular_headers)
        
        # Check response - should be forbidden for non-admin users
        self.assertIn(response.status_code, [401, 403])

    def test_system_status_endpoint_missing_auth(self):
        """Test GET /api/v1/system/status without authentication."""
        response = self.client.get('/api/v1/system/status')
        self.assertIn(response.status_code, [401, 403])

    @patch('api.system_resources.get_system_logs')
    @patch('api.rbac.UserRoleManager.get_user_roles')
    @patch('api.system_resources.get_user_id')
    def test_system_logs_endpoint_admin(self, mock_get_user_id, mock_get_user_roles, mock_get_system_logs):
        """Test GET /api/v1/system/logs with admin user."""
        # Mock admin user
        mock_get_user_id.return_value = self.admin_user_id
        mock_get_user_roles.return_value = [{"name": "admin"}]
        
        # Mock system logs
        mock_get_system_logs.return_value = {
            "logs": [
                {"timestamp": "2025-04-20T12:00:00", "level": "INFO", "message": "System started"},
                {"timestamp": "2025-04-20T12:01:00", "level": "INFO", "message": "Database connected"},
                {"timestamp": "2025-04-20T12:30:00", "level": "WARNING", "message": "High memory usage"}
            ],
            "total_count": 3
        }
        
        # Make request
        response = self.client.get('/api/v1/system/logs', headers=self.admin_headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn("logs", data)
        self.assertEqual(len(data["logs"]), 3)
        self.assertEqual(data["total_count"], 3)

    @patch('api.system_resources.get_system_logs')
    @patch('api.rbac.UserRoleManager.get_user_roles')
    @patch('api.system_resources.get_user_id')
    def test_system_logs_endpoint_with_filters(self, mock_get_user_id, mock_get_user_roles, mock_get_system_logs):
        """Test GET /api/v1/system/logs with filters."""
        # Mock admin user
        mock_get_user_id.return_value = self.admin_user_id
        mock_get_user_roles.return_value = [{"name": "admin"}]
        
        # Mock system logs
        mock_get_system_logs.return_value = {
            "logs": [
                {"timestamp": "2025-04-20T12:30:00", "level": "WARNING", "message": "High memory usage"}
            ],
            "total_count": 1
        }
        
        # Make request with filters
        response = self.client.get('/api/v1/system/logs?level=WARNING&limit=10&start_date=2025-04-20', 
                                  headers=self.admin_headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn("logs", data)
        self.assertEqual(len(data["logs"]), 1)
        self.assertEqual(data["logs"][0]["level"], "WARNING")
        
        # Verify that the function was called with the correct parameters
        mock_get_system_logs.assert_called_once()
        args, kwargs = mock_get_system_logs.call_args
        self.assertEqual(kwargs.get("level"), "WARNING")
        self.assertEqual(kwargs.get("limit"), "10")
        self.assertEqual(kwargs.get("start_date"), "2025-04-20")

    @patch('api.system_resources.get_system_metrics')
    @patch('api.rbac.UserRoleManager.get_user_roles')
    @patch('api.system_resources.get_user_id')
    def test_system_metrics_endpoint_admin(self, mock_get_user_id, mock_get_user_roles, mock_get_system_metrics):
        """Test GET /api/v1/system/metrics with admin user."""
        # Mock admin user
        mock_get_user_id.return_value = self.admin_user_id
        mock_get_user_roles.return_value = [{"name": "admin"}]
        
        # Mock system metrics
        mock_get_system_metrics.return_value = {
            "system": {
                "cpu_usage": [
                    {"timestamp": "2025-04-20T12:00:00", "value": 25},
                    {"timestamp": "2025-04-20T12:15:00", "value": 32}
                ],
                "memory_usage": [
                    {"timestamp": "2025-04-20T12:00:00", "value": 40},
                    {"timestamp": "2025-04-20T12:15:00", "value": 45}
                ]
            },
            "database": {
                "connections": [
                    {"timestamp": "2025-04-20T12:00:00", "value": 5},
                    {"timestamp": "2025-04-20T12:15:00", "value": 8}
                ],
                "query_time": [
                    {"timestamp": "2025-04-20T12:00:00", "value": 0.05},
                    {"timestamp": "2025-04-20T12:15:00", "value": 0.06}
                ]
            },
            "api": {
                "requests_per_minute": [
                    {"timestamp": "2025-04-20T12:00:00", "value": 12},
                    {"timestamp": "2025-04-20T12:15:00", "value": 18}
                ],
                "average_response_time": [
                    {"timestamp": "2025-04-20T12:00:00", "value": 0.12},
                    {"timestamp": "2025-04-20T12:15:00", "value": 0.15}
                ]
            }
        }
        
        # Make request
        response = self.client.get('/api/v1/system/metrics', headers=self.admin_headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn("system", data)
        self.assertIn("database", data)
        self.assertIn("api", data)
        self.assertIn("cpu_usage", data["system"])
        self.assertEqual(len(data["system"]["cpu_usage"]), 2)

    @patch('api.system_resources.get_system_metrics')
    @patch('api.rbac.UserRoleManager.get_user_roles')
    @patch('api.system_resources.get_user_id')
    def test_system_metrics_endpoint_with_timeframe(self, mock_get_user_id, mock_get_user_roles, mock_get_system_metrics):
        """Test GET /api/v1/system/metrics with timeframe parameter."""
        # Mock admin user
        mock_get_user_id.return_value = self.admin_user_id
        mock_get_user_roles.return_value = [{"name": "admin"}]
        
        # Mock system metrics
        mock_get_system_metrics.return_value = {
            "system": {
                "cpu_usage": [
                    {"timestamp": "2025-04-20T12:00:00", "value": 25}
                ]
            },
            "database": {
                "connections": [
                    {"timestamp": "2025-04-20T12:00:00", "value": 5}
                ]
            },
            "api": {
                "requests_per_minute": [
                    {"timestamp": "2025-04-20T12:00:00", "value": 12}
                ]
            }
        }
        
        # Make request with timeframe
        response = self.client.get('/api/v1/system/metrics?timeframe=hour', headers=self.admin_headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        
        # Verify that the function was called with the correct parameters
        mock_get_system_metrics.assert_called_once()
        args, kwargs = mock_get_system_metrics.call_args
        self.assertEqual(kwargs.get("timeframe"), "hour")

if __name__ == '__main__':
    unittest.main()