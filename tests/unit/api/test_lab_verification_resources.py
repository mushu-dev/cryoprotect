"""
Unit tests for the LabVerificationResource API endpoints.

This module contains tests for the LabVerificationResource API endpoints, including:
- GET /experiments/{experiment_id}/lab_verification
- POST /experiments/{experiment_id}/lab_verification
- PUT /verifications/{verification_id}
"""

import json
import unittest
from unittest.mock import patch, MagicMock
from flask import url_for
from api.models import LabVerification
from tests.base_test_case import MockSupabaseBaseTestCase


class TestLabVerificationResource(MockSupabaseBaseTestCase):
    """Test case for LabVerificationResource API endpoints."""

    def setUp(self):
        """Set up test case."""
        super().setUp()
        
        # Sample data for testing
        self.experiment_id = "12345678-1234-5678-1234-567812345678"
        self.verification_id = "87654321-4321-8765-4321-876543210987"
        self.verifier = "Lab Technician 1"
        self.equipment_used = "Microscope Model XYZ-123"
        self.comments = "Sample appears to match expected results."
        
        # Mock verification data
        self.verification_data = {
            "id": self.verification_id,
            "experiment_id": self.experiment_id,
            "verification_status": LabVerification.VERIFIED,
            "verifier": self.verifier,
            "equipment_used": self.equipment_used,
            "comments": self.comments,
            "created_at": "2025-04-23T10:00:00Z",
            "updated_at": "2025-04-23T10:00:00Z"
        }
        
        # Mock auth token
        self.auth_headers = {
            'Authorization': 'Bearer test_token'
        }

    @patch('api.lab_verification_resources.LabVerification.get_verification')
    def test_get_verification(self, mock_get_verification):
        """Test GET /experiments/{experiment_id}/lab_verification endpoint."""
        # Mock the model method
        mock_get_verification.return_value = self.verification_data
        
        # Make the request
        response = self.client.get(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['id'], self.verification_id)
        self.assertEqual(data['experiment_id'], self.experiment_id)
        self.assertEqual(data['verification_status'], LabVerification.VERIFIED)
        self.assertEqual(data['verifier'], self.verifier)
        self.assertEqual(data['equipment_used'], self.equipment_used)
        self.assertEqual(data['comments'], self.comments)
        
        # Verify the model method was called correctly
        mock_get_verification.assert_called_once_with(self.experiment_id)

    @patch('api.lab_verification_resources.LabVerification.get_verification')
    def test_get_verification_not_found(self, mock_get_verification):
        """Test GET endpoint when verification is not found."""
        # Mock the model method to return None (not found)
        mock_get_verification.return_value = None
        
        # Make the request
        response = self.client.get(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 404)
        
        # Verify the model method was called correctly
        mock_get_verification.assert_called_once_with(self.experiment_id)

    @patch('api.lab_verification_resources.LabVerification.record_verification')
    def test_post_verification(self, mock_record_verification):
        """Test POST /experiments/{experiment_id}/lab_verification endpoint."""
        # Mock the model method
        mock_record_verification.return_value = self.verification_data
        
        # Request payload
        payload = {
            'verification_status': LabVerification.VERIFIED,
            'verifier': self.verifier,
            'equipment_used': self.equipment_used,
            'comments': self.comments
        }
        
        # Make the request
        response = self.client.post(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers,
            json=payload
        )
        
        # Assertions
        self.assertEqual(response.status_code, 201)
        data = json.loads(response.data)
        self.assertEqual(data['id'], self.verification_id)
        self.assertEqual(data['experiment_id'], self.experiment_id)
        self.assertEqual(data['verification_status'], LabVerification.VERIFIED)
        self.assertEqual(data['verifier'], self.verifier)
        self.assertEqual(data['equipment_used'], self.equipment_used)
        self.assertEqual(data['comments'], self.comments)
        
        # Verify the model method was called correctly
        mock_record_verification.assert_called_once_with(
            experiment_id=self.experiment_id,
            verification_status=LabVerification.VERIFIED,
            verifier=self.verifier,
            equipment_used=self.equipment_used,
            comments=self.comments
        )

    @patch('api.lab_verification_resources.LabVerification.record_verification')
    def test_post_verification_missing_required_fields(self, mock_record_verification):
        """Test POST endpoint with missing required fields."""
        # Mock the model method to raise ValueError
        mock_record_verification.side_effect = ValueError("Required fields missing")
        
        # Request payload with missing fields
        payload = {
            'verification_status': LabVerification.VERIFIED,
            # Missing 'verifier'
            'equipment_used': self.equipment_used
        }
        
        # Make the request
        response = self.client.post(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers,
            json=payload
        )
        
        # Assertions
        self.assertEqual(response.status_code, 400)

    @patch('api.lab_verification_resources.LabVerification.update_verification_status')
    def test_put_verification_status(self, mock_update_verification_status):
        """Test PUT /verifications/{verification_id} endpoint."""
        # Updated data
        updated_status = LabVerification.REJECTED
        updated_comments = "Sample does not meet quality standards."
        
        updated_data = self.verification_data.copy()
        updated_data["verification_status"] = updated_status
        updated_data["comments"] = updated_comments
        
        # Mock the model method
        mock_update_verification_status.return_value = updated_data
        
        # Request payload
        payload = {
            'verification_status': updated_status,
            'comments': updated_comments
        }
        
        # Make the request
        response = self.client.put(
            f'/api/verifications/{self.verification_id}',
            headers=self.auth_headers,
            json=payload
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['id'], self.verification_id)
        self.assertEqual(data['verification_status'], updated_status)
        self.assertEqual(data['comments'], updated_comments)
        
        # Verify the model method was called correctly
        mock_update_verification_status.assert_called_once_with(
            verification_id=self.verification_id,
            new_status=updated_status,
            comments=updated_comments
        )

    @patch('api.lab_verification_resources.LabVerification.update_verification_status')
    def test_put_verification_status_invalid_status(self, mock_update_verification_status):
        """Test PUT endpoint with invalid status."""
        # Mock the model method to raise ValueError
        mock_update_verification_status.side_effect = ValueError("Invalid verification status")
        
        # Request payload with invalid status
        payload = {
            'verification_status': 'invalid_status',
            'comments': self.comments
        }
        
        # Make the request
        response = self.client.put(
            f'/api/verifications/{self.verification_id}',
            headers=self.auth_headers,
            json=payload
        )
        
        # Assertions
        self.assertEqual(response.status_code, 400)


if __name__ == "__main__":
    unittest.main()