"""
Integration tests for the lab verification workflow.

This module contains tests that verify the complete lab verification workflow,
including database operations, API endpoints, and their interactions.
"""

import json
import unittest
from unittest.mock import patch, MagicMock
from api.models import LabVerification, Experiment
from tests.base_test_case import MockSupabaseBaseTestCase


class TestLabVerificationWorkflow(MockSupabaseBaseTestCase):
    """Test case for the lab verification workflow."""

    def setUp(self):
        """Set up test case."""
        super().setUp()
        
        # Sample data for testing
        self.experiment_id = "12345678-1234-5678-1234-567812345678"
        self.verification_id = "87654321-4321-8765-4321-876543210987"
        self.verifier = "Lab Technician 1"
        self.equipment_used = "Microscope Model XYZ-123"
        self.comments = "Sample appears to match expected results."
        
        # Mock experiment data
        self.experiment_data = {
            "id": self.experiment_id,
            "mixture_id": "98765432-9876-5432-9876-987654321098",
            "property_name": "Cell Viability",
            "numeric_value": 85.5,
            "experimental_conditions": "Standard laboratory conditions",
            "date_performed": "2025-04-20",
            "created_at": "2025-04-20T14:30:00Z"
        }
        
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

    @patch('api.models.Experiment.get')
    def test_experiment_exists_before_verification(self, mock_experiment_get):
        """Test that an experiment exists before attempting verification."""
        # Mock the Experiment.get method
        mock_experiment_get.return_value = self.experiment_data
        
        # Verify the experiment exists
        experiment = Experiment.get(self.experiment_id)
        self.assertIsNotNone(experiment)
        self.assertEqual(experiment["id"], self.experiment_id)
        
        # Verify the mock was called correctly
        mock_experiment_get.assert_called_once_with(self.experiment_id)

    @patch('api.models.LabVerification.record_verification')
    @patch('api.models.LabVerification.get_verification')
    def test_full_verification_workflow(self, mock_get_verification, mock_record_verification):
        """Test the complete lab verification workflow."""
        # 1. Initially, no verification exists
        mock_get_verification.return_value = None
        
        # Make GET request to check if verification exists
        response = self.client.get(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers
        )
        
        # Verify response (should be 404 Not Found)
        self.assertEqual(response.status_code, 404)
        
        # 2. Create a new verification
        mock_record_verification.return_value = self.verification_data
        
        # Request payload
        payload = {
            'verification_status': LabVerification.VERIFIED,
            'verifier': self.verifier,
            'equipment_used': self.equipment_used,
            'comments': self.comments
        }
        
        # Make POST request to create verification
        response = self.client.post(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers,
            json=payload
        )
        
        # Verify response
        self.assertEqual(response.status_code, 201)
        data = json.loads(response.data)
        self.assertEqual(data['id'], self.verification_id)
        self.assertEqual(data['verification_status'], LabVerification.VERIFIED)
        
        # 3. Now verification exists
        mock_get_verification.return_value = self.verification_data
        
        # Make GET request again
        response = self.client.get(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers
        )
        
        # Verify response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['id'], self.verification_id)
        self.assertEqual(data['verification_status'], LabVerification.VERIFIED)

    @patch('api.models.LabVerification.update_verification_status')
    @patch('api.models.LabVerification.get_verification')
    @patch('api.models.LabVerification.record_verification')
    def test_verification_status_update_workflow(self, mock_record_verification, 
                                               mock_get_verification, 
                                               mock_update_verification_status):
        """Test the workflow for updating verification status."""
        # 1. Create initial verification with PENDING status
        initial_data = self.verification_data.copy()
        initial_data["verification_status"] = LabVerification.PENDING
        mock_record_verification.return_value = initial_data
        
        # Request payload
        payload = {
            'verification_status': LabVerification.PENDING,
            'verifier': self.verifier,
            'equipment_used': self.equipment_used,
            'comments': "Initial verification, pending final review."
        }
        
        # Make POST request to create verification
        response = self.client.post(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers,
            json=payload
        )
        
        # Verify response
        self.assertEqual(response.status_code, 201)
        data = json.loads(response.data)
        self.assertEqual(data['verification_status'], LabVerification.PENDING)
        
        # 2. Verification exists with PENDING status
        mock_get_verification.return_value = initial_data
        
        # Make GET request
        response = self.client.get(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers
        )
        
        # Verify response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['verification_status'], LabVerification.PENDING)
        
        # 3. Update verification status to VERIFIED
        updated_data = initial_data.copy()
        updated_data["verification_status"] = LabVerification.VERIFIED
        updated_data["comments"] = "Verification complete, all criteria met."
        mock_update_verification_status.return_value = updated_data
        
        # Request payload for update
        update_payload = {
            'verification_status': LabVerification.VERIFIED,
            'comments': "Verification complete, all criteria met."
        }
        
        # Make PUT request to update verification
        response = self.client.put(
            f'/api/verifications/{self.verification_id}',
            headers=self.auth_headers,
            json=update_payload
        )
        
        # Verify response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['verification_status'], LabVerification.VERIFIED)
        self.assertEqual(data['comments'], "Verification complete, all criteria met.")
        
        # Verify the model method was called correctly
        mock_update_verification_status.assert_called_once_with(
            verification_id=self.verification_id,
            new_status=LabVerification.VERIFIED,
            comments="Verification complete, all criteria met."
        )

    @patch('api.models.LabVerification.record_verification')
    def test_verification_with_invalid_data(self, mock_record_verification):
        """Test verification with invalid data."""
        # Mock the model method to raise ValueError
        mock_record_verification.side_effect = ValueError("Invalid verification status")
        
        # Request payload with invalid status
        payload = {
            'verification_status': 'invalid_status',
            'verifier': self.verifier,
            'equipment_used': self.equipment_used
        }
        
        # Make POST request
        response = self.client.post(
            f'/api/experiments/{self.experiment_id}/lab_verification',
            headers=self.auth_headers,
            json=payload
        )
        
        # Verify response
        self.assertEqual(response.status_code, 400)
        data = json.loads(response.data)
        self.assertIn('error', data)


if __name__ == "__main__":
    unittest.main()