"""
Unit tests for the LabVerification model.

This module contains tests for the LabVerification model methods, including:
- record_verification
- get_verification
- update_verification_status
"""

import unittest
from unittest.mock import patch, MagicMock
from api.models import LabVerification
from tests.base_test_case import MockSupabaseBaseTestCase


class TestLabVerificationModel(MockSupabaseBaseTestCase):
    """Test case for LabVerification model."""

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

    def test_record_verification_success(self):
        """Test successful verification recording."""
        # Mock the Supabase response
        mock_response = MagicMock()
        mock_response.data = [self.verification_data]
        mock_response.error = None
        
        # Set up the mock for the insert method
        mock_table = self.mock_supabase.table.return_value
        mock_table.insert.return_value.execute.return_value = mock_response
        
        # Call the method
        result = LabVerification.record_verification(
            experiment_id=self.experiment_id,
            verification_status=LabVerification.VERIFIED,
            verifier=self.verifier,
            equipment_used=self.equipment_used,
            comments=self.comments
        )
        
        # Assertions
        self.assertEqual(result, self.verification_data)
        mock_table.insert.assert_called_once()
        
        # Verify the data passed to insert
        call_args = mock_table.insert.call_args[0][0]
        self.assertEqual(call_args["experiment_id"], self.experiment_id)
        self.assertEqual(call_args["verification_status"], LabVerification.VERIFIED)
        self.assertEqual(call_args["verifier"], self.verifier)
        self.assertEqual(call_args["equipment_used"], self.equipment_used)
        self.assertEqual(call_args["comments"], self.comments)

    def test_record_verification_invalid_status(self):
        """Test verification recording with invalid status."""
        # Test with invalid status
        with self.assertRaises(ValueError) as context:
            LabVerification.record_verification(
                experiment_id=self.experiment_id,
                verification_status="invalid_status",
                verifier=self.verifier,
                equipment_used=self.equipment_used
            )
        
        self.assertIn("Invalid verification status", str(context.exception))

    def test_record_verification_missing_required_fields(self):
        """Test verification recording with missing required fields."""
        # Test with missing experiment_id
        with self.assertRaises(ValueError) as context:
            LabVerification.record_verification(
                experiment_id="",
                verification_status=LabVerification.VERIFIED,
                verifier=self.verifier,
                equipment_used=self.equipment_used
            )
        
        self.assertIn("experiment_id, verifier, and equipment_used are required", str(context.exception))
        
        # Test with missing verifier
        with self.assertRaises(ValueError) as context:
            LabVerification.record_verification(
                experiment_id=self.experiment_id,
                verification_status=LabVerification.VERIFIED,
                verifier="",
                equipment_used=self.equipment_used
            )
        
        self.assertIn("experiment_id, verifier, and equipment_used are required", str(context.exception))
        
        # Test with missing equipment_used
        with self.assertRaises(ValueError) as context:
            LabVerification.record_verification(
                experiment_id=self.experiment_id,
                verification_status=LabVerification.VERIFIED,
                verifier=self.verifier,
                equipment_used=""
            )
        
        self.assertIn("experiment_id, verifier, and equipment_used are required", str(context.exception))

    def test_get_verification_success(self):
        """Test successful retrieval of verification."""
        # Mock the Supabase response
        mock_response = MagicMock()
        mock_response.data = [self.verification_data]
        mock_response.error = None
        
        # Set up the mock for the select method
        mock_table = self.mock_supabase.table.return_value
        mock_select = mock_table.select.return_value
        mock_eq = mock_select.eq.return_value
        mock_eq.execute.return_value = mock_response
        
        # Call the method
        result = LabVerification.get_verification(experiment_id=self.experiment_id)
        
        # Assertions
        self.assertEqual(result, self.verification_data)
        mock_table.select.assert_called_once_with("*")
        mock_select.eq.assert_called_once_with("experiment_id", self.experiment_id)
        mock_eq.execute.assert_called_once()

    def test_get_verification_not_found(self):
        """Test retrieval of non-existent verification."""
        # Mock the Supabase response for not found
        mock_response = MagicMock()
        mock_response.data = []
        mock_response.error = None
        
        # Set up the mock for the select method
        mock_table = self.mock_supabase.table.return_value
        mock_select = mock_table.select.return_value
        mock_eq = mock_select.eq.return_value
        mock_eq.execute.return_value = mock_response
        
        # Call the method
        result = LabVerification.get_verification(experiment_id=self.experiment_id)
        
        # Assertions
        self.assertIsNone(result)
        mock_table.select.assert_called_once_with("*")
        mock_select.eq.assert_called_once_with("experiment_id", self.experiment_id)
        mock_eq.execute.assert_called_once()

    def test_get_verification_missing_id(self):
        """Test retrieval with missing experiment_id."""
        # Test with missing experiment_id
        with self.assertRaises(ValueError) as context:
            LabVerification.get_verification(experiment_id="")
        
        self.assertIn("experiment_id is required", str(context.exception))

    def test_update_verification_status_success(self):
        """Test successful update of verification status."""
        # Updated data
        updated_status = LabVerification.REJECTED
        updated_comments = "Sample does not meet quality standards."
        
        updated_data = self.verification_data.copy()
        updated_data["verification_status"] = updated_status
        updated_data["comments"] = updated_comments
        
        # Mock the Supabase response
        mock_response = MagicMock()
        mock_response.data = [updated_data]
        mock_response.error = None
        
        # Set up the mock for the update method
        mock_table = self.mock_supabase.table.return_value
        mock_update = mock_table.update.return_value
        mock_eq = mock_update.eq.return_value
        mock_eq.execute.return_value = mock_response
        
        # Call the method
        result = LabVerification.update_verification_status(
            verification_id=self.verification_id,
            new_status=updated_status,
            comments=updated_comments
        )
        
        # Assertions
        self.assertEqual(result, updated_data)
        mock_table.update.assert_called_once()
        
        # Verify the data passed to update
        call_args = mock_table.update.call_args[0][0]
        self.assertEqual(call_args["verification_status"], updated_status)
        self.assertEqual(call_args["comments"], updated_comments)
        
        mock_update.eq.assert_called_once_with("id", self.verification_id)
        mock_eq.execute.assert_called_once()

    def test_update_verification_status_invalid_status(self):
        """Test update with invalid status."""
        # Test with invalid status
        with self.assertRaises(ValueError) as context:
            LabVerification.update_verification_status(
                verification_id=self.verification_id,
                new_status="invalid_status"
            )
        
        self.assertIn("Invalid new_status value", str(context.exception))

    def test_update_verification_status_missing_id(self):
        """Test update with missing verification_id."""
        # Test with missing verification_id
        with self.assertRaises(ValueError) as context:
            LabVerification.update_verification_status(
                verification_id="",
                new_status=LabVerification.PENDING
            )
        
        self.assertIn("verification_id is required", str(context.exception))


if __name__ == "__main__":
    unittest.main()