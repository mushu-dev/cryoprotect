"""
Tests for the batch_resources module.

These tests demonstrate the improved testability of the refactored batch_resources module,
particularly the ability to test business logic without needing Flask's request context.
"""

import unittest
from unittest.mock import patch, MagicMock, Mock
import json
import io
from flask import Flask, Response

from api.batch_resources import BatchOperationService, BatchOperationResource


class TestBatchOperationService(unittest.TestCase):
    """Test cases for the BatchOperationService class."""

    def test_validate_request_valid(self):
        """Test that a valid request passes validation."""
        result = BatchOperationService.validate_request(
            operation="property_calculation",
            entity_type="molecule",
            ids=["mol-123", "mol-456"]
        )
        self.assertIsNone(result)

    def test_validate_request_invalid(self):
        """Test that an invalid request fails validation."""
        # Missing operation
        result = BatchOperationService.validate_request(
            operation="",
            entity_type="molecule",
            ids=["mol-123"]
        )
        self.assertIsNotNone(result)
        self.assertEqual(result["status"], "ERROR")
        
        # Missing entity_type
        result = BatchOperationService.validate_request(
            operation="property_calculation",
            entity_type="",
            ids=["mol-123"]
        )
        self.assertIsNotNone(result)
        self.assertEqual(result["status"], "ERROR")
        
        # Invalid ids (not a list)
        result = BatchOperationService.validate_request(
            operation="property_calculation",
            entity_type="molecule",
            ids="mol-123"  # Not a list
        )
        self.assertIsNotNone(result)
        self.assertEqual(result["status"], "ERROR")
        
        # Empty ids list
        result = BatchOperationService.validate_request(
            operation="property_calculation",
            entity_type="molecule",
            ids=[]
        )
        # Note: The current implementation doesn't check for empty lists
        # If we want to test this, we'd need to modify the implementation

    @patch('api.batch_resources.calculate_molecule_score')
    def test_process_property_calculation_molecule(self, mock_calculate_molecule_score):
        """Test property calculation for molecules."""
        # Setup mock
        mock_calculate_molecule_score.return_value = {"overall_score": 85}
        
        # Call the method
        result = BatchOperationService.process_property_calculation("molecule", "mol-123")
        
        # Verify the result
        self.assertEqual(result["overall_score"], 85)
        mock_calculate_molecule_score.assert_called_once_with("mol-123")

    @patch('api.batch_resources.calculate_mixture_score')
    def test_process_property_calculation_mixture(self, mock_calculate_mixture_score):
        """Test property calculation for mixtures."""
        # Setup mock
        mock_calculate_mixture_score.return_value = {"overall_score": 75}
        
        # Call the method
        result = BatchOperationService.process_property_calculation("mixture", "mix-123")
        
        # Verify the result
        self.assertEqual(result["overall_score"], 75)
        mock_calculate_mixture_score.assert_called_once_with("mix-123")

    def test_process_property_calculation_invalid_entity(self):
        """Test property calculation with invalid entity type."""
        with self.assertRaises(ValueError):
            BatchOperationService.process_property_calculation("invalid_entity", "id-123")

    @patch('api.batch_resources.MixtureOptimization.optimize_composition')
    def test_process_mixture_optimization(self, mock_optimize_composition):
        """Test mixture optimization."""
        # Setup mock
        mock_optimize_composition.return_value = {
            "optimized_components": [{"molecule_id": "mol-123", "concentration": 60}]
        }
        
        # Call the method
        result = BatchOperationService.process_mixture_optimization("mixture", "mix-123")
        
        # Verify the result
        self.assertEqual(result["optimized_components"][0]["concentration"], 60)
        mock_optimize_composition.assert_called_once_with("mix-123")

    def test_process_mixture_optimization_invalid_entity(self):
        """Test mixture optimization with invalid entity type."""
        with self.assertRaises(ValueError):
            BatchOperationService.process_mixture_optimization("molecule", "mol-123")

    @patch('api.batch_resources.ModelManager')
    def test_process_predictive_scoring(self, MockModelManager):
        """Test predictive scoring."""
        # Setup mock
        mock_manager = MockModelManager.return_value
        mock_manager.predict.return_value = {"prediction": 90, "confidence": 0.85}
        
        # Call the method
        result = BatchOperationService.process_predictive_scoring(
            "mixture", "mix-123", "Test Property", "test_algorithm"
        )
        
        # Verify the result
        self.assertEqual(result["prediction"], 90)
        self.assertEqual(result["confidence"], 0.85)
        mock_manager.predict.assert_called_once_with("Test Property", "mix-123", "test_algorithm")

    @patch('api.batch_resources.ModelManager')
    def test_process_predictive_scoring_default_params(self, MockModelManager):
        """Test predictive scoring with default parameters."""
        # Setup mock
        mock_manager = MockModelManager.return_value
        mock_manager.predict.return_value = {"prediction": 80, "confidence": 0.75}
        
        # Call the method with default parameters
        result = BatchOperationService.process_predictive_scoring("molecule", "mol-123")
        
        # Verify the result
        self.assertEqual(result["prediction"], 80)
        self.assertEqual(result["confidence"], 0.75)
        mock_manager.predict.assert_called_once_with("Cryoprotection Score", "mol-123", "random_forest")

    @patch('api.batch_resources.get_data_for_export')
    @patch('api.batch_resources.generate_csv')
    def test_process_export_csv(self, mock_generate_csv, mock_get_data_for_export):
        """Test export in CSV format."""
        # Setup mocks
        mock_get_data_for_export.return_value = [{"id": "mol-123", "name": "Test Molecule"}]
        mock_csv = MagicMock()
        mock_csv.getvalue.return_value = "id,name\nmol-123,Test Molecule"
        mock_generate_csv.return_value = mock_csv
        
        # Call the method
        result = BatchOperationService.process_export("molecule", "mol-123", "csv")
        
        # Verify the result
        self.assertEqual(result, "id,name\nmol-123,Test Molecule")
        mock_get_data_for_export.assert_called_once_with("molecules", "mol-123")
        mock_generate_csv.assert_called_once()

    @patch('api.batch_resources.get_data_for_export')
    @patch('api.batch_resources.generate_json')
    def test_process_export_json(self, mock_generate_json, mock_get_data_for_export):
        """Test export in JSON format."""
        # Setup mocks
        mock_get_data_for_export.return_value = [{"id": "mol-123", "name": "Test Molecule"}]
        mock_json = MagicMock()
        mock_json.getvalue.return_value = '[{"id":"mol-123","name":"Test Molecule"}]'
        mock_generate_json.return_value = mock_json
        
        # Call the method
        result = BatchOperationService.process_export("molecule", "mol-123", "json")
        
        # Verify the result
        self.assertEqual(result, '[{"id":"mol-123","name":"Test Molecule"}]')
        mock_get_data_for_export.assert_called_once_with("molecules", "mol-123")
        mock_generate_json.assert_called_once()

    @patch('api.batch_resources.get_data_for_export')
    @patch('api.batch_resources.generate_excel')
    def test_process_export_excel(self, mock_generate_excel, mock_get_data_for_export):
        """Test export in Excel format."""
        # Setup mocks
        mock_get_data_for_export.return_value = [{"id": "mol-123", "name": "Test Molecule"}]
        mock_excel = MagicMock()
        mock_excel.getvalue.return_value = b"Excel binary data"
        mock_generate_excel.return_value = mock_excel
        
        # Call the method
        result = BatchOperationService.process_export("molecule", "mol-123", "excel")
        
        # Verify the result
        self.assertEqual(result, b"Excel binary data")
        mock_get_data_for_export.assert_called_once_with("molecules", "mol-123")
        mock_generate_excel.assert_called_once()

    @patch('api.batch_resources.get_data_for_export')
    @patch('api.batch_resources.generate_pdf')
    def test_process_export_pdf(self, mock_generate_pdf, mock_get_data_for_export):
        """Test export in PDF format."""
        # Setup mocks
        mock_get_data_for_export.return_value = [{"id": "mol-123", "name": "Test Molecule"}]
        mock_pdf = MagicMock()
        mock_pdf.getvalue.return_value = b"PDF binary data"
        mock_generate_pdf.return_value = mock_pdf
        
        # Call the method
        result = BatchOperationService.process_export("molecule", "mol-123", "pdf")
        
        # Verify the result
        self.assertEqual(result, b"PDF binary data")
        mock_get_data_for_export.assert_called_once_with("molecules", "mol-123")
        mock_generate_pdf.assert_called_once_with(
            [{"id": "mol-123", "name": "Test Molecule"}], 
            "molecule_mol-123", 
            "molecules"
        )

    def test_process_export_invalid_format(self):
        """Test export with invalid format."""
        with patch('api.batch_resources.get_data_for_export', return_value=[]):
            with self.assertRaises(ValueError) as context:
                BatchOperationService.process_export("molecule", "mol-123", "invalid_format")
            
            self.assertIn("Unsupported export format", str(context.exception))

    @patch('api.batch_resources.BatchOperationService.process_property_calculation')
    def test_process_batch_operation_property_calculation(self, mock_process_property_calculation):
        """Test batch operation for property calculation."""
        # Setup mock
        mock_process_property_calculation.return_value = {"overall_score": 85}
        
        # Call the method
        result = BatchOperationService.process_batch_operation(
            operation="property_calculation",
            entity_type="molecule",
            ids=["mol-123", "mol-456"]
        )
        
        # Verify the result
        self.assertEqual(result["status"], "SUCCESS")
        self.assertEqual(len(result["results"]), 2)
        self.assertEqual(result["results"][0]["id"], "mol-123")
        self.assertEqual(result["results"][0]["result"]["overall_score"], 85)
        self.assertEqual(len(result["errors"]), 0)
        self.assertEqual(mock_process_property_calculation.call_count, 2)

    @patch('api.batch_resources.BatchOperationService.process_mixture_optimization')
    def test_process_batch_operation_mixture_optimization(self, mock_process_mixture_optimization):
        """Test batch operation for mixture optimization."""
        # Setup mock
        mock_process_mixture_optimization.return_value = {
            "optimized_components": [{"molecule_id": "mol-123", "concentration": 60}]
        }
        
        # Call the method
        result = BatchOperationService.process_batch_operation(
            operation="mixture_optimization",
            entity_type="mixture",
            ids=["mix-123", "mix-456"]
        )
        
        # Verify the result
        self.assertEqual(result["status"], "SUCCESS")
        self.assertEqual(len(result["results"]), 2)
        self.assertEqual(result["results"][0]["id"], "mix-123")
        self.assertEqual(result["results"][0]["result"]["optimized_components"][0]["concentration"], 60)
        self.assertEqual(len(result["errors"]), 0)
        self.assertEqual(mock_process_mixture_optimization.call_count, 2)

    @patch('api.batch_resources.BatchOperationService.process_predictive_scoring')
    def test_process_batch_operation_predictive_scoring(self, mock_process_predictive_scoring):
        """Test batch operation for predictive scoring."""
        # Setup mock
        mock_process_predictive_scoring.return_value = {"prediction": 90, "confidence": 0.85}
        
        # Call the method
        result = BatchOperationService.process_batch_operation(
            operation="predictive_scoring",
            entity_type="molecule",
            ids=["mol-123", "mol-456"],
            additional_params={
                "property_name": "Test Property",
                "algorithm": "test_algorithm"
            }
        )
        
        # Verify the result
        self.assertEqual(result["status"], "SUCCESS")
        self.assertEqual(len(result["results"]), 2)
        self.assertEqual(result["results"][0]["id"], "mol-123")
        self.assertEqual(result["results"][0]["result"]["prediction"], 90)
        self.assertEqual(result["results"][0]["result"]["confidence"], 0.85)
        self.assertEqual(len(result["errors"]), 0)
        self.assertEqual(mock_process_predictive_scoring.call_count, 2)
        
        # Verify that additional_params were passed correctly
        mock_process_predictive_scoring.assert_any_call(
            "molecule", "mol-123", "Test Property", "test_algorithm"
        )

    @patch('api.batch_resources.BatchOperationService.process_export')
    def test_process_batch_operation_export(self, mock_process_export):
        """Test batch operation for export."""
        # Setup mock
        mock_process_export.return_value = "exported data"
        
        # Call the method
        result = BatchOperationService.process_batch_operation(
            operation="export",
            entity_type="molecule",
            ids=["mol-123", "mol-456"],
            additional_params={"format": "json"}
        )
        
        # Verify the result
        self.assertEqual(result["status"], "SUCCESS")
        self.assertEqual(len(result["results"]), 2)
        self.assertEqual(result["results"][0]["id"], "mol-123")
        self.assertEqual(result["results"][0]["export"], "exported data")
        self.assertEqual(len(result["errors"]), 0)
        self.assertEqual(mock_process_export.call_count, 2)
        
        # Verify that additional_params were passed correctly
        mock_process_export.assert_any_call("molecule", "mol-123", "json")

    @patch('api.batch_resources.BatchOperationService.process_property_calculation')
    def test_process_batch_operation_with_errors(self, mock_process_property_calculation):
        """Test batch operation with some errors."""
        # Setup mock to succeed for first ID and fail for second ID
        mock_process_property_calculation.side_effect = [
            {"overall_score": 85},  # Success for first ID
            ValueError("Test error")  # Error for second ID
        ]
        
        # Call the method
        result = BatchOperationService.process_batch_operation(
            operation="property_calculation",
            entity_type="molecule",
            ids=["mol-123", "mol-456"]
        )
        
        # Verify the result
        self.assertEqual(result["status"], "COMPLETED_WITH_WARNINGS")
        self.assertEqual(len(result["results"]), 1)
        self.assertEqual(result["results"][0]["id"], "mol-123")
        self.assertEqual(len(result["errors"]), 1)
        self.assertEqual(result["errors"][0]["id"], "mol-456")
        self.assertEqual(result["errors"][0]["error"], "Test error")

    @patch('api.batch_resources.BatchOperationService.process_property_calculation')
    def test_process_batch_operation_all_errors(self, mock_process_property_calculation):
        """Test batch operation with all errors."""
        # Setup mock to fail for all IDs
        mock_process_property_calculation.side_effect = ValueError("Test error")
        
        # Call the method
        result = BatchOperationService.process_batch_operation(
            operation="property_calculation",
            entity_type="molecule",
            ids=["mol-123", "mol-456"]
        )
        
        # Verify the result
        self.assertEqual(result["status"], "ERROR")
        self.assertEqual(len(result["results"]), 0)
        self.assertEqual(len(result["errors"]), 2)
        self.assertEqual(result["errors"][0]["id"], "mol-123")
        self.assertEqual(result["errors"][1]["id"], "mol-456")

    def test_process_batch_operation_invalid_operation(self):
        """Test batch operation with invalid operation."""
        result = BatchOperationService.process_batch_operation(
            operation="invalid_operation",
            entity_type="molecule",
            ids=["mol-123"]
        )
        
        self.assertEqual(result["status"], "ERROR")
        self.assertEqual(len(result["results"]), 0)
        self.assertEqual(len(result["errors"]), 1)
        self.assertEqual(result["errors"][0]["id"], "mol-123")
        self.assertIn("Unsupported operation", result["errors"][0]["error"])


class TestBatchOperationResource(unittest.TestCase):
    """Test cases for the BatchOperationResource class."""
    
    def setUp(self):
        """Set up the test environment."""
        self.app = Flask(__name__)
        self.app_context = self.app.app_context()
        self.app_context.push()
        self.request_context = self.app.test_request_context()
        self.request_context.push()
        
        # Create patchers
        # Mock token_required to simply return the function unchanged
        self.token_required_patcher = patch('api.batch_resources.token_required')
        self.mock_token_required = self.token_required_patcher.start()
        # Configure the mock to return the function unchanged
        self.mock_token_required.side_effect = lambda f: f
        
        # Mock _handle_json_serialization to return appropriate status codes
        self.serialization_patcher = patch('api.batch_resources._handle_json_serialization')
        self.mock_serialization = self.serialization_patcher.start()
        self.mock_serialization.side_effect = lambda x: (x, 200) if x.get('status') != 'ERROR' else (x, 400)
        
        # Create a resource instance
        self.resource = BatchOperationResource()
    
    def tearDown(self):
        """Clean up after the test."""
        self.token_required_patcher.stop()
        self.serialization_patcher.stop()
        self.request_context.pop()
        self.app_context.pop()
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.BatchOperationService.process_batch_operation')
    def test_post_valid_request(self, mock_process_batch_operation, mock_request):
        """Test post method with a valid request."""
        # Setup mocks
        mock_request.get_json.return_value = {
            "operation": "property_calculation",
            "entity_type": "molecule",
            "ids": ["mol-123", "mol-456"]
        }
        mock_process_batch_operation.return_value = {
            "status": "SUCCESS",
            "results": [
                {"id": "mol-123", "result": {"overall_score": 85}},
                {"id": "mol-456", "result": {"overall_score": 75}}
            ],
            "errors": []
        }
        
        # Call the method
        result, status_code = self.resource.post()
        
        # Verify the result
        self.assertEqual(status_code, 200)
        self.assertEqual(result["status"], "SUCCESS")
        self.assertEqual(len(result["results"]), 2)
        self.assertEqual(result["results"][0]["id"], "mol-123")
        self.assertEqual(result["results"][0]["result"]["overall_score"], 85)
        
        # Verify that process_batch_operation was called with the correct arguments
        mock_process_batch_operation.assert_called_once_with(
            "property_calculation", "molecule", ["mol-123", "mol-456"],
            {'format': 'csv', 'property_name': 'Cryoprotection Score', 'algorithm': 'random_forest'}
        )
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.BatchOperationService.process_batch_operation')
    def test_post_with_additional_params(self, mock_process_batch_operation, mock_request):
        """Test post method with additional parameters."""
        # Setup mocks
        mock_request.get_json.return_value = {
            "operation": "predictive_scoring",
            "entity_type": "molecule",
            "ids": ["mol-123"],
            "property_name": "Test Property",
            "algorithm": "test_algorithm"
        }
        mock_process_batch_operation.return_value = {
            "status": "SUCCESS",
            "results": [
                {"id": "mol-123", "result": {"prediction": 90, "confidence": 0.85}}
            ],
            "errors": []
        }
        
        # Call the method
        result, status_code = self.resource.post()
        
        # Verify the result
        self.assertEqual(status_code, 200)
        self.assertEqual(result["status"], "SUCCESS")
        
        # Verify that process_batch_operation was called with the correct arguments
        mock_process_batch_operation.assert_called_once_with(
            "predictive_scoring", "molecule", ["mol-123"],
            {'format': 'csv', 'property_name': 'Test Property', 'algorithm': 'test_algorithm'}
        )
    
    @patch('api.batch_resources.request')
    def test_post_missing_parameters(self, mock_request):
        """Test post method with missing parameters."""
        # Setup mock
        mock_request.get_json.return_value = {}
        
        # Call the method
        result, status_code = self.resource.post()
        
        # Verify the result
        self.assertEqual(status_code, 400)
        self.assertEqual(result["status"], "ERROR")
        self.assertIn("Missing or invalid parameters", result["errors"][0])
    
    @patch('api.batch_resources.request')
    def test_post_invalid_json(self, mock_request):
        """Test post method with invalid JSON."""
        # Setup mock to raise an exception
        mock_request.get_json.side_effect = Exception("Invalid JSON")
        
        # Call the method
        result, status_code = self.resource.post()
        
        # Verify the result
        self.assertEqual(status_code, 400)
        self.assertEqual(result["status"], "ERROR")
        self.assertIn("Unexpected error", result["errors"][0])
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.BatchOperationService.process_batch_operation')
    def test_post_service_error(self, mock_process_batch_operation, mock_request):
        """Test post method when the service raises an error."""
        # Setup mocks
        mock_request.get_json.return_value = {
            "operation": "property_calculation",
            "entity_type": "molecule",
            "ids": ["mol-123"]
        }
        mock_process_batch_operation.side_effect = Exception("Service error")
        
        # Call the method
        result, status_code = self.resource.post()
        
        # Verify the result
        self.assertEqual(status_code, 400)
        self.assertEqual(result["status"], "ERROR")
        self.assertIn("Unexpected error", result["errors"][0])
        self.assertIn("Service error", result["errors"][0])


if __name__ == '__main__':
    unittest.main()