import unittest
from unittest.mock import patch, MagicMock
import sys
import os
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flask import Flask, request, current_app
from api.batch_resources import BatchOperationResource

class TestBatchOperationResourceDirect(unittest.TestCase):
    """
    Direct tests for BatchOperationResource that bypass Flask test client
    to avoid serialization issues and measure code coverage.
    """
    
    @patch('api.batch_resources.token_required', lambda f: f)
    @patch('api.batch_resources._handle_json_serialization', lambda d: d)
    def setUp(self):
        self.app = Flask(__name__)
        self.ctx = self.app.app_context()
        self.ctx.push()
        self.resource = BatchOperationResource()
        
    def tearDown(self):
        self.ctx.pop()
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    def test_missing_params(self, mock_app, mock_request):
        # Mock request.get_json to return empty dict
        mock_request.get_json.return_value = {}
        
        # Call post method directly
        result = self.resource.post()
        
        # Check result
        self.assertEqual(result[0]['status'], 'ERROR')
        self.assertIn('Missing or invalid parameters', result[0]['errors'][0])
        self.assertEqual(result[1], 400)  # Status code
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    @patch('api.batch_resources.calculate_molecule_score', return_value=42)
    def test_property_calculation_molecule(self, mock_calc, mock_app, mock_request):
        # Mock request.get_json
        mock_request.get_json.return_value = {
            "operation": "property_calculation",
            "entity_type": "molecule",
            "ids": ["id1"]
        }
        
        # Call post method directly
        result = self.resource.post()
        
        # Check result
        self.assertEqual(result[0]['status'], 'SUCCESS')
        self.assertEqual(len(result[0]['results']), 1)
        self.assertEqual(result[0]['results'][0]['result'], 42)
        self.assertEqual(result[1], 200)  # Status code
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    @patch('api.batch_resources.MixtureProperty')
    def test_property_calculation_mixture(self, mock_mixture_prop, mock_app, mock_request):
        # Mock request.get_json
        mock_request.get_json.return_value = {
            "operation": "property_calculation",
            "entity_type": "mixture",
            "ids": ["id2"]
        }
        
        # Mock MixtureProperty.predict_mixture_properties
        mock_mixture_prop.predict_mixture_properties.return_value = 99
        
        # Call post method directly
        result = self.resource.post()
        
        # Check result
        self.assertEqual(result[0]['status'], 'SUCCESS')
        self.assertEqual(len(result[0]['results']), 1)
        self.assertEqual(result[0]['results'][0]['result'], 99)
        self.assertEqual(result[1], 200)  # Status code
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    @patch('api.batch_resources.MixtureOptimization')
    def test_mixture_optimization(self, mock_opt, mock_app, mock_request):
        # Mock request.get_json
        mock_request.get_json.return_value = {
            "operation": "mixture_optimization",
            "entity_type": "mixture",
            "ids": ["id3"]
        }
        
        # Mock MixtureOptimization.optimize_composition
        mock_opt.optimize_composition.return_value = 123
        
        # Call post method directly
        result = self.resource.post()
        
        # Check result
        self.assertEqual(result[0]['status'], 'SUCCESS')
        self.assertEqual(len(result[0]['results']), 1)
        self.assertEqual(result[0]['results'][0]['result'], 123)
        self.assertEqual(result[1], 200)  # Status code
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    @patch('api.batch_resources.ModelManager')
    def test_predictive_scoring(self, mock_model_mgr, mock_app, mock_request):
        # Mock request.get_json
        mock_request.get_json.return_value = {
            "operation": "predictive_scoring",
            "entity_type": "molecule",
            "ids": ["id4"]
        }
        
        # Mock ModelManager.predict
        mock_model_mgr.return_value.predict.return_value = 0.5
        
        # Call post method directly
        result = self.resource.post()
        
        # Check result
        self.assertEqual(result[0]['status'], 'SUCCESS')
        self.assertEqual(len(result[0]['results']), 1)
        self.assertEqual(result[0]['results'][0]['result'], 0.5)
        self.assertEqual(result[1], 200)  # Status code
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    @patch('api.batch_resources.get_data_for_export', return_value=[{"foo": "bar"}])
    @patch('api.batch_resources.generate_csv')
    def test_export(self, mock_csv, mock_export, mock_app, mock_request):
        # Mock request.get_json
        mock_request.get_json.return_value = {
            "operation": "export",
            "entity_type": "molecule",
            "ids": ["id5"]
        }
        
        # Mock generate_csv
        mock_csv_result = MagicMock()
        mock_csv_result.getvalue.return_value = "csvdata"
        mock_csv.return_value = mock_csv_result
        
        # Call post method directly
        result = self.resource.post()
        
        # Check result
        self.assertEqual(result[0]['status'], 'SUCCESS')
        self.assertEqual(len(result[0]['results']), 1)
        self.assertEqual(result[0]['results'][0]['export'], "csvdata")
        self.assertEqual(result[1], 200)  # Status code
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    def test_unsupported_operation(self, mock_app, mock_request):
        # Mock request.get_json
        mock_request.get_json.return_value = {
            "operation": "not_supported",
            "entity_type": "molecule",
            "ids": ["id6"]
        }
        
        # Call post method directly
        result = self.resource.post()
        
        # Check result
        self.assertEqual(result[0]['status'], 'ERROR')
        self.assertEqual(len(result[0]['errors']), 1)
        self.assertIn('Unsupported operation', result[0]['errors'][0]['error'])
        self.assertEqual(result[1], 400)  # Status code
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    def test_property_calculation_unsupported_entity(self, mock_app, mock_request):
        # Mock request.get_json
        mock_request.get_json.return_value = {
            "operation": "property_calculation",
            "entity_type": "unsupported",
            "ids": ["id7"]
        }
        
        # Call post method directly
        result = self.resource.post()
        
        # Check result
        self.assertEqual(result[0]['status'], 'ERROR')
        self.assertEqual(len(result[0]['errors']), 1)
        self.assertIn('Unsupported entity_type', result[0]['errors'][0]['error'])
        self.assertEqual(result[1], 400)  # Status code
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    def test_mixture_optimization_unsupported_entity(self, mock_app, mock_request):
        # Mock request.get_json
        mock_request.get_json.return_value = {
            "operation": "mixture_optimization",
            "entity_type": "unsupported",
            "ids": ["id8"]
        }
        
        # Call post method directly
        result = self.resource.post()
        
        # Check result
        self.assertEqual(result[0]['status'], 'ERROR')
        self.assertEqual(len(result[0]['errors']), 1)
        self.assertIn('Unsupported entity_type', result[0]['errors'][0]['error'])
        self.assertEqual(result[1], 400)  # Status code
    
    @patch('api.batch_resources.request')
    @patch('api.batch_resources.current_app')
    def test_batch_with_partial_errors(self, mock_app, mock_request):
        # Mock request.get_json
        mock_request.get_json.return_value = {
            "operation": "property_calculation",
            "entity_type": "molecule",
            "ids": ["id1", "id_error", "id2"]
        }
        
        # Mock calculate_molecule_score to succeed for id1 and id2, but fail for id_error
        def mock_calculate(item_id):
            if item_id == "id_error":
                raise ValueError("Test error")
            return 42
        
        with patch('api.batch_resources.calculate_molecule_score', side_effect=mock_calculate):
            # Call post method directly
            result = self.resource.post()
            
            # Check result
            self.assertEqual(result[0]['status'], 'COMPLETED_WITH_WARNINGS')
            self.assertEqual(len(result[0]['results']), 2)  # Two successful results
            self.assertEqual(len(result[0]['errors']), 1)   # One error
            self.assertEqual(result[1], 200)  # Status code

if __name__ == '__main__':
    unittest.main()