"""
API Tests for Predictive Models Batch Endpoints in CryoProtect v2

Covers:
- /api/v1/predictive-models/train-batch
- /api/v1/predictive-models/evaluate
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

class TestPredictiveModelsBatchEndpoints(MockSupabaseBaseTestCase):
    """Integration tests for predictive models batch API endpoints."""

    def setUp(self):
        super().setUp()
        self.client = self.app.test_client()
        self.auth_token = str(uuid.uuid4())
        self.headers = {'Authorization': f'Bearer {self.auth_token}'}

    @patch('api.predictive_models_resources.batch_train_models')
    def test_train_batch_endpoint_success(self, mock_batch_train):
        """Test POST /api/v1/predictive-models/train-batch with valid data."""
        # Mock the batch training function
        mock_batch_train.return_value = {
            "models": [
                {
                    "property_name": "Freezing Point",
                    "algorithm": "random_forest",
                    "status": "success",
                    "metrics": {
                        "r2": 0.85,
                        "rmse": 2.3
                    },
                    "model_id": "model-123"
                },
                {
                    "property_name": "Vitrification",
                    "algorithm": "gradient_boosting",
                    "status": "success",
                    "metrics": {
                        "r2": 0.78,
                        "rmse": 3.1
                    },
                    "model_id": "model-456"
                }
            ],
            "training_time": 45.2
        }
        
        # Test payload
        payload = {
            "models": [
                {
                    "property_name": "Freezing Point",
                    "algorithm": "random_forest",
                    "hyperparameters": {
                        "n_estimators": 100,
                        "max_depth": 10
                    }
                },
                {
                    "property_name": "Vitrification",
                    "algorithm": "gradient_boosting",
                    "hyperparameters": {
                        "n_estimators": 150,
                        "learning_rate": 0.1
                    }
                }
            ],
            "training_data_source": "experiments",
            "cross_validation": True
        }
        
        # Make request
        response = self.client.post('/api/v1/predictive-models/train-batch', 
                                   json=payload, 
                                   headers=self.headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn("models", data)
        self.assertEqual(len(data["models"]), 2)
        self.assertIn("training_time", data)
        
        # Check first model details
        model1 = data["models"][0]
        self.assertEqual(model1["property_name"], "Freezing Point")
        self.assertEqual(model1["algorithm"], "random_forest")
        self.assertEqual(model1["status"], "success")
        self.assertIn("metrics", model1)
        self.assertIn("model_id", model1)

    def test_train_batch_endpoint_missing_auth(self):
        """Test POST /api/v1/predictive-models/train-batch without authentication."""
        payload = {
            "models": [
                {
                    "property_name": "Freezing Point",
                    "algorithm": "random_forest"
                }
            ]
        }
        
        response = self.client.post('/api/v1/predictive-models/train-batch', json=payload)
        self.assertIn(response.status_code, [401, 403])

    def test_train_batch_endpoint_invalid_payload(self):
        """Test POST /api/v1/predictive-models/train-batch with invalid payload."""
        # Empty payload
        response = self.client.post('/api/v1/predictive-models/train-batch', 
                                   json={}, 
                                   headers=self.headers)
        self.assertIn(response.status_code, [400, 422])
        
        # Missing required fields
        response = self.client.post('/api/v1/predictive-models/train-batch', 
                                   json={"training_data_source": "experiments"}, 
                                   headers=self.headers)
        self.assertIn(response.status_code, [400, 422])

    @patch('api.predictive_models_resources.batch_train_models')
    def test_train_batch_endpoint_with_errors(self, mock_batch_train):
        """Test POST /api/v1/predictive-models/train-batch with partial success."""
        # Mock the batch training function with partial success
        mock_batch_train.return_value = {
            "models": [
                {
                    "property_name": "Freezing Point",
                    "algorithm": "random_forest",
                    "status": "success",
                    "metrics": {
                        "r2": 0.85,
                        "rmse": 2.3
                    },
                    "model_id": "model-123"
                },
                {
                    "property_name": "Vitrification",
                    "algorithm": "invalid_algorithm",
                    "status": "error",
                    "error": "Invalid algorithm specified"
                }
            ],
            "training_time": 25.7
        }
        
        # Test payload
        payload = {
            "models": [
                {
                    "property_name": "Freezing Point",
                    "algorithm": "random_forest"
                },
                {
                    "property_name": "Vitrification",
                    "algorithm": "invalid_algorithm"
                }
            ],
            "training_data_source": "experiments"
        }
        
        # Make request
        response = self.client.post('/api/v1/predictive-models/train-batch', 
                                   json=payload, 
                                   headers=self.headers)
        
        # Check response - should still be 200 with partial success
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn("models", data)
        self.assertEqual(len(data["models"]), 2)
        
        # Check error model
        error_model = data["models"][1]
        self.assertEqual(error_model["status"], "error")
        self.assertIn("error", error_model)

    @patch('api.predictive_models_resources.evaluate_model')
    def test_evaluate_model_endpoint_success(self, mock_evaluate):
        """Test POST /api/v1/predictive-models/evaluate with valid data."""
        # Mock the model evaluation function
        mock_evaluate.return_value = {
            "model_id": "model-123",
            "property_name": "Freezing Point",
            "algorithm": "random_forest",
            "evaluation_method": "test_split",
            "metrics": {
                "r2": 0.82,
                "rmse": 2.5,
                "mae": 1.8,
                "mse": 6.25
            },
            "feature_importance": {
                "h_bond_donors": 0.35,
                "logp": 0.25,
                "molecular_weight": 0.20,
                "tpsa": 0.15,
                "rotatable_bonds": 0.05
            },
            "evaluation_time": 3.2
        }
        
        # Test payload
        payload = {
            "model_id": "model-123",
            "evaluation_method": "test_split",
            "test_size": 0.2,
            "random_seed": 42
        }
        
        # Make request
        response = self.client.post('/api/v1/predictive-models/evaluate', 
                                   json=payload, 
                                   headers=self.headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["model_id"], "model-123")
        self.assertEqual(data["evaluation_method"], "test_split")
        self.assertIn("metrics", data)
        self.assertIn("feature_importance", data)
        
        # Check metrics
        metrics = data["metrics"]
        self.assertIn("r2", metrics)
        self.assertIn("rmse", metrics)
        self.assertIn("mae", metrics)
        self.assertIn("mse", metrics)

    @patch('api.predictive_models_resources.evaluate_model')
    def test_evaluate_model_endpoint_cross_validation(self, mock_evaluate):
        """Test POST /api/v1/predictive-models/evaluate with cross-validation."""
        # Mock the model evaluation function
        mock_evaluate.return_value = {
            "model_id": "model-123",
            "property_name": "Freezing Point",
            "algorithm": "random_forest",
            "evaluation_method": "cross_validation",
            "metrics": {
                "cv_r2_mean": 0.80,
                "cv_r2_std": 0.05,
                "cv_rmse_mean": 2.7,
                "cv_rmse_std": 0.3
            },
            "evaluation_time": 8.5
        }
        
        # Test payload
        payload = {
            "model_id": "model-123",
            "evaluation_method": "cross_validation",
            "cv_folds": 5,
            "random_seed": 42
        }
        
        # Make request
        response = self.client.post('/api/v1/predictive-models/evaluate', 
                                   json=payload, 
                                   headers=self.headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["evaluation_method"], "cross_validation")
        
        # Check CV metrics
        metrics = data["metrics"]
        self.assertIn("cv_r2_mean", metrics)
        self.assertIn("cv_r2_std", metrics)
        self.assertIn("cv_rmse_mean", metrics)
        self.assertIn("cv_rmse_std", metrics)
        
        # Verify that the function was called with the correct parameters
        mock_evaluate.assert_called_once()
        args, kwargs = mock_evaluate.call_args
        self.assertEqual(kwargs.get("model_id"), "model-123")
        self.assertEqual(kwargs.get("evaluation_method"), "cross_validation")
        self.assertEqual(kwargs.get("cv_folds"), 5)
        self.assertEqual(kwargs.get("random_seed"), 42)

    def test_evaluate_model_endpoint_missing_auth(self):
        """Test POST /api/v1/predictive-models/evaluate without authentication."""
        payload = {
            "model_id": "model-123",
            "evaluation_method": "test_split"
        }
        
        response = self.client.post('/api/v1/predictive-models/evaluate', json=payload)
        self.assertIn(response.status_code, [401, 403])

    def test_evaluate_model_endpoint_invalid_payload(self):
        """Test POST /api/v1/predictive-models/evaluate with invalid payload."""
        # Missing model_id
        response = self.client.post('/api/v1/predictive-models/evaluate', 
                                   json={"evaluation_method": "test_split"}, 
                                   headers=self.headers)
        self.assertIn(response.status_code, [400, 422])
        
        # Invalid evaluation method
        response = self.client.post('/api/v1/predictive-models/evaluate', 
                                   json={"model_id": "model-123", "evaluation_method": "invalid_method"}, 
                                   headers=self.headers)
        self.assertIn(response.status_code, [400, 422])

if __name__ == '__main__':
    unittest.main()