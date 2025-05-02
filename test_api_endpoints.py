#!/usr/bin/env python3
"""
CryoProtect v2 - API Endpoint Tests

This script tests the API endpoints of the CryoProtect v2 application.
It verifies that all endpoints are functioning correctly, including request validation,
response formatting, error handling, and authentication.
"""

import os
import sys
import json
import uuid
import logging
import requests
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple

# Add the parent directory to the path so we can import the app modules
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("api_endpoint_test.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Import the Supabase client for test data setup
try:
    from config import SUPABASE_URL, SUPABASE_KEY, SUPABASE_SERVICE_KEY
    from supabase import create_client, Client
except ImportError as e:
    logger.error(f"Failed to import required modules: {str(e)}")
    logger.error("Please ensure that the required modules are installed.")
    sys.exit(1)

class APIEndpointTest:
    """Test class for API endpoints."""

    def __init__(self, base_url: str = "http://localhost:5000"):
        """Initialize the test class."""
        self.base_url = base_url
        self.supabase: Optional[Client] = None
        self.auth_token: Optional[str] = None
        self.test_results: Dict[str, Any] = {
            "status": "Not Started",
            "total_tests": 0,
            "passed_tests": 0,
            "failed_tests": 0,
            "skipped_tests": 0,
            "test_cases": []
        }
        self.test_data = {}

    def connect_to_database(self) -> bool:
        """Connect to the Supabase database."""
        try:
            self.supabase = create_client(SUPABASE_URL, SUPABASE_SERVICE_KEY)
            logger.info("Connected to Supabase database")
            return True
        except Exception as e:
            logger.error(f"Failed to connect to Supabase database: {str(e)}")
            return False

    def run_tests(self) -> Dict[str, Any]:
        """Run all API endpoint tests."""
        self.test_results["status"] = "Running"
        self.test_results["start_time"] = datetime.now().isoformat()

        # Connect to the database for test data setup
        if not self.connect_to_database():
            self.test_results["status"] = "Failed"
            self.test_results["end_time"] = datetime.now().isoformat()
            return self.test_results

        # Setup test data
        self.setup_test_data()

        # Authenticate
        self.authenticate()

        # Run the tests
        self.test_health_endpoint()
        self.test_molecules_endpoints()
        self.test_mixtures_endpoints()
        self.test_predictions_endpoints()
        self.test_experiments_endpoints()
        self.test_comparison_endpoint()
        self.test_rdkit_endpoints()
        self.test_error_handling()

        # Cleanup test data
        self.cleanup_test_data()

        # Calculate test results
        self.test_results["total_tests"] = len(self.test_results["test_cases"])
        self.test_results["passed_tests"] = sum(1 for tc in self.test_results["test_cases"] if tc["status"] == "Passed")
        self.test_results["failed_tests"] = sum(1 for tc in self.test_results["test_cases"] if tc["status"] == "Failed")
        self.test_results["skipped_tests"] = sum(1 for tc in self.test_results["test_cases"] if tc["status"] == "Skipped")
        
        if self.test_results["failed_tests"] == 0:
            self.test_results["status"] = "Passed"
        else:
            self.test_results["status"] = "Failed"
        
        self.test_results["end_time"] = datetime.now().isoformat()
        
        return self.test_results

    def add_test_result(self, test_id: str, test_name: str, status: str, message: str, response: Optional[requests.Response] = None) -> None:
        """Add a test result to the test results."""
        result = {
            "id": test_id,
            "name": test_name,
            "status": status,
            "message": message
        }
        
        if response:
            try:
                result["response"] = {
                    "status_code": response.status_code,
                    "content": response.json() if response.headers.get("content-type") == "application/json" else response.text
                }
            except Exception:
                result["response"] = {
                    "status_code": response.status_code,
                    "content": response.text
                }
        
        self.test_results["test_cases"].append(result)
        
        if status == "Passed":
            logger.info(f"Test {test_id} - {test_name}: PASSED")
        elif status == "Failed":
            logger.error(f"Test {test_id} - {test_name}: FAILED - {message}")
        else:
            logger.warning(f"Test {test_id} - {test_name}: SKIPPED - {message}")

    def setup_test_data(self) -> None:
        """Set up test data for API endpoint tests."""
        try:
            logger.info("Setting up test data for API endpoint tests")
            
            # Create a test user
            self.create_test_user()
            
            # Create test molecules
            self.create_test_molecules()
            
            # Create test mixtures
            self.create_test_mixtures()
            
            # Create test predictions
            self.create_test_predictions()
            
            # Create test experiments
            self.create_test_experiments()
            
            logger.info("Test data setup complete")
        except Exception as e:
            logger.error(f"Failed to set up test data: {str(e)}")
            raise

    def create_test_user(self) -> None:
        """Create a test user for API endpoint tests."""
        try:
            # Check if test user already exists
            test_email = "api_test_user@example.com"
            response = self.supabase.table("auth.users").select("id, email").eq("email", test_email).execute()
            
            if response.data:
                user_id = response.data[0]["id"]
                logger.info(f"Test user {test_email} already exists with ID {user_id}")
            else:
                # Create user with auth admin API
                password = "Test123456!"
                user_data = {
                    "email": test_email,
                    "password": password,
                    "email_confirm": True
                }
                
                response = self.supabase.auth.admin.create_user(user_data)
                user_id = response.user.id
                logger.info(f"Created test user {test_email} with ID {user_id}")
            
            # Store the user ID for later use
            self.test_data["user_id"] = user_id
            self.test_data["user_email"] = test_email
            self.test_data["user_password"] = "Test123456!"
        except Exception as e:
            logger.error(f"Failed to create test user: {str(e)}")
            raise

    def create_test_molecules(self) -> None:
        """Create test molecules for API endpoint tests."""
        try:
            # Create test molecules
            molecules = [
                {
                    "id": str(uuid.uuid4()),
                    "name": "Test Ethanol",
                    "formula": "C2H6O",
                    "smiles": "CCO",
                    "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                    "created_by": self.test_data["user_id"]
                },
                {
                    "id": str(uuid.uuid4()),
                    "name": "Test Glycerol",
                    "formula": "C3H8O3",
                    "smiles": "C(C(CO)O)O",
                    "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
                    "created_by": self.test_data["user_id"]
                }
            ]
            
            # Insert the molecules
            for molecule in molecules:
                self.supabase.table("molecules").insert(molecule).execute()
                logger.info(f"Created test molecule {molecule['id']} - {molecule['name']}")
            
            # Store the molecule IDs for later use
            self.test_data["molecules"] = molecules
        except Exception as e:
            logger.error(f"Failed to create test molecules: {str(e)}")
            raise

    def create_test_mixtures(self) -> None:
        """Create test mixtures for API endpoint tests."""
        try:
            # Create test mixtures
            mixtures = [
                {
                    "id": str(uuid.uuid4()),
                    "name": "Test Mixture 1",
                    "description": "A test mixture of ethanol and glycerol",
                    "created_by": self.test_data["user_id"]
                }
            ]
            
            # Insert the mixtures
            for mixture in mixtures:
                self.supabase.table("mixtures").insert(mixture).execute()
                logger.info(f"Created test mixture {mixture['id']} - {mixture['name']}")
            
            # Create mixture components
            components = [
                {
                    "id": str(uuid.uuid4()),
                    "mixture_id": mixtures[0]["id"],
                    "molecule_id": self.test_data["molecules"][0]["id"],
                    "concentration": 70,
                    "concentration_unit": "%v/v",
                    "created_by": self.test_data["user_id"]
                },
                {
                    "id": str(uuid.uuid4()),
                    "mixture_id": mixtures[0]["id"],
                    "molecule_id": self.test_data["molecules"][1]["id"],
                    "concentration": 30,
                    "concentration_unit": "%v/v",
                    "created_by": self.test_data["user_id"]
                }
            ]
            
            # Insert the components
            for component in components:
                self.supabase.table("mixture_components").insert(component).execute()
                logger.info(f"Created test mixture component {component['id']}")
            
            # Store the mixture IDs for later use
            self.test_data["mixtures"] = mixtures
            self.test_data["mixture_components"] = components
        except Exception as e:
            logger.error(f"Failed to create test mixtures: {str(e)}")
            raise

    def create_test_predictions(self) -> None:
        """Create test predictions for API endpoint tests."""
        try:
            # Get property type ID for Glass Transition Temperature
            response = self.supabase.table("property_types").select("id").eq("name", "Glass Transition Temperature").execute()
            if not response.data:
                # Create the property type
                property_type = {
                    "id": str(uuid.uuid4()),
                    "name": "Glass Transition Temperature",
                    "data_type": "numeric",
                    "description": "The temperature at which a material transitions from a hard, glassy state to a soft, rubbery state",
                    "units": "°C",
                    "created_by": self.test_data["user_id"]
                }
                self.supabase.table("property_types").insert(property_type).execute()
                property_type_id = property_type["id"]
                logger.info(f"Created test property type {property_type_id} - {property_type['name']}")
            else:
                property_type_id = response.data[0]["id"]
            
            # Get calculation method ID for ML Model v1
            response = self.supabase.table("calculation_methods").select("id").eq("name", "ML Model v1").execute()
            if not response.data:
                # Create the calculation method
                calculation_method = {
                    "id": str(uuid.uuid4()),
                    "name": "ML Model v1",
                    "description": "Machine learning model for property prediction",
                    "version": "1.0",
                    "created_by": self.test_data["user_id"]
                }
                self.supabase.table("calculation_methods").insert(calculation_method).execute()
                calculation_method_id = calculation_method["id"]
                logger.info(f"Created test calculation method {calculation_method_id} - {calculation_method['name']}")
            else:
                calculation_method_id = response.data[0]["id"]
            
            # Create test predictions
            predictions = [
                {
                    "id": str(uuid.uuid4()),
                    "mixture_id": self.test_data["mixtures"][0]["id"],
                    "molecule_id": None,
                    "property_type_id": property_type_id,
                    "calculation_method_id": calculation_method_id,
                    "numeric_value": -45.2,
                    "confidence": 0.85,
                    "created_by": self.test_data["user_id"]
                }
            ]
            
            # Insert the predictions
            for prediction in predictions:
                self.supabase.table("predictions").insert(prediction).execute()
                logger.info(f"Created test prediction {prediction['id']}")
            
            # Store the prediction IDs for later use
            self.test_data["predictions"] = predictions
            self.test_data["property_type_id"] = property_type_id
            self.test_data["calculation_method_id"] = calculation_method_id
        except Exception as e:
            logger.error(f"Failed to create test predictions: {str(e)}")
            raise

    def create_test_experiments(self) -> None:
        """Create test experiments for API endpoint tests."""
        try:
            # Create test experiments
            experiments = [
                {
                    "id": str(uuid.uuid4()),
                    "mixture_id": self.test_data["mixtures"][0]["id"],
                    "molecule_id": None,
                    "property_type_id": self.test_data["property_type_id"],
                    "numeric_value": -47.3,
                    "experimental_conditions": "Measured using DSC at 10°C/min",
                    "date_performed": "2025-04-15",
                    "created_by": self.test_data["user_id"]
                }
            ]
            
            # Insert the experiments
            for experiment in experiments:
                self.supabase.table("experiments").insert(experiment).execute()
                logger.info(f"Created test experiment {experiment['id']}")
            
            # Store the experiment IDs for later use
            self.test_data["experiments"] = experiments
        except Exception as e:
            logger.error(f"Failed to create test experiments: {str(e)}")
            raise

    def cleanup_test_data(self) -> None:
        """Clean up test data after tests are complete."""
        try:
            logger.info("Cleaning up test data")
            
            # Delete test experiments
            if "experiments" in self.test_data:
                for experiment in self.test_data["experiments"]:
                    self.supabase.table("experiments").delete().eq("id", experiment["id"]).execute()
                    logger.info(f"Deleted test experiment {experiment['id']}")
            
            # Delete test predictions
            if "predictions" in self.test_data:
                for prediction in self.test_data["predictions"]:
                    self.supabase.table("predictions").delete().eq("id", prediction["id"]).execute()
                    logger.info(f"Deleted test prediction {prediction['id']}")
            
            # Delete test mixture components
            if "mixture_components" in self.test_data:
                for component in self.test_data["mixture_components"]:
                    self.supabase.table("mixture_components").delete().eq("id", component["id"]).execute()
                    logger.info(f"Deleted test mixture component {component['id']}")
            
            # Delete test mixtures
            if "mixtures" in self.test_data:
                for mixture in self.test_data["mixtures"]:
                    self.supabase.table("mixtures").delete().eq("id", mixture["id"]).execute()
                    logger.info(f"Deleted test mixture {mixture['id']}")
            
            # Delete test molecules
            if "molecules" in self.test_data:
                for molecule in self.test_data["molecules"]:
                    self.supabase.table("molecules").delete().eq("id", molecule["id"]).execute()
                    logger.info(f"Deleted test molecule {molecule['id']}")
            
            # Note: We don't delete the test user, property types, or calculation methods
            
            logger.info("Test data cleanup complete")
        except Exception as e:
            logger.error(f"Failed to clean up test data: {str(e)}")
            # Don't raise the exception, as we want to continue with the test results

    def authenticate(self) -> None:
        """Authenticate with the API."""
        try:
            # Sign in with the test user
            response = self.supabase.auth.sign_in_with_password({
                "email": self.test_data["user_email"],
                "password": self.test_data["user_password"]
            })
            
            # Get the access token
            self.auth_token = response.session.access_token
            logger.info("Authentication successful")
        except Exception as e:
            logger.error(f"Failed to authenticate: {str(e)}")
            self.auth_token = None

    def get_headers(self) -> Dict[str, str]:
        """Get headers for API requests."""
        headers = {
            "Content-Type": "application/json"
        }
        
        if self.auth_token:
            headers["Authorization"] = f"Bearer {self.auth_token}"
        
        return headers

    def test_health_endpoint(self) -> None:
        """Test the health endpoint."""
        try:
            # Send request to the health endpoint
            response = requests.get(f"{self.base_url}/health")
            
            # Check the response
            if response.status_code == 200:
                data = response.json()
                if data.get("status") == "ok":
                    self.add_test_result("API-1", "Health Endpoint", "Passed", "Health endpoint returned status ok", response)
                else:
                    self.add_test_result("API-1", "Health Endpoint", "Failed", f"Health endpoint returned unexpected status: {data.get('status')}", response)
            else:
                self.add_test_result("API-1", "Health Endpoint", "Failed", f"Health endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-1", "Health Endpoint", "Failed", f"Error testing health endpoint: {str(e)}")

    def test_molecules_endpoints(self) -> None:
        """Test the molecules endpoints."""
        # Test GET /molecules
        try:
            response = requests.get(f"{self.base_url}/api/v1/molecules", headers=self.get_headers())
            
            if response.status_code == 200:
                data = response.json()
                if isinstance(data, list):
                    self.add_test_result("API-2.1", "GET /molecules", "Passed", "Endpoint returned a list of molecules", response)
                else:
                    self.add_test_result("API-2.1", "GET /molecules", "Failed", "Endpoint did not return a list", response)
            else:
                self.add_test_result("API-2.1", "GET /molecules", "Failed", f"Endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-2.1", "GET /molecules", "Failed", f"Error testing GET /molecules: {str(e)}")
        
        # Test GET /molecules/{id}
        try:
            molecule_id = self.test_data["molecules"][0]["id"]
            response = requests.get(f"{self.base_url}/api/v1/molecules/{molecule_id}", headers=self.get_headers())
            
            if response.status_code == 200:
                data = response.json()
                if data.get("id") == molecule_id:
                    self.add_test_result("API-2.2", "GET /molecules/{id}", "Passed", "Endpoint returned the correct molecule", response)
                else:
                    self.add_test_result("API-2.2", "GET /molecules/{id}", "Failed", "Endpoint returned the wrong molecule", response)
            else:
                self.add_test_result("API-2.2", "GET /molecules/{id}", "Failed", f"Endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-2.2", "GET /molecules/{id}", "Failed", f"Error testing GET /molecules/{{id}}: {str(e)}")
        
        # Test POST /molecules
        try:
            # Create a new molecule
            new_molecule = {
                "name": "Test DMSO",
                "formula": "C2H6OS",
                "smiles": "CS(=O)C",
                "inchikey": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N"
            }
            
            response = requests.post(f"{self.base_url}/api/v1/molecules", headers=self.get_headers(), json=new_molecule)
            
            if response.status_code == 201:
                data = response.json()
                if data.get("name") == new_molecule["name"]:
                    self.add_test_result("API-2.3", "POST /molecules", "Passed", "Endpoint created a new molecule", response)
                    
                    # Store the new molecule ID for cleanup
                    if "created_molecules" not in self.test_data:
                        self.test_data["created_molecules"] = []
                    self.test_data["created_molecules"].append(data)
                else:
                    self.add_test_result("API-2.3", "POST /molecules", "Failed", "Endpoint created a molecule with the wrong name", response)
            else:
                self.add_test_result("API-2.3", "POST /molecules", "Failed", f"Endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-2.3", "POST /molecules", "Failed", f"Error testing POST /molecules: {str(e)}")

    def test_mixtures_endpoints(self) -> None:
        """Test the mixtures endpoints."""
        # Test GET /mixtures
        try:
            response = requests.get(f"{self.base_url}/api/v1/mixtures", headers=self.get_headers())
            
            if response.status_code == 200:
                data = response.json()
                if isinstance(data, list):
                    self.add_test_result("API-3.1", "GET /mixtures", "Passed", "Endpoint returned a list of mixtures", response)
                else:
                    self.add_test_result("API-3.1", "GET /mixtures", "Failed", "Endpoint did not return a list", response)
            else:
                self.add_test_result("API-3.1", "GET /mixtures", "Failed", f"Endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-3.1", "GET /mixtures", "Failed", f"Error testing GET /mixtures: {str(e)}")
        
        # Test GET /mixtures/{id}
        try:
            mixture_id = self.test_data["mixtures"][0]["id"]
            response = requests.get(f"{self.base_url}/api/v1/mixtures/{mixture_id}", headers=self.get_headers())
            
            if response.status_code == 200:
                data = response.json()
                if data.get("id") == mixture_id:
                    self.add_test_result("API-3.2", "GET /mixtures/{id}", "Passed", "Endpoint returned the correct mixture", response)
                else:
                    self.add_test_result("API-3.2", "GET /mixtures/{id}", "Failed", "Endpoint returned the wrong mixture", response)
            else:
                self.add_test_result("API-3.2", "GET /mixtures/{id}", "Failed", f"Endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-3.2", "GET /mixtures/{id}", "Failed", f"Error testing GET /mixtures/{{id}}: {str(e)}")
        
        # Test POST /mixtures
        try:
            # Create a new mixture
            new_mixture = {
                "name": "Test Mixture 2",
                "description": "A test mixture of DMSO and glycerol",
                "components": [
                    {
                        "molecule_id": self.test_data["molecules"][1]["id"],
                        "concentration": 50,
                        "concentration_unit": "%v/v"
                    }
                ]
            }
            
            # Add the created molecule if available
            if "created_molecules" in self.test_data and self.test_data["created_molecules"]:
                new_mixture["components"].append({
                    "molecule_id": self.test_data["created_molecules"][0]["id"],
                    "concentration": 50,
                    "concentration_unit": "%v/v"
                })
            
            response = requests.post(f"{self.base_url}/api/v1/mixtures", headers=self.get_headers(), json=new_mixture)
            
            if response.status_code == 201:
                data = response.json()
                if data.get("name") == new_mixture["name"]:
                    self.add_test_result("API-3.3", "POST /mixtures", "Passed", "Endpoint created a new mixture", response)
                    
                    # Store the new mixture ID for cleanup
                    if "created_mixtures" not in self.test_data:
                        self.test_data["created_mixtures"] = []
                    self.test_data["created_mixtures"].append(data)
                else:
                    self.add_test_result("API-3.3", "POST /mixtures", "Failed", "Endpoint created a mixture with the wrong name", response)
            else:
                self.add_test_result("API-3.3", "POST /mixtures", "Failed", f"Endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-3.3", "POST /mixtures", "Failed", f"Error testing POST /mixtures: {str(e)}")

    def test_predictions_endpoints(self) -> None:
        """Test the predictions endpoints."""
        # Test GET /mixtures/{id}/predictions
        try:
            mixture_id = self.test_data["mixtures"][0]["id"]
            response = requests.get(f"{self.base_url}/api/v1/mixtures/{mixture_id}/predictions", headers=self.get_headers())
            
            if response.status_code == 200:
                data = response.json()
                if isinstance(data, list):
                    self.add_test_result("API-4.1", "GET /mixtures/{id}/predictions", "Passed", "Endpoint returned a list of predictions", response)
                else:
                    self.add_test_result("API-4.1", "GET /mixtures/{id}/predictions", "Failed", "Endpoint did not return a list", response)
            else:
                self.add_test_result("API-4.1", "GET /mixtures/{id}/predictions", "Failed", f"Endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-4.1", "GET /mixtures/{id}/predictions", "Failed", f"Error testing GET /mixtures/{{id}}/predictions: {str(e)}")
        
        # Test POST /mixtures/{id}/predictions
        try:
            # Create a new prediction
            mixture_id = self.test_data["mixtures"][0]["id"]
            new_prediction = {
                "property_name": "Freezing Point",
                "value": -15.3,
                "confidence": 0.9,
                "calculation_method": "CryoProtect Scoring"
            }
            
            response = requests.post(f"{self.base_url}/api/v1/mixtures/{mixture_id}/predictions", headers=self.get_headers(), json=new_prediction)
            
            if response.status_code == 201:
                data = response.json()
                if data.get("property_name") == new_prediction["property_name"]:
                    self.add_test_result("API-4.2", "POST /mixtures/{id}/predictions", "Passed", "Endpoint created a new prediction", response)
                    
                    # Store the new prediction ID for cleanup
                    if "created_predictions" not in self.test_data:
                        self.test_data["created_predictions"] = []
                    self.test_data["created_predictions"].append(data)
                else:
                    self.add_test_result("API-4.2", "POST /mixtures/{id}/predictions", "Failed", "Endpoint created a prediction with the wrong property name", response)
            else:
                self.add_test_result("API-4.2", "POST /mixtures/{id}/predictions", "Failed", f"Endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-4.2", "POST /mixtures/{id}/predictions", "Failed", f"Error testing POST /mixtures/{{id}}/predictions: {str(e)}")

    def test_experiments_endpoints(self) -> None:
        """Test the experiments endpoints."""
        # Test GET /mixtures/{id}/experiments
        try:
            mixture_id = self.test_data["mixtures"][0]["id"]
            response = requests.get(f"{self.base_url}/api/v1/mixtures/{mixture_id}/experiments", headers=self.get_headers())
            
            if response.status_code == 200:
                data = response.json()
                if isinstance(data, list):
                    self.add_test_result("API-5.1", "GET /mixtures/{id}/experiments", "Passed", "Endpoint returned a list of experiments", response)
                else:
                    self.add_test_result("API-5.1", "GET /mixtures/{id}/experiments", "Failed", "Endpoint did not return a list", response)
            else:
                self.add_test_result("API-5.1", "GET /mixtures/{id}/experiments", "Failed", f"Endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-5.1", "GET /mixtures/{id}/experiments", "Failed", f"Error testing GET /mixtures/{{id}}/experiments: {str(e)}")
        
        # Test POST /mixtures/{id}/experiments
        try:
            # Create a new experiment
