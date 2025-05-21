#!/usr/bin/env python3
"""
CryoProtect v2 - End-to-End Workflow Tests

This script tests the end-to-end workflows of the CryoProtect v2 application.
It verifies that complete user workflows function correctly.
"""

import os
import sys
import json
import uuid
import logging
import requests
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("end_to_end_test.log"),
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

class EndToEndWorkflowTest:
    """Test class for end-to-end workflows."""

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
        """Run all end-to-end workflow tests."""
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
        self.test_molecule_creation_and_analysis()
        self.test_mixture_creation_and_analysis()
        self.test_prediction_and_experiment_comparison()

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
        """Set up test data for end-to-end workflow tests."""
        try:
            logger.info("Setting up test data for end-to-end workflow tests")
            
            # Create a test user
            self.create_test_user()
            
            logger.info("Test data setup complete")
        except Exception as e:
            logger.error(f"Failed to set up test data: {str(e)}")
            raise

    def create_test_user(self) -> None:
        """Create a test user for end-to-end workflow tests."""
        try:
            # Check if test user already exists
            test_email = "e2e_test_user@example.com"
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
            
            # Note: We don't delete the test user
            
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

    def test_molecule_creation_and_analysis(self) -> None:
        """Test the workflow of creating and analyzing a molecule."""
        try:
            # Step 1: Create a new molecule
            new_molecule = {
                "name": "Test Ethanol",
                "formula": "C2H6O",
                "smiles": "CCO",
                "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
            }
            
            response = requests.post(f"{self.base_url}/api/v1/molecules", headers=self.get_headers(), json=new_molecule)
            
            if response.status_code == 201:
                molecule_data = response.json()
                molecule_id = molecule_data["id"]
                self.add_test_result("E2E-1.1", "Create Molecule", "Passed", f"Created molecule with ID {molecule_id}", response)
                
                # Store the molecule for cleanup
                if "molecules" not in self.test_data:
                    self.test_data["molecules"] = []
                self.test_data["molecules"].append(molecule_data)
            else:
                self.add_test_result("E2E-1.1", "Create Molecule", "Failed", f"Failed to create molecule: {response.status_code}", response)
                return  # Skip the rest of the test if molecule creation fails
            
            # Step 2: Calculate molecular properties
            response = requests.post(f"{self.base_url}/api/v1/rdkit/properties", headers=self.get_headers(), json={
                "molecule_data": new_molecule["smiles"],
                "input_format": "smiles"
            })
            
            if response.status_code == 200:
                properties_data = response.json()
                if "hydrogen_bonding" in properties_data and "logp" in properties_data:
                    self.add_test_result("E2E-1.2", "Calculate Properties", "Passed", "Calculated molecular properties", response)
                else:
                    self.add_test_result("E2E-1.2", "Calculate Properties", "Failed", "Missing expected properties in response", response)
            else:
                self.add_test_result("E2E-1.2", "Calculate Properties", "Failed", f"Failed to calculate properties: {response.status_code}", response)
            
            # Step 3: Generate visualization
            response = requests.post(f"{self.base_url}/api/v1/rdkit/visualize", headers=self.get_headers(), json={
                "molecule_data": new_molecule["smiles"],
                "input_format": "smiles",
                "width": 400,
                "height": 300
            })
            
            if response.status_code == 200:
                visualization_data = response.json()
                if "svg" in visualization_data:
                    self.add_test_result("E2E-1.3", "Generate Visualization", "Passed", "Generated molecular visualization", response)
                else:
                    self.add_test_result("E2E-1.3", "Generate Visualization", "Failed", "Missing SVG in response", response)
            else:
                self.add_test_result("E2E-1.3", "Generate Visualization", "Failed", f"Failed to generate visualization: {response.status_code}", response)
            
            # Step 4: Perform substructure search
            response = requests.post(f"{self.base_url}/api/v1/rdkit/substructure", headers=self.get_headers(), json={
                "query_mol_data": "[OH]",
                "target_mol_data": new_molecule["smiles"],
                "query_format": "smarts",
                "target_format": "smiles"
            })
            
            if response.status_code == 200:
                substructure_data = response.json()
                if "match" in substructure_data and substructure_data["match"]:
                    self.add_test_result("E2E-1.4", "Substructure Search", "Passed", "Found hydroxyl group in ethanol", response)
                else:
                    self.add_test_result("E2E-1.4", "Substructure Search", "Failed", "Did not find hydroxyl group in ethanol", response)
            else:
                self.add_test_result("E2E-1.4", "Substructure Search", "Failed", f"Failed to perform substructure search: {response.status_code}", response)
            
            # Overall workflow status
            passed_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("E2E-1.") and tc["status"] == "Passed")
            total_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("E2E-1."))
            
            if passed_steps == total_steps:
                self.add_test_result("E2E-1", "Molecule Creation and Analysis Workflow", "Passed", f"All {passed_steps}/{total_steps} steps passed")
            else:
                self.add_test_result("E2E-1", "Molecule Creation and Analysis Workflow", "Failed", f"Only {passed_steps}/{total_steps} steps passed")
        except Exception as e:
            self.add_test_result("E2E-1", "Molecule Creation and Analysis Workflow", "Failed", f"Error testing workflow: {str(e)}")

    def test_mixture_creation_and_analysis(self) -> None:
        """Test the workflow of creating and analyzing a mixture."""
        try:
            # Step 1: Create molecules for the mixture
            molecules = [
                {
                    "name": "Test Ethanol for Mixture",
                    "formula": "C2H6O",
                    "smiles": "CCO",
                    "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N-1"  # Modified to avoid unique constraint
                },
                {
                    "name": "Test Glycerol for Mixture",
                    "formula": "C3H8O3",
                    "smiles": "C(C(CO)O)O",
                    "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N-1"  # Modified to avoid unique constraint
                }
            ]
            
            molecule_ids = []
            for molecule in molecules:
                response = requests.post(f"{self.base_url}/api/v1/molecules", headers=self.get_headers(), json=molecule)
                
                if response.status_code == 201:
                    molecule_data = response.json()
                    molecule_ids.append(molecule_data["id"])
                    
                    # Store the molecule for cleanup
                    if "molecules" not in self.test_data:
                        self.test_data["molecules"] = []
                    self.test_data["molecules"].append(molecule_data)
                else:
                    self.add_test_result("E2E-2.1", "Create Molecules for Mixture", "Failed", f"Failed to create molecule: {response.status_code}", response)
                    return  # Skip the rest of the test if molecule creation fails
            
            if len(molecule_ids) == 2:
                self.add_test_result("E2E-2.1", "Create Molecules for Mixture", "Passed", f"Created {len(molecule_ids)} molecules for mixture")
            else:
                self.add_test_result("E2E-2.1", "Create Molecules for Mixture", "Failed", f"Only created {len(molecule_ids)}/2 molecules")
                return  # Skip the rest of the test if not all molecules were created
            
            # Step 2: Create a mixture
            new_mixture = {
                "name": "Test Ethanol-Glycerol Mixture",
                "description": "A test mixture of ethanol and glycerol",
                "components": [
                    {
                        "molecule_id": molecule_ids[0],
                        "concentration": 70,
                        "concentration_unit": "%v/v"
                    },
                    {
                        "molecule_id": molecule_ids[1],
                        "concentration": 30,
                        "concentration_unit": "%v/v"
                    }
                ]
            }
            
            response = requests.post(f"{self.base_url}/api/v1/mixtures", headers=self.get_headers(), json=new_mixture)
            
            if response.status_code == 201:
                mixture_data = response.json()
                mixture_id = mixture_data["id"]
                self.add_test_result("E2E-2.2", "Create Mixture", "Passed", f"Created mixture with ID {mixture_id}", response)
                
                # Store the mixture for cleanup
                if "mixtures" not in self.test_data:
                    self.test_data["mixtures"] = []
                self.test_data["mixtures"].append(mixture_data)
            else:
                self.add_test_result("E2E-2.2", "Create Mixture", "Failed", f"Failed to create mixture: {response.status_code}", response)
                return  # Skip the rest of the test if mixture creation fails
            
            # Step 3: Generate predictions for the mixture
            # First, check if we have a property type and calculation method
            property_type_id = None
            calculation_method_id = None
            
            # Get property type ID for Glass Transition Temperature
            response = self.supabase.table("property_types").select("id").eq("name", "Glass Transition Temperature").execute()
            if response.data:
                property_type_id = response.data[0]["id"]
            else:
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
            
            # Get calculation method ID for ML Model v1
            response = self.supabase.table("calculation_methods").select("id").eq("name", "ML Model v1").execute()
            if response.data:
                calculation_method_id = response.data[0]["id"]
            else:
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
            
            # Now create a prediction
            new_prediction = {
                "property_name": "Glass Transition Temperature",
                "value": -45.2,
                "confidence": 0.85,
                "calculation_method": "ML Model v1"
            }
            
            response = requests.post(f"{self.base_url}/api/v1/mixtures/{mixture_id}/predictions", headers=self.get_headers(), json=new_prediction)
            
            if response.status_code == 201:
                prediction_data = response.json()
                prediction_id = prediction_data["id"]
                self.add_test_result("E2E-2.3", "Create Prediction", "Passed", f"Created prediction with ID {prediction_id}", response)
                
                # Store the prediction for cleanup
                if "predictions" not in self.test_data:
                    self.test_data["predictions"] = []
                self.test_data["predictions"].append(prediction_data)
            else:
                self.add_test_result("E2E-2.3", "Create Prediction", "Failed", f"Failed to create prediction: {response.status_code}", response)
            
            # Overall workflow status
            passed_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("E2E-2.") and tc["status"] == "Passed")
            total_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("E2E-2."))
            
            if passed_steps == total_steps:
                self.add_test_result("E2E-2", "Mixture Creation and Analysis Workflow", "Passed", f"All {passed_steps}/{total_steps} steps passed")
            else:
                self.add_test_result("E2E-2", "Mixture Creation and Analysis Workflow", "Failed", f"Only {passed_steps}/{total_steps} steps passed")
        except Exception as e:
            self.add_test_result("E2E-2", "Mixture Creation and Analysis Workflow", "Failed", f"Error testing workflow: {str(e)}")

    def test_prediction_and_experiment_comparison(self) -> None:
        """Test the workflow of comparing predictions with experiments."""
        try:
            # Check if we have a mixture from the previous test
            if "mixtures" not in self.test_data or not self.test_data["mixtures"]:
                # Create a new mixture for this test
                self.test_mixture_creation_and_analysis()
                
                # Check again if we have a mixture
                if "mixtures" not in self.test_data or not self.test_data["mixtures"]:
                    self.add_test_result("E2E-3.1", "Get Mixture", "Failed", "No mixture available for testing")
                    return  # Skip the rest of the test if no mixture is available
            
            mixture_id = self.test_data["mixtures"][0]["id"]
            self.add_test_result("E2E-3.1", "Get Mixture", "Passed", f"Using mixture with ID {mixture_id}")
            
            # Step 1: Check if we have a prediction from the previous test
            if "predictions" not in self.test_data or not self.test_data["predictions"]:
                # Create a new prediction for this test
                new_prediction = {
                    "property_name": "Glass Transition Temperature",
                    "value": -45.2,
                    "confidence": 0.85,
                    "calculation_method": "ML Model v1"
                }
                
                response = requests.post(f"{self.base_url}/api/v1/mixtures/{mixture_id}/predictions", headers=self.get_headers(), json=new_prediction)
                
                if response.status_code == 201:
                    prediction_data = response.json()
                    prediction_id = prediction_data["id"]
                    self.add_test_result("E2E-3.2", "Create Prediction", "Passed", f"Created prediction with ID {prediction_id}", response)
                    
                    # Store the prediction for cleanup
                    if "predictions" not in self.test_data:
                        self.test_data["predictions"] = []
                    self.test_data["predictions"].append(prediction_data)
                else:
                    self.add_test_result("E2E-3.2", "Create Prediction", "Failed", f"Failed to create prediction: {response.status_code}", response)
                    return  # Skip the rest of the test if prediction creation fails
            else:
                self.add_test_result("E2E-3.2", "Create Prediction", "Passed", "Using existing prediction")
            
            # Step 2: Record an experiment
            new_experiment = {
                "property_name": "Glass Transition Temperature",
                "value": -47.3,
                "experimental_conditions": "Measured using DSC at 10°C/min",
                "date_performed": "2025-04-15"
            }
            
            response = requests.post(f"{self.base_url}/api/v1/mixtures/{mixture_id}/experiments", headers=self.get_headers(), json=new_experiment)
            
            if response.status_code == 201:
                experiment_data = response.json()
                experiment_id = experiment_data["id"]
                self.add_test_result("E2E-3.3", "Record Experiment", "Passed", f"Recorded experiment with ID {experiment_id}", response)
                
                # Store the experiment for cleanup
                if "experiments" not in self.test_data:
                    self.test_data["experiments"] = []
                self.test_data["experiments"].append(experiment_data)
            else:
                self.add_test_result("E2E-3.3", "Record Experiment", "Failed", f"Failed to record experiment: {response.status_code}", response)
                return  # Skip the rest of the test if experiment recording fails
            
            # Step 3: Compare prediction with experiment
            response = requests.get(f"{self.base_url}/api/v1/mixtures/{mixture_id}/compare?property_name=Glass Transition Temperature", headers=self.get_headers())
            
            if response.status_code == 200:
                comparison_data = response.json()
                if "prediction" in comparison_data and "experiment" in comparison_data and "difference" in comparison_data:
                    self.add_test_result("E2E-3.4", "Compare Prediction with Experiment", "Passed", "Successfully compared prediction with experiment", response)
                else:
                    self.add_test_result("E2E-3.4", "Compare Prediction with Experiment", "Failed", "Missing expected data in comparison response", response)
            else:
                self.add_test_result("E2E-3.4", "Compare Prediction with Experiment", "Failed", f"Failed to compare prediction with experiment: {response.status_code}", response)
            
            # Overall workflow status
            passed_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("E2E-3.") and tc["status"] == "Passed")
            total_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("E2E-3."))
            
            if passed_steps == total_steps:
                self.add_test_result("E2E-3", "Prediction and Experiment Comparison Workflow", "Passed", f"All {passed_steps}/{total_steps} steps passed")
            else:
                self.add_test_result("E2E-3", "Prediction and Experiment Comparison Workflow", "Failed", f"Only {passed_steps}/{total_steps} steps passed")
        except Exception as e:
            self.add_test_result("E2E-3", "Prediction and Experiment Comparison Workflow", "Failed", f"Error testing workflow: {str(e)}")

    def save_results(self, filename: str) -> None:
        """Save the test results to a file."""
        try:
            with open(filename, 'w') as f:
                json.dump(self.test_results, f, indent=2)
            logger.info(f"Test results saved to {filename}")
        except Exception as e:
            logger.error(f"Failed to save test results: {str(e)}")

def main():
    """Main function."""
    logger.info("Starting end-to-end workflow tests")
    
    # Create and run the tests
    test = EndToEndWorkflowTest()
    results = test.run_tests()
    
    # Save the results
    test.save_results("end_to_end_test_results.json")
    
    # Print summary
    logger.info(f"Test Status: {results['status']}")
    logger.info(f"Total Tests: {results['total_tests']}")
    logger.info(f"Passed Tests: {results['passed_tests']}")
    logger.info(f"Failed Tests: {results['failed_tests']}")
    logger.info(f"Skipped Tests: {results['skipped_tests']}")
    
    # Exit with appropriate status code
    if results["status"] == "Passed":
        logger.info("All tests passed!")
        sys.exit(0)
    else:
        logger.error("Some tests failed. See log for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()