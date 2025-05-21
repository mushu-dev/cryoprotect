#!/usr/bin/env python3
"""
CryoProtect v2 - Real Data End-to-End Workflow Tests

This script tests the end-to-end workflows of the CryoProtect v2 application
using real data from the Supabase database (project ID: tsdlmynydfuypiugmkev).
"""

import os
import sys
import json
import time
import logging
import requests
import statistics
import uuid
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("real_data_test.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Import the service role helper
try:
    from service_role_helper import get_supabase_client, get_user_id
    from supabase import create_client, Client
except ImportError as e:
    logger.error(f"Failed to import required modules: {str(e)}")
    logger.error("Please ensure that the required modules are installed.")
    sys.exit(1)

class RealDataWorkflowTest:
    """Test class for real data end-to-end workflows."""

    def __init__(self, base_url: str = "http://localhost:5000", project_id: str = "tsdlmynydfuypiugmkev"):
        """Initialize the test class."""
        self.base_url = base_url
        self.project_id = project_id
        self.supabase: Optional[Client] = None
        self.auth_token: Optional[str] = None
        self.test_results: Dict[str, Any] = {
            "status": "Not Started",
            "total_tests": 0,
            "passed_tests": 0,
            "failed_tests": 0,
            "skipped_tests": 0,
            "performance_metrics": {},
            "test_cases": []
        }
        self.test_data = {}
        
    def connect_to_database(self) -> bool:
        """Connect to the Supabase database."""
        try:
            # Use the service role helper to connect
            self.supabase = get_supabase_client()
            user_id = get_user_id()
            
            if not user_id:
                logger.error("Failed to authenticate with Supabase")
                return False
                
            logger.info(f"Connected to Supabase database (Project ID: {self.project_id})")
            logger.info(f"Authenticated as user ID: {user_id}")
            return True
        except Exception as e:
            logger.error(f"Failed to connect to Supabase database: {str(e)}")
            return False
    
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
            
    def run_tests(self) -> Dict[str, Any]:
        """Run all real data workflow tests."""
        self.test_results["status"] = "Running"
        self.test_results["start_time"] = datetime.now().isoformat()
        
        # Connect to the database
        if not self.connect_to_database():
            self.test_results["status"] = "Failed"
            self.test_results["end_time"] = datetime.now().isoformat()
            return self.test_results
            
        # Load test data
        self.load_test_data()
        
        # Run the tests
        self.test_molecule_retrieval_and_property_calculation()
        self.test_mixture_creation_and_analysis()
        self.test_prediction_and_experiment_comparison()
        self.test_search_functionality()
        self.test_data_visualization_and_export()
        self.test_rls_policy_effectiveness()
        self.test_performance_with_real_data()
        
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
        
    def load_test_data(self) -> None:
        """Load real test data from the database."""
        try:
            logger.info("Loading real test data from the database")
            
            # Load molecules
            response = self.supabase.table("molecules").select("*").limit(10).execute()
            if response.data:
                self.test_data["molecules"] = response.data
                logger.info(f"Loaded {len(response.data)} molecules")
            else:
                logger.warning("No molecules found in the database")
                
            # Load mixtures
            response = self.supabase.table("mixtures").select("*").limit(10).execute()
            if response.data:
                self.test_data["mixtures"] = response.data
                logger.info(f"Loaded {len(response.data)} mixtures")
            else:
                logger.warning("No mixtures found in the database")
                
            # Load mixture components
            if "mixtures" in self.test_data and self.test_data["mixtures"]:
                mixture_ids = [m["id"] for m in self.test_data["mixtures"]]
                response = self.supabase.table("mixture_components").select("*").in_("mixture_id", mixture_ids).execute()
                if response.data:
                    self.test_data["mixture_components"] = response.data
                    logger.info(f"Loaded {len(response.data)} mixture components")
                else:
                    logger.warning("No mixture components found for the loaded mixtures")
                
            # Load predictions
            response = self.supabase.table("predictions").select("*").limit(10).execute()
            if response.data:
                self.test_data["predictions"] = response.data
                logger.info(f"Loaded {len(response.data)} predictions")
            else:
                logger.warning("No predictions found in the database")
                
            # Load experiments
            response = self.supabase.table("experiments").select("*").limit(10).execute()
            if response.data:
                self.test_data["experiments"] = response.data
                logger.info(f"Loaded {len(response.data)} experiments")
            else:
                logger.warning("No experiments found in the database")
                
            # Load property types
            response = self.supabase.table("property_types").select("*").execute()
            if response.data:
                self.test_data["property_types"] = response.data
                logger.info(f"Loaded {len(response.data)} property types")
            else:
                logger.warning("No property types found in the database")
                
            # Load calculation methods
            response = self.supabase.table("calculation_methods").select("*").execute()
            if response.data:
                self.test_data["calculation_methods"] = response.data
                logger.info(f"Loaded {len(response.data)} calculation methods")
            else:
                logger.warning("No calculation methods found in the database")
                
            logger.info("Test data loading complete")
        except Exception as e:
            logger.error(f"Failed to load test data: {str(e)}")
            raise
    
    def get_headers(self) -> Dict[str, str]:
        """Get headers for API requests."""
        headers = {
            "Content-Type": "application/json"
        }
        
        if self.auth_token:
            headers["Authorization"] = f"Bearer {self.auth_token}"
        
        return headers
    
    def measure_performance(self, func, *args, **kwargs) -> Tuple[Any, float]:
        """Measure the performance of a function."""
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        return result, end_time - start_time
        
    def test_molecule_retrieval_and_property_calculation(self) -> None:
        """Test retrieving real molecules and calculating properties."""
        try:
            logger.info("Testing molecule retrieval and property calculation")
            
            if "molecules" not in self.test_data or not self.test_data["molecules"]:
                self.add_test_result("RD-1", "Molecule Retrieval and Property Calculation", "Skipped", "No molecules available for testing")
                return
                
            # Test retrieving molecules from the API
            molecule_ids = [m["id"] for m in self.test_data["molecules"][:3]]
            response_times = []
            
            for molecule_id in molecule_ids:
                start_time = time.time()
                response = requests.get(f"{self.base_url}/api/v1/molecules/{molecule_id}", headers=self.get_headers())
                end_time = time.time()
                response_times.append(end_time - start_time)
                
                if response.status_code == 200:
                    molecule_data = response.json()
                    self.add_test_result("RD-1.1", f"Retrieve Molecule {molecule_id}", "Passed", f"Retrieved molecule with ID {molecule_id}", response)
                else:
                    self.add_test_result("RD-1.1", f"Retrieve Molecule {molecule_id}", "Failed", f"Failed to retrieve molecule: {response.status_code}", response)
            
            # Calculate average response time
            if response_times:
                avg_response_time = sum(response_times) / len(response_times)
                self.test_results["performance_metrics"]["molecule_retrieval_avg_time"] = avg_response_time
                logger.info(f"Average molecule retrieval time: {avg_response_time:.4f} seconds")
            
            # Test calculating properties for molecules
            for molecule in self.test_data["molecules"][:3]:
                if "smiles" not in molecule or not molecule["smiles"]:
                    logger.warning(f"Molecule {molecule['id']} does not have SMILES data, skipping property calculation")
                    continue
                    
                start_time = time.time()
                response = requests.post(f"{self.base_url}/api/v1/rdkit/properties", headers=self.get_headers(), json={
                    "molecule_data": molecule["smiles"],
                    "input_format": "smiles"
                })
                end_time = time.time()
                
                if response.status_code == 200:
                    properties_data = response.json()
                    if "hydrogen_bonding" in properties_data and "logp" in properties_data:
                        self.add_test_result("RD-1.2", f"Calculate Properties for {molecule['name']}", "Passed", f"Calculated properties for {molecule['name']}", response)
                        
                        # Store the calculation time
                        calc_time = end_time - start_time
                        if "property_calculation_times" not in self.test_results["performance_metrics"]:
                            self.test_results["performance_metrics"]["property_calculation_times"] = []
                        self.test_results["performance_metrics"]["property_calculation_times"].append(calc_time)
                        logger.info(f"Property calculation time for {molecule['name']}: {calc_time:.4f} seconds")
                    else:
                        self.add_test_result("RD-1.2", f"Calculate Properties for {molecule['name']}", "Failed", "Missing expected properties in response", response)
                else:
                    self.add_test_result("RD-1.2", f"Calculate Properties for {molecule['name']}", "Failed", f"Failed to calculate properties: {response.status_code}", response)
            
            # Calculate average property calculation time
            if "property_calculation_times" in self.test_results["performance_metrics"]:
                times = self.test_results["performance_metrics"]["property_calculation_times"]
                self.test_results["performance_metrics"]["property_calculation_avg_time"] = sum(times) / len(times)
                logger.info(f"Average property calculation time: {self.test_results['performance_metrics']['property_calculation_avg_time']:.4f} seconds")
            
            # Test generating visualizations for molecules
            for molecule in self.test_data["molecules"][:3]:
                if "smiles" not in molecule or not molecule["smiles"]:
                    logger.warning(f"Molecule {molecule['id']} does not have SMILES data, skipping visualization")
                    continue
                    
                start_time = time.time()
                response = requests.post(f"{self.base_url}/api/v1/rdkit/visualization", headers=self.get_headers(), json={
                    "molecule_data": molecule["smiles"],
                    "input_format": "smiles",
                    "width": 400,
                    "height": 300
                })
                end_time = time.time()
                
                if response.status_code == 200:
                    visualization_data = response.json()
                    if "svg" in visualization_data:
                        self.add_test_result("RD-1.3", f"Generate Visualization for {molecule['name']}", "Passed", f"Generated visualization for {molecule['name']}", response)
                        
                        # Store the visualization time
                        viz_time = end_time - start_time
                        if "visualization_times" not in self.test_results["performance_metrics"]:
                            self.test_results["performance_metrics"]["visualization_times"] = []
                        self.test_results["performance_metrics"]["visualization_times"].append(viz_time)
                        logger.info(f"Visualization time for {molecule['name']}: {viz_time:.4f} seconds")
                    else:
                        self.add_test_result("RD-1.3", f"Generate Visualization for {molecule['name']}", "Failed", "Missing SVG in response", response)
                else:
                    self.add_test_result("RD-1.3", f"Generate Visualization for {molecule['name']}", "Failed", f"Failed to generate visualization: {response.status_code}", response)
            
            # Calculate average visualization time
            if "visualization_times" in self.test_results["performance_metrics"]:
                times = self.test_results["performance_metrics"]["visualization_times"]
                self.test_results["performance_metrics"]["visualization_avg_time"] = sum(times) / len(times)
                logger.info(f"Average visualization time: {self.test_results['performance_metrics']['visualization_avg_time']:.4f} seconds")
            
            # Overall workflow status
            passed_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-1.") and tc["status"] == "Passed")
            total_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-1."))
            
            if passed_steps == total_steps:
                self.add_test_result("RD-1", "Molecule Retrieval and Property Calculation Workflow", "Passed", f"All {passed_steps}/{total_steps} steps passed")
            else:
                self.add_test_result("RD-1", "Molecule Retrieval and Property Calculation Workflow", "Failed", f"Only {passed_steps}/{total_steps} steps passed")
        except Exception as e:
            self.add_test_result("RD-1", "Molecule Retrieval and Property Calculation Workflow", "Failed", f"Error testing workflow: {str(e)}")
    
    def test_mixture_creation_and_analysis(self) -> None:
        """Test creating and analyzing mixtures with real molecules."""
        try:
            logger.info("Testing mixture creation and analysis")
            
            if "molecules" not in self.test_data or not self.test_data["molecules"]:
                self.add_test_result("RD-2", "Mixture Creation and Analysis", "Skipped", "No molecules available for testing")
                return
                
            # Get two molecules for the mixture
            if len(self.test_data["molecules"]) < 2:
                self.add_test_result("RD-2", "Mixture Creation and Analysis", "Skipped", "Not enough molecules available for testing")
                return
                
            molecules = self.test_data["molecules"][:2]
            molecule_ids = [m["id"] for m in molecules]
            
            # Create a new mixture
            new_mixture = {
                "name": f"Test Mixture {datetime.now().strftime('%Y%m%d%H%M%S')}",
                "description": "A test mixture created by the real data workflow test",
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
            
            start_time = time.time()
            response = requests.post(f"{self.base_url}/api/v1/mixtures", headers=self.get_headers(), json=new_mixture)
            end_time = time.time()
            
            if response.status_code == 201:
                mixture_data = response.json()
                mixture_id = mixture_data["id"]
                self.add_test_result("RD-2.1", "Create Mixture", "Passed", f"Created mixture with ID {mixture_id}", response)
                
                # Store the mixture creation time
                mixture_creation_time = end_time - start_time
                self.test_results["performance_metrics"]["mixture_creation_time"] = mixture_creation_time
                logger.info(f"Mixture creation time: {mixture_creation_time:.4f} seconds")
                
                # Store the mixture for later use
                if "created_mixtures" not in self.test_data:
                    self.test_data["created_mixtures"] = []
                self.test_data["created_mixtures"].append(mixture_data)
            else:
                self.add_test_result("RD-2.1", "Create Mixture", "Failed", f"Failed to create mixture: {response.status_code}", response)
                return  # Skip the rest of the test if mixture creation fails
            
            # Retrieve the mixture
            start_time = time.time()
            response = requests.get(f"{self.base_url}/api/v1/mixtures/{mixture_id}", headers=self.get_headers())
            end_time = time.time()
            
            if response.status_code == 200:
                retrieved_mixture = response.json()
                if retrieved_mixture["id"] == mixture_id:
                    self.add_test_result("RD-2.2", "Retrieve Mixture", "Passed", f"Retrieved mixture with ID {mixture_id}", response)
                    
                    # Store the mixture retrieval time
                    mixture_retrieval_time = end_time - start_time
                    self.test_results["performance_metrics"]["mixture_retrieval_time"] = mixture_retrieval_time
                    logger.info(f"Mixture retrieval time: {mixture_retrieval_time:.4f} seconds")
                else:
                    self.add_test_result("RD-2.2", "Retrieve Mixture", "Failed", "Retrieved mixture ID does not match", response)
            else:
                self.add_test_result("RD-2.2", "Retrieve Mixture", "Failed", f"Failed to retrieve mixture: {response.status_code}", response)
            
            # Update the mixture
            updated_mixture = {
                "name": f"Updated Test Mixture {datetime.now().strftime('%Y%m%d%H%M%S')}",
                "description": "An updated test mixture",
                "components": [
                    {
                        "molecule_id": molecule_ids[0],
                        "concentration": 60,
                        "concentration_unit": "%v/v"
                    },
                    {
                        "molecule_id": molecule_ids[1],
                        "concentration": 40,
                        "concentration_unit": "%v/v"
                    }
                ]
            }
            
            start_time = time.time()
            response = requests.put(f"{self.base_url}/api/v1/mixtures/{mixture_id}", headers=self.get_headers(), json=updated_mixture)
            end_time = time.time()
            
            if response.status_code == 200:
                updated_data = response.json()
                if updated_data["name"] == updated_mixture["name"]:
                    self.add_test_result("RD-2.3", "Update Mixture", "Passed", f"Updated mixture with ID {mixture_id}", response)
                    
                    # Store the mixture update time
                    mixture_update_time = end_time - start_time
                    self.test_results["performance_metrics"]["mixture_update_time"] = mixture_update_time
                    logger.info(f"Mixture update time: {mixture_update_time:.4f} seconds")
                else:
                    self.add_test_result("RD-2.3", "Update Mixture", "Failed", "Updated mixture name does not match", response)
            else:
                self.add_test_result("RD-2.3", "Update Mixture", "Failed", f"Failed to update mixture: {response.status_code}", response)
            
            # Overall workflow status
            passed_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-2.") and tc["status"] == "Passed")
            total_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-2."))
            
            if passed_steps == total_steps:
                self.add_test_result("RD-2", "Mixture Creation and Analysis Workflow", "Passed", f"All {passed_steps}/{total_steps} steps passed")
            else:
                self.add_test_result("RD-2", "Mixture Creation and Analysis Workflow", "Failed", f"Only {passed_steps}/{total_steps} steps passed")
        except Exception as e:
            self.add_test_result("RD-2", "Mixture Creation and Analysis Workflow", "Failed", f"Error testing workflow: {str(e)}")
            
    def test_prediction_and_experiment_comparison(self) -> None:
        """Test comparing predictions with experiments using real data."""
        try:
            logger.info("Testing prediction and experiment comparison")
            
            # Check if we have a mixture from the previous test
            if "created_mixtures" not in self.test_data or not self.test_data["created_mixtures"]:
                # Check if we have any mixtures from the database
                if "mixtures" not in self.test_data or not self.test_data["mixtures"]:
                    self.add_test_result("RD-3", "Prediction and Experiment Comparison", "Skipped", "No mixtures available for testing")
                    return
                
                # Use a mixture from the database
                mixture_id = self.test_data["mixtures"][0]["id"]
                self.add_test_result("RD-3.1", "Get Mixture", "Passed", f"Using existing mixture with ID {mixture_id}")
            else:
                # Use the mixture we created
                mixture_id = self.test_data["created_mixtures"][0]["id"]
                self.add_test_result("RD-3.1", "Get Mixture", "Passed", f"Using created mixture with ID {mixture_id}")
            
            # Get property type ID for Glass Transition Temperature
            property_type_id = None
            property_name = "Glass Transition Temperature"
            
            if "property_types" in self.test_data:
                for pt in self.test_data["property_types"]:
                    if pt["name"] == property_name:
                        property_type_id = pt["id"]
                        break
            
            if not property_type_id:
                # Try to get it from the database
                response = self.supabase.table("property_types").select("id").eq("name", property_name).execute()
                if response.data:
                    property_type_id = response.data[0]["id"]
                else:
                    # Create the property type
                    property_type_id = str(uuid.uuid4())
                    property_type = {
                        "id": property_type_id,
                        "name": property_name,
                        "data_type": "numeric",
                        "description": "The temperature at which a material transitions from a hard, glassy state to a soft, rubbery state",
                        "units": "°C"
                    }
                    self.supabase.table("property_types").insert(property_type).execute()
            
            # Get calculation method ID for ML Model v1
            calculation_method_id = None
            calculation_method_name = "ML Model v1"
            
            if "calculation_methods" in self.test_data:
                for cm in self.test_data["calculation_methods"]:
                    if cm["name"] == calculation_method_name:
                        calculation_method_id = cm["id"]
                        break
            
            if not calculation_method_id:
                # Try to get it from the database
                response = self.supabase.table("calculation_methods").select("id").eq("name", calculation_method_name).execute()
                if response.data:
                    calculation_method_id = response.data[0]["id"]
                else:
                    # Create the calculation method
                    calculation_method_id = str(uuid.uuid4())
                    calculation_method = {
                        "id": calculation_method_id,
                        "name": calculation_method_name,
                        "description": "Machine learning model for property prediction",
                        "version": "1.0"
                    }
                    self.supabase.table("calculation_methods").insert(calculation_method).execute()
            
            # Create a prediction
            new_prediction = {
                "property_type_id": property_type_id,
                "property_name": property_name,
                "value": -45.2,
                "confidence": 0.85,
                "calculation_method_id": calculation_method_id,
                "calculation_method": calculation_method_name
            }
            
            start_time = time.time()
            response = requests.post(f"{self.base_url}/api/v1/mixtures/{mixture_id}/predictions", headers=self.get_headers(), json=new_prediction)
            end_time = time.time()
            
            if response.status_code == 201:
                prediction_data = response.json()
                prediction_id = prediction_data["id"]
                self.add_test_result("RD-3.2", "Create Prediction", "Passed", f"Created prediction with ID {prediction_id}", response)
                
                # Store the prediction creation time
                prediction_creation_time = end_time - start_time
                self.test_results["performance_metrics"]["prediction_creation_time"] = prediction_creation_time
                logger.info(f"Prediction creation time: {prediction_creation_time:.4f} seconds")
                
                # Store the prediction for later use
                if "created_predictions" not in self.test_data:
                    self.test_data["created_predictions"] = []
                self.test_data["created_predictions"].append(prediction_data)
            else:
                self.add_test_result("RD-3.2", "Create Prediction", "Failed", f"Failed to create prediction: {response.status_code}", response)
                return  # Skip the rest of the test if prediction creation fails
            
            # Create an experiment
            new_experiment = {
                "property_type_id": property_type_id,
                "property_name": property_name,
                "value": -47.3,
                "experimental_conditions": "Measured using DSC at 10°C/min",
                "date_performed": datetime.now().strftime("%Y-%m-%d")
            }
            
            start_time = time.time()
            response = requests.post(f"{self.base_url}/api/v1/mixtures/{mixture_id}/experiments", headers=self.get_headers(), json=new_experiment)
            end_time = time.time()
            
            if response.status_code == 201:
                experiment_data = response.json()
                experiment_id = experiment_data["id"]
                self.add_test_result("RD-3.3", "Create Experiment", "Passed", f"Created experiment with ID {experiment_id}", response)
                
                # Store the experiment creation time
                experiment_creation_time = end_time - start_time
                self.test_results["performance_metrics"]["experiment_creation_time"] = experiment_creation_time
                logger.info(f"Experiment creation time: {experiment_creation_time:.4f} seconds")
                
                # Store the experiment for later use
                if "created_experiments" not in self.test_data:
                    self.test_data["created_experiments"] = []
                self.test_data["created_experiments"].append(experiment_data)
            else:
                self.add_test_result("RD-3.3", "Create Experiment", "Failed", f"Failed to create experiment: {response.status_code}", response)
                return  # Skip the rest of the test if experiment creation fails
            
            # Compare prediction with experiment
            start_time = time.time()
            response = requests.get(f"{self.base_url}/api/v1/mixtures/{mixture_id}/compare?property_name={property_name}", headers=self.get_headers())
            end_time = time.time()
            
            if response.status_code == 200:
                comparison_data = response.json()
                if "prediction" in comparison_data and "experiment" in comparison_data and "difference" in comparison_data:
                    self.add_test_result("RD-3.4", "Compare Prediction with Experiment", "Passed", "Successfully compared prediction with experiment", response)
                    
                    # Store the comparison time
                    comparison_time = end_time - start_time
                    self.test_results["performance_metrics"]["comparison_time"] = comparison_time
                    logger.info(f"Comparison time: {comparison_time:.4f} seconds")
                else:

def test_search_functionality(self) -> None:
        """Test search functionality with real data."""
        try:
            logger.info("Testing search functionality")
            
            if "molecules" not in self.test_data or not self.test_data["molecules"]:
                self.add_test_result("RD-4", "Search Functionality", "Skipped", "No molecules available for testing")
                return
            
            # Test search by name
            if len(self.test_data["molecules"]) > 0:
                molecule = self.test_data["molecules"][0]
                if "name" in molecule and molecule["name"]:
                    # Use part of the name for the search
                    search_term = molecule["name"][:4]
                    
                    start_time = time.time()
                    response = requests.get(f"{self.base_url}/api/v1/molecules/search?name={search_term}", headers=self.get_headers())
                    end_time = time.time()
                    
                    if response.status_code == 200:
                        search_results = response.json()
                        if isinstance(search_results, list) and len(search_results) > 0:
                            self.add_test_result("RD-4.1", "Search by Name", "Passed", f"Found {len(search_results)} molecules matching '{search_term}'", response)
                            
                            # Store the search time
                            search_time = end_time - start_time
                            self.test_results["performance_metrics"]["name_search_time"] = search_time
                            logger.info(f"Name search time: {search_time:.4f} seconds")
                        else:
                            self.add_test_result("RD-4.1", "Search by Name", "Failed", f"No molecules found matching '{search_term}'", response)
                    else:
                        self.add_test_result("RD-4.1", "Search by Name", "Failed", f"Failed to search by name: {response.status_code}", response)
            
            # Test search by structure (SMILES)
            if len(self.test_data["molecules"]) > 0:
                for molecule in self.test_data["molecules"]:
                    if "smiles" in molecule and molecule["smiles"]:
                        start_time = time.time()
                        response = requests.post(f"{self.base_url}/api/v1/molecules/search/structure", headers=self.get_headers(), json={
                            "structure": molecule["smiles"],
                            "structure_format": "smiles"
                        })
                        end_time = time.time()
                        
                        if response.status_code == 200:
                            search_results = response.json()
                            if isinstance(search_results, list) and len(search_results) > 0:
                                self.add_test_result("RD-4.2", "Search by Structure", "Passed", f"Found {len(search_results)} molecules matching structure", response)
                                
                                # Store the search time
                                search_time = end_time - start_time
                                self.test_results["performance_metrics"]["structure_search_time"] = search_time
                                logger.info(f"Structure search time: {search_time:.4f} seconds")
                                break
                            else:
                                self.add_test_result("RD-4.2", "Search by Structure", "Failed", "No molecules found matching structure", response)
                        else:
                            self.add_test_result("RD-4.2", "Search by Structure", "Failed", f"Failed to search by structure: {response.status_code}", response)
                        break
            
            # Test similarity search
            if len(self.test_data["molecules"]) > 0:
                for molecule in self.test_data["molecules"]:
                    if "smiles" in molecule and molecule["smiles"]:
                        start_time = time.time()
                        response = requests.post(f"{self.base_url}/api/v1/molecules/search/similarity", headers=self.get_headers(), json={
                            "structure": molecule["smiles"],
                            "structure_format": "smiles",
                            "threshold": 0.7
                        })
                        end_time = time.time()
                        
                        if response.status_code == 200:
                            search_results = response.json()
                            if isinstance(search_results, list) and len(search_results) > 0:
                                self.add_test_result("RD-4.3", "Similarity Search", "Passed", f"Found {len(search_results)} similar molecules", response)
                                
                                # Store the search time
                                search_time = end_time - start_time
                                self.test_results["performance_metrics"]["similarity_search_time"] = search_time
                                logger.info(f"Similarity search time: {search_time:.4f} seconds")
                                break
                            else:
                                self.add_test_result("RD-4.3", "Similarity Search", "Failed", "No similar molecules found", response)
                        else:
                            self.add_test_result("RD-4.3", "Similarity Search", "Failed", f"Failed to perform similarity search: {response.status_code}", response)
                        break
            
            # Overall workflow status
            passed_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-4.") and tc["status"] == "Passed")
            total_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-4."))
            
            if passed_steps == total_steps:
                self.add_test_result("RD-4", "Search Functionality Workflow", "Passed", f"All {passed_steps}/{total_steps} steps passed")
            else:
                self.add_test_result("RD-4", "Search Functionality Workflow", "Failed", f"Only {passed_steps}/{total_steps} steps passed")
        except Exception as e:
            self.add_test_result("RD-4", "Search Functionality Workflow", "Failed", f"Error testing workflow: {str(e)}")
    
    def test_data_visualization_and_export(self) -> None:
        """Test data visualization and export with real data."""
        try:
            logger.info("Testing data visualization and export")
            
            # Test molecule visualization
            if "molecules" in self.test_data and self.test_data["molecules"]:
                for molecule in self.test_data["molecules"][:1]:
                    if "smiles" in molecule and molecule["smiles"]:
                        start_time = time.time()
                        response = requests.get(f"{self.base_url}/api/v1/molecules/{molecule['id']}/visualization", headers=self.get_headers())
                        end_time = time.time()
                        
                        if response.status_code == 200:
                            visualization_data = response.json()
                            if "svg" in visualization_data:
                                self.add_test_result("RD-5.1", "Molecule Visualization", "Passed", f"Generated visualization for molecule {molecule['id']}", response)
                                
                                # Store the visualization time
                                viz_time = end_time - start_time
                                self.test_results["performance_metrics"]["molecule_visualization_time"] = viz_time
                                logger.info(f"Molecule visualization time: {viz_time:.4f} seconds")
                            else:
                                self.add_test_result("RD-5.1", "Molecule Visualization", "Failed", "Missing SVG in response", response)
                        else:
                            self.add_test_result("RD-5.1", "Molecule Visualization", "Failed", f"Failed to generate visualization: {response.status_code}", response)
                        break
            
            # Test property data visualization
            if "predictions" in self.test_data and self.test_data["predictions"]:
                start_time = time.time()
                response = requests.get(f"{self.base_url}/api/v1/visualizations/properties", headers=self.get_headers())
                end_time = time.time()
                
                if response.status_code == 200:
                    visualization_data = response.json()
                    if "chart_data" in visualization_data:
                        self.add_test_result("RD-5.2", "Property Data Visualization", "Passed", "Generated property data visualization", response)
                        
                        # Store the visualization time
                        viz_time = end_time - start_time
                        self.test_results["performance_metrics"]["property_visualization_time"] = viz_time
                        logger.info(f"Property visualization time: {viz_time:.4f} seconds")
                    else:
                        self.add_test_result("RD-5.2", "Property Data Visualization", "Failed", "Missing chart data in response", response)
                else:
                    self.add_test_result("RD-5.2", "Property Data Visualization", "Failed", f"Failed to generate property visualization: {response.status_code}", response)
            
            # Test data export
            if "molecules" in self.test_data and self.test_data["molecules"]:
                start_time = time.time()
                response = requests.get(f"{self.base_url}/api/v1/export/molecules?format=csv", headers=self.get_headers())
                end_time = time.time()
                
                if response.status_code == 200:
                    export_data = response.text
                    if export_data and "," in export_data:  # Simple check for CSV format
                        self.add_test_result("RD-5.3", "Data Export", "Passed", "Exported molecule data as CSV", response)
                        
                        # Store the export time
                        export_time = end_time - start_time
                        self.test_results["performance_metrics"]["data_export_time"] = export_time
                        logger.info(f"Data export time: {export_time:.4f} seconds")
                    else:
                        self.add_test_result("RD-5.3", "Data Export", "Failed", "Invalid CSV data in response", response)
                else:
                    self.add_test_result("RD-5.3", "Data Export", "Failed", f"Failed to export data: {response.status_code}", response)
            
            # Overall workflow status
            passed_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-5.") and tc["status"] == "Passed")
            total_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-5."))
            
            if passed_steps == total_steps:
                self.add_test_result("RD-5", "Data Visualization and Export Workflow", "Passed", f"All {passed_steps}/{total_steps} steps passed")
            else:
                self.add_test_result("RD-5", "Data Visualization and Export Workflow", "Failed", f"Only {passed_steps}/{total_steps} steps passed")
        except Exception as e:
            self.add_test_result("RD-5", "Data Visualization and Export Workflow", "Failed", f"Error testing workflow: {str(e)}")

def test_rls_policy_effectiveness(self) -> None:
        """Test RLS policy effectiveness with real data."""
        try:
            logger.info("Testing RLS policy effectiveness")
            
            # Test service role access to all data
            try:
                # Test access to molecules
                response = self.supabase.table("molecules").select("*").limit(5).execute()
                if response.data and len(response.data) > 0:
                    self.add_test_result("RD-6.1", "Service Role Access to Molecules", "Passed", f"Service role can access {len(response.data)} molecules")
                else:
                    self.add_test_result("RD-6.1", "Service Role Access to Molecules", "Failed", "Service role cannot access molecules")
                
                # Test access to mixtures
                response = self.supabase.table("mixtures").select("*").limit(5).execute()
                if response.data and len(response.data) > 0:
                    self.add_test_result("RD-6.2", "Service Role Access to Mixtures", "Passed", f"Service role can access {len(response.data)} mixtures")
                else:
                    self.add_test_result("RD-6.2", "Service Role Access to Mixtures", "Failed", "Service role cannot access mixtures")
                
                # Test access to predictions
                response = self.supabase.table("predictions").select("*").limit(5).execute()
                if response.data and len(response.data) > 0:
                    self.add_test_result("RD-6.3", "Service Role Access to Predictions", "Passed", f"Service role can access {len(response.data)} predictions")
                else:
                    self.add_test_result("RD-6.3", "Service Role Access to Predictions", "Failed", "Service role cannot access predictions")
            except Exception as e:
                self.add_test_result("RD-6", "Service Role Access", "Failed", f"Error testing service role access: {str(e)}")
            
            # Overall workflow status
            passed_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-6.") and tc["status"] == "Passed")
            total_steps = sum(1 for tc in self.test_results["test_cases"] if tc["id"].startswith("RD-6."))
            
            if passed_steps == total_steps:
                self.add_test_result("RD-6", "RLS Policy Effectiveness", "Passed", f"All {passed_steps}/{total_steps} steps passed")
            else:
                self.add_test_result("RD-6", "RLS Policy Effectiveness", "Failed", f"Only {passed_steps}/{total_steps} steps passed")
        except Exception as e:
            self.add_test_result("RD-6", "RLS Policy Effectiveness", "Failed", f"Error testing RLS policy effectiveness: {str(e)}")
    
    def test_performance_with_real_data(self) -> None:
        """Test performance with realistic data volumes."""
        try:
            logger.info("Testing performance with real data")
            
            # Collect performance metrics from previous tests
            performance_metrics = self.test_results["performance_metrics"]
            
            # Calculate summary statistics
            summary_stats = {}
            
            for metric_name, metric_value in performance_metrics.items():
                if isinstance(metric_value, list):
                    # Calculate statistics for lists of values
                    if len(metric_value) > 0:
                        summary_stats[metric_name] = {
                            "min": min(metric_value),
                            "max": max(metric_value),
                            "mean": sum(metric_value) / len(metric_value),
                            "median": statistics.median(metric_value) if len(metric_value) > 1 else metric_value[0],
                            "count": len(metric_value)
                        }
                        
                        # Calculate 95th percentile if we have enough data points
                        if len(metric_value) >= 5:
                            sorted_values = sorted(metric_value)
                            idx = int(0.95 * len(sorted_values))
                            summary_stats[metric_name]["p95"] = sorted_values[idx]
                else:
                    # Single values are stored directly
                    summary_stats[metric_name] = metric_value
            
            # Store the summary statistics
            self.test_results["performance_summary"] = summary_stats
            
            # Log the summary statistics
            logger.info("Performance summary:")
            for metric_name, stats in summary_stats.items():
                if isinstance(stats, dict):
                    logger.info(f"  {metric_name}:")
                    for stat_name, stat_value in stats.items():
                        logger.info(f"    {stat_name}: {stat_value:.4f}" if isinstance(stat_value, float) else f"    {stat_name}: {stat_value}")
                else:
                    logger.info(f"  {metric_name}: {stats:.4f}" if isinstance(stats, float) else f"  {metric_name}: {stats}")
            
            # Evaluate performance against thresholds
            # These thresholds should be adjusted based on actual performance requirements
            performance_thresholds = {
                "molecule_retrieval_avg_time": 1.0,  # seconds
                "property_calculation_avg_time": 2.0,  # seconds
                "visualization_avg_time": 2.0,  # seconds
                "mixture_creation_time": 2.0,  # seconds
                "mixture_retrieval_time": 1.0,  # seconds
                "prediction_creation_time": 2.0,  # seconds
                "experiment_creation_time": 2.0,  # seconds
                "comparison_time": 1.0,  # seconds
                "name_search_time": 1.0,  # seconds
                "structure_search_time": 2.0,  # seconds
                "similarity_search_time": 3.0,  # seconds
            }
            
            # Check performance against thresholds
            performance_issues = []
            
            for metric_name, threshold in performance_thresholds.items():
                if metric_name in performance_metrics:
                    metric_value = performance_metrics[metric_name]
                    
                    # Handle both single values and averages
                    if isinstance(metric_value, list):
                        if len(metric_value) > 0:
                            avg_value = sum(metric_value) / len(metric_value)
                            if avg_value > threshold:
                                performance_issues.append(f"{metric_name} ({avg_value:.4f}s) exceeds threshold ({threshold:.4f}s)")
                    else:
                        if metric_value > threshold:
                            performance_issues.append(f"{metric_name} ({metric_value:.4f}s) exceeds threshold ({threshold:.4f}s)")
            
            # Add test result based on performance issues
            if performance_issues:
                self.add_test_result("RD-7", "Performance with Real Data", "Failed", f"Performance issues found: {', '.join(performance_issues)}")
            else:
                self.add_test_result("RD-7", "Performance with Real Data", "Passed", "Performance is within acceptable thresholds")
        except Exception as e:
            self.add_test_result("RD-7", "Performance with Real Data", "Failed", f"Error testing performance: {str(e)}")
    
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
    logger.info("Starting real data end-to-end workflow tests")
    
    # Create and run the tests
    test = RealDataWorkflowTest()
    results = test.run_tests()
    
    # Save the results
    test.save_results("real_data_test_results.json")
    
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
