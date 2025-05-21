#!/usr/bin/env python3
"""
CryoProtect Production Workflow Test

This script runs production-grade tests for CryoProtect, focusing on real data workflows.
It uses the test data loaded from core_cryoprotectants.json, mixtures.json, and
edge_cases.json to verify proper functioning of the system in a production-like environment.
"""

import os
import sys
import json
import time
import logging
import requests
import uuid
import argparse
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("production_workflow_test.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class ProductionWorkflowTest:
    """Production workflow test for CryoProtect."""
    
    def __init__(self, app_url: str = "http://localhost:5001", rdkit_url: str = "http://localhost:5002"):
        """
        Initialize the test class.
        
        Args:
            app_url: URL of the main application API
            rdkit_url: URL of the RDKit service API
        """
        self.app_url = app_url
        self.rdkit_url = rdkit_url
        self.test_data = {}
        self.test_results = {
            "status": "Not Started",
            "start_time": None,
            "end_time": None,
            "total_tests": 0,
            "passed_tests": 0,
            "failed_tests": 0,
            "skipped_tests": 0,
            "performance_metrics": {},
            "test_cases": []
        }
    
    def load_test_data(self) -> None:
        """Load test data from JSON files."""
        logger.info("Loading test data...")
        
        # Define the test data files
        test_data_dir = os.path.join("tests", "test_data")
        files = {
            "core_cryoprotectants": os.path.join(test_data_dir, "core_cryoprotectants.json"),
            "mixtures": os.path.join(test_data_dir, "mixtures.json"),
            "edge_cases": os.path.join(test_data_dir, "edge_cases.json")
        }
        
        # Load each file
        for data_type, file_path in files.items():
            try:
                if os.path.exists(file_path):
                    with open(file_path, 'r') as f:
                        self.test_data[data_type] = json.load(f)
                    logger.info(f"Loaded {data_type} test data from {file_path}")
                else:
                    logger.warning(f"Test data file not found: {file_path}")
            except Exception as e:
                logger.error(f"Error loading test data from {file_path}: {str(e)}")
        
        # Prepare molecules list for testing
        self.test_data["molecules"] = []
        if "core_cryoprotectants" in self.test_data and "molecules" in self.test_data["core_cryoprotectants"]:
            self.test_data["molecules"].extend(self.test_data["core_cryoprotectants"]["molecules"])
        if "edge_cases" in self.test_data and "molecules" in self.test_data["edge_cases"]:
            self.test_data["molecules"].extend(self.test_data["edge_cases"]["molecules"])
        
        logger.info(f"Loaded {len(self.test_data.get('molecules', []))} molecules for testing")
    
    def add_test_result(self, test_id: str, test_name: str, status: str, message: str, response=None, duration=None) -> None:
        """
        Add a test result.
        
        Args:
            test_id: Identifier for the test
            test_name: Human-readable name of the test
            status: "Passed", "Failed", or "Skipped"
            message: Description of the test result
            response: Optional response object
            duration: Optional test duration in seconds
        """
        result = {
            "id": test_id,
            "name": test_name,
            "status": status,
            "message": message,
            "timestamp": datetime.now().isoformat()
        }
        
        if duration is not None:
            result["duration"] = duration
        
        if response:
            try:
                if hasattr(response, 'json'):
                    result["response"] = {
                        "status_code": response.status_code,
                        "content_type": response.headers.get("content-type", ""),
                        "content": response.json() if response.headers.get("content-type") == "application/json" else None
                    }
                else:
                    result["response"] = {
                        "content": response
                    }
            except Exception as e:
                result["response"] = {
                    "error": str(e)
                }
        
        self.test_results["test_cases"].append(result)
        
        # Update counters
        if status == "Passed":
            self.test_results["passed_tests"] += 1
            logger.info(f"Test {test_id} - {test_name}: PASSED {f'in {duration:.2f}s' if duration else ''}")
        elif status == "Failed":
            self.test_results["failed_tests"] += 1
            logger.error(f"Test {test_id} - {test_name}: FAILED - {message} {f'in {duration:.2f}s' if duration else ''}")
        else:  # Skipped
            self.test_results["skipped_tests"] += 1
            logger.warning(f"Test {test_id} - {test_name}: SKIPPED - {message}")
        
        self.test_results["total_tests"] = len(self.test_results["test_cases"])
    
    def measure_execution_time(self, func, *args, **kwargs) -> Tuple[Any, float]:
        """
        Measure execution time of a function.
        
        Args:
            func: Function to measure
            *args: Arguments to pass to the function
            **kwargs: Keyword arguments to pass to the function
            
        Returns:
            Tuple of (function result, execution time in seconds)
        """
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        return result, end_time - start_time
    
    def check_service_health(self) -> bool:
        """
        Check if the application and RDKit services are healthy.
        
        Returns:
            True if both services are healthy, False otherwise
        """
        try:
            # Check app health
            response = requests.get(f"{self.app_url}/health", timeout=10)
            app_healthy = response.status_code == 200
            
            if app_healthy:
                self.add_test_result("HEALTH-1", "App Health Check", "Passed", "App service is healthy", response)
            else:
                self.add_test_result("HEALTH-1", "App Health Check", "Failed", f"App service returned status code {response.status_code}", response)
            
            # Check RDKit health
            response = requests.get(f"{self.rdkit_url}/health", timeout=10)
            rdkit_healthy = response.status_code == 200
            
            if rdkit_healthy:
                self.add_test_result("HEALTH-2", "RDKit Health Check", "Passed", "RDKit service is healthy", response)
            else:
                self.add_test_result("HEALTH-2", "RDKit Health Check", "Failed", f"RDKit service returned status code {response.status_code}", response)
            
            # Check RDKit integration
            response = requests.get(f"{self.app_url}/rdkit/check", timeout=10)
            integration_healthy = response.status_code == 200
            
            if integration_healthy:
                self.add_test_result("HEALTH-3", "RDKit Integration Check", "Passed", "RDKit integration is working", response)
            else:
                self.add_test_result("HEALTH-3", "RDKit Integration Check", "Failed", f"RDKit integration check returned status code {response.status_code}", response)
            
            return app_healthy and rdkit_healthy and integration_healthy
            
        except requests.exceptions.RequestException as e:
            self.add_test_result("HEALTH", "Service Health Check", "Failed", f"Error connecting to services: {str(e)}")
            return False
        except Exception as e:
            self.add_test_result("HEALTH", "Service Health Check", "Failed", f"Unexpected error: {str(e)}")
            return False
    
    def test_molecule_property_calculation(self) -> None:
        """Test molecule property calculation."""
        logger.info("Testing molecule property calculation...")
        
        if not self.test_data.get("molecules"):
            self.add_test_result("PROP-1", "Molecule Property Calculation", "Skipped", "No test molecules available")
            return
        
        # Test molecule property calculation with a few molecules
        response_times = []
        for i, molecule in enumerate(self.test_data["molecules"][:5]):
            try:
                if "smiles" not in molecule:
                    logger.warning(f"Molecule {i} has no SMILES data, skipping")
                    continue
                
                # Call the property calculation endpoint
                response, duration = self.measure_execution_time(
                    requests.get,
                    f"{self.app_url}/molecule/{molecule['smiles']}",
                    timeout=30
                )
                
                response_times.append(duration)
                
                if response.status_code == 200:
                    response_json = response.json()
                    if "properties" in response_json:
                        properties = response_json["properties"]
                        property_keys = properties.keys()
                        expected_keys = ["logp", "molecular_weight", "h_acceptors", "h_donors"]
                        
                        # Check if all expected properties are present
                        if all(key in property_keys for key in expected_keys):
                            self.add_test_result(
                                f"PROP-1.{i+1}",
                                f"Calculate Properties for {molecule.get('name', 'Unknown Molecule')}",
                                "Passed",
                                f"Successfully calculated properties in {duration:.2f}s",
                                response,
                                duration
                            )
                        else:
                            missing_keys = [key for key in expected_keys if key not in property_keys]
                            self.add_test_result(
                                f"PROP-1.{i+1}",
                                f"Calculate Properties for {molecule.get('name', 'Unknown Molecule')}",
                                "Failed",
                                f"Missing expected properties: {', '.join(missing_keys)}",
                                response,
                                duration
                            )
                    else:
                        self.add_test_result(
                            f"PROP-1.{i+1}",
                            f"Calculate Properties for {molecule.get('name', 'Unknown Molecule')}",
                            "Failed",
                            "Response does not contain properties",
                            response,
                            duration
                        )
                else:
                    self.add_test_result(
                        f"PROP-1.{i+1}",
                        f"Calculate Properties for {molecule.get('name', 'Unknown Molecule')}",
                        "Failed",
                        f"API returned status code {response.status_code}",
                        response,
                        duration
                    )
            except Exception as e:
                self.add_test_result(
                    f"PROP-1.{i+1}",
                    f"Calculate Properties for {molecule.get('name', 'Unknown Molecule')}",
                    "Failed",
                    f"Error: {str(e)}",
                    None,
                    None
                )
        
        # Record performance metrics
        if response_times:
            avg_time = sum(response_times) / len(response_times)
            max_time = max(response_times)
            min_time = min(response_times)
            
            self.test_results["performance_metrics"]["property_calculation"] = {
                "average": avg_time,
                "min": min_time,
                "max": max_time,
                "count": len(response_times)
            }
            
            logger.info(f"Property calculation performance: avg={avg_time:.2f}s, min={min_time:.2f}s, max={max_time:.2f}s")
    
    def test_molecule_visualization(self) -> None:
        """Test molecule visualization."""
        logger.info("Testing molecule visualization...")
        
        if not self.test_data.get("molecules"):
            self.add_test_result("VIZ-1", "Molecule Visualization", "Skipped", "No test molecules available")
            return
        
        # Test molecule visualization with a few molecules
        response_times = []
        for i, molecule in enumerate(self.test_data["molecules"][:3]):
            try:
                if "smiles" not in molecule:
                    logger.warning(f"Molecule {i} has no SMILES data, skipping")
                    continue
                
                # Call the visualization endpoint
                response, duration = self.measure_execution_time(
                    requests.post,
                    f"{self.app_url}/api/v1/rdkit/visualization",
                    json={
                        "molecule_data": molecule["smiles"],
                        "input_format": "smiles",
                        "width": 400,
                        "height": 300
                    },
                    timeout=30
                )
                
                response_times.append(duration)
                
                if response.status_code == 200:
                    response_json = response.json()
                    if "svg" in response_json and response_json["svg"]:
                        self.add_test_result(
                            f"VIZ-1.{i+1}",
                            f"Visualize {molecule.get('name', 'Unknown Molecule')}",
                            "Passed",
                            f"Successfully generated visualization in {duration:.2f}s",
                            response,
                            duration
                        )
                    else:
                        self.add_test_result(
                            f"VIZ-1.{i+1}",
                            f"Visualize {molecule.get('name', 'Unknown Molecule')}",
                            "Failed",
                            "Response does not contain SVG data",
                            response,
                            duration
                        )
                else:
                    self.add_test_result(
                        f"VIZ-1.{i+1}",
                        f"Visualize {molecule.get('name', 'Unknown Molecule')}",
                        "Failed",
                        f"API returned status code {response.status_code}",
                        response,
                        duration
                    )
            except Exception as e:
                self.add_test_result(
                    f"VIZ-1.{i+1}",
                    f"Visualize {molecule.get('name', 'Unknown Molecule')}",
                    "Failed",
                    f"Error: {str(e)}",
                    None,
                    None
                )
        
        # Record performance metrics
        if response_times:
            avg_time = sum(response_times) / len(response_times)
            max_time = max(response_times)
            min_time = min(response_times)
            
            self.test_results["performance_metrics"]["visualization"] = {
                "average": avg_time,
                "min": min_time,
                "max": max_time,
                "count": len(response_times)
            }
            
            logger.info(f"Visualization performance: avg={avg_time:.2f}s, min={min_time:.2f}s, max={max_time:.2f}s")
    
    def test_mixture_analysis(self) -> None:
        """Test mixture analysis."""
        logger.info("Testing mixture analysis...")
        
        if not self.test_data.get("molecules") or len(self.test_data.get("molecules", [])) < 2:
            self.add_test_result("MIX-1", "Mixture Analysis", "Skipped", "Not enough test molecules available")
            return
        
        # Create a test mixture using the first two molecules
        molecules = self.test_data["molecules"][:2]
        mixture = {
            "name": f"Test Mixture {datetime.now().strftime('%Y%m%d%H%M%S')}",
            "description": "A test mixture created for production testing",
            "components": [
                {
                    "molecule_smiles": molecules[0]["smiles"],
                    "concentration": 70,
                    "concentration_unit": "%v/v"
                },
                {
                    "molecule_smiles": molecules[1]["smiles"],
                    "concentration": 30,
                    "concentration_unit": "%v/v"
                }
            ]
        }
        
        try:
            # Create the mixture
            response, duration = self.measure_execution_time(
                requests.post,
                f"{self.app_url}/api/v1/mixtures",
                json=mixture,
                timeout=30
            )
            
            if response.status_code in (200, 201):
                mixture_data = response.json()
                mixture_id = mixture_data.get("id")
                
                self.add_test_result(
                    "MIX-1.1",
                    "Create Mixture",
                    "Passed",
                    f"Successfully created mixture in {duration:.2f}s",
                    response,
                    duration
                )
                
                self.test_results["performance_metrics"]["mixture_creation"] = duration
                
                # Analyze the mixture properties
                response, duration = self.measure_execution_time(
                    requests.get,
                    f"{self.app_url}/api/v1/mixtures/{mixture_id}/properties",
                    timeout=30
                )
                
                if response.status_code == 200:
                    property_data = response.json()
                    
                    if "properties" in property_data and property_data["properties"]:
                        self.add_test_result(
                            "MIX-1.2",
                            "Analyze Mixture Properties",
                            "Passed",
                            f"Successfully analyzed mixture properties in {duration:.2f}s",
                            response,
                            duration
                        )
                    else:
                        self.add_test_result(
                            "MIX-1.2",
                            "Analyze Mixture Properties",
                            "Failed",
                            "Response does not contain properties data",
                            response,
                            duration
                        )
                else:
                    self.add_test_result(
                        "MIX-1.2",
                        "Analyze Mixture Properties",
                        "Failed",
                        f"API returned status code {response.status_code}",
                        response,
                        duration
                    )
                
                self.test_results["performance_metrics"]["mixture_analysis"] = duration
            else:
                self.add_test_result(
                    "MIX-1.1",
                    "Create Mixture",
                    "Failed",
                    f"API returned status code {response.status_code}",
                    response,
                    duration
                )
        except Exception as e:
            self.add_test_result(
                "MIX-1",
                "Mixture Analysis",
                "Failed",
                f"Error: {str(e)}",
                None,
                None
            )
    
    def test_search_functionality(self) -> None:
        """Test search functionality."""
        logger.info("Testing search functionality...")
        
        if not self.test_data.get("molecules"):
            self.add_test_result("SEARCH-1", "Search Functionality", "Skipped", "No test molecules available")
            return
        
        # Test keyword search using parts of molecule names
        response_times = []
        for i, molecule in enumerate(self.test_data["molecules"][:3]):
            try:
                if "name" not in molecule:
                    logger.warning(f"Molecule {i} has no name, skipping")
                    continue
                
                # Extract search term from molecule name (first word)
                search_term = molecule["name"].split()[0]
                
                # Call the search endpoint
                response, duration = self.measure_execution_time(
                    requests.get,
                    f"{self.app_url}/api/v1/molecules/search?name={search_term}",
                    timeout=30
                )
                
                response_times.append(duration)
                
                if response.status_code == 200:
                    search_results = response.json()
                    
                    if isinstance(search_results, list) and len(search_results) > 0:
                        self.add_test_result(
                            f"SEARCH-1.{i+1}",
                            f"Search for '{search_term}'",
                            "Passed",
                            f"Found {len(search_results)} results in {duration:.2f}s",
                            response,
                            duration
                        )
                    else:
                        self.add_test_result(
                            f"SEARCH-1.{i+1}",
                            f"Search for '{search_term}'",
                            "Failed",
                            "No search results found",
                            response,
                            duration
                        )
                else:
                    self.add_test_result(
                        f"SEARCH-1.{i+1}",
                        f"Search for '{search_term}'",
                        "Failed",
                        f"API returned status code {response.status_code}",
                        response,
                        duration
                    )
            except Exception as e:
                self.add_test_result(
                    f"SEARCH-1.{i+1}",
                    f"Search for molecule",
                    "Failed",
                    f"Error: {str(e)}",
                    None,
                    None
                )
        
        # Record performance metrics
        if response_times:
            avg_time = sum(response_times) / len(response_times)
            max_time = max(response_times)
            min_time = min(response_times)
            
            self.test_results["performance_metrics"]["keyword_search"] = {
                "average": avg_time,
                "min": min_time,
                "max": max_time,
                "count": len(response_times)
            }
            
            logger.info(f"Keyword search performance: avg={avg_time:.2f}s, min={min_time:.2f}s, max={max_time:.2f}s")
        
        # Test structure search using SMILES
        structure_response_times = []
        for i, molecule in enumerate(self.test_data["molecules"][:2]):
            try:
                if "smiles" not in molecule:
                    logger.warning(f"Molecule {i} has no SMILES data, skipping")
                    continue
                
                # Call the structure search endpoint
                response, duration = self.measure_execution_time(
                    requests.post,
                    f"{self.app_url}/api/v1/molecules/search/structure",
                    json={
                        "structure": molecule["smiles"],
                        "structure_format": "smiles"
                    },
                    timeout=30
                )
                
                structure_response_times.append(duration)
                
                if response.status_code == 200:
                    search_results = response.json()
                    
                    if isinstance(search_results, list) and len(search_results) > 0:
                        self.add_test_result(
                            f"SEARCH-2.{i+1}",
                            f"Structure Search",
                            "Passed",
                            f"Found {len(search_results)} results in {duration:.2f}s",
                            response,
                            duration
                        )
                    else:
                        self.add_test_result(
                            f"SEARCH-2.{i+1}",
                            f"Structure Search",
                            "Failed",
                            "No search results found",
                            response,
                            duration
                        )
                else:
                    self.add_test_result(
                        f"SEARCH-2.{i+1}",
                        f"Structure Search",
                        "Failed",
                        f"API returned status code {response.status_code}",
                        response,
                        duration
                    )
            except Exception as e:
                self.add_test_result(
                    f"SEARCH-2.{i+1}",
                    f"Structure Search",
                    "Failed",
                    f"Error: {str(e)}",
                    None,
                    None
                )
        
        # Record performance metrics
        if structure_response_times:
            avg_time = sum(structure_response_times) / len(structure_response_times)
            max_time = max(structure_response_times)
            min_time = min(structure_response_times)
            
            self.test_results["performance_metrics"]["structure_search"] = {
                "average": avg_time,
                "min": min_time,
                "max": max_time,
                "count": len(structure_response_times)
            }
            
            logger.info(f"Structure search performance: avg={avg_time:.2f}s, min={min_time:.2f}s, max={max_time:.2f}s")
    
    def test_batch_processing(self) -> None:
        """Test batch processing functionality."""
        logger.info("Testing batch processing...")
        
        if not self.test_data.get("molecules") or len(self.test_data.get("molecules", [])) < 5:
            self.add_test_result("BATCH-1", "Batch Processing", "Skipped", "Not enough test molecules available")
            return
        
        try:
            # Prepare batch request data
            molecules = self.test_data["molecules"][:5]
            batch_data = {
                "molecules": [
                    {"smiles": m["smiles"], "name": m.get("name", f"Molecule-{i}")}
                    for i, m in enumerate(molecules)
                ],
                "calculations": ["logp", "molecular_weight", "tpsa", "h_bond_donors", "h_bond_acceptors"]
            }
            
            # Call the batch processing endpoint
            response, duration = self.measure_execution_time(
                requests.post,
                f"{self.app_url}/api/v1/batch/properties",
                json=batch_data,
                timeout=60  # Higher timeout for batch processing
            )
            
            if response.status_code == 200:
                batch_results = response.json()
                
                if "results" in batch_results and len(batch_results["results"]) == len(batch_data["molecules"]):
                    self.add_test_result(
                        "BATCH-1",
                        "Batch Property Calculation",
                        "Passed",
                        f"Successfully processed {len(batch_data['molecules'])} molecules in {duration:.2f}s",
                        response,
                        duration
                    )
                    
                    # Verify individual results
                    all_properties_present = True
                    for result in batch_results["results"]:
                        if "properties" not in result or not all(calc in result["properties"] for calc in batch_data["calculations"]):
                            all_properties_present = False
                            break
                    
                    if not all_properties_present:
                        self.add_test_result(
                            "BATCH-1.1",
                            "Batch Result Validation",
                            "Failed",
                            "Some properties are missing in the batch results",
                            batch_results,
                            None
                        )
                    else:
                        self.add_test_result(
                            "BATCH-1.1",
                            "Batch Result Validation",
                            "Passed",
                            "All properties are present in the batch results",
                            None,
                            None
                        )
                else:
                    self.add_test_result(
                        "BATCH-1",
                        "Batch Property Calculation",
                        "Failed",
                        "Response does not contain expected results",
                        response,
                        duration
                    )
            else:
                self.add_test_result(
                    "BATCH-1",
                    "Batch Property Calculation",
                    "Failed",
                    f"API returned status code {response.status_code}",
                    response,
                    duration
                )
            
            self.test_results["performance_metrics"]["batch_processing"] = duration
            logger.info(f"Batch processing performance: {duration:.2f}s for {len(batch_data['molecules'])} molecules")
                
        except Exception as e:
            self.add_test_result(
                "BATCH-1",
                "Batch Processing",
                "Failed",
                f"Error: {str(e)}",
                None,
                None
            )
    
    def test_error_handling(self) -> None:
        """Test error handling."""
        logger.info("Testing error handling...")
        
        try:
            # Test with invalid SMILES
            invalid_smiles = "NOT_A_VALID_SMILES"
            response, duration = self.measure_execution_time(
                requests.get,
                f"{self.app_url}/molecule/{invalid_smiles}",
                timeout=30
            )
            
            # For error handling, we expect a 4xx status code with an error message
            if 400 <= response.status_code < 500:
                try:
                    error_data = response.json()
                    if "error" in error_data or "message" in error_data:
                        self.add_test_result(
                            "ERROR-1",
                            "Invalid SMILES Error Handling",
                            "Passed",
                            f"API correctly returned error status {response.status_code} with message",
                            response,
                            duration
                        )
                    else:
                        self.add_test_result(
                            "ERROR-1",
                            "Invalid SMILES Error Handling",
                            "Failed",
                            f"API returned status {response.status_code} but without error message",
                            response,
                            duration
                        )
                except:
                    self.add_test_result(
                        "ERROR-1",
                        "Invalid SMILES Error Handling",
                        "Failed",
                        f"API returned status {response.status_code} but response is not valid JSON",
                        response,
                        duration
                    )
            else:
                self.add_test_result(
                    "ERROR-1",
                    "Invalid SMILES Error Handling",
                    "Failed",
                    f"API returned unexpected status code {response.status_code}",
                    response,
                    duration
                )
            
            # Test with non-existent endpoint
            response, duration = self.measure_execution_time(
                requests.get,
                f"{self.app_url}/nonexistent_endpoint",
                timeout=30
            )
            
            if response.status_code == 404:
                self.add_test_result(
                    "ERROR-2",
                    "Non-existent Endpoint Error Handling",
                    "Passed",
                    "API correctly returned 404 for non-existent endpoint",
                    response,
                    duration
                )
            else:
                self.add_test_result(
                    "ERROR-2",
                    "Non-existent Endpoint Error Handling",
                    "Failed",
                    f"API returned unexpected status code {response.status_code}",
                    response,
                    duration
                )
                
        except Exception as e:
            self.add_test_result(
                "ERROR",
                "Error Handling",
                "Failed",
                f"Error: {str(e)}",
                None,
                None
            )
    
    def analyze_performance(self) -> None:
        """Analyze performance metrics and report findings."""
        logger.info("Analyzing performance metrics...")
        
        # Define performance thresholds (adjust based on requirements)
        thresholds = {
            "property_calculation.average": 2.0,  # seconds
            "visualization.average": 2.0,  # seconds
            "mixture_creation": 2.0,  # seconds
            "mixture_analysis": 2.0,  # seconds
            "keyword_search.average": 1.0,  # seconds
            "structure_search.average": 2.0,  # seconds
            "batch_processing": 5.0,  # seconds
        }
        
        # Check performance against thresholds
        performance_issues = []
        
        for metric_path, threshold in thresholds.items():
            parts = metric_path.split('.')
            metric_value = None
            
            # Navigate the metrics dictionary
            current = self.test_results["performance_metrics"]
            for part in parts:
                if part in current:
                    current = current[part]
                else:
                    break
            else:
                # We found the metric
                metric_value = current
            
            if metric_value is not None:
                if metric_value > threshold:
                    performance_issues.append(f"{metric_path} ({metric_value:.2f}s) exceeds threshold ({threshold:.2f}s)")
        
        # Add performance analysis result
        if performance_issues:
            self.add_test_result(
                "PERF",
                "Performance Analysis",
                "Failed",
                f"Performance issues found: {'; '.join(performance_issues)}",
                self.test_results["performance_metrics"],
                None
            )
        else:
            self.add_test_result(
                "PERF",
                "Performance Analysis",
                "Passed",
                "All performance metrics are within acceptable thresholds",
                self.test_results["performance_metrics"],
                None
            )
    
    def run_tests(self) -> Dict[str, Any]:
        """
        Run all tests and return results.
        
        Returns:
            Dictionary of test results
        """
        self.test_results["start_time"] = datetime.now().isoformat()
        self.test_results["status"] = "Running"
        
        logger.info("Starting production workflow tests...")
        
        # Load test data
        self.load_test_data()
        
        # Check service health first
        if not self.check_service_health():
            self.test_results["status"] = "Failed"
            self.test_results["end_time"] = datetime.now().isoformat()
            logger.error("Service health check failed. Aborting tests.")
            return self.test_results
        
        # Run functional tests
        self.test_molecule_property_calculation()
        self.test_molecule_visualization()
        self.test_mixture_analysis()
        self.test_search_functionality()
        self.test_batch_processing()
        self.test_error_handling()
        
        # Analyze performance
        self.analyze_performance()
        
        # Set final status
        if self.test_results["failed_tests"] > 0:
            self.test_results["status"] = "Failed"
        else:
            self.test_results["status"] = "Passed"
        
        self.test_results["end_time"] = datetime.now().isoformat()
        
        # Log summary
        logger.info(f"Test Status: {self.test_results['status']}")
        logger.info(f"Total Tests: {self.test_results['total_tests']}")
        logger.info(f"Passed Tests: {self.test_results['passed_tests']}")
        logger.info(f"Failed Tests: {self.test_results['failed_tests']}")
        logger.info(f"Skipped Tests: {self.test_results['skipped_tests']}")
        
        return self.test_results
    
    def save_results(self, filename: str) -> None:
        """
        Save test results to a file.
        
        Args:
            filename: Path to save the results to
        """
        try:
            with open(filename, 'w') as f:
                json.dump(self.test_results, f, indent=2)
            logger.info(f"Test results saved to {filename}")
        except Exception as e:
            logger.error(f"Failed to save test results: {str(e)}")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="CryoProtect Production Workflow Test")
    parser.add_argument("--app", default="http://localhost:5001", help="URL of the main application API")
    parser.add_argument("--rdkit", default="http://localhost:5002", help="URL of the RDKit service API")
    parser.add_argument("--output", default="production_workflow_results.json", help="Output file for test results")
    
    args = parser.parse_args()
    
    # Run tests
    test = ProductionWorkflowTest(args.app, args.rdkit)
    results = test.run_tests()
    
    # Save results
    test.save_results(args.output)
    
    # Exit with appropriate status code
    if results["status"] == "Passed":
        logger.info("All tests passed!")
        sys.exit(0)
    else:
        logger.error("Some tests failed. See log for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()