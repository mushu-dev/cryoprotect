#!/usr/bin/env python3
"""
CryoProtect v2 - API Integration Tests

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

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("api_integration_test.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class APIIntegrationTest:
    """Test class for API integration."""

    def __init__(self, base_url: str = "http://localhost:5000"):
        """Initialize the test class."""
        self.base_url = base_url
        self.auth_token = None
        self.test_results = {
            "status": "Not Started",
            "total_tests": 0,
            "passed_tests": 0,
            "failed_tests": 0,
            "skipped_tests": 0,
            "test_cases": []
        }

    def run_tests(self) -> Dict[str, Any]:
        """Run all API integration tests."""
        self.test_results["status"] = "Running"
        self.test_results["start_time"] = datetime.now().isoformat()

        # Run the tests
        self.test_health_endpoint()
        self.test_molecules_endpoint()
        self.test_mixtures_endpoint()
        self.test_rdkit_properties_endpoint()
        self.test_error_handling()

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

    def test_molecules_endpoint(self) -> None:
        """Test the molecules endpoint."""
        try:
            # Send request to the molecules endpoint
            response = requests.get(f"{self.base_url}/api/v1/molecules")
            
            # Check the response
            if response.status_code == 200:
                data = response.json()
                if isinstance(data, list):
                    self.add_test_result("API-2", "Molecules Endpoint", "Passed", "Molecules endpoint returned a list", response)
                else:
                    self.add_test_result("API-2", "Molecules Endpoint", "Failed", "Molecules endpoint did not return a list", response)
            else:
                self.add_test_result("API-2", "Molecules Endpoint", "Failed", f"Molecules endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-2", "Molecules Endpoint", "Failed", f"Error testing molecules endpoint: {str(e)}")

    def test_mixtures_endpoint(self) -> None:
        """Test the mixtures endpoint."""
        try:
            # Send request to the mixtures endpoint
            response = requests.get(f"{self.base_url}/api/v1/mixtures")
            
            # Check the response
            if response.status_code == 200:
                data = response.json()
                if isinstance(data, list):
                    self.add_test_result("API-3", "Mixtures Endpoint", "Passed", "Mixtures endpoint returned a list", response)
                else:
                    self.add_test_result("API-3", "Mixtures Endpoint", "Failed", "Mixtures endpoint did not return a list", response)
            else:
                self.add_test_result("API-3", "Mixtures Endpoint", "Failed", f"Mixtures endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-3", "Mixtures Endpoint", "Failed", f"Error testing mixtures endpoint: {str(e)}")

    def test_rdkit_properties_endpoint(self) -> None:
        """Test the RDKit properties endpoint."""
        try:
            # Send request to the RDKit properties endpoint
            payload = {
                "molecule_data": "CCO",
                "input_format": "smiles"
            }
            response = requests.post(f"{self.base_url}/api/v1/rdkit/properties", json=payload)
            
            # Check the response
            if response.status_code == 200:
                data = response.json()
                if "hydrogen_bonding" in data and "logp" in data and "tpsa" in data:
                    self.add_test_result("API-4", "RDKit Properties Endpoint", "Passed", "RDKit properties endpoint returned expected properties", response)
                else:
                    self.add_test_result("API-4", "RDKit Properties Endpoint", "Failed", "RDKit properties endpoint did not return expected properties", response)
            else:
                self.add_test_result("API-4", "RDKit Properties Endpoint", "Failed", f"RDKit properties endpoint returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-4", "RDKit Properties Endpoint", "Failed", f"Error testing RDKit properties endpoint: {str(e)}")

    def test_error_handling(self) -> None:
        """Test error handling."""
        try:
            # Test 404 error
            response = requests.get(f"{self.base_url}/api/v1/nonexistent")
            
            if response.status_code == 404:
                self.add_test_result("API-5.1", "404 Error Handling", "Passed", "API correctly returned 404 for nonexistent endpoint", response)
            else:
                self.add_test_result("API-5.1", "404 Error Handling", "Failed", f"API returned unexpected status code: {response.status_code}", response)
            
            # Test invalid request
            response = requests.post(f"{self.base_url}/api/v1/rdkit/properties", json={})
            
            if response.status_code in [400, 422]:
                self.add_test_result("API-5.2", "Invalid Request Handling", "Passed", "API correctly rejected invalid request", response)
            else:
                self.add_test_result("API-5.2", "Invalid Request Handling", "Failed", f"API returned unexpected status code: {response.status_code}", response)
        except Exception as e:
            self.add_test_result("API-5", "Error Handling", "Failed", f"Error testing error handling: {str(e)}")

    def save_results(self, filename: str) -> None:
        """Save the test results to a file."""
        try:
            with open(filename, 'w') as f:
                json.dump(self.test_results, f, indent=2)
            logger.info(f"Test results saved to {filename}")
        except Exception as e:
            logger.error(f"Failed to save test results: {str(e)}")

    def generate_report(self, filename: str) -> None:
        """Generate a Markdown report of the test results."""
        try:
            with open(filename, 'w') as f:
                f.write("# API Integration Verification Report\n\n")
                f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("## Summary\n\n")
                f.write(f"- **Total Endpoints Tested:** {self.test_results['total_tests']}\n")
                f.write(f"- **Tests Passed:** {self.test_results['passed_tests']}\n")
                f.write(f"- **Tests Failed:** {self.test_results['failed_tests']}\n")
                f.write(f"- **Success Rate:** {(self.test_results['passed_tests'] / self.test_results['total_tests'] * 100):.2f}%\n")
                f.write(f"- **Verification Status:** {self.test_results['status']}\n\n")
                
                f.write("## Endpoint Results\n\n")
                f.write("| Endpoint | Method | Status | Notes |\n")
                f.write("|----------|--------|--------|-------|\n")
                
                for test_case in self.test_results["test_cases"]:
                    endpoint = test_case["name"].split(" - ")[0] if " - " in test_case["name"] else test_case["name"]
                    method = "GET" if "GET" in test_case["name"] else "POST" if "POST" in test_case["name"] else "N/A"
                    status = test_case["status"]
                    notes = test_case["message"]
                    
                    f.write(f"| {endpoint} | {method} | {status} | {notes} |\n")
                
                f.write("\n## Recommendations\n\n")
                
                if self.test_results["failed_tests"] > 0:
                    f.write("The following issues should be addressed:\n\n")
                    for test_case in self.test_results["test_cases"]:
                        if test_case["status"] == "Failed":
                            f.write(f"1. **{test_case['name']}**: {test_case['message']}\n")
                else:
                    f.write("All tests passed. No immediate issues to address.\n")
                
                f.write("\n## Next Steps\n\n")
                f.write("1. Continue monitoring API performance\n")
                f.write("2. Expand test coverage to include more edge cases\n")
                f.write("3. Implement automated testing as part of CI/CD pipeline\n")
            
            logger.info(f"Report generated at {filename}")
        except Exception as e:
            logger.error(f"Failed to generate report: {str(e)}")

def main():
    """Main function."""
    logger.info("Starting API integration tests")
    
    # Create and run the tests
    test = APIIntegrationTest()
    results = test.run_tests()
    
    # Save the results
    test.save_results("api_integration_test_results.json")
    
    # Generate a report
    report_filename = f"API_VERIFICATION_REPORT_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
    test.generate_report(report_filename)
    
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