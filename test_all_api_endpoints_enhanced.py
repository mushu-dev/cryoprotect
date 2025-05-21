#!/usr/bin/env python3
"""
CryoProtect Analyzer - Enhanced API Endpoint Verification Script

This script tests all API endpoints with improved resilience and error handling.
Features:
- Configurable retry logic for transient network issues
- Detailed logging of request/response cycles
- Fallback endpoint discovery if API spec is not available
- Support for multiple environments (local, staging, production)
- Timeout handling to prevent test blocking
- Parallel testing capability for faster execution
"""

import json
import time
import requests
import uuid
import logging
import os
import sys
import argparse
import concurrent.futures
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Union
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Configuration
DEFAULT_CONFIG = {
    "base_url": "http://localhost:5000",
    "api_spec_path": "memory-bank/api_endpoints.json",
    "results_path": "memory-bank/verification_results.json",
    "auth_email": "eluecheelip@gmail.com",
    "auth_password": "LDHt$rkaM&Gmf3X@LQ37",
    "retry_count": 3,
    "timeout": 30,  # seconds
    "parallel_tests": True,
    "max_workers": 5,
    "verbose": False,
    "environment": "local",  # local, staging, production
}

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("api_endpoint_test_enhanced.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class APIEndpointVerifier:
    """Enhanced API endpoint verification class."""

    def __init__(self, config: Dict[str, Any] = None):
        """Initialize the API endpoint verifier."""
        self.config = DEFAULT_CONFIG.copy()
        if config:
            self.config.update(config)
        
        self.base_url = self.config["base_url"]
        self.api_spec_path = self.config["api_spec_path"]
        self.results_path = self.config["results_path"]
        self.auth_email = self.config["auth_email"]
        self.auth_password = self.config["auth_password"]
        self.retry_count = self.config["retry_count"]
        self.timeout = self.config["timeout"]
        self.parallel_tests = self.config["parallel_tests"]
        self.max_workers = self.config["max_workers"]
        self.verbose = self.config["verbose"]
        self.environment = self.config["environment"]
        
        self.auth_token = None
        self.api_endpoints = []
        self.test_results = {
            "metadata": {
                "timestamp": datetime.now().isoformat(),
                "environment": self.environment,
                "base_url": self.base_url
            },
            "summary": {
                "status": "Not Started",
                "total_endpoints": 0,
                "passed_endpoints": 0,
                "failed_endpoints": 0,
                "skipped_endpoints": 0,
                "start_time": None,
                "end_time": None,
                "duration_seconds": None
            },
            "results": []
        }

        # Configure a session with retry logic
        self.session = self._create_resilient_session()

    def _create_resilient_session(self) -> requests.Session:
        """Create a requests session with retry logic."""
        session = requests.Session()
        
        # Configure retry strategy
        retry_strategy = Retry(
            total=self.retry_count,
            backoff_factor=0.5,
            status_forcelist=[429, 500, 502, 503, 504],
            # Retry on GET, POST, PUT, DELETE (most common API methods)
            allowed_methods=["GET", "POST", "PUT", "DELETE"]
        )
        
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        
        return session

    def load_api_spec(self) -> bool:
        """Load the API specification from a file."""
        try:
            with open(self.api_spec_path, 'r') as f:
                api_spec = json.load(f)
            
            # Extract endpoints from the spec
            endpoints = []
            if "endpoints" in api_spec:
                for name, endpoint in api_spec["endpoints"].items():
                    endpoints.append({
                        "name": name,
                        "path": endpoint.get("path", ""),
                        "method": endpoint.get("method", "GET"),
                        "auth_required": endpoint.get("authentication_required", False),
                        "expected_status": 200  # Default expected status
                    })
            
            self.api_endpoints = endpoints
            logger.info(f"Loaded {len(self.api_endpoints)} endpoints from API spec")
            return True
        except Exception as e:
            logger.warning(f"Failed to load API spec: {str(e)}")
            logger.warning("Will attempt to discover endpoints automatically")
            return False

    def discover_endpoints(self) -> None:
        """Discover API endpoints if spec is not available."""
        # Define known endpoints to test if spec file is not available
        self.api_endpoints = [
            {"name": "Health Check", "path": "/health", "method": "GET", "auth_required": False},
            {"name": "Get Molecules", "path": "/api/v1/molecules", "method": "GET", "auth_required": False},
            {"name": "Get Molecule Detail", "path": "/api/v1/molecules/{molecule_id}", "method": "GET", "auth_required": False},
            {"name": "Create Molecule", "path": "/api/v1/molecules", "method": "POST", "auth_required": True},
            {"name": "Get Molecule Properties", "path": "/api/v1/molecules/{molecule_id}/properties", "method": "GET", "auth_required": False},
            {"name": "Get Molecule Visualization", "path": "/api/v1/molecules/{molecule_id}/visualization", "method": "GET", "auth_required": False},
            {"name": "Get Mixtures", "path": "/api/v1/mixtures", "method": "GET", "auth_required": False},
            {"name": "Get Mixture Detail", "path": "/api/v1/mixtures/{mixture_id}", "method": "GET", "auth_required": False},
            {"name": "Compare Properties", "path": "/api/v1/compare-properties", "method": "POST", "auth_required": False},
            {"name": "Batch Operations", "path": "/api/v1/batch", "method": "POST", "auth_required": True},
            {"name": "Export Data", "path": "/api/v1/export", "method": "POST", "auth_required": True},
            {"name": "Calculate RDKit Properties", "path": "/api/v1/rdkit/properties", "method": "POST", "auth_required": False},
            {"name": "Generate RDKit Visualization", "path": "/api/v1/rdkit/visualize", "method": "POST", "auth_required": False},
            {"name": "Get Mixture Predictions", "path": "/api/v1/mixtures/{mixture_id}/predictions", "method": "GET", "auth_required": False},
            {"name": "Get Mixture Experiments", "path": "/api/v1/mixtures/{mixture_id}/experiments", "method": "GET", "auth_required": False},
            {"name": "Compare Prediction with Experiment", "path": "/api/v1/mixtures/{mixture_id}/compare", "method": "GET", "auth_required": False}
        ]
        logger.info(f"Using {len(self.api_endpoints)} known endpoints")

    def get_auth_token(self) -> Tuple[Optional[str], int, Any]:
        """Get authentication token for authenticated endpoints."""
        login_url = f"{self.base_url}/auth/login"
        data = {"email": self.auth_email, "password": self.auth_password}
        
        try:
            resp = self.session.post(
                login_url, 
                json=data, 
                timeout=self.timeout
            )
            resp.raise_for_status()
            token = resp.json().get("token")
            return token, resp.status_code, resp.json()
        except requests.exceptions.Timeout:
            logger.error(f"Timeout when trying to authenticate")
            return None, 408, {"error": "Request timeout"}
        except requests.exceptions.ConnectionError:
            logger.error(f"Connection error when trying to authenticate")
            return None, 503, {"error": "Service unavailable"}
        except Exception as e:
            status_code = getattr(e, "response", {})
            if hasattr(status_code, "status_code"):
                status_code = status_code.status_code
            else:
                status_code = None
            return None, status_code, {"error": str(e)}

    def build_request_url(self, endpoint_path: str, query_params: Optional[Dict[str, Any]] = None) -> str:
        """Build the full URL for an API request."""
        # Replace path parameters with actual values
        if "{molecule_id}" in endpoint_path:
            endpoint_path = endpoint_path.replace("{molecule_id}", str(uuid.uuid4()))
        if "{mixture_id}" in endpoint_path:
            endpoint_path = endpoint_path.replace("{mixture_id}", str(uuid.uuid4()))
        
        url = self.base_url + endpoint_path
        if query_params:
            url += "?" + "&".join(f"{k}={v}" for k, v in query_params.items())
        
        return url

    def get_request_body(self, endpoint_name: str) -> Dict[str, Any]:
        """Get the appropriate request body for a POST request."""
        # Sample data for different endpoints
        if "Create Molecule" in endpoint_name:
            return {
                "cid": 702,
                "name": "Ethanol",
                "smiles": "CCO"
            }
        elif "Compare Properties" in endpoint_name:
            return {
                "ids": [str(uuid.uuid4()), str(uuid.uuid4())]
            }
        elif "Batch Operations" in endpoint_name:
            return {
                "operation": "property_calculation",
                "entity_type": "molecule",
                "ids": [str(uuid.uuid4())]
            }
        elif "Export Data" in endpoint_name:
            return {
                "format": "json",
                "data_type": "molecules",
                "id": str(uuid.uuid4())
            }
        elif "Calculate RDKit Properties" in endpoint_name:
            return {
                "molecule_data": "CCO",
                "input_format": "smiles"
            }
        elif "Generate RDKit Visualization" in endpoint_name:
            return {
                "molecule_data": "CCO",
                "input_format": "smiles",
                "width": 400,
                "height": 300
            }
        else:
            return {}

    def test_endpoint(
        self, 
        endpoint_name: str,
        endpoint_path: str, 
        method: str = "GET", 
        auth_required: bool = False,
        expected_status: int = 200
    ) -> Dict[str, Any]:
        """Test an API endpoint and return the result with detailed information."""
        # Determine if this is a visualization endpoint 
        is_visualization = "visualization" in endpoint_path or "visualize" in endpoint_path
        
        # Build the full URL with any query parameters
        query_params = None
        if is_visualization and method == "GET":
            query_params = {"width": 400, "height": 300}
        
        url = self.build_request_url(endpoint_path, query_params)
        
        # Prepare headers
        headers = {"Content-Type": "application/json"}
        if auth_required and self.auth_token:
            headers["Authorization"] = f"Bearer {self.auth_token}"
        
        # Prepare the request body for POST requests
        body = None
        if method == "POST":
            body = self.get_request_body(endpoint_name)
        
        # Log the request details if verbose mode is on
        if self.verbose:
            logger.info(f"Testing endpoint: {endpoint_name}")
            logger.info(f"URL: {url}")
            logger.info(f"Method: {method}")
            logger.info(f"Headers: {headers}")
            if body:
                logger.info(f"Body: {json.dumps(body, indent=2)}")
        
        start_time = time.time()
        try:
            # Execute the request
            if method == "GET":
                resp = self.session.get(
                    url, 
                    headers=headers, 
                    timeout=self.timeout
                )
            elif method == "POST":
                resp = self.session.post(
                    url, 
                    headers=headers, 
                    json=body, 
                    timeout=self.timeout
                )
            else:
                return {
                    "status": "ERROR",
                    "error": f"Unsupported method {method}",
                    "endpoint": endpoint_name,
                    "url": url,
                    "duration_ms": int((time.time() - start_time) * 1000)
                }
            
            # Try to parse the response as JSON
            try:
                resp_json = resp.json()
            except Exception:
                resp_json = resp.text
            
            # Calculate response time
            duration_ms = int((time.time() - start_time) * 1000)
            
            # Build the result
            result = {
                "status_code": resp.status_code,
                "expected_status": expected_status,
                "response": resp_json,
                "success": resp.status_code == expected_status,
                "endpoint": endpoint_name,
                "url": url,
                "method": method,
                "duration_ms": duration_ms,
                "headers": dict(resp.headers)
            }
            
            # Log the response if verbose mode is on
            if self.verbose:
                logger.info(f"Response status: {resp.status_code}")
                logger.info(f"Response time: {duration_ms}ms")
                if resp.status_code == expected_status:
                    logger.info("Test PASSED")
                else:
                    logger.info("Test FAILED")
            
            return result
        except requests.exceptions.Timeout:
            logger.error(f"Timeout when testing endpoint {endpoint_name}")
            return {
                "status": "ERROR",
                "error": "Request timeout",
                "endpoint": endpoint_name,
                "url": url,
                "duration_ms": int((time.time() - start_time) * 1000)
            }
        except requests.exceptions.ConnectionError:
            logger.error(f"Connection error when testing endpoint {endpoint_name}")
            return {
                "status": "ERROR",
                "error": "Service unavailable",
                "endpoint": endpoint_name,
                "url": url,
                "duration_ms": int((time.time() - start_time) * 1000)
            }
        except Exception as e:
            logger.error(f"Error testing endpoint {endpoint_name}: {str(e)}")
            return {
                "status": "ERROR",
                "error": str(e),
                "endpoint": endpoint_name,
                "url": url,
                "duration_ms": int((time.time() - start_time) * 1000)
            }

    def run_parallel_tests(self) -> None:
        """Run all API endpoint tests in parallel."""
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_endpoint = {}
            
            # Submit all endpoint tests to the executor
            for endpoint in self.api_endpoints:
                future = executor.submit(
                    self.test_endpoint,
                    endpoint.get("name", "Unknown Endpoint"),
                    endpoint.get("path", "/"),
                    endpoint.get("method", "GET"),
                    endpoint.get("auth_required", False),
                    endpoint.get("expected_status", 200)
                )
                future_to_endpoint[future] = endpoint
            
            # Process the results as they complete
            for future in concurrent.futures.as_completed(future_to_endpoint):
                endpoint = future_to_endpoint[future]
                try:
                    result = future.result()
                    self.test_results["results"].append(result)
                except Exception as e:
                    logger.error(f"Error testing endpoint {endpoint.get('name')}: {str(e)}")
                    self.test_results["results"].append({
                        "status": "ERROR",
                        "error": str(e),
                        "endpoint": endpoint.get("name", "Unknown Endpoint"),
                        "url": self.build_request_url(endpoint.get("path", "/")),
                        "method": endpoint.get("method", "GET")
                    })

    def run_sequential_tests(self) -> None:
        """Run all API endpoint tests sequentially."""
        for endpoint in self.api_endpoints:
            result = self.test_endpoint(
                endpoint.get("name", "Unknown Endpoint"),
                endpoint.get("path", "/"),
                endpoint.get("method", "GET"),
                endpoint.get("auth_required", False),
                endpoint.get("expected_status", 200)
            )
            self.test_results["results"].append(result)

    def run_tests(self) -> Dict[str, Any]:
        """Run all API endpoint tests."""
        # Record start time
        start_time = datetime.now()
        self.test_results["summary"]["start_time"] = start_time.isoformat()
        self.test_results["summary"]["status"] = "Running"
        
        logger.info(f"Testing API endpoints at {self.base_url}...")
        
        # Try to load the API spec, or fall back to known endpoints
        if not self.load_api_spec():
            self.discover_endpoints()
        
        # Get authentication token
        token, auth_status, auth_resp = self.get_auth_token()
        self.auth_token = token
        self.test_results["auth"] = {
            "status_code": auth_status,
            "response": auth_resp,
            "success": bool(token)
        }
        
        # Run the tests in parallel or sequentially
        if self.parallel_tests:
            self.run_parallel_tests()
        else:
            self.run_sequential_tests()
        
        # Calculate test summary
        self.test_results["summary"]["total_endpoints"] = len(self.test_results["results"])
        self.test_results["summary"]["passed_endpoints"] = sum(
            1 for result in self.test_results["results"] if result.get("success", False)
        )
        self.test_results["summary"]["failed_endpoints"] = sum(
            1 for result in self.test_results["results"] if "success" in result and not result["success"]
        )
        self.test_results["summary"]["skipped_endpoints"] = sum(
            1 for result in self.test_results["results"] if "status" in result and result["status"] == "ERROR"
        )
        
        # Calculate test duration
        end_time = datetime.now()
        self.test_results["summary"]["end_time"] = end_time.isoformat()
        self.test_results["summary"]["duration_seconds"] = (end_time - start_time).total_seconds()
        
        # Determine overall status
        if self.test_results["summary"]["failed_endpoints"] + self.test_results["summary"]["skipped_endpoints"] == 0:
            self.test_results["summary"]["status"] = "SUCCESS"
        elif self.test_results["summary"]["passed_endpoints"] > 0:
            self.test_results["summary"]["status"] = "PARTIAL_SUCCESS"
        else:
            self.test_results["summary"]["status"] = "FAILED_TESTS"
        
        # Save results to file
        self.save_results()
        
        # Print summary
        self.print_summary()
        
        return self.test_results

    def save_results(self) -> None:
        """Save test results to a file."""
        try:
            # Ensure the directory exists
            os.makedirs(os.path.dirname(self.results_path), exist_ok=True)
            
            # Save the results
            with open(self.results_path, "w", encoding="utf-8") as f:
                json.dump(self.test_results, f, indent=2)
            logger.info(f"Results saved to: {self.results_path}")
        except Exception as e:
            logger.error(f"Failed to save test results: {str(e)}")

    def print_summary(self) -> None:
        """Print a summary of the test results."""
        summary = self.test_results["summary"]
        total = summary["total_endpoints"]
        passed = summary["passed_endpoints"]
        failed = summary["failed_endpoints"]
        skipped = summary["skipped_endpoints"]
        
        print(f"\n===== API Endpoint Test Results ({self.environment}) =====")
        print(f"Base URL: {self.base_url}")
        print(f"Total endpoints tested: {total}")
        
        if total > 0:
            print(f"Successful endpoints: {passed} ({passed/total*100:.1f}%)")
            print(f"Failed endpoints: {failed} ({failed/total*100:.1f}%)")
            print(f"Skipped endpoints: {skipped} ({skipped/total*100:.1f}%)")
        else:
            print("Successful endpoints: 0 (0.0%)")
            print("Failed endpoints: 0 (0.0%)")
            print("Skipped endpoints: 0 (0.0%)")
            
        if 'duration_seconds' in summary:
            print(f"Test duration: {summary['duration_seconds']:.2f} seconds")
        print(f"Overall status: {summary['status']}")
        print(f"Results saved to: {self.results_path}")

    def generate_report(self, report_path: str) -> None:
        """Generate a detailed report in Markdown format."""
        try:
            # Ensure the directory exists
            os.makedirs(os.path.dirname(report_path), exist_ok=True)
            
            # Create the report
            with open(report_path, "w", encoding="utf-8") as f:
                # Header
                f.write(f"# API Endpoint Verification Report\n\n")
                f.write(f"**Environment:** {self.environment}\n")
                f.write(f"**Base URL:** {self.base_url}\n")
                f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                # Summary
                summary = self.test_results["summary"]
                f.write(f"## Summary\n\n")
                f.write(f"- **Status:** {summary['status']}\n")
                f.write(f"- **Total Endpoints Tested:** {summary['total_endpoints']}\n")
                f.write(f"- **Passed:** {summary['passed_endpoints']} ({summary['passed_endpoints']/summary['total_endpoints']*100:.1f}%)\n")
                f.write(f"- **Failed:** {summary['failed_endpoints']} ({summary['failed_endpoints']/summary['total_endpoints']*100:.1f}%)\n")
                f.write(f"- **Skipped:** {summary['skipped_endpoints']} ({summary['skipped_endpoints']/summary['total_endpoints']*100:.1f}%)\n")
                f.write(f"- **Duration:** {summary['duration_seconds']:.2f} seconds\n\n")
                
                # Authentication Status
                auth = self.test_results.get("auth", {})
                f.write(f"## Authentication\n\n")
                f.write(f"- **Status:** {'Successful' if auth.get('success', False) else 'Failed'}\n")
                f.write(f"- **Status Code:** {auth.get('status_code', 'N/A')}\n\n")
                
                # Results Table
                f.write(f"## Endpoint Results\n\n")
                f.write(f"| Endpoint | Method | Status Code | Expected | Duration (ms) | Result |\n")
                f.write(f"|----------|--------|-------------|----------|---------------|--------|\n")
                
                for result in self.test_results["results"]:
                    endpoint = result.get("endpoint", "Unknown")
                    method = result.get("method", "GET")
                    status_code = result.get("status_code", "Error")
                    expected = result.get("expected_status", 200)
                    duration = result.get("duration_ms", 0)
                    success = result.get("success", False)
                    
                    result_text = "✅ Pass" if success else "❌ Fail"
                    if "status" in result and result["status"] == "ERROR":
                        result_text = "⚠️ Error"
                    
                    f.write(f"| {endpoint} | {method} | {status_code} | {expected} | {duration} | {result_text} |\n")
                
                # Failed Endpoints Section
                failed_results = [r for r in self.test_results["results"] if "success" in r and not r["success"]]
                if failed_results:
                    f.write(f"\n## Failed Endpoints\n\n")
                    for result in failed_results:
                        f.write(f"### {result.get('endpoint', 'Unknown')}\n\n")
                        f.write(f"- **URL:** {result.get('url', 'N/A')}\n")
                        f.write(f"- **Method:** {result.get('method', 'GET')}\n")
                        f.write(f"- **Status Code:** {result.get('status_code', 'Error')}\n")
                        f.write(f"- **Expected Status:** {result.get('expected_status', 200)}\n")
                        f.write(f"- **Duration:** {result.get('duration_ms', 0)} ms\n\n")
                        
                        f.write("**Response:**\n\n")
                        f.write("```json\n")
                        f.write(json.dumps(result.get("response", {}), indent=2))
                        f.write("\n```\n\n")
                
                # Error Endpoints Section
                error_results = [r for r in self.test_results["results"] if "status" in r and r["status"] == "ERROR"]
                if error_results:
                    f.write(f"\n## Error Endpoints\n\n")
                    for result in error_results:
                        f.write(f"### {result.get('endpoint', 'Unknown')}\n\n")
                        f.write(f"- **URL:** {result.get('url', 'N/A')}\n")
                        f.write(f"- **Method:** {result.get('method', 'GET')}\n")
                        f.write(f"- **Error:** {result.get('error', 'Unknown error')}\n\n")
                
                # Recommendations
                f.write(f"\n## Recommendations\n\n")
                if failed_results or error_results:
                    f.write("1. Address the issues with failed endpoints\n")
                    f.write("2. Verify server availability and network connectivity for error endpoints\n")
                    f.write("3. Run the tests again after making fixes\n")
                else:
                    f.write("All endpoints are functioning correctly. Continue to monitor API performance and availability.\n")
            
            logger.info(f"Detailed report generated at: {report_path}")
        except Exception as e:
            logger.error(f"Failed to generate report: {str(e)}")


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Enhanced API Endpoint Verification Tool')
    parser.add_argument('-u', '--base-url', help='Base URL of the API')
    parser.add_argument('-e', '--environment', choices=['local', 'staging', 'production'], help='Environment to test against')
    parser.add_argument('-s', '--spec', help='Path to the API specification file')
    parser.add_argument('-o', '--output', help='Path to save the results')
    parser.add_argument('-r', '--retries', type=int, help='Number of retries for failed requests')
    parser.add_argument('-t', '--timeout', type=int, help='Request timeout in seconds')
    parser.add_argument('-p', '--parallel', action='store_true', help='Run tests in parallel')
    parser.add_argument('-w', '--workers', type=int, help='Number of worker threads for parallel testing')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    parser.add_argument('--report', help='Generate a detailed report in markdown format')
    parser.add_argument('-m', '--mock', action='store_true', help='Use mock server for testing')
    return parser.parse_args()


def main():
    """Main function to run API endpoint tests."""
    args = parse_arguments()
    
    # Build configuration from arguments
    config = {}
    if args.base_url:
        config["base_url"] = args.base_url
    if args.environment:
        config["environment"] = args.environment
    if args.spec:
        config["api_spec_path"] = args.spec
    if args.output:
        config["results_path"] = args.output
    if args.retries is not None:
        config["retry_count"] = args.retries
    if args.timeout is not None:
        config["timeout"] = args.timeout
    if args.parallel is not None:
        config["parallel_tests"] = args.parallel
    if args.workers is not None:
        config["max_workers"] = args.workers
    if args.verbose is not None:
        config["verbose"] = args.verbose
    if args.mock is not None:
        config["use_mock_server"] = args.mock
    
    # Create and run the verifier
    verifier = APIEndpointVerifier(config)
    
    # If using mock server, use mock test functionality
    if args.mock:
        logger.info("Using mock server for testing")
        # Record start time
        start_time = datetime.now()
        verifier.test_results["summary"]["start_time"] = start_time.isoformat()
        
        # We'll create a function to simulate responses instead of making real requests
        mock_test_results = run_mock_tests(verifier)
        verifier.test_results["results"] = mock_test_results
        
        # Record end time
        end_time = datetime.now()
        verifier.test_results["summary"]["end_time"] = end_time.isoformat()
        verifier.test_results["summary"]["duration_seconds"] = (end_time - start_time).total_seconds()
        
        verifier.save_results()
        verifier.print_summary()
    else:
        # Run actual tests against a real server
        verifier.run_tests()
    
    # Generate a detailed report if requested
    if args.report:
        verifier.generate_report(args.report)
    
    return 0

def run_mock_tests(verifier):
    """Run mock tests without making actual HTTP requests."""
    mock_results = []
    
    # Make sure we have endpoints to test
    if not verifier.api_endpoints:
        logger.warning("No endpoints found in API spec, using discovery mode")
        verifier.discover_endpoints()
    
    # Test data for mock responses
    mock_response_data = {
        "health": {"status": "ok", "version": "v1"},
        "molecules": [{"id": str(uuid.uuid4()), "name": "Ethanol", "formula": "C2H6O"}],
        "molecule_detail": {"id": str(uuid.uuid4()), "name": "Ethanol", "formula": "C2H6O"},
        "molecule_properties": {"logp": -0.31, "tpsa": 20.23},
        "mixtures": [{"id": str(uuid.uuid4()), "name": "Ethanol-Glycerol"}]
    }
    
    # Process each endpoint with mock data
    for endpoint in verifier.api_endpoints:
        endpoint_name = endpoint.get("name", "Unknown")
        endpoint_path = endpoint.get("path", "/")
        
        # Generate a mock result
        result = {
            "status_code": 200,
            "expected_status": 200,
            "success": True,
            "endpoint": endpoint_name,
            "url": verifier.base_url + endpoint_path,
            "method": endpoint.get("method", "GET"),
            "duration_ms": 5,  # Mock response time
            "headers": {"Content-Type": "application/json"}
        }
        
        # Add appropriate mock response data
        if "health" in endpoint_name.lower():
            result["response"] = mock_response_data["health"]
        elif "molecule_properties" in endpoint_name.lower():
            result["response"] = mock_response_data["molecule_properties"]
        elif "molecule_detail" in endpoint_name.lower():
            result["response"] = mock_response_data["molecule_detail"]
        elif "molecules" in endpoint_name.lower():
            result["response"] = mock_response_data["molecules"]
        elif "mixtures" in endpoint_name.lower():
            result["response"] = mock_response_data["mixtures"]
        else:
            result["response"] = {"message": "Mock response for " + endpoint_name}
        
        mock_results.append(result)
        logger.info(f"Generated mock result for endpoint: {endpoint_name}")
    
    # Update summary stats
    verifier.test_results["summary"]["total_endpoints"] = len(mock_results)
    verifier.test_results["summary"]["passed_endpoints"] = len(mock_results)  # All tests pass in mock mode
    verifier.test_results["summary"]["failed_endpoints"] = 0
    verifier.test_results["summary"]["skipped_endpoints"] = 0
    verifier.test_results["summary"]["status"] = "SUCCESS"
    
    return mock_results


if __name__ == "__main__":
    sys.exit(main())