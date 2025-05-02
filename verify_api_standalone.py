#!/usr/bin/env python3
"""
CryoProtect v2 - Standalone API Endpoint Verification

This script independently verifies the API endpoints without running the full fix process.
It can be used to diagnose issues with specific endpoints.

Usage:
    python verify_api_standalone.py [--endpoint <endpoint_name>]
"""

import os
import sys
import json
import uuid
import argparse
import requests
import logging
from datetime import datetime
from tabulate import tabulate

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("api_verification_standalone.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# Import Flask app for local testing
try:
    from app import create_app
    logger.info("Successfully imported Flask app")
except Exception as e:
    logger.error(f"Error importing Flask app: {e}")
    print(f"Error importing Flask app: {e}")
    sys.exit(1)

class APIEndpointVerifier:
    """Verifies API endpoints for CryoProtect v2."""
    
    def __init__(self, base_url="http://localhost:5000"):
        """Initialize the verifier with the base URL."""
        self.base_url = base_url
        self.results = []
        
        try:
            self.app = create_app(testing=True)
            self.client = self.app.test_client()
            logger.info("Test client created successfully")
        except Exception as e:
            logger.error(f"Error creating test client: {e}")
            raise
        
        # Sample data for testing
        self.sample_molecule_id = str(uuid.uuid4())
        self.sample_mixture_id = str(uuid.uuid4())
    
    def verify_endpoint(self, endpoint, method="GET", data=None, expected_status=200, description=""):
        """Verify an API endpoint."""
        url = f"{self.base_url}{endpoint}"
        
        try:
            logger.info(f"Testing endpoint: {method} {endpoint}")
            with self.app.app_context():
                if method == "GET":
                    response = self.client.get(url)
                elif method == "POST":
                    response = self.client.post(url, json=data)
                elif method == "PUT":
                    response = self.client.put(url, json=data)
                elif method == "DELETE":
                    response = self.client.delete(url)
                else:
                    raise ValueError(f"Unsupported method: {method}")
                
                status = response.status_code
                
                # Try to parse response as JSON
                try:
                    response_data = response.get_json()
                except:
                    response_data = response.data.decode('utf-8')
                
                # Check if the endpoint exists
                if status == 404 and "404 Not Found" in str(response_data):
                    implemented = "Not Implemented"
                    functional = "N/A"
                    error_handling = "N/A"
                    notes = "Endpoint not found"
                elif status == expected_status:
                    implemented = "Implemented"
                    functional = "Functional"
                    error_handling = "Proper"
                    notes = "Returns expected response"
                else:
                    implemented = "Implemented"
                    functional = "Not Functional"
                    error_handling = "Improper" if status >= 500 else "Proper"
                    notes = f"Unexpected status code: {status}"
                
                self.results.append({
                    "endpoint": endpoint,
                    "method": method,
                    "description": description,
                    "implemented": implemented,
                    "functional": functional,
                    "error_handling": error_handling,
                    "status_code": status,
                    "notes": notes,
                    "response": str(response_data)[:200] + "..." if len(str(response_data)) > 200 else str(response_data)
                })
                
                logger.info(f"Endpoint {method} {endpoint}: {implemented}, {functional}, Status: {status}")
                return status, response_data
        except Exception as e:
            logger.error(f"Error testing endpoint {method} {endpoint}: {str(e)}")
            self.results.append({
                "endpoint": endpoint,
                "method": method,
                "description": description,
                "implemented": "Error",
                "functional": "Error",
                "error_handling": "Error",
                "status_code": "Error",
                "notes": str(e),
                "response": str(e)
            })
            return None, str(e)
    
    def get_endpoints(self):
        """Get a list of all endpoints to test."""
        return [
            {
                "endpoint": "/health",
                "method": "GET",
                "description": "Health check endpoint"
            },
            {
                "endpoint": "/api/v1/molecules",
                "method": "GET",
                "description": "Get a list of all molecules"
            },
            {
                "endpoint": f"/api/v1/molecules/{self.sample_molecule_id}",
                "method": "GET",
                "description": "Get detailed information about a specific molecule"
            },
            {
                "endpoint": "/api/v1/molecules",
                "method": "POST",
                "data": {
                    "cid": 702,
                    "name": "Ethanol",
                    "smiles": "CCO"
                },
                "description": "Create a new molecule"
            },
            {
                "endpoint": f"/api/v1/molecules/{self.sample_molecule_id}/properties",
                "method": "GET",
                "description": "Get properties of a specific molecule"
            },
            {
                "endpoint": f"/api/v1/molecules/{self.sample_molecule_id}/visualization",
                "method": "GET",
                "description": "Get SVG visualization of a specific molecule"
            },
            {
                "endpoint": "/api/v1/mixtures",
                "method": "GET",
                "description": "Get a list of all mixtures"
            },
            {
                "endpoint": f"/api/v1/mixtures/{self.sample_mixture_id}",
                "method": "GET",
                "description": "Get detailed information about a specific mixture"
            },
            {
                "endpoint": "/api/v1/compare-properties",
                "method": "POST",
                "data": {
                    "ids": [self.sample_molecule_id, self.sample_mixture_id]
                },
                "description": "Compare molecular and mixture properties"
            },
            {
                "endpoint": "/api/v1/batch",
                "method": "POST",
                "data": {
                    "operation": "property_calculation",
                    "entity_type": "molecule",
                    "ids": [self.sample_molecule_id]
                },
                "description": "Perform batch operations"
            },
            {
                "endpoint": "/api/v1/export",
                "method": "POST",
                "data": {
                    "format": "json",
                    "data_type": "molecules",
                    "id": self.sample_molecule_id
                },
                "description": "Export protocols and results"
            },
            {
                "endpoint": "/api/v1/rdkit/properties",
                "method": "POST",
                "data": {
                    "molecule_data": "CCO",
                    "input_format": "smiles"
                },
                "description": "Calculate molecular properties using RDKit"
            },
            {
                "endpoint": "/api/v1/rdkit/visualize",
                "method": "POST",
                "data": {
                    "molecule_data": "CCO",
                    "input_format": "smiles",
                    "width": 400,
                    "height": 300
                },
                "description": "Generate molecule visualization using RDKit"
            },
            {
                "endpoint": f"/api/v1/mixtures/{self.sample_mixture_id}/predictions",
                "method": "GET",
                "description": "Get predictions for a mixture"
            },
            {
                "endpoint": f"/api/v1/mixtures/{self.sample_mixture_id}/experiments",
                "method": "GET",
                "description": "Get experiments for a mixture"
            },
            {
                "endpoint": f"/api/v1/mixtures/{self.sample_mixture_id}/compare",
                "method": "GET",
                "description": "Compare prediction with experiment for a mixture"
            }
        ]
    
    def verify_all_endpoints(self):
        """Verify all required API endpoints."""
        print("Verifying API endpoints...")
        
        endpoints = self.get_endpoints()
        for endpoint_info in endpoints:
            self.verify_endpoint(
                endpoint_info["endpoint"],
                method=endpoint_info.get("method", "GET"),
                data=endpoint_info.get("data"),
                description=endpoint_info.get("description", "")
            )
    
    def verify_specific_endpoint(self, endpoint_name):
        """Verify a specific endpoint by name."""
        print(f"Verifying endpoint: {endpoint_name}")
        
        endpoints = self.get_endpoints()
        found = False
        
        for endpoint_info in endpoints:
            if endpoint_name in endpoint_info["endpoint"]:
                found = True
                self.verify_endpoint(
                    endpoint_info["endpoint"],
                    method=endpoint_info.get("method", "GET"),
                    data=endpoint_info.get("data"),
                    description=endpoint_info.get("description", "")
                )
        
        if not found:
            print(f"Endpoint '{endpoint_name}' not found in the list of endpoints to test.")
            return False
        
        return True
    
    def generate_report(self):
        """Generate a report of the verification results."""
        print("\n===== API Endpoint Verification Report =====")
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Base URL: {self.base_url}")
        print("\n")
        
        # Prepare table data
        table_data = []
        for result in self.results:
            table_data.append([
                f"{result['method']} {result['endpoint']}",
                result['description'],
                result['implemented'],
                result['functional'],
                result['error_handling'],
                result['status_code'],
                result['notes']
            ])
        
        # Print table
        headers = ["Endpoint", "Description", "Implemented", "Functional", "Error Handling", "Status Code", "Notes"]
        print(tabulate(table_data, headers=headers, tablefmt="grid"))
        
        # Summary
        implemented_count = sum(1 for r in self.results if r['implemented'] == "Implemented")
        functional_count = sum(1 for r in self.results if r['functional'] == "Functional")
        total_count = len(self.results)
        
        print("\n===== Summary =====")
        print(f"Total Endpoints: {total_count}")
        print(f"Implemented: {implemented_count} ({implemented_count/total_count*100:.1f}%)")
        print(f"Functional: {functional_count} ({functional_count/total_count*100:.1f}%)")
        
        # Discrepancies between documentation and implementation
        print("\n===== Discrepancies =====")
        discrepancies = [r for r in self.results if r['implemented'] != "Implemented" or r['functional'] != "Functional"]
        if discrepancies:
            for d in discrepancies:
                print(f"- {d['method']} {d['endpoint']}: {d['notes']}")
        else:
            print("No discrepancies found.")
        
        # Recommendations
        print("\n===== Recommendations =====")
        if discrepancies:
            print("1. Implement missing endpoints")
            print("2. Fix non-functional endpoints")
            print("3. Ensure proper error handling for all endpoints")
        else:
            print("All endpoints are implemented and functional.")
        
        # Save detailed results
        detailed_results = []
        for result in self.results:
            detailed_results.append({
                "endpoint": f"{result['method']} {result['endpoint']}",
                "description": result['description'],
                "implemented": result['implemented'],
                "functional": result['functional'],
                "error_handling": result['error_handling'],
                "status_code": result['status_code'],
                "notes": result['notes'],
                "response": result.get('response', '')
            })
        
        report = {
            "status": "SUCCESS" if not discrepancies else "COMPLETED_WITH_WARNINGS",
            "summary": {
                "total": total_count,
                "implemented": implemented_count,
                "functional": functional_count,
                "discrepancies": len(discrepancies)
            },
            "results": detailed_results
        }
        
        # Save report to file
        report_dir = os.path.join(os.path.dirname(__file__), 'reports')
        os.makedirs(report_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        report_path = os.path.join(report_dir, f'api_verification_standalone_{timestamp}.json')
        
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"\nDetailed report saved to: {report_path}")
        
        return report

def main():
    """Main function to run the verification."""
    parser = argparse.ArgumentParser(description='Verify API endpoints')
    parser.add_argument('--endpoint', help='Specific endpoint to verify (partial match)')
    args = parser.parse_args()
    
    try:
        verifier = APIEndpointVerifier()
        
        if args.endpoint:
            verifier.verify_specific_endpoint(args.endpoint)
        else:
            verifier.verify_all_endpoints()
        
        report = verifier.generate_report()
        
        # Return exit code based on verification results
        return 0 if report["status"] == "SUCCESS" else 1
    
    except Exception as e:
        logger.error(f"Error in verification: {e}")
        print(f"Error in verification: {e}")
        return 1

if __name__ == '__main__':
    sys.exit(main())