#!/usr/bin/env python
"""
CryoProtect Analyzer - API Endpoint Verification Script

This script verifies all required API endpoints for CryoProtect v2.
It checks if endpoints exist, return expected responses, and handle errors appropriately.
"""

import os
import sys
import json
import uuid
import requests
from datetime import datetime
from tabulate import tabulate

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import Flask app for local testing
from app import create_app

class APIEndpointVerifier:
    """Verifies API endpoints for CryoProtect v2."""
    
    def __init__(self, base_url="http://localhost:5000"):
        """Initialize the verifier with the base URL."""
        self.base_url = base_url
        self.results = []
        self.app = create_app(testing=True)
        self.client = self.app.test_client()
        
        # Sample data for testing
        self.sample_molecule_id = str(uuid.uuid4())
        self.sample_mixture_id = str(uuid.uuid4())
        
    def verify_endpoint(self, endpoint, method="GET", data=None, expected_status=200, description=""):
        """Verify an API endpoint."""
        url = f"{self.base_url}{endpoint}"
        
        try:
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
                    "notes": notes
                })
                
                return status, response_data
        except Exception as e:
            self.results.append({
                "endpoint": endpoint,
                "method": method,
                "description": description,
                "implemented": "Error",
                "functional": "Error",
                "error_handling": "Error",
                "status_code": "Error",
                "notes": str(e)
            })
            return None, str(e)
    
    def verify_all_endpoints(self):
        """Verify all required API endpoints."""
        print("Verifying API endpoints...")
        
        # Health check endpoint
        self.verify_endpoint(
            "/health",
            method="GET",
            description="Health check endpoint"
        )
        
        # Molecule endpoints
        self.verify_endpoint(
            "/api/v1/molecules",
            method="GET",
            description="Get a list of all molecules"
        )
        
        self.verify_endpoint(
            f"/api/v1/molecules/{self.sample_molecule_id}",
            method="GET",
            description="Get detailed information about a specific molecule"
        )
        
        self.verify_endpoint(
            "/api/v1/molecules",
            method="POST",
            data={
                "cid": 702,
                "name": "Ethanol",
                "smiles": "CCO"
            },
            description="Create a new molecule"
        )
        
        self.verify_endpoint(
            f"/api/v1/molecules/{self.sample_molecule_id}/properties",
            method="GET",
            description="Get properties of a specific molecule"
        )
        
        self.verify_endpoint(
            f"/api/v1/molecules/{self.sample_molecule_id}/visualization",
            method="GET",
            description="Get SVG visualization of a specific molecule"
        )
        
        # Mixture endpoints
        self.verify_endpoint(
            "/api/v1/mixtures",
            method="GET",
            description="Get a list of all mixtures"
        )
        
        self.verify_endpoint(
            f"/api/v1/mixtures/{self.sample_mixture_id}",
            method="GET",
            description="Get detailed information about a specific mixture"
        )
        
        # Compare properties endpoint
        self.verify_endpoint(
            "/api/v1/compare-properties",
            method="POST",
            data={
                "ids": [self.sample_molecule_id, self.sample_mixture_id]
            },
            description="Compare molecular and mixture properties"
        )
        
        # Batch operations endpoint
        self.verify_endpoint(
            "/api/v1/batch",
            method="POST",
            data={
                "operation": "property_calculation",
                "entity_type": "molecule",
                "ids": [self.sample_molecule_id]
            },
            description="Perform batch operations"
        )
        
        # Export endpoint
        self.verify_endpoint(
            "/api/v1/export",
            method="POST",
            data={
                "format": "json",
                "data_type": "molecules",
                "id": self.sample_molecule_id
            },
            description="Export protocols and results"
        )
        
        # RDKit endpoints
        self.verify_endpoint(
            "/api/v1/rdkit/properties",
            method="POST",
            data={
                "molecule_data": "CCO",
                "input_format": "smiles"
            },
            description="Calculate molecular properties using RDKit"
        )
        
        self.verify_endpoint(
            "/api/v1/rdkit/visualize",
            method="POST",
            data={
                "molecule_data": "CCO",
                "input_format": "smiles",
                "width": 400,
                "height": 300
            },
            description="Generate molecule visualization using RDKit"
        )
        
        # Prediction and experiment endpoints
        self.verify_endpoint(
            f"/api/v1/mixtures/{self.sample_mixture_id}/predictions",
            method="GET",
            description="Get predictions for a mixture"
        )
        
        self.verify_endpoint(
            f"/api/v1/mixtures/{self.sample_mixture_id}/experiments",
            method="GET",
            description="Get experiments for a mixture"
        )
        
        self.verify_endpoint(
            f"/api/v1/mixtures/{self.sample_mixture_id}/compare",
            method="GET",
            description="Compare prediction with experiment for a mixture"
        )
    
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
        
        return {
            "status": "SUCCESS" if not discrepancies else "COMPLETED_WITH_WARNINGS",
            "summary": {
                "total": total_count,
                "implemented": implemented_count,
                "functional": functional_count,
                "discrepancies": len(discrepancies)
            },
            "results": self.results
        }

def main():
    """Main function to run the verification."""
    verifier = APIEndpointVerifier()
    verifier.verify_all_endpoints()
    report = verifier.generate_report()
    
    # Save report to file
    report_dir = os.path.join(os.path.dirname(__file__), 'reports')
    os.makedirs(report_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_path = os.path.join(report_dir, f'api_verification_report_{timestamp}.json')
    
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\nReport saved to: {report_path}")
    
    # Return exit code based on verification results
    return 0 if report["status"] == "SUCCESS" else 1

if __name__ == '__main__':
    sys.exit(main())