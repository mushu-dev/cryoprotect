#!/usr/bin/env python3
"""
CryoProtect Analyzer - API Endpoint Verification Script

This script tests all 16 API endpoints to verify they are functioning correctly
now that the database has been populated with test data.
"""

import json
import requests
import uuid
from datetime import datetime

# Configuration
BASE_URL = "http://localhost:5000"
API_SPEC_PATH = "memory-bank/api_endpoints.json"
RESULTS_PATH = "memory-bank/verification_results.json"
AUTH_EMAIL = "eluecheelip@gmail.com"
AUTH_PASSWORD = "LDHt$rkaM&Gmf3X@LQ37"

def get_auth_token():
    """Get authentication token for authenticated endpoints."""
    login_url = f"{BASE_URL}/auth/login"
    data = {"email": AUTH_EMAIL, "password": AUTH_PASSWORD}
    try:
        resp = requests.post(login_url, json=data)
        resp.raise_for_status()
        token = resp.json().get("token")
        return token, resp.status_code, resp.json()
    except Exception as e:
        return None, getattr(e, "response", {}).get("status_code", None), str(e)

def test_endpoint(endpoint, method="GET", auth_token=None, query_params=None, body=None, expected_status=200, headers=None):
    """Test an API endpoint and return the result."""
    url = BASE_URL + endpoint
    if query_params:
        url += "?" + "&".join(f"{k}={v}" for k, v in query_params.items())
    
    req_headers = headers.copy() if headers else {}
    if auth_token:
        req_headers["Authorization"] = f"Bearer {auth_token}"
    
    try:
        if method == "GET":
            resp = requests.get(url, headers=req_headers)
        elif method == "POST":
            resp = requests.post(url, headers=req_headers, json=body)
        else:
            return {"status": "ERROR", "error": f"Unsupported method {method}"}
        
        try:
            resp_json = resp.json()
        except Exception:
            resp_json = resp.text
        
        return {
            "status_code": resp.status_code,
            "expected_status": expected_status,
            "response": resp_json,
            "success": resp.status_code == expected_status
        }
    except Exception as e:
        return {"status": "ERROR", "error": str(e)}

def main():
    """Main function to test all API endpoints."""
    print(f"Testing API endpoints at {BASE_URL}...")
    
    # Get authentication token
    token, auth_status, auth_resp = get_auth_token()
    results = {
        "auth": {
            "status_code": auth_status,
            "response": auth_resp,
            "success": bool(token)
        }
    }
    
    # Generate sample IDs for testing
    sample_molecule_id = str(uuid.uuid4())
    sample_mixture_id = str(uuid.uuid4())
    
    # Test all 16 endpoints
    
    # 1. Health check endpoint
    results["health"] = test_endpoint(
        "/health",
        method="GET"
    )
    
    # 2. Get a list of all molecules
    results["molecules"] = test_endpoint(
        "/api/v1/molecules",
        method="GET"
    )
    
    # 3. Get detailed information about a specific molecule
    # First, get a real molecule ID from the list if available
    molecule_id = sample_molecule_id
    if results["molecules"]["success"]:
        try:
            molecules = results["molecules"]["response"].get("molecules", [])
            if molecules:
                molecule_id = molecules[0]["id"]
        except Exception:
            pass
    
    results["molecule_detail"] = test_endpoint(
        f"/api/v1/molecules/{molecule_id}",
        method="GET"
    )
    
    # 4. Create a new molecule
    results["molecule_create"] = test_endpoint(
        "/api/v1/molecules",
        method="POST",
        auth_token=token,
        body={
            "cid": 702,
            "name": "Ethanol",
            "smiles": "CCO"
        }
    )
    
    # 5. Get properties of a specific molecule
    results["molecule_properties"] = test_endpoint(
        f"/api/v1/molecules/{molecule_id}/properties",
        method="GET"
    )
    
    # 6. Get SVG visualization of a specific molecule
    results["molecule_visualization"] = test_endpoint(
        f"/api/v1/molecules/{molecule_id}/visualization",
        method="GET",
        query_params={"width": 400, "height": 300}
    )
    
    # 7. Get a list of all mixtures
    results["mixtures"] = test_endpoint(
        "/api/v1/mixtures",
        method="GET"
    )
    
    # 8. Get detailed information about a specific mixture
    # First, get a real mixture ID from the list if available
    mixture_id = sample_mixture_id
    if results["mixtures"]["success"]:
        try:
            mixtures = results["mixtures"]["response"].get("mixtures", [])
            if mixtures:
                mixture_id = mixtures[0]["id"]
        except Exception:
            pass
    
    results["mixture_detail"] = test_endpoint(
        f"/api/v1/mixtures/{mixture_id}",
        method="GET"
    )
    
    # 9. Compare molecular and mixture properties
    results["compare_properties"] = test_endpoint(
        "/api/v1/compare-properties",
        method="POST",
        body={
            "ids": [molecule_id, mixture_id]
        }
    )
    
    # 10. Perform batch operations
    results["batch_operations"] = test_endpoint(
        "/api/v1/batch",
        method="POST",
        auth_token=token,
        body={
            "operation": "property_calculation",
            "entity_type": "molecule",
            "ids": [molecule_id]
        }
    )
    
    # 11. Export protocols and results
    results["export"] = test_endpoint(
        "/api/v1/export",
        method="POST",
        auth_token=token,
        body={
            "format": "json",
            "data_type": "molecules",
            "id": molecule_id
        }
    )
    
    # 12. Calculate molecular properties using RDKit
    results["rdkit_properties"] = test_endpoint(
        "/api/v1/rdkit/properties",
        method="POST",
        body={
            "molecule_data": "CCO",
            "input_format": "smiles"
        }
    )
    
    # 13. Generate molecule visualization using RDKit
    results["rdkit_visualize"] = test_endpoint(
        "/api/v1/rdkit/visualize",
        method="POST",
        body={
            "molecule_data": "CCO",
            "input_format": "smiles",
            "width": 400,
            "height": 300
        }
    )
    
    # 14. Get predictions for a mixture
    results["mixture_predictions"] = test_endpoint(
        f"/api/v1/mixtures/{mixture_id}/predictions",
        method="GET"
    )
    
    # 15. Get experiments for a mixture
    results["mixture_experiments"] = test_endpoint(
        f"/api/v1/mixtures/{mixture_id}/experiments",
        method="GET"
    )
    
    # 16. Compare prediction with experiment for a mixture
    results["mixture_compare"] = test_endpoint(
        f"/api/v1/mixtures/{mixture_id}/compare",
        method="GET"
    )
    
    # Save results to file
    with open(RESULTS_PATH, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    
    # Print summary
    success_count = sum(1 for result in results.values() if result.get("success", False))
    total_count = len(results)
    
    print(f"\n===== API Endpoint Test Results =====")
    print(f"Total endpoints tested: {total_count}")
    print(f"Successful endpoints: {success_count} ({success_count/total_count*100:.1f}%)")
    print(f"Failed endpoints: {total_count - success_count} ({(total_count - success_count)/total_count*100:.1f}%)")
    print(f"Results saved to: {RESULTS_PATH}")
    
    return {
        "status": "SUCCESS" if success_count == total_count else "PARTIAL_SUCCESS" if success_count > 0 else "FAILED_TESTS",
        "summary": {
            "total": total_count,
            "successful": success_count,
            "failed": total_count - success_count
        },
        "results_file": RESULTS_PATH
    }

if __name__ == "__main__":
    main()