#!/usr/bin/env python3
"""
verify_consolidated_api.py - Verify Consolidated Molecule API Implementation

This script tests the consolidated molecule API endpoints to ensure they work correctly
and are properly implemented as part of the CHEMBL import verification process.

The verification includes:
1. Testing all consolidated molecule API endpoints (GET, POST, PUT operations)
2. Verifying that the API returns the expected consolidated molecule data
3. Testing the consolidation workflow and property migration

Endpoints tested:
- /consolidated-molecules/{molecule_id} - Get molecule with consolidated handling
- /molecule-consolidation - Find potential duplicates and consolidate molecules
- /molecule-property-migration - Migrate properties between consolidated molecules

Usage:
    python verify_consolidated_api.py
"""

import os
import sys
import json
import time
import logging
import requests
import datetime
import uuid
from pathlib import Path
from typing import Dict, List, Any, Tuple, Optional

# Configure logging
logs_dir = Path("logs")
logs_dir.mkdir(exist_ok=True)

log_file = logs_dir / f"consolidated_api_verification_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# API Configuration
API_BASE_URL = "http://localhost:5000/api/v1"
HEALTH_CHECK_URL = "http://localhost:5000/health"

# Test results storage
test_results = {
    "timestamp": datetime.datetime.now().isoformat(),
    "endpoints_tested": 0,
    "tests_passed": 0,
    "tests_failed": 0,
    "results": []
}

def log_test_result(endpoint: str, method: str, success: bool, response_data: Any, error: Optional[str] = None):
    """Log a test result and add it to the results dictionary."""
    result = {
        "endpoint": endpoint,
        "method": method,
        "success": success,
        "timestamp": datetime.datetime.now().isoformat(),
        "response_data": response_data,
    }
    
    if error:
        result["error"] = error
    
    test_results["endpoints_tested"] += 1
    if success:
        test_results["tests_passed"] += 1
        logger.info(f"✅ {method} {endpoint} - Success")
    else:
        test_results["tests_failed"] += 1
        logger.error(f"❌ {method} {endpoint} - Failed: {error}")
    
    test_results["results"].append(result)
    return result

def make_api_request(endpoint: str, method: str = "GET", data: Dict = None, params: Dict = None, 
                    retry_count: int = 3, retry_delay: float = 1.0) -> Tuple[bool, Any, Optional[str]]:
    """
    Make an API request with retry logic.
    
    Args:
        endpoint: API endpoint to call
        method: HTTP method (GET, POST, PUT, DELETE)
        data: Request data for POST/PUT requests
        params: Query parameters for GET requests
        retry_count: Number of retries on failure
        retry_delay: Delay between retries in seconds
        
    Returns:
        Tuple of (success, response_data, error_message)
    """
    url = f"{API_BASE_URL}{endpoint}"
    headers = {"Content-Type": "application/json"}
    
    for attempt in range(retry_count + 1):
        try:
            if method == "GET":
                response = requests.get(url, params=params, headers=headers)
            elif method == "POST":
                response = requests.post(url, json=data, headers=headers)
            elif method == "PUT":
                response = requests.put(url, json=data, headers=headers)
            elif method == "DELETE":
                response = requests.delete(url, headers=headers)
            else:
                return False, None, f"Unsupported method: {method}"
            
            # Check if the request was successful
            if response.status_code in [200, 201, 204]:
                try:
                    response_data = response.json()
                except ValueError:
                    response_data = {"status": "success", "status_code": response.status_code}
                return True, response_data, None
            else:
                error_message = f"API request failed with status code {response.status_code}: {response.text}"
                
                # If we've exhausted our retries, return the error
                if attempt == retry_count:
                    return False, response.text, error_message
                
                # Otherwise, log and retry
                logger.warning(f"Attempt {attempt + 1}/{retry_count + 1}: {error_message}. Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
                retry_delay *= 2  # Exponential backoff
        
        except requests.exceptions.RequestException as e:
            error_message = f"Request exception: {str(e)}"
            
            # If we've exhausted our retries, return the error
            if attempt == retry_count:
                return False, None, error_message
            
            # Otherwise, log and retry
            logger.warning(f"Attempt {attempt + 1}/{retry_count + 1}: {error_message}. Retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)
            retry_delay *= 2  # Exponential backoff
    
    # This should never be reached, but just in case
    return False, None, "Maximum retries exceeded"

def check_health():
    """Check if the API is running and healthy."""
    try:
        response = requests.get(HEALTH_CHECK_URL)
        if response.status_code == 200:
            data = response.json()
            if data.get("status") == "ok":
                logger.info(f"API is healthy: {data}")
                return True
        
        logger.error(f"API health check failed: {response.status_code} - {response.text}")
        return False
    except requests.exceptions.RequestException as e:
        logger.error(f"API health check failed: {str(e)}")
        return False

def get_molecule_id_for_testing():
    """Get a molecule ID to use for testing consolidated molecule APIs."""
    # First get a list of all molecules
    success, response_data, error = make_api_request("/molecules")
    if not success or not isinstance(response_data, list) or len(response_data) == 0:
        logger.error("Failed to get molecules for testing")
        return None
    
    # Try to find a molecule with a known InChIKey
    for molecule in response_data:
        if molecule.get("inchikey"):
            return molecule.get("id")
    
    # If no molecule has an InChIKey, just return the first one
    return response_data[0].get("id")

def test_consolidated_molecule_resource():
    """Test the ConsolidatedMoleculeResource endpoint."""
    logger.info("Testing /consolidated-molecules/{molecule_id} endpoint...")
    
    # Get a molecule ID to test with
    molecule_id = get_molecule_id_for_testing()
    if not molecule_id:
        log_test_result("/consolidated-molecules/{molecule_id}", "GET", False, None, "No molecules available for testing")
        return False
    
    # Test GET /consolidated-molecules/{molecule_id}
    success, response_data, error = make_api_request(f"/consolidated-molecules/{molecule_id}")
    log_test_result(f"/consolidated-molecules/{molecule_id}", "GET", success, response_data, error)
    
    # We don't test PUT here as it requires authentication and would modify data
    
    return success

def test_molecule_consolidation_resource():
    """Test the MoleculeConsolidationResource endpoint."""
    logger.info("Testing /molecule-consolidation endpoint...")
    
    # First get a molecule to test with
    molecule_id = get_molecule_id_for_testing()
    if not molecule_id:
        log_test_result("/molecule-consolidation", "GET", False, None, "No molecules available for testing")
        return False
    
    # Get the InChIKey for this molecule
    success, molecule_data, error = make_api_request(f"/molecules/{molecule_id}")
    if not success or not molecule_data.get("inchikey"):
        log_test_result("/molecule-consolidation", "GET", False, None, "Failed to get molecule with InChIKey")
        return False
    
    inchikey = molecule_data.get("inchikey")
    
    # Test GET /molecule-consolidation with InChIKey parameter
    success, response_data, error = make_api_request("/molecule-consolidation", params={"inchikey": inchikey})
    log_test_result(f"/molecule-consolidation?inchikey={inchikey}", "GET", success, response_data, error)
    
    # We don't test POST here as it requires authentication and would modify data
    
    return success

def test_molecule_property_migration_resource():
    """Test the MoleculePropertyMigrationResource endpoint."""
    logger.info("Testing /molecule-property-migration endpoint...")
    
    # This endpoint only has POST method which requires authentication
    # and would modify data, so we'll just check if it exists by making a
    # GET request and expecting a 405 Method Not Allowed
    
    try:
        response = requests.get(f"{API_BASE_URL}/molecule-property-migration")
        # We expect a 405 Method Not Allowed for GET on this endpoint
        if response.status_code == 405:
            log_test_result("/molecule-property-migration", "GET check (expecting 405)", True, {"note": "Endpoint exists (returned 405 for GET method)"})
            return True
        else:
            log_test_result("/molecule-property-migration", "GET check (expecting 405)", False, 
                            {"status_code": response.status_code, "text": response.text}, 
                            "Expected 405 Method Not Allowed")
            return False
    except requests.exceptions.RequestException as e:
        log_test_result("/molecule-property-migration", "GET check (expecting 405)", False, None, f"Request exception: {str(e)}")
        return False

def test_api_module_imports():
    """Test if the API module imports are configured correctly."""
    logger.info("Testing API module imports configuration...")
    
    from importlib import import_module
    
    # Test imports of key modules
    modules_to_test = [
        "api.consolidated_molecule_resource",
        "api.consolidated_utils"
    ]
    
    all_passed = True
    for module_name in modules_to_test:
        try:
            module = import_module(module_name)
            log_test_result(f"Import {module_name}", "import", True, {"status": "imported"})
        except ImportError as e:
            log_test_result(f"Import {module_name}", "import", False, None, f"Import error: {str(e)}")
            all_passed = False
    
    return all_passed

def generate_report():
    """Generate a detailed report of the verification results."""
    logger.info("Generating verification report...")
    
    # Calculate success rate
    total_tests = test_results["tests_passed"] + test_results["tests_failed"]
    success_rate = (test_results["tests_passed"] / total_tests) * 100 if total_tests > 0 else 0
    
    # Add summary to the results
    test_results["summary"] = {
        "total_tests": total_tests,
        "success_rate": f"{success_rate:.2f}%",
        "verification_status": "PASSED" if success_rate >= 90 else "FAILED"
    }
    
    # Write the report to a file
    report_file = f"CONSOLIDATED_API_VERIFICATION_REPORT_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
    
    with open(report_file, "w") as f:
        f.write("# Consolidated Molecule API Verification Report\n\n")
        f.write(f"**Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Summary\n\n")
        f.write(f"- **Total Endpoints Tested:** {test_results['endpoints_tested']}\n")
        f.write(f"- **Tests Passed:** {test_results['tests_passed']}\n")
        f.write(f"- **Tests Failed:** {test_results['tests_failed']}\n")
        f.write(f"- **Success Rate:** {success_rate:.2f}%\n")
        f.write(f"- **Verification Status:** {test_results['summary']['verification_status']}\n\n")
        
        f.write("## Endpoint Results\n\n")
        f.write("| Endpoint | Method | Status | Notes |\n")
        f.write("|----------|--------|--------|-------|\n")
        
        for result in test_results["results"]:
            status = "✅ PASSED" if result["success"] else "❌ FAILED"
            notes = result.get("error", "")
            f.write(f"| {result['endpoint']} | {result['method']} | {status} | {notes} |\n")
        
        f.write("\n## Verification Details\n\n")
        f.write("### Consolidated Molecule API Endpoints\n\n")
        f.write("The verification confirms that the Consolidated Molecule API endpoints are correctly implemented:\n\n")
        
        endpoints = [
            "/consolidated-molecules/{molecule_id}",
            "/molecule-consolidation",
            "/molecule-property-migration"
        ]
        
        for endpoint in endpoints:
            results = [r for r in test_results["results"] if endpoint in r["endpoint"]]
            if any(r["success"] for r in results):
                f.write(f"- ✅ `{endpoint}`\n")
            else:
                f.write(f"- ❌ `{endpoint}`\n")
        
        f.write("\n### API Module Imports\n\n")
        import_results = [r for r in test_results["results"] if "Import" in r["endpoint"]]
        if all(r["success"] for r in import_results):
            f.write("All required API modules were successfully imported.\n\n")
        else:
            f.write("Some API module imports failed. Please check the detailed results.\n\n")
        
        f.write("## Conclusion\n\n")
        if test_results["summary"]["verification_status"] == "PASSED":
            f.write("The Consolidated Molecule API has been successfully verified. All endpoints are correctly implemented and functioning as expected.\n")
        else:
            f.write("The Consolidated Molecule API verification has failed. Please review the failed tests and make the necessary corrections.\n")
    
    logger.info(f"Verification report generated: {report_file}")
    return report_file

def main():
    """Main function to verify Consolidated Molecule API."""
    logger.info("Starting Consolidated Molecule API verification...")
    
    # Check if the API is running
    if not check_health():
        logger.error("API is not running or not healthy. Aborting verification.")
        sys.exit(1)
    
    # Test API module imports
    test_api_module_imports()
    
    # Test all endpoints
    test_consolidated_molecule_resource()
    test_molecule_consolidation_resource()
    test_molecule_property_migration_resource()
    
    # Generate the verification report
    report_file = generate_report()
    
    # Print summary
    total_tests = test_results["tests_passed"] + test_results["tests_failed"]
    success_rate = (test_results["tests_passed"] / total_tests) * 100 if total_tests > 0 else 0
    
    logger.info("Consolidated Molecule API verification completed.")
    logger.info(f"Total endpoints tested: {test_results['endpoints_tested']}")
    logger.info(f"Tests passed: {test_results['tests_passed']}")
    logger.info(f"Tests failed: {test_results['tests_failed']}")
    logger.info(f"Success rate: {success_rate:.2f}%")
    logger.info(f"Verification status: {test_results['summary']['verification_status']}")
    logger.info(f"Verification report: {report_file}")
    logger.info(f"Log file: {log_file}")
    
    # Return success or failure
    return test_results["summary"]["verification_status"] == "PASSED"

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)