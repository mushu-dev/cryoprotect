#!/usr/bin/env python3
"""
ChEMBL Data API Verification Test

This script tests that the website/API endpoints return real ChEMBL data
as specified in .specs/chembl_data_verification.md. It verifies that:
1. API endpoints return non-empty responses with real ChEMBL data
2. Returned data includes canonical fields and matches DB content
3. Both authenticated and anonymous access are tested (if allowed)
4. Results are logged and reported

Related Task: task-imp-wv-1-1-api-check
"""

import os
import sys
import json
import pytest
import requests
import datetime
from typing import Dict, List, Any, Optional
import logging

# Add the parent directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/chembl_api_verification.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Configuration
API_BASE_URL = os.getenv('API_BASE_URL', 'http://localhost:5000/api')
SUPABASE_URL = os.getenv('SUPABASE_URL')
SUPABASE_KEY = os.getenv('SUPABASE_KEY')
REPORT_DIR = os.getenv('REPORT_DIR', 'reports')

# Ensure report directory exists
os.makedirs(REPORT_DIR, exist_ok=True)

# Import Supabase client if available
try:
    from supabase import create_client, Client
    supabase_available = True
except ImportError:
    logger.warning("Supabase client not available. Database verification will be skipped.")
    supabase_available = False

class ChEMBLAPIVerifier:
    """Verifies that API endpoints return real ChEMBL data."""
    
    def __init__(self, api_base_url: str, supabase_url: Optional[str] = None, supabase_key: Optional[str] = None):
        """
        Initialize the verifier.
        
        Args:
            api_base_url: Base URL for the API
            supabase_url: Supabase URL for database verification
            supabase_key: Supabase key for database verification
        """
        self.api_base_url = api_base_url
        self.supabase_url = supabase_url
        self.supabase_key = supabase_key
        self.supabase_client = None
        self.results = {
            "timestamp": datetime.datetime.now().isoformat(),
            "api_base_url": api_base_url,
            "endpoints_tested": [],
            "endpoints_passed": [],
            "endpoints_failed": [],
            "chembl_data_found": False,
            "details": {}
        }
        
        # Initialize Supabase client if credentials are provided
        if supabase_available and supabase_url and supabase_key:
            try:
                self.supabase_client = create_client(supabase_url, supabase_key)
                logger.info("Supabase client initialized successfully.")
            except Exception as e:
                logger.error(f"Failed to initialize Supabase client: {e}")
    
    def verify_endpoint(self, endpoint: str, auth_token: Optional[str] = None, 
                       expected_fields: List[str] = None, params: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Verify that an API endpoint returns valid ChEMBL data.
        
        Args:
            endpoint: API endpoint to test
            auth_token: Optional authentication token
            expected_fields: List of fields expected in the response
            params: Query parameters
            
        Returns:
            Dictionary with verification results
        """
        full_url = f"{self.api_base_url}/{endpoint}"
        headers = {}
        if auth_token:
            headers['Authorization'] = f'Bearer {auth_token}'
        
        logger.info(f"Testing endpoint: {full_url}")
        
        try:
            response = requests.get(full_url, headers=headers, params=params)
            
            # Check if response is successful
            if response.status_code != 200:
                result = {
                    "endpoint": endpoint,
                    "status": "failed",
                    "status_code": response.status_code,
                    "error": f"Unexpected status code: {response.status_code}",
                    "response": response.text[:200] + "..." if len(response.text) > 200 else response.text
                }
                logger.error(f"Endpoint {endpoint} failed with status code {response.status_code}")
                self.results["endpoints_failed"].append(endpoint)
                return result
            
            # Parse response
            try:
                data = response.json()
            except json.JSONDecodeError:
                result = {
                    "endpoint": endpoint,
                    "status": "failed",
                    "status_code": response.status_code,
                    "error": "Response is not valid JSON",
                    "response": response.text[:200] + "..." if len(response.text) > 200 else response.text
                }
                logger.error(f"Endpoint {endpoint} returned invalid JSON")
                self.results["endpoints_failed"].append(endpoint)
                return result
            
            # Check if response is empty
            if isinstance(data, list) and len(data) == 0:
                result = {
                    "endpoint": endpoint,
                    "status": "failed",
                    "status_code": response.status_code,
                    "error": "Response is empty",
                    "response": data
                }
                logger.error(f"Endpoint {endpoint} returned empty list")
                self.results["endpoints_failed"].append(endpoint)
                return result
            
            # Extract the first item if response is a list
            first_item = data[0] if isinstance(data, list) else data
            
            # Check for expected fields
            if expected_fields:
                missing_fields = [field for field in expected_fields if field not in first_item]
                if missing_fields:
                    result = {
                        "endpoint": endpoint,
                        "status": "failed",
                        "status_code": response.status_code,
                        "error": f"Missing expected fields: {missing_fields}",
                        "response": first_item
                    }
                    logger.error(f"Endpoint {endpoint} missing fields: {missing_fields}")
                    self.results["endpoints_failed"].append(endpoint)
                    return result
            
            # Check for ChEMBL data
            chembl_data_found = False
            chembl_id = None
            
            # Look for ChEMBL ID in various possible fields
            for field in ['chembl_id', 'data_source', 'source_id', 'external_id']:
                if field in first_item and first_item[field] and 'CHEMBL' in str(first_item[field]):
                    chembl_data_found = True
                    chembl_id = first_item[field]
                    break
            
            # If no direct ChEMBL ID found, check if there's a reference to ChEMBL in any field
            if not chembl_data_found:
                for key, value in first_item.items():
                    if isinstance(value, str) and 'CHEMBL' in value:
                        chembl_data_found = True
                        chembl_id = value
                        break
            
            # Verify with database if Supabase client is available
            db_verification = None
            if self.supabase_client and chembl_id:
                db_verification = self.verify_with_database(chembl_id, first_item)
            
            # Prepare result
            result = {
                "endpoint": endpoint,
                "status": "passed" if chembl_data_found else "failed",
                "status_code": response.status_code,
                "chembl_data_found": chembl_data_found,
                "chembl_id": chembl_id,
                "sample_data": first_item,
                "db_verification": db_verification
            }
            
            if chembl_data_found:
                logger.info(f"Endpoint {endpoint} passed with ChEMBL ID {chembl_id}")
                self.results["endpoints_passed"].append(endpoint)
                self.results["chembl_data_found"] = True
            else:
                logger.error(f"Endpoint {endpoint} failed: No ChEMBL data found")
                self.results["endpoints_failed"].append(endpoint)
            
            return result
            
        except Exception as e:
            result = {
                "endpoint": endpoint,
                "status": "failed",
                "error": str(e),
                "exception_type": type(e).__name__
            }
            logger.error(f"Exception testing endpoint {endpoint}: {e}")
            self.results["endpoints_failed"].append(endpoint)
            return result
    
    def verify_with_database(self, chembl_id: str, api_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Verify that API data matches database content.
        
        Args:
            chembl_id: ChEMBL ID to verify
            api_data: API response data
            
        Returns:
            Dictionary with verification results
        """
        try:
            # Query the database for the molecule with the given ChEMBL ID
            # First try looking in the data_source column
            response = self.supabase_client.table('molecules').select('*').eq('data_source', chembl_id).execute()
            
            if not response.data:
                # Try other possible columns
                for column in ['chembl_id', 'source_id', 'external_id']:
                    response = self.supabase_client.table('molecules').select('*').eq(column, chembl_id).execute()
                    if response.data:
                        break
            
            if not response.data:
                # Try a more flexible search
                response = self.supabase_client.table('molecules').select('*').ilike('data_source', f'%{chembl_id}%').execute()
            
            if not response.data:
                return {
                    "status": "failed",
                    "error": f"No matching record found in database for ChEMBL ID {chembl_id}"
                }
            
            db_data = response.data[0]
            
            # Compare API data with database data
            matching_fields = []
            mismatched_fields = []
            
            for key in api_data:
                if key in db_data:
                    if api_data[key] == db_data[key]:
                        matching_fields.append(key)
                    else:
                        mismatched_fields.append({
                            "field": key,
                            "api_value": api_data[key],
                            "db_value": db_data[key]
                        })
            
            return {
                "status": "passed" if len(matching_fields) > len(mismatched_fields) else "failed",
                "matching_fields": matching_fields,
                "mismatched_fields": mismatched_fields,
                "db_record_found": True,
                "db_data": db_data
            }
            
        except Exception as e:
            logger.error(f"Database verification failed: {e}")
            return {
                "status": "failed",
                "error": str(e),
                "exception_type": type(e).__name__
            }
    
    def run_verification(self, auth_token: Optional[str] = None) -> Dict[str, Any]:
        """
        Run verification on all relevant endpoints.
        
        Args:
            auth_token: Optional authentication token
            
        Returns:
            Dictionary with verification results
        """
        # Define endpoints to test
        endpoints = [
            {"path": "molecules", "expected_fields": ["name", "smiles", "formula"]},
            {"path": "molecular_properties", "expected_fields": ["molecule_id", "property_type_id", "value"]},
            {"path": "property_types", "expected_fields": ["name", "description", "units"]}
        ]
        
        # Test each endpoint
        for endpoint_info in endpoints:
            endpoint = endpoint_info["path"]
            expected_fields = endpoint_info.get("expected_fields", [])
            params = endpoint_info.get("params", {})
            
            result = self.verify_endpoint(endpoint, auth_token, expected_fields, params)
            self.results["endpoints_tested"].append(endpoint)
            self.results["details"][endpoint] = result
        
        # Generate timestamp for report filename
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        report_path = os.path.join(REPORT_DIR, f"chembl_api_verification_{timestamp}.json")
        
        # Write results to report file
        with open(report_path, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        logger.info(f"Verification report saved to {report_path}")
        
        return self.results


@pytest.mark.integration
def test_api_returns_real_chembl_data():
    """Test that API endpoints return real ChEMBL data."""
    verifier = ChEMBLAPIVerifier(API_BASE_URL, SUPABASE_URL, SUPABASE_KEY)
    
    # Test with anonymous access
    logger.info("Testing API endpoints with anonymous access")
    results = verifier.run_verification()
    
    # Check if any endpoints passed
    assert len(results["endpoints_passed"]) > 0, "No endpoints passed verification"
    
    # Check if ChEMBL data was found
    assert results["chembl_data_found"], "No ChEMBL data found in API responses"
    
    # Log summary
    logger.info(f"API verification completed. {len(results['endpoints_passed'])}/{len(results['endpoints_tested'])} endpoints passed.")
    
    # Return results for potential further use
    return results


if __name__ == "__main__":
    # Run the test directly if script is executed
    logger.info("Starting ChEMBL API verification")
    verifier = ChEMBLAPIVerifier(API_BASE_URL, SUPABASE_URL, SUPABASE_KEY)
    results = verifier.run_verification()
    
    # Print summary
    print("\n=== ChEMBL API Verification Summary ===")
    print(f"Endpoints tested: {len(results['endpoints_tested'])}")
    print(f"Endpoints passed: {len(results['endpoints_passed'])}")
    print(f"Endpoints failed: {len(results['endpoints_failed'])}")
    print(f"ChEMBL data found: {results['chembl_data_found']}")
    print(f"Report saved to: {REPORT_DIR}/chembl_api_verification_*.json")
    
    # Exit with appropriate status code
    sys.exit(0 if results["chembl_data_found"] and len(results["endpoints_failed"]) == 0 else 1)
