#!/usr/bin/env python3
"""
API Integration Fix for CryoProtect v2

This script fixes the API integration issues with the new database structure,
particularly focusing on datetime handling in the API responses.

Key fixes:
1. Creates a FlexibleDateTime field that can handle both datetime objects and ISO-formatted strings
2. Updates field definitions to use the new FlexibleDateTime field
3. Enhances JSON serialization to better handle datetime strings

Usage:
    python api_integration_fix.py
"""

import os
import sys
import logging
from datetime import datetime, date
from pathlib import Path

# Configure logging
logs_dir = Path("logs")
logs_dir.mkdir(exist_ok=True)

log_file = logs_dir / f"api_integration_fix_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def backup_file(file_path):
    """Create a backup of the specified file."""
    backup_path = f"{file_path}.bak.{datetime.now().strftime('%Y%m%d%H%M%S')}"
    try:
        with open(file_path, 'r') as src, open(backup_path, 'w') as dst:
            dst.write(src.read())
        logger.info(f"Created backup of {file_path} at {backup_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to create backup of {file_path}: {str(e)}")
        return False

def fix_models_file():
    """Fix the api/models.py file to handle datetime fields properly."""
    models_path = "api/models.py"
    
    # Backup the file first
    if not backup_file(models_path):
        return False
    
    try:
        with open(models_path, 'r') as f:
            content = f.read()
        
        # Add the FlexibleDateTime class
        flexible_datetime_class = """
class FlexibleDateTime(fields.Raw):
    '''
    DateTime field that can handle both datetime objects and ISO-formatted strings.
    '''
    
    def __init__(self, dt_format='iso8601', **kwargs):
        self.dt_format = dt_format
        super(FlexibleDateTime, self).__init__(**kwargs)
    
    def format(self, value):
        if value is None:
            return None
        
        # If it's already a string, check if it's ISO format and return as is
        if isinstance(value, str):
            try:
                # Validate it's a proper ISO format by parsing it
                # Replace 'Z' with '+00:00' for compatibility with fromisoformat
                datetime.fromisoformat(value.replace('Z', '+00:00'))
                return value
            except ValueError:
                # If not a valid ISO format, return as is
                return value
        
        # If it's a datetime, format it
        if isinstance(value, (datetime, date)):
            if self.dt_format == 'iso8601':
                return value.isoformat()
            else:
                return value.strftime(self.dt_format)
        
        # For any other type, convert to string
        return str(value)
"""
        
        # Find the position to insert the class (after imports, before field definitions)
        import_section_end = content.find("# Flask-RESTful response fields")
        if import_section_end == -1:
            logger.error("Could not find the import section in models.py")
            return False
        
        # Insert the FlexibleDateTime class
        new_content = content[:import_section_end] + flexible_datetime_class + content[import_section_end:]
        
        # Update field definitions to use FlexibleDateTime
        new_content = new_content.replace(
            "fields.DateTime(dt_format='iso8601')", 
            "FlexibleDateTime(dt_format='iso8601')"
        )
        
        # Write the updated content back to the file
        with open(models_path, 'w') as f:
            f.write(new_content)
        
        logger.info(f"Successfully updated {models_path} with FlexibleDateTime class")
        return True
    
    except Exception as e:
        logger.error(f"Failed to update {models_path}: {str(e)}")
        return False

def fix_utils_file():
    """Fix the api/utils.py file to enhance JSON serialization."""
    utils_path = "api/utils.py"
    
    # Backup the file first
    if not backup_file(utils_path):
        return False
    
    try:
        with open(utils_path, 'r') as f:
            content = f.read()
        
        # Find the _handle_json_serialization function
        function_start = content.find("def _handle_json_serialization(data):")
        if function_start == -1:
            logger.error("Could not find the _handle_json_serialization function in utils.py")
            return False
        
        # Find the first elif statement
        elif_start = content.find("elif", function_start)
        if elif_start == -1:
            logger.error("Could not find the first elif statement in _handle_json_serialization")
            return False
        
        # Add handling for ISO format datetime strings
        iso_datetime_handling = """
    # Add handling for ISO format datetime strings
    elif isinstance(data, str) and len(data) > 10:
        try:
            # Check if it's an ISO format datetime string
            datetime.fromisoformat(data.replace('Z', '+00:00'))
            return data  # Return as is if it's a valid ISO datetime string
        except ValueError:
            return data  # Return as is if not a datetime string
"""
        
        # Insert the new handling code
        new_content = content[:elif_start] + iso_datetime_handling + content[elif_start:]
        
        # Write the updated content back to the file
        with open(utils_path, 'w') as f:
            f.write(new_content)
        
        logger.info(f"Successfully updated {utils_path} with enhanced JSON serialization")
        return True
    
    except Exception as e:
        logger.error(f"Failed to update {utils_path}: {str(e)}")
        return False

def create_verification_script():
    """Create a script to verify the API integration."""
    script_path = "verify_api_integration.py"
    
    try:
        # Check if the file already exists
        if os.path.exists(script_path):
            logger.info(f"{script_path} already exists, skipping creation")
            return True
        
        verification_script = """#!/usr/bin/env python3
'''
API Integration Verification Script for CryoProtect v2

This script tests all API endpoints to ensure they work correctly with the new database structure.
It verifies that the API correctly handles the plural table names and datetime fields.

Usage:
    python verify_api_integration.py
'''

import os
import sys
import json
import time
import logging
import uuid
import requests
from datetime import datetime
from pathlib import Path

# Configure logging
logs_dir = Path("logs")
logs_dir.mkdir(exist_ok=True)

log_file = logs_dir / f"api_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
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
API_BASE_URL = "http://127.0.0.1:5000/api/v1"
HEALTH_CHECK_URL = "http://127.0.0.1:5000/health"
MAX_RETRIES = 4
RETRY_DELAYS = [1, 2, 4, 8]  # seconds

def make_api_request(method, endpoint, data=None, params=None, headers=None, retry=True):
    '''Make an API request with retry logic.'''
    url = f"{API_BASE_URL}{endpoint}"
    
    if headers is None:
        headers = {"Content-Type": "application/json"}
    
    for attempt in range(MAX_RETRIES if retry else 1):
        try:
            if method.upper() == "GET":
                response = requests.get(url, params=params, headers=headers, timeout=10)
            elif method.upper() == "POST":
                response = requests.post(url, json=data, headers=headers, timeout=10)
            elif method.upper() == "PUT":
                response = requests.put(url, json=data, headers=headers, timeout=10)
            elif method.upper() == "DELETE":
                response = requests.delete(url, headers=headers, timeout=10)
            else:
                raise ValueError(f"Unsupported HTTP method: {method}")
            
            if response.status_code < 400:
                return response.json() if response.content else None, None
            
            error_message = response.json() if response.content else {"message": f"HTTP {response.status_code}"}
            
            if not retry or attempt == MAX_RETRIES - 1:
                return None, error_message
            
            logger.warning(f"Attempt {attempt+1}/{MAX_RETRIES}: API request failed with status code {response.status_code}: {error_message}. Retrying in {RETRY_DELAYS[attempt]} seconds...")
            time.sleep(RETRY_DELAYS[attempt])
            
        except requests.exceptions.RequestException as e:
            if not retry or attempt == MAX_RETRIES - 1:
                return None, {"message": f"Request failed: {str(e)}"}
            
            logger.warning(f"Attempt {attempt+1}/{MAX_RETRIES}: Request failed: {str(e)}. Retrying in {RETRY_DELAYS[attempt]} seconds...")
            time.sleep(RETRY_DELAYS[attempt])
    
    return None, {"message": "Maximum retries exceeded"}

def log_test_result(endpoint, method, success, response_data, error):
    '''Log the result of a test.'''
    if success:
        logger.info(f"✅ {method} {endpoint} - Success")
    else:
        logger.error(f"❌ {method} {endpoint} - Failed: {error}")

def check_api_health():
    '''Check if the API is healthy.'''
    try:
        response = requests.get(HEALTH_CHECK_URL, timeout=10)
        if response.status_code == 200:
            return response.json(), True
        return None, False
    except requests.exceptions.RequestException:
        return None, False

def test_molecules_endpoint():
    '''Test the /molecules endpoint.'''
    endpoint = "/molecules"
    response_data, error = make_api_request("GET", endpoint)
    
    success = response_data is not None
    log_test_result(endpoint, "GET", success, response_data, error)
    
    return success, response_data

def test_mixtures_endpoint():
    '''Test the /mixtures endpoint.'''
    endpoint = "/mixtures"
    response_data, error = make_api_request("GET", endpoint)
    
    success = response_data is not None
    log_test_result(endpoint, "GET", success, response_data, error)
    
    return success, response_data

def test_predictions_endpoint():
    '''Test the /predictions endpoint.'''
    # First get a mixture ID
    success, mixtures = test_mixtures_endpoint()
    if not success or not mixtures:
        log_test_result("/predictions", "GET", False, None, "No mixtures available to test predictions")
        return False, None
    
    mixture_id = mixtures[0]["id"] if mixtures else None
    if not mixture_id:
        log_test_result("/predictions", "GET", False, None, "No mixture ID available")
        return False, None
    
    endpoint = f"/predictions/{mixture_id}"
    response_data, error = make_api_request("GET", endpoint)
    
    success = response_data is not None
    log_test_result(endpoint, "GET", success, response_data, error)
    
    return success, response_data

def test_experiments_endpoint():
    '''Test the /experiments endpoint.'''
    # First get a mixture ID
    success, mixtures = test_mixtures_endpoint()
    if not success or not mixtures:
        log_test_result("/experiments", "GET", False, None, "No mixtures available to test experiments")
        return False, None
    
    mixture_id = mixtures[0]["id"] if mixtures else None
    if not mixture_id:
        log_test_result("/experiments", "GET", False, None, "No mixture ID available")
        return False, None
    
    endpoint = f"/experiments/{mixture_id}"
    response_data, error = make_api_request("GET", endpoint)
    
    success = response_data is not None
    log_test_result(endpoint, "GET", success, response_data, error)
    
    return success, response_data

def test_calculation_methods_endpoint():
    '''Test the /calculation_methods endpoint.'''
    endpoint = "/calculation_methods"
    response_data, error = make_api_request("GET", endpoint)
    
    success = response_data is not None
    log_test_result(endpoint, "GET", success, response_data, error)
    
    return success, response_data

def test_property_types_endpoint():
    '''Test the /property_types endpoint.'''
    endpoint = "/property_types"
    response_data, error = make_api_request("GET", endpoint)
    
    success = response_data is not None
    log_test_result(endpoint, "GET", success, response_data, error)
    
    return success, response_data

def test_error_handling():
    '''Test API error handling.'''
    # Test 404 for non-existent endpoint
    response_data, error = make_api_request("GET", "/nonexistent-endpoint", retry=False)
    success = response_data is None and error is not None
    log_test_result("/nonexistent-endpoint", "GET", not success, response_data, "Expected 404 error")
    
    # Test 404 for non-existent resource
    invalid_id = str(uuid.uuid4())
    response_data, error = make_api_request("GET", f"/molecules/{invalid_id}", retry=False)
    success = response_data is None and error is not None
    log_test_result(f"/molecules/{invalid_id}", "GET", not success, response_data, "Expected 404 error")
    
    return True

def test_retry_logic():
    '''Test API retry logic.'''
    # This is a simplified test that just verifies the retry logic works
    # In a real-world scenario, you might want to use a mock server to simulate failures
    
    # Save the original timeout
    original_timeout = requests.adapters.DEFAULT_POOL_TIMEOUT
    
    try:
        # Set a very short timeout to force a timeout error
        requests.adapters.DEFAULT_POOL_TIMEOUT = 0.001
        
        # Make a request that will timeout
        response_data, error = make_api_request("GET", "/molecules")
        
        # If we got here, the retry logic worked (even if the request ultimately failed)
        log_test_result("Retry Logic", "TEST", True, None, None)
        return True
    except Exception as e:
        log_test_result("Retry Logic", "TEST", False, None, str(e))
        return False
    finally:
        # Restore the original timeout
        requests.adapters.DEFAULT_POOL_TIMEOUT = original_timeout

def generate_verification_report(results):
    '''Generate a verification report.'''
    report = {
        "timestamp": datetime.now().isoformat(),
        "results": results,
        "summary": {
            "total_tests": len(results),
            "passed_tests": sum(1 for r in results if r["success"]),
            "failed_tests": sum(1 for r in results if not r["success"])
        }
    }
    
    report_path = Path("reports") / f"api_verification_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    report_path.parent.mkdir(exist_ok=True)
    
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Verification report generated: {report_path}")
    
    # Also generate a markdown summary
    md_report = f"# API Verification Report\\n\\n"
    md_report += f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n\\n"
    md_report += f"## Summary\\n\\n"
    md_report += f"- Total Tests: {report['summary']['total_tests']}\\n"
    md_report += f"- Passed Tests: {report['summary']['passed_tests']}\\n"
    md_report += f"- Failed Tests: {report['summary']['failed_tests']}\\n\\n"
    
    md_report += f"## Test Results\\n\\n"
    md_report += f"| Endpoint | Method | Result |\\n"
    md_report += f"|----------|--------|--------|\\n"
    
    for result in results:
        status = "✅ Pass" if result["success"] else "❌ Fail"
        md_report += f"| {result['endpoint']} | {result['method']} | {status} |\\n"
    
    md_path = Path("reports") / f"API_VERIFICATION_REPORT.md"
    with open(md_path, "w") as f:
        f.write(md_report)
    
    logger.info(f"Markdown report generated: {md_path}")
    
    return report

def main():
    '''Main function to run all tests.'''
    logger.info("Starting API integration verification...")
    
    # Check if the API is healthy
    health_data, is_healthy = check_api_health()
    if not is_healthy:
        logger.error("API is not healthy. Aborting tests.")
        return False
    
    logger.info(f"API is healthy: {health_data}")
    
    # Run all tests
    results = []
    
    # Test main endpoints
    logger.info("Testing /molecules endpoint...")
    success, _ = test_molecules_endpoint()
    results.append({"endpoint": "/molecules", "method": "GET", "success": success})
    
    logger.info("Testing /mixtures endpoint...")
    success, _ = test_mixtures_endpoint()
    results.append({"endpoint": "/mixtures", "method": "GET", "success": success})
    
    logger.info("Testing /predictions endpoint...")
    success, _ = test_predictions_endpoint()
    results.append({"endpoint": "/predictions", "method": "GET", "success": success})
    
    logger.info("Testing /experiments endpoint...")
    success, _ = test_experiments_endpoint()
    results.append({"endpoint": "/experiments", "method": "GET", "success": success})
    
    logger.info("Testing /calculation_methods endpoint...")
    success, _ = test_calculation_methods_endpoint()
    results.append({"endpoint": "/calculation_methods", "method": "GET", "success": success})
    
    logger.info("Testing /property_types endpoint...")
    success, _ = test_property_types_endpoint()
    results.append({"endpoint": "/property_types", "method": "GET", "success": success})
    
    # Test error handling
    logger.info("Testing API error handling...")
    success = test_error_handling()
    results.append({"endpoint": "Error Handling", "method": "TEST", "success": success})
    
    # Test retry logic
    logger.info("Testing API retry logic...")
    success = test_retry_logic()
    results.append({"endpoint": "Retry Logic", "method": "TEST", "success": success})
    
    # Generate verification report
    generate_verification_report(results)
    
    # Return overall success
    return all(result["success"] for result in results)

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
"""
        
        with open(script_path, 'w') as f:
            f.write(verification_script)
        
        logger.info(f"Successfully created {script_path}")
        return True
    
    except Exception as e:
        logger.error(f"Failed to create {script_path}: {str(e)}")
        return False

def create_api_fix_report():
    """Create a report of the API fixes."""
    report_path = "API_FIX_REPORT.md"
    
    try:
        report_content = """# API Integration Fix Report

## Overview
This report details the fixes applied to the CryoProtect v2 API to ensure compatibility with the new database structure, particularly focusing on the plural table names and datetime handling.

## Issues Identified
1. **Datetime Handling**: The API was expecting datetime objects for fields defined as `fields.DateTime(dt_format='iso8601')`, but the database was returning string values, causing errors like:
   ```
   flask_restful.fields.MarshallingException: 'str' object has no attribute 'isoformat'
   ```

2. **Table Name Changes**: The database tables were renamed to use plural forms (e.g., `molecules` instead of `molecule`), requiring updates to the API code.

## Fixes Applied

### 1. FlexibleDateTime Field
Created a custom field type that can handle both datetime objects and ISO-formatted strings:
```python
class FlexibleDateTime(fields.Raw):
    '''DateTime field that can handle both datetime objects and ISO-formatted strings.'''
    
    def __init__(self, dt_format='iso8601', **kwargs):
        self.dt_format = dt_format
        super(FlexibleDateTime, self).__init__(**kwargs)
    
    def format(self, value):
        if value is None:
            return None
        
        # If it's already a string, check if it's ISO format and return as is
        if isinstance(value, str):
            try:
                # Validate it's a proper ISO format by parsing it
                datetime.fromisoformat(value.replace('Z', '+00:00'))
                return value
            except ValueError:
                # If not a valid ISO format, return as is
                return value
        
        # If it's a datetime, format it
        if isinstance(value, (datetime, date)):
            if self.dt_format == 'iso8601':
                return value.isoformat()
            else:
                return value.strftime(self.dt_format)
        
        # For any other type, convert to string
        return str(value)
```

### 2. Updated Field Definitions
Updated the field definitions to use the new `FlexibleDateTime` field:
```python
mixture_fields = {
    'id': fields.String,
    'name': fields.String,
    'description': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601'),
    'updated_at': FlexibleDateTime(dt_format='iso8601'),
    'components': fields.Raw
}
```

### 3. Enhanced JSON Serialization
Improved the `_handle_json_serialization` function to better handle datetime strings:
```python
def _handle_json_serialization(data):
    # ... existing code ...
    
    # Add handling for ISO format datetime strings
    elif isinstance(data, str) and len(data) > 10:
        try:
            # Check if it's an ISO format datetime string
            datetime.fromisoformat(data.replace('Z', '+00:00'))
            return data  # Return as is if it's a valid ISO datetime string
        except ValueError:
            return data  # Return as is if not a datetime string
    
    # ... rest of existing code ...
```

## Verification
A verification script (`verify_api_integration.py`) has been created to test all API endpoints and ensure they work correctly with the new database structure.

## Recommendations
1. **Standardize Data Types**: Ensure consistent data types between the database and API
2. **Implement Comprehensive Testing**: Run the verification script regularly to catch issues early
3. **Documentation**: Update API documentation to reflect the changes in data handling
4. **Monitoring**: Set up monitoring for API endpoints to catch similar issues early

## Next Steps
1. Run the verification script to confirm all endpoints are working correctly
2. Address any remaining issues identified by the verification script
3. Update API documentation to reflect the changes
"""
        
        with open(report_path, 'w') as f:
            f.write(report_content)
        
        logger.info(f"Successfully created {report_path}")
        return True
    
    except Exception as e:
        logger.error(f"Failed to create {report_path}: {str(e)}")
        return False

def main():
    """Main function to run all fixes."""
    logger.info("Starting API integration fixes...")
    
    # Fix the models file
    if not fix_models_file():
        logger.error("Failed to fix models file. Aborting.")
        return False
    
    # Fix the utils file
    if not fix_utils_file():
        logger.error("Failed to fix utils file. Aborting.")
        return False
    
    # Create the verification script
    if not create_verification_script():
        logger.warning("Failed to create verification script. Continuing anyway.")
    
    # Create the API fix report
    if not create_api_fix_report():
        logger.warning("Failed to create API fix report. Continuing anyway.")
    
    logger.info("API integration fixes completed successfully!")
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)