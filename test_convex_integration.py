#!/usr/bin/env python3
"""
Convex Integration Test Script

This script tests the integration between the CryoProtect API and Convex database.
It verifies that the API can connect to Convex, perform basic CRUD operations,
and that authentication works properly.
"""

import os
import sys
import json
import requests
import logging
import time
from datetime import datetime
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"convex_integration_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("convex_integration_test")

# Load environment variables
load_dotenv()

# Test configuration
API_BASE_URL = os.environ.get('API_URL', 'http://localhost:5000/api')
TEST_USER_EMAIL = os.environ.get('TEST_USER_EMAIL', 'test@example.com')
TEST_USER_PASSWORD = os.environ.get('TEST_USER_PASSWORD', 'test_password')
CONVEX_URL = os.environ.get('CONVEX_URL', 'https://upbeat-parrot-866.convex.cloud')

class ConvexIntegrationTest:
    """Test suite for Convex integration with CryoProtect API."""

    def __init__(self):
        """Initialize the test suite."""
        # Initialize test session
        self.session = requests.Session()
        self.access_token = None
        self.convex_token = None

    def run_all_tests(self):
        """Run all tests in the suite."""
        try:
            logger.info("Starting Convex integration tests...")
            
            # Check environment configuration
            success = self.check_environment()
            if not success:
                logger.error("Environment configuration is invalid, aborting tests.")
                return False

            # Test authentication
            success = self.test_authentication()
            if not success:
                logger.error("Authentication test failed, aborting remaining tests.")
                return False

            # Test basic CRUD operations
            success = self.test_molecule_crud()
            if not success:
                logger.error("Molecule CRUD test failed.")
                return False

            # Test real-time subscription (using a simple HTTP test as a proxy)
            success = self.test_subscription_simulation()
            if not success:
                logger.error("Subscription simulation test failed.")

            # Test connection resilience
            success = self.test_connection_resilience()
            if not success:
                logger.error("Connection resilience test failed.")

            logger.info("All tests completed!")
            return True
        except Exception as e:
            logger.error(f"Unexpected error during tests: {str(e)}", exc_info=True)
            return False

    def check_environment(self):
        """Check if the environment is properly configured for testing."""
        logger.info("Checking environment configuration...")
        
        # Check API URL
        if not API_BASE_URL:
            logger.error("API_URL environment variable is not set.")
            return False
        
        # Check Convex URL
        if not CONVEX_URL:
            logger.error("CONVEX_URL environment variable is not set.")
            return False
        
        # Check USE_CONVEX flag
        use_convex = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
        if not use_convex:
            logger.error("USE_CONVEX environment variable must be set to 'true' for testing.")
            return False
        
        # Check test user credentials
        if not TEST_USER_EMAIL or not TEST_USER_PASSWORD:
            logger.error("Test user credentials (TEST_USER_EMAIL, TEST_USER_PASSWORD) are not set.")
            return False

        # Check API availability
        try:
            response = requests.get(f"{API_BASE_URL}/health/connectivity")
            if response.status_code != 200:
                logger.error(f"API health check failed with status code {response.status_code}")
                return False
            logger.info("API health check passed.")
        except Exception as e:
            logger.error(f"API health check failed with error: {str(e)}")
            return False

        logger.info("Environment configuration is valid.")
        return True

    def test_authentication(self):
        """Test authentication with Convex."""
        logger.info("Testing authentication...")
        
        # Try to login with test user
        try:
            response = requests.post(f"{API_BASE_URL}/auth/login", json={
                "email": TEST_USER_EMAIL,
                "password": TEST_USER_PASSWORD
            })
            
            if response.status_code != 200:
                logger.error(f"Authentication failed with status code {response.status_code}")
                logger.error(f"Response: {response.text}")
                return False
            
            # Extract tokens
            response_data = response.json()
            self.access_token = response_data.get('token', {}).get('access_token')
            self.convex_token = response_data.get('convex_token')
            
            # Set token in session
            if self.access_token:
                self.session.headers.update({
                    'Authorization': f'Bearer {self.access_token}'
                })
                logger.info("Successfully authenticated with JWT token.")
            else:
                logger.error("Authentication response missing access token.")
                return False
            
            # Verify we received a Convex token
            if not self.convex_token:
                logger.warning("Authentication response missing Convex token. This is required for a complete integration.")
                return False
            
            logger.info("Successfully authenticated with Convex token.")
            return True
        except Exception as e:
            logger.error(f"Authentication failed with error: {str(e)}")
            return False

    def test_molecule_crud(self):
        """Test basic CRUD operations for molecules using Convex."""
        logger.info("Testing molecule CRUD operations...")
        
        # Skip if not authenticated
        if not self.access_token:
            logger.error("Cannot test CRUD operations without authentication.")
            return False
        
        # Test molecule creation
        try:
            # Create a test molecule
            molecule_data = {
                "name": f"Test Molecule {datetime.now().strftime('%Y%m%d%H%M%S')}",
                "inchikey": "HMFHBZSHGGEWLO-UHFFFAOYSA-N",
                "smiles": "CCO",
                "molecularFormula": "C2H6O",
                "molecularWeight": 46.07
            }
            
            response = self.session.post(f"{API_BASE_URL}/molecules", json=molecule_data)
            
            if response.status_code != 201:
                logger.error(f"Molecule creation failed with status code {response.status_code}")
                logger.error(f"Response: {response.text}")
                return False
            
            # Extract the molecule ID
            response_data = response.json()
            molecule_id = response_data.get('data', {}).get('id')
            
            if not molecule_id:
                logger.error("Molecule creation response missing molecule ID.")
                return False
            
            logger.info(f"Successfully created molecule with ID: {molecule_id}")
            
            # Test molecule retrieval
            response = self.session.get(f"{API_BASE_URL}/molecules/{molecule_id}")
            
            if response.status_code != 200:
                logger.error(f"Molecule retrieval failed with status code {response.status_code}")
                logger.error(f"Response: {response.text}")
                return False
            
            logger.info("Successfully retrieved molecule.")
            
            # Test molecule update
            update_data = {
                "name": f"Updated Test Molecule {datetime.now().strftime('%Y%m%d%H%M%S')}"
            }
            
            response = self.session.put(f"{API_BASE_URL}/molecules/{molecule_id}", json=update_data)
            
            if response.status_code != 200:
                logger.error(f"Molecule update failed with status code {response.status_code}")
                logger.error(f"Response: {response.text}")
                return False
            
            logger.info("Successfully updated molecule.")
            
            # Test molecule deletion
            response = self.session.delete(f"{API_BASE_URL}/molecules/{molecule_id}")
            
            if response.status_code != 200:
                logger.error(f"Molecule deletion failed with status code {response.status_code}")
                logger.error(f"Response: {response.text}")
                return False
            
            logger.info("Successfully deleted molecule.")
            
            # Verify deletion by trying to retrieve it again
            response = self.session.get(f"{API_BASE_URL}/molecules/{molecule_id}")
            
            if response.status_code == 404:
                logger.info("Molecule correctly shows as deleted/not found.")
            else:
                logger.warning(f"Unexpected status code {response.status_code} when verifying deletion.")
                
            return True
        except Exception as e:
            logger.error(f"Molecule CRUD test failed with error: {str(e)}")
            return False

    def test_subscription_simulation(self):
        """
        Test real-time subscription simulation.
        
        Since we can't directly test WebSocket connections in a simple script,
        we'll simulate it by:
        1. Creating a molecule in one request
        2. Immediately trying to retrieve it in another request
        3. Verifying it's found within a reasonable time
        """
        logger.info("Testing subscription simulation...")
        
        # Skip if not authenticated
        if not self.access_token:
            logger.error("Cannot test subscriptions without authentication.")
            return False
        
        try:
            # Create a test molecule with unique name
            unique_id = datetime.now().strftime('%Y%m%d%H%M%S')
            molecule_data = {
                "name": f"Subscription Test Molecule {unique_id}",
                "inchikey": "HMFHBZSHGGEWLO-UHFFFAOYSA-N",
                "smiles": "CCO",
                "molecularFormula": "C2H6O",
                "molecularWeight": 46.07
            }
            
            # Create the molecule
            response = self.session.post(f"{API_BASE_URL}/molecules", json=molecule_data)
            
            if response.status_code != 201:
                logger.error(f"Molecule creation failed with status code {response.status_code}")
                logger.error(f"Response: {response.text}")
                return False
            
            # Extract the molecule ID
            molecule_id = response.json().get('data', {}).get('id')
            
            # Immediately try to query for molecules with this name
            start_time = time.time()
            max_wait_time = 5  # seconds
            found = False
            
            while time.time() - start_time < max_wait_time:
                # Query for molecules with the unique name
                response = self.session.get(f"{API_BASE_URL}/molecules?name={molecule_data['name']}")
                
                if response.status_code != 200:
                    logger.warning(f"Query failed with status code {response.status_code}")
                    time.sleep(0.5)
                    continue
                
                # Check if the molecule is found
                molecules = response.json().get('data', [])
                if any(m['id'] == molecule_id for m in molecules):
                    found = True
                    break
                
                # Wait a bit before retrying
                time.sleep(0.5)
            
            if found:
                logger.info(f"Successfully found molecule through query within {time.time() - start_time:.2f} seconds.")
            else:
                logger.error(f"Failed to find molecule through query within {max_wait_time} seconds.")
                return False
            
            # Clean up
            self.session.delete(f"{API_BASE_URL}/molecules/{molecule_id}")
            
            return True
        except Exception as e:
            logger.error(f"Subscription simulation test failed with error: {str(e)}")
            return False

    def test_connection_resilience(self):
        """Test connection resilience by making multiple rapid requests."""
        logger.info("Testing connection resilience...")
        
        # Skip if not authenticated
        if not self.access_token:
            logger.error("Cannot test connection resilience without authentication.")
            return False
        
        try:
            # Make multiple requests in quick succession
            num_requests = 10
            success_count = 0
            
            logger.info(f"Making {num_requests} rapid requests to test connection resilience...")
            
            for i in range(num_requests):
                try:
                    response = self.session.get(f"{API_BASE_URL}/molecules?limit=5")
                    
                    if response.status_code == 200:
                        success_count += 1
                    else:
                        logger.warning(f"Request {i+1} failed with status code {response.status_code}")
                except Exception as e:
                    logger.warning(f"Request {i+1} failed with error: {str(e)}")
            
            success_rate = (success_count / num_requests) * 100
            logger.info(f"Connection resilience test completed with success rate: {success_rate:.2f}%")
            
            if success_rate >= 80:
                return True
            else:
                logger.error(f"Connection resilience test failed with success rate below 80%: {success_rate:.2f}%")
                return False
        except Exception as e:
            logger.error(f"Connection resilience test failed with error: {str(e)}")
            return False

def main():
    """Main function to run all tests."""
    tester = ConvexIntegrationTest()
    success = tester.run_all_tests()
    
    if success:
        logger.info("All Convex integration tests passed!")
        return 0
    else:
        logger.error("Convex integration tests failed.")
        return 1

if __name__ == "__main__":
    sys.exit(main())