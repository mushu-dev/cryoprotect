#!/usr/bin/env python3
"""
CryoProtect - Toxicity Optimization Tests

This script performs comprehensive testing of the toxicity data optimization,
including:
1. Functional testing of all endpoints
2. Data integrity verification
3. Caching mechanism validation
4. Performance comparison

Usage:
    python test_toxicity_optimization.py [--api-url http://localhost:5000]
"""

import os
import sys
import time
import json
import argparse
import logging
import unittest
import statistics
import requests
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f'toxicity_optimization_test_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)

class ToxicityOptimizationTest(unittest.TestCase):
    """Test suite for toxicity optimization."""
    
    api_url = "http://localhost:5000"
    test_molecule_ids = []
    
    @classmethod
    def setUpClass(cls):
        """Initialize test data."""
        # Get sample molecule IDs
        try:
            response = requests.get(f"{cls.api_url}/api/molecules?limit=5")
            data = response.json()
            cls.test_molecule_ids = [m['id'] for m in data.get('molecules', [])]
            
            if not cls.test_molecule_ids:
                logger.error("No test molecule IDs found. Tests will fail.")
                cls.test_molecule_ids = ["00000000-0000-0000-0000-000000000000"]  # Dummy ID
            
            logger.info(f"Using test molecule IDs: {cls.test_molecule_ids}")
        except Exception as e:
            logger.error(f"Failed to get test molecule IDs: {str(e)}")
            cls.test_molecule_ids = ["00000000-0000-0000-0000-000000000000"]  # Dummy ID
    
    def test_01_endpoints_accessible(self):
        """Test that all optimized endpoints are accessible."""
        # Define endpoints to test
        endpoints = [
            '/api/toxicity/summary/molecule/{id}',
            '/api/toxicity/ld50/molecule/{id}',
            '/api/toxicity/tox21/molecule/{id}',
            '/api/toxicity/classifications/molecule/{id}',
            '/api/toxicity/scores/molecule/{id}',
            '/api/toxicity/similar/{id}'
        ]
        
        # Test each endpoint with first test molecule
        if not self.test_molecule_ids:
            self.skipTest("No test molecule IDs available")
            
        molecule_id = self.test_molecule_ids[0]
        session = requests.Session()
        
        for endpoint in endpoints:
            url = self.api_url + endpoint.format(id=molecule_id)
            response = session.get(url)
            
            # Some molecules might not have all types of data, so accept 404
            self.assertIn(response.status_code, [200, 404], 
                         f"Endpoint {endpoint} returned status {response.status_code}")
            
            # If successful, validate response format
            if response.status_code == 200:
                data = response.json()
                self.assertIsInstance(data, dict, f"Response from {endpoint} is not a dictionary")
                self.assertIn('molecule_id', data, f"Response from {endpoint} missing molecule_id field")
    
    def test_02_data_integrity(self):
        """Test data integrity between original and optimized endpoints."""
        if not self.test_molecule_ids:
            self.skipTest("No test molecule IDs available")
            
        # Test molecule
        molecule_id = self.test_molecule_ids[0]
        session = requests.Session()
        
        # Get data from original endpoint
        original_url = f"{self.api_url}/api/toxicity/molecule/{molecule_id}"
        original_response = session.get(original_url)
        
        # Skip if no data in original endpoint
        if original_response.status_code == 404:
            self.skipTest(f"No toxicity data for molecule {molecule_id}")
        
        # Get data from optimized endpoint
        optimized_url = f"{self.api_url}/api/toxicity/summary/molecule/{molecule_id}"
        optimized_response = session.get(optimized_url)
        
        # Both should be successful
        self.assertEqual(original_response.status_code, 200, 
                         f"Original endpoint returned {original_response.status_code}")
        self.assertEqual(optimized_response.status_code, 200,
                         f"Optimized endpoint returned {optimized_response.status_code}")
        
        # Check data structure
        original_data = original_response.json()
        optimized_data = optimized_response.json()
        
        self.assertIn('molecule_id', original_data, "Original data missing molecule_id")
        self.assertIn('molecule_id', optimized_data, "Optimized data missing molecule_id")
        
        # Verify molecule ID matches
        self.assertEqual(original_data['molecule_id'], optimized_data['molecule_id'],
                         "Molecule IDs don't match between endpoints")
    
    def test_03_caching_mechanism(self):
        """Test ETag caching mechanism."""
        if not self.test_molecule_ids:
            self.skipTest("No test molecule IDs available")
            
        # Test molecule
        molecule_id = self.test_molecule_ids[0]
        session = requests.Session()
        
        # Test URLs
        urls = [
            f"{self.api_url}/api/toxicity/summary/molecule/{molecule_id}",
            f"{self.api_url}/api/toxicity/ld50/molecule/{molecule_id}",
            f"{self.api_url}/api/toxicity/scores/molecule/{molecule_id}"
        ]
        
        for url in urls:
            # Make initial request
            response = session.get(url)
            
            # Skip if data not found
            if response.status_code == 404:
                logger.info(f"No data found for {url}, skipping cache test")
                continue
                
            # Verify response has ETag
            self.assertIn('ETag', response.headers, f"ETag header missing for {url}")
            etag = response.headers['ETag']
            
            # Test If-None-Match with correct ETag
            response2 = session.get(url, headers={'If-None-Match': etag})
            self.assertEqual(response2.status_code, 304, 
                            f"Expected 304 Not Modified for {url} with matching ETag")
            
            # Test If-None-Match with incorrect ETag
            response3 = session.get(url, headers={'If-None-Match': 'invalid-etag'})
            self.assertEqual(response3.status_code, 200,
                            f"Expected 200 OK for {url} with non-matching ETag")
            
            # Verify Cache-Control header exists
            self.assertIn('Cache-Control', response.headers, 
                         f"Cache-Control header missing for {url}")
    
    def test_04_bulk_endpoint(self):
        """Test bulk toxicity endpoint."""
        if len(self.test_molecule_ids) < 2:
            self.skipTest("Not enough test molecule IDs for bulk test")
        
        session = requests.Session()
        
        # Test bulk endpoint with 2 molecules
        test_ids = self.test_molecule_ids[:2]
        bulk_url = f"{self.api_url}/api/toxicity/bulk/molecules"
        
        payload = {
            'molecule_ids': test_ids,
            'data_type': 'summary'
        }
        
        response = session.post(bulk_url, json=payload)
        self.assertEqual(response.status_code, 200, 
                         f"Bulk endpoint returned {response.status_code}")
        
        # Verify response structure
        data = response.json()
        self.assertIn('molecule_count', data, "Bulk response missing molecule_count")
        self.assertIn('results', data, "Bulk response missing results")
        self.assertEqual(data['molecule_count'], len(test_ids), 
                         "Molecule count mismatch in bulk response")
    
    def test_05_performance_comparison(self):
        """Compare performance between original and optimized endpoints."""
        if not self.test_molecule_ids:
            self.skipTest("No test molecule IDs available")
        
        session = requests.Session()
        results = {'original': [], 'optimized': []}
        
        # Test each molecule
        for molecule_id in self.test_molecule_ids:
            # Original endpoint
            original_url = f"{self.api_url}/api/toxicity/molecule/{molecule_id}"
            start_time = time.time()
            original_response = session.get(original_url)
            original_time = time.time() - start_time
            
            # Skip if no data
            if original_response.status_code == 404:
                continue
            
            # Optimized endpoint
            optimized_url = f"{self.api_url}/api/toxicity/summary/molecule/{molecule_id}"
            start_time = time.time()
            optimized_response = session.get(optimized_url)
            optimized_time = time.time() - start_time
            
            # Record times
            results['original'].append(original_time)
            results['optimized'].append(optimized_time)
        
        # Skip if no successful requests
        if not results['original'] or not results['optimized']:
            self.skipTest("No successful requests for performance comparison")
        
        # Calculate averages
        avg_original = statistics.mean(results['original'])
        avg_optimized = statistics.mean(results['optimized'])
        
        # Log results
        logger.info(f"Average response time (original): {avg_original:.4f}s")
        logger.info(f"Average response time (optimized): {avg_optimized:.4f}s")
        
        # Calculate improvement percentage
        improvement = (avg_original - avg_optimized) / avg_original * 100
        logger.info(f"Performance improvement: {improvement:.2f}%")
        
        # Assert some improvement (at least 30%)
        self.assertGreater(improvement, 30, 
                          f"Performance improvement ({improvement:.2f}%) below 30% threshold")

def main():
    parser = argparse.ArgumentParser(description="Test toxicity optimization")
    parser.add_argument("--api-url", default="http://localhost:5000", help="Base URL for the API")
    args = parser.parse_args()
    
    # Set API URL for tests
    ToxicityOptimizationTest.api_url = args.api_url
    
    # Create test suite
    test_suite = unittest.TestLoader().loadTestsFromTestCase(ToxicityOptimizationTest)
    
    # Run tests
    test_result = unittest.TextTestRunner(verbosity=2).run(test_suite)
    
    # Return result code
    return 0 if test_result.wasSuccessful() else 1

if __name__ == "__main__":
    sys.exit(main())