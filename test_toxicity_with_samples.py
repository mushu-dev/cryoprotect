#!/usr/bin/env python3
"""
CryoProtect - Toxicity Optimization Test with Sample Data

This script tests the toxicity optimization using sample data,
allowing for more reliable testing even without a full database.

Usage:
    python test_toxicity_with_samples.py
"""

import time
import json
import logging
import requests
import uuid
import random
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Try to load real molecule IDs if available
def get_sample_molecules(count=5):
    """Get sample molecule IDs for testing."""
    try:
        # Try to get real molecules from the database
        import os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from database.db import get_db

        try:
            db = get_db()
            result = db.execute("SELECT id FROM molecules ORDER BY id LIMIT %s", (count,))
            molecules = [str(row[0]) for row in result]
            if molecules:
                logger.info(f"Using {len(molecules)} real molecule IDs from database")
                return molecules
        except Exception as db_err:
            logger.warning(f"Could not retrieve real molecule IDs: {str(db_err)}")

        # Try to load from ChEMBL test file
        if os.path.exists("chembl_test_ids.json"):
            try:
                with open("chembl_test_ids.json", "r") as f:
                    chembl_ids = json.load(f)
                if chembl_ids and len(chembl_ids) >= count:
                    molecules = chembl_ids[:count]
                    logger.info(f"Using {len(molecules)} molecule IDs from ChEMBL test file")
                    return molecules
            except Exception as file_err:
                logger.warning(f"Could not load ChEMBL test IDs: {str(file_err)}")
    except ImportError:
        pass

    # Generate random UUIDs as fallback
    molecules = [str(uuid.uuid4()) for _ in range(count)]
    logger.info(f"Using {len(molecules)} generated UUIDs for testing")
    return molecules

# Generate sample molecules
SAMPLE_MOLECULES = get_sample_molecules(10)

class ToxicityPerformanceTest:
    """Test toxicity optimization with sample data."""
    
    def __init__(self, api_url="http://localhost:5000"):
        """Initialize test with API URL."""
        self.api_url = api_url
        self.session = requests.Session()
        self.results = {
            'original': defaultdict(list),
            'optimized': defaultdict(list)
        }
    
    def time_request(self, url, etag=None):
        """Time a request and return elapsed time and response."""
        headers = {}
        if etag:
            headers['If-None-Match'] = etag
            
        start_time = time.time()
        response = self.session.get(url, headers=headers)
        elapsed = time.time() - start_time
        
        return response, elapsed
    
    def test_original_endpoints(self):
        """Test original toxicity endpoints."""
        logger.info("Testing original endpoints")
        
        for molecule_id in SAMPLE_MOLECULES:
            # Test toxicity data endpoint
            url = f"{self.api_url}/api/toxicity/molecule/{molecule_id}"
            response, elapsed = self.time_request(url)
            
            self.results['original']['toxicity'].append({
                'url': url,
                'status_code': response.status_code,
                'time': elapsed,
                'response_size': len(response.content) if response.content else 0
            })
            
            # Test toxicity score endpoint
            url = f"{self.api_url}/api/toxicity/scores/molecule/{molecule_id}"
            response, elapsed = self.time_request(url)
            
            self.results['original']['scores'].append({
                'url': url,
                'status_code': response.status_code,
                'time': elapsed,
                'response_size': len(response.content) if response.content else 0
            })
    
    def test_optimized_endpoints(self):
        """Test optimized toxicity endpoints."""
        logger.info("Testing optimized endpoints")
        
        for molecule_id in SAMPLE_MOLECULES:
            # Test toxicity summary endpoint
            url = f"{self.api_url}/api/toxicity/summary/molecule/{molecule_id}"
            response, elapsed = self.time_request(url)
            
            result = {
                'url': url,
                'status_code': response.status_code,
                'time': elapsed,
                'response_size': len(response.content) if response.content else 0
            }
            
            # Test ETag caching if available
            if 'ETag' in response.headers:
                etag = response.headers['ETag']
                cached_response, cached_elapsed = self.time_request(url, etag)
                
                result['etag'] = etag
                result['cached_status_code'] = cached_response.status_code
                result['cached_time'] = cached_elapsed
            
            self.results['optimized']['summary'].append(result)
            
            # Test LD50 endpoint
            url = f"{self.api_url}/api/toxicity/ld50/molecule/{molecule_id}"
            response, elapsed = self.time_request(url)
            
            self.results['optimized']['ld50'].append({
                'url': url,
                'status_code': response.status_code,
                'time': elapsed,
                'response_size': len(response.content) if response.content else 0
            })
            
            # Test Tox21 endpoint
            url = f"{self.api_url}/api/toxicity/tox21/molecule/{molecule_id}"
            response, elapsed = self.time_request(url)
            
            self.results['optimized']['tox21'].append({
                'url': url,
                'status_code': response.status_code,
                'time': elapsed,
                'response_size': len(response.content) if response.content else 0
            })
            
            # Test classification endpoint
            url = f"{self.api_url}/api/toxicity/classifications/molecule/{molecule_id}"
            response, elapsed = self.time_request(url)
            
            self.results['optimized']['classifications'].append({
                'url': url,
                'status_code': response.status_code,
                'time': elapsed,
                'response_size': len(response.content) if response.content else 0
            })
    
    def test_bulk_endpoints(self):
        """Test bulk endpoints."""
        logger.info("Testing bulk endpoints")
        
        # Test bulk endpoint
        url = f"{self.api_url}/api/toxicity/bulk/molecules"
        payload = {
            'molecule_ids': SAMPLE_MOLECULES[:3],
            'data_type': 'summary'
        }
        
        start_time = time.time()
        response = self.session.post(url, json=payload)
        elapsed = time.time() - start_time
        
        self.results['optimized']['bulk'] = {
            'url': url,
            'status_code': response.status_code,
            'time': elapsed,
            'response_size': len(response.content) if response.content else 0,
            'molecules_count': len(payload['molecule_ids'])
        }
    
    def run_tests(self):
        """Run all tests."""
        logger.info("Starting toxicity performance tests")
        
        # Test original endpoints
        self.test_original_endpoints()
        
        # Test optimized endpoints
        self.test_optimized_endpoints()
        
        # Test bulk endpoints
        self.test_bulk_endpoints()
        
        # Print results
        self.print_results()
        
        logger.info("Tests completed")
        return self.results
    
    def print_results(self):
        """Print test results."""
        print("\n" + "=" * 80)
        print("TOXICITY OPTIMIZATION TEST RESULTS")
        print("=" * 80)
        
        print("\nOriginal Endpoints:")
        print("-" * 80)
        for endpoint, timings in self.results['original'].items():
            if timings:
                avg_time = sum(t['time'] for t in timings) / len(timings)
                avg_size = sum(t['response_size'] for t in timings) / len(timings) / 1024  # KB
                success_rate = sum(1 for t in timings if t['status_code'] == 200) / len(timings) * 100
                print(f"Endpoint: {endpoint}")
                print(f"  Average Time: {avg_time:.4f}s")
                print(f"  Average Size: {avg_size:.2f} KB")
                print(f"  Success Rate: {success_rate:.0f}%")
        
        print("\nOptimized Endpoints:")
        print("-" * 80)
        for endpoint, timings in self.results['optimized'].items():
            if endpoint != 'bulk' and timings:
                avg_time = sum(t['time'] for t in timings) / len(timings)
                avg_size = sum(t['response_size'] for t in timings) / len(timings) / 1024  # KB
                success_rate = sum(1 for t in timings if t['status_code'] == 200) / len(timings) * 100
                
                # Check for ETag results
                has_etag = any('etag' in t for t in timings)
                if has_etag:
                    etag_timings = [t for t in timings if 'etag' in t]
                    avg_cached_time = sum(t['cached_time'] for t in etag_timings) / len(etag_timings)
                    cache_hit_rate = sum(1 for t in etag_timings if t['cached_status_code'] == 304) / len(etag_timings) * 100
                    
                    print(f"Endpoint: {endpoint}")
                    print(f"  Average Time: {avg_time:.4f}s")
                    print(f"  Average Cached Time: {avg_cached_time:.4f}s")
                    print(f"  Cache Improvement: {((avg_time - avg_cached_time) / avg_time * 100):.1f}%")
                    print(f"  Cache Hit Rate: {cache_hit_rate:.0f}%")
                    print(f"  Average Size: {avg_size:.2f} KB")
                    print(f"  Success Rate: {success_rate:.0f}%")
                else:
                    print(f"Endpoint: {endpoint}")
                    print(f"  Average Time: {avg_time:.4f}s")
                    print(f"  Average Size: {avg_size:.2f} KB")
                    print(f"  Success Rate: {success_rate:.0f}%")
        
        # Print bulk endpoint results
        if 'bulk' in self.results['optimized']:
            bulk = self.results['optimized']['bulk']
            print("\nBulk Endpoint:")
            print("-" * 80)
            print(f"Molecules: {bulk['molecules_count']}")
            print(f"Time: {bulk['time']:.4f}s")
            print(f"Size: {bulk['response_size'] / 1024:.2f} KB")
            print(f"Status: {'Success' if bulk['status_code'] == 200 else 'Failed'}")
        
        # Comparison (if we have both original and optimized results)
        if self.results['original'] and self.results['optimized']:
            print("\n" + "=" * 80)
            print("PERFORMANCE COMPARISON")
            print("=" * 80)
            
            # Compare toxicity data vs. summary
            if 'toxicity' in self.results['original'] and 'summary' in self.results['optimized']:
                orig_times = [t['time'] for t in self.results['original']['toxicity']]
                opt_times = [t['time'] for t in self.results['optimized']['summary']]
                
                if orig_times and opt_times:
                    orig_avg = sum(orig_times) / len(orig_times)
                    opt_avg = sum(opt_times) / len(opt_times)
                    improvement = (orig_avg - opt_avg) / orig_avg * 100
                    
                    print(f"Basic toxicity data:")
                    print(f"  Original average time: {orig_avg:.4f}s")
                    print(f"  Optimized average time: {opt_avg:.4f}s")
                    print(f"  Performance improvement: {improvement:.1f}%")
            
            # Compare score endpoints
            if 'scores' in self.results['original'] and 'summary' in self.results['optimized']:
                orig_times = [t['time'] for t in self.results['original']['scores']]
                opt_times = [t['time'] for t in self.results['optimized']['summary']]
                
                if orig_times and opt_times:
                    orig_avg = sum(orig_times) / len(orig_times)
                    opt_avg = sum(opt_times) / len(opt_times)
                    improvement = (orig_avg - opt_avg) / orig_avg * 100
                    
                    print(f"\nToxicity scores:")
                    print(f"  Original average time: {orig_avg:.4f}s")
                    print(f"  Optimized average time: {opt_avg:.4f}s")
                    print(f"  Performance improvement: {improvement:.1f}%")
    
    def save_results(self, filename):
        """Save results to a file."""
        try:
            with open(filename, 'w') as f:
                json.dump(self.results, f, indent=2)
            logger.info(f"Results saved to {filename}")
        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")

def main():
    # Run tests
    test = ToxicityPerformanceTest()
    results = test.run_tests()
    
    # Save results
    test.save_results("toxicity_performance_results.json")
    
    return 0

if __name__ == "__main__":
    main()