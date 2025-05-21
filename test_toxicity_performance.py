#!/usr/bin/env python3
"""
CryoProtect - Test Toxicity Data Performance

This script benchmarks the performance of the original and optimized toxicity
data endpoints to measure the improvement in response times and resource usage.

Usage:
    python test_toxicity_performance.py
"""

import os
import time
import random
import argparse
import logging
import json
import statistics
import requests
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f'toxicity_performance_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger()

# Global variables
results = {
    'original': {},
    'optimized': {}
}
results_lock = Lock()

def time_request(url, session, label, test_etags=False):
    """
    Time a request to the given URL and record the results.
    
    Args:
        url: The URL to request
        session: The requests session to use
        label: The label for this test
        test_etags: Whether to test ETag caching
    
    Returns:
        dict: Request timing information
    """
    # Initial request
    headers = {}
    
    # Measure the time
    start_time = time.time()
    response = session.get(url, headers=headers)
    end_time = time.time()
    
    # Check if request was successful
    if response.status_code != 200:
        logger.warning(f"Request failed with status {response.status_code}: {url}")
        return None
        
    # Calculate timings
    timing = {
        'url': url,
        'status_code': response.status_code,
        'time': end_time - start_time,
        'response_size': len(response.content)
    }
    
    # Test ETag caching if enabled
    if test_etags and 'ETag' in response.headers:
        etag = response.headers['ETag']
        
        # Make request with If-None-Match header
        etag_headers = {'If-None-Match': etag}
        etag_start_time = time.time()
        etag_response = session.get(url, headers=etag_headers)
        etag_end_time = time.time()
        
        # Record ETag response
        timing['etag_status_code'] = etag_response.status_code
        timing['etag_time'] = etag_end_time - etag_start_time
    
    # Record the result
    with results_lock:
        if label not in results:
            results[label] = {}
        
        endpoint = url.split('/')[-2] + '/' + url.split('/')[-1]
        if endpoint not in results[label]:
            results[label][endpoint] = []
        
        results[label][endpoint].append(timing)
    
    return timing

def get_molecule_ids(api_url, limit=10):
    """Get a list of molecule IDs from the database."""
    try:
        url = f"{api_url}/api/molecules?limit={limit}"
        response = requests.get(url)
        data = response.json()
        return [m['id'] for m in data.get('molecules', [])]
    except Exception as e:
        logger.error(f"Failed to get molecule IDs: {str(e)}")
        return []

def run_performance_test(api_url, use_optimized=False, num_requests=10, concurrent=2, test_etags=True):
    """
    Run performance tests on the toxicity endpoints.
    
    Args:
        api_url: Base URL for the API
        use_optimized: Whether to use the optimized endpoints
        num_requests: Number of requests to make per endpoint
        concurrent: Number of concurrent requests
        test_etags: Whether to test ETag caching
    """
    # Get molecule IDs
    molecule_ids = get_molecule_ids(api_url, limit=20)
    
    if not molecule_ids:
        logger.error("No molecule IDs found")
        return
    
    logger.info(f"Retrieved {len(molecule_ids)} molecule IDs")
    
    # Define endpoints to test
    if use_optimized:
        logger.info("Testing optimized endpoints")
        endpoints = [
            '/api/toxicity/summary/molecule/{id}',
            '/api/toxicity/ld50/molecule/{id}',
            '/api/toxicity/tox21/molecule/{id}',
            '/api/toxicity/classifications/molecule/{id}',
            '/api/toxicity/scores/molecule/{id}',
            '/api/toxicity/similar/{id}'
        ]
        label = 'optimized'
    else:
        logger.info("Testing original endpoints")
        endpoints = [
            '/api/toxicity/molecule/{id}',
            '/api/toxicity/scores/molecule/{id}'
        ]
        label = 'original'
    
    # Create a session for requests
    session = requests.Session()
    
    # Generate URLs for testing
    urls = []
    for _ in range(num_requests):
        for endpoint in endpoints:
            # Choose a random molecule ID
            molecule_id = random.choice(molecule_ids)
            # Construct the URL
            url = api_url + endpoint.format(id=molecule_id)
            urls.append((url, label))
    
    # Shuffle URLs to randomize the order
    random.shuffle(urls)
    
    # Execute requests with thread pool
    with ThreadPoolExecutor(max_workers=concurrent) as executor:
        futures = []
        for url, lbl in urls:
            futures.append(executor.submit(time_request, url, session, lbl, test_etags))
        
        # Wait for all requests to complete
        for future in futures:
            future.result()

def print_results():
    """Print the test results in a formatted table."""
    print("\n" + "=" * 80)
    print("PERFORMANCE TEST RESULTS")
    print("=" * 80)
    
    # Print summary for each test type
    for label, endpoints in results.items():
        print(f"\n{label.upper()} ENDPOINTS:")
        print("-" * 80)
        print(f"{'Endpoint':<40} {'Avg Time (s)':<15} {'Med Time (s)':<15} {'Resp Size (KB)':<15}")
        print("-" * 80)
        
        for endpoint, timings in endpoints.items():
            # Calculate stats
            times = [t['time'] for t in timings]
            sizes = [t['response_size'] for t in timings]
            
            avg_time = statistics.mean(times) if times else 0
            med_time = statistics.median(times) if times else 0
            avg_size = statistics.mean(sizes) / 1024 if sizes else 0
            
            # Print row
            print(f"{endpoint:<40} {avg_time:<15.4f} {med_time:<15.4f} {avg_size:<15.2f}")
        
        # Print ETag stats if available
        etag_timings = []
        for endpoint_timings in endpoints.values():
            for timing in endpoint_timings:
                if 'etag_time' in timing:
                    etag_timings.append(timing['etag_time'])
        
        if etag_timings:
            etag_avg = statistics.mean(etag_timings)
            etag_med = statistics.median(etag_timings)
            print("\nETag Response Time:")
            print(f"Average: {etag_avg:.4f}s, Median: {etag_med:.4f}s")
    
    # Print comparison if both original and optimized results exist
    if 'original' in results and 'optimized' in results:
        print("\n" + "=" * 80)
        print("PERFORMANCE COMPARISON")
        print("=" * 80)
        
        # Get all times for each type
        original_times = []
        for timings in results['original'].values():
            original_times.extend([t['time'] for t in timings])
        
        optimized_times = []
        for timings in results['optimized'].values():
            optimized_times.extend([t['time'] for t in timings])
        
        # Calculate averages
        if original_times and optimized_times:
            orig_avg = statistics.mean(original_times)
            opt_avg = statistics.mean(optimized_times)
            improvement = (orig_avg - opt_avg) / orig_avg * 100
            
            print(f"Original average response time: {orig_avg:.4f}s")
            print(f"Optimized average response time: {opt_avg:.4f}s")
            print(f"Performance improvement: {improvement:.2f}%")

def save_results(output_file):
    """Save the results to a JSON file."""
    try:
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Results saved to {output_file}")
    except Exception as e:
        logger.error(f"Failed to save results: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description="Test toxicity data performance")
    parser.add_argument("--api-url", default="http://localhost:5000", help="Base URL for the API")
    parser.add_argument("--num-requests", type=int, default=10, help="Number of requests per endpoint")
    parser.add_argument("--concurrent", type=int, default=2, help="Number of concurrent requests")
    parser.add_argument("--output", default="toxicity_performance_results.json", help="Output file for results")
    parser.add_argument("--test-etags", action="store_true", help="Test ETag caching")
    args = parser.parse_args()
    
    try:
        logger.info("Starting performance test")
        
        # Run tests for original endpoints
        run_performance_test(
            args.api_url, 
            use_optimized=False, 
            num_requests=args.num_requests,
            concurrent=args.concurrent,
            test_etags=args.test_etags
        )
        
        # Run tests for optimized endpoints
        run_performance_test(
            args.api_url, 
            use_optimized=True, 
            num_requests=args.num_requests,
            concurrent=args.concurrent,
            test_etags=args.test_etags
        )
        
        # Print and save results
        print_results()
        save_results(args.output)
        
        logger.info("Performance test completed")
    except Exception as e:
        logger.error(f"Performance test failed: {str(e)}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())