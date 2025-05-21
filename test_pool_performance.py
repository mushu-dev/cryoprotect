#!/usr/bin/env python3
"""
Test script for connection pool performance.
"""

import requests
import time
import threading
from concurrent.futures import ThreadPoolExecutor
import argparse
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

def make_request(i):
    """Make a request to the API"""
    start = time.time()
    url = os.environ.get("TEST_URL", "http://localhost:5000/api/molecules")
    
    try:
        response = requests.get(f'{url}?limit=1', timeout=10)
        status_code = response.status_code
        duration = time.time() - start
        
        return {
            'request_id': i,
            'status_code': status_code,
            'duration': duration,
            'success': 200 <= status_code < 300
        }
    except Exception as e:
        duration = time.time() - start
        return {
            'request_id': i,
            'status_code': 0,
            'duration': duration,
            'success': False,
            'error': str(e)
        }

def run_test(num_requests=100, max_workers=10):
    """Run a load test with multiple concurrent requests"""
    print(f"Running load test with {num_requests} requests using {max_workers} concurrent workers...")
    print("----------------------------------------")
    
    results = []
    start_time = time.time()
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(make_request, i) for i in range(num_requests)]
        
        for future in futures:
            results.append(future.result())
    
    total_time = time.time() - start_time
    
    # Calculate statistics
    durations = [r['duration'] for r in results]
    successes = [r for r in results if r.get('success', False)]
    failures = [r for r in results if not r.get('success', False)]
    
    avg_duration = sum(durations) / len(durations) if durations else 0
    min_duration = min(durations) if durations else 0
    max_duration = max(durations) if durations else 0
    
    # Print results
    print(f"Test completed in {total_time:.3f} seconds")
    print(f"Results for {num_requests} requests with {max_workers} workers:")
    print(f"  Success rate: {len(successes)}/{num_requests} ({len(successes)/num_requests*100:.1f}%)")
    print(f"  Average request time: {avg_duration:.3f}s")
    print(f"  Min request time: {min_duration:.3f}s")
    print(f"  Max request time: {max_duration:.3f}s")
    print(f"  Requests per second: {num_requests/total_time:.1f}")
    
    # Print failure information if any
    if failures:
        print(f"\nFailures: {len(failures)}")
        for f in failures[:5]:  # Show first 5 failures
            print(f"  Request {f['request_id']}: {f.get('error', 'Unknown error')}")
        
        if len(failures) > 5:
            print(f"  ... and {len(failures) - 5} more")
    
    # Get pool stats if available
    try:
        stats_url = os.environ.get("STATS_URL", "http://localhost:5000/pool/stats")
        response = requests.get(stats_url, timeout=5)
        if response.status_code == 200:
            stats = response.json()
            if 'stats' in stats:
                print("\nPool stats:")
                for k, v in stats['stats'].items():
                    print(f"  {k}: {v}")
    except Exception as e:
        print(f"\nFailed to get pool stats: {str(e)}")
    
    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test connection pool performance')
    parser.add_argument('--requests', type=int, default=100, help='Number of requests to make')
    parser.add_argument('--workers', type=int, default=10, help='Number of concurrent workers')
    parser.add_argument('--url', type=str, help='URL to test (defaults to http://localhost:5000/api/molecules)')
    parser.add_argument('--stats-url', type=str, help='URL for pool stats (defaults to http://localhost:5000/pool/stats)')
    
    args = parser.parse_args()
    
    if args.url:
        os.environ["TEST_URL"] = args.url
    
    if args.stats_url:
        os.environ["STATS_URL"] = args.stats_url
    
    run_test(num_requests=args.requests, max_workers=args.workers)