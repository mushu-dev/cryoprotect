#!/usr/bin/env python
"""
Example usage of the ResilientPubChemClient.

This script demonstrates how to use the ResilientPubChemClient
for various PubChem API operations with resilience features.
"""

import os
import logging
import time
from pprint import pprint

from client import ResilientPubChemClient

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def main():
    """Main function demonstrating ResilientPubChemClient usage."""
    print("Initializing ResilientPubChemClient...")
    
    # Create client with custom settings
    client = ResilientPubChemClient(
        cache_dir="example_cache",
        weekday_requests_per_second=2.0,
        weekend_requests_per_second=5.0,
        max_retries=3,
        failure_threshold=3,
        recovery_timeout=30,
        cache_ttl=86400,  # 1 day
        memory_cache_size=100,
        enable_scheduler=True
    )
    
    # Example 1: Get molecule properties
    print("\nExample 1: Get molecule properties for ethanol (CID 702)")
    properties = client.get_molecule_properties(702)
    print("Ethanol properties:")
    pprint(properties)
    
    # Example 2: Get compound synonyms
    print("\nExample 2: Get synonyms for glucose (CID 5793)")
    synonyms = client.get_compound_synonyms(5793)
    print("Glucose synonyms:")
    pprint(synonyms.get("Synonyms", [])[:10])  # Show first 10 synonyms
    
    # Example 3: Search for compounds
    print("\nExample 3: Search for compounds with 'aspirin'")
    results = client.search_compounds("aspirin")
    print(f"Found {len(results.get('CIDs', []))} compounds for 'aspirin'")
    print(f"First 5 CIDs: {results.get('CIDs', [])[:5]}")
    
    # Example 4: Demonstrate caching
    print("\nExample 4: Demonstrate caching")
    print("First call (should hit API):")
    start_time = time.time()
    client.get_molecule_properties(5288)  # Ibuprofen
    first_call_time = time.time() - start_time
    
    print("Second call (should hit cache):")
    start_time = time.time()
    client.get_molecule_properties(5288)  # Ibuprofen again
    second_call_time = time.time() - start_time
    
    print(f"First call time: {first_call_time:.4f} seconds")
    print(f"Second call time: {second_call_time:.4f} seconds")
    print(f"Speed improvement: {first_call_time / second_call_time:.1f}x faster")
    
    # Example 5: Schedule weekend job
    print("\nExample 5: Schedule weekend job for bulk operations")
    cids = [2244, 5793, 5288, 702, 1983]  # Acetaminophen, Glucose, Ibuprofen, Ethanol, Caffeine
    job_id = client.prefetch_molecule_properties(cids)
    print(f"Scheduled job with ID: {job_id}")
    
    # Check job status
    status = client.get_job_status(job_id)
    print(f"Job status: {status}")
    
    # Example 6: Get statistics
    print("\nExample 6: Get client statistics")
    cache_stats = client.get_cache_stats()
    rate_limiter_stats = client.get_rate_limiter_stats()
    circuit_breaker_stats = client.get_circuit_breaker_stats()
    
    print("Cache statistics:")
    pprint(cache_stats)
    
    print("\nRate limiter statistics:")
    pprint(rate_limiter_stats)
    
    print("\nCircuit breaker statistics:")
    pprint(circuit_breaker_stats)
    
    # Clean up
    print("\nCleaning up...")
    client.scheduler.stop()
    
    print("\nExample completed successfully!")

if __name__ == "__main__":
    main()