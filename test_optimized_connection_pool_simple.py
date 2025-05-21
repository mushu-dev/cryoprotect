#!/usr/bin/env python3
"""
Simple test for the optimized connection pool implementation.
This test verifies basic functionality and resilience features.
"""

import os
import sys
import time
import logging
import threading
import psycopg2
from typing import List, Dict, Any
from optimized_connection_pool import OptimizedConnectionPool, CircuitBreaker

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configuration for the connection pool
POOL_CONFIG = {
    'host': os.getenv('SUPABASE_DB_HOST', 'localhost'),
    'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
    'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
    'user': os.getenv('SUPABASE_DB_USER', 'postgres'),
    'password': os.getenv('SUPABASE_DB_PASSWORD', ''),
    'min_connections': 2,
    'max_connections': 5,
    'connection_timeout': 10,
    'connection_lifetime': 60,  # Short for testing
    'idle_timeout': 30,  # Short for testing
    'health_check_interval': 5,  # Short for testing
    'retry_attempts': 3,
    'initial_retry_delay': 0.1,
    'max_retry_delay': 1.0,
    'retry_jitter_factor': 0.1,
    'circuit_breaker_threshold': 3,
    'circuit_breaker_timeout': 10,  # Short for testing
    'circuit_breaker_reset': 2,
    'metrics_reporting_interval': 10  # Short for testing
}

def test_basic_functionality():
    """Test basic connection pool functionality."""
    logger.info("Testing basic connection pool functionality")
    
    # Create pool
    pool = OptimizedConnectionPool(POOL_CONFIG)
    
    try:
        # Perform a simple query
        with pool.get_connection() as conn:
            with conn.cursor() as cursor:
                cursor.execute("SELECT 1 as test")
                result = cursor.fetchone()
                assert result[0] == 1, f"Expected 1, got {result[0]}"
                logger.info("Basic query test passed")
        
        # Wait for metrics to be collected
        time.sleep(2)
        
        # Get metrics
        metrics = pool.get_metrics()
        logger.info(f"Connection pool metrics: {metrics}")
        
        # Verify metrics
        assert metrics['created'] >= 1, "Expected at least 1 connection created"
        assert metrics['peak_connections'] >= 1, "Expected peak connections >= 1"
        
        logger.info("Basic functionality test passed")
    finally:
        # Close the pool
        pool.close()

def test_concurrent_connections():
    """Test concurrent connections."""
    logger.info("Testing concurrent connections")
    
    # Create pool
    pool = OptimizedConnectionPool(POOL_CONFIG)
    results = []
    errors = []
    
    def worker(worker_id):
        """Worker thread function."""
        try:
            with pool.get_connection() as conn:
                with conn.cursor() as cursor:
                    # Simulate some work
                    cursor.execute("SELECT pg_sleep(0.5), %s as worker_id", (worker_id,))
                    result = cursor.fetchone()
                    results.append(result[1])
        except Exception as e:
            errors.append(str(e))
    
    try:
        # Create and start threads
        threads = []
        for i in range(10):  # Try to get 10 connections from a pool with max=5
            thread = threading.Thread(target=worker, args=(i,))
            threads.append(thread)
            thread.start()
        
        # Wait for threads to finish
        for thread in threads:
            thread.join()
        
        # Wait for metrics to be collected
        time.sleep(2)
        
        # Get metrics
        metrics = pool.get_metrics()
        logger.info(f"Connection pool metrics: {metrics}")
        
        # Verify results
        logger.info(f"Successful workers: {len(results)}, Errors: {len(errors)}")
        assert len(results) + len(errors) == 10, "Expected 10 results or errors"
        
        # Verify pool metrics
        assert metrics['peak_connections'] <= POOL_CONFIG['max_connections'], \
            f"Peak connections ({metrics['peak_connections']}) exceeded max ({POOL_CONFIG['max_connections']})"
        
        logger.info("Concurrent connections test passed")
    finally:
        # Close the pool
        pool.close()

def test_circuit_breaker():
    """Test circuit breaker functionality."""
    logger.info("Testing circuit breaker functionality")
    
    # Create a circuit breaker
    breaker = CircuitBreaker(threshold=3, timeout=5, reset_count=2)
    
    # Initially should be closed
    assert breaker.get_state() == CircuitBreaker.CLOSED, "Circuit breaker should start closed"
    
    # Record some failures
    for i in range(3):
        breaker.record_failure()
    
    # Should be open now
    assert breaker.get_state() == CircuitBreaker.OPEN, "Circuit breaker should be open after 3 failures"
    
    # Request should be rejected
    assert not breaker.allow_request(), "Circuit breaker should reject requests when open"
    
    # Wait for timeout
    logger.info("Waiting for circuit breaker timeout (5 seconds)...")
    time.sleep(6)
    
    # Should be half-open now
    assert breaker.allow_request(), "Circuit breaker should allow requests after timeout"
    assert breaker.get_state() == CircuitBreaker.HALF_OPEN, "Circuit breaker should be half-open after timeout"
    
    # Record successes
    for i in range(2):
        breaker.record_success()
    
    # Should be closed now
    assert breaker.get_state() == CircuitBreaker.CLOSED, "Circuit breaker should be closed after 2 successes"
    
    logger.info("Circuit breaker test passed")

def main():
    """Run all tests."""
    try:
        test_basic_functionality()
        test_concurrent_connections()
        test_circuit_breaker()
        
        logger.info("All tests passed!")
        return 0
    except Exception as e:
        logger.error(f"Test failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())