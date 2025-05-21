#!/usr/bin/env python3
"""
Stress Test for Connection Pool

This script puts the connection pool through rigorous testing to verify:
1. Performance under high concurrent load
2. Resilience to connection failures
3. Circuit breaker effectiveness
4. Connection pooling efficiency
5. Connection lifecycle management
6. Dynamic pool scaling capabilities

Used as part of the database verification process to ensure connection pooling is
properly optimized before moving to production.
"""

import os
import sys
import time
import random
import logging
import threading
import statistics
import argparse
import json
import concurrent.futures
from datetime import datetime, timedelta
from typing import Dict, List, Any, Tuple

# Conditionally import the appropriate connection pool
try:
    from optimized_connection_pool import (
        OptimizedConnectionPool, ConnectionManager, 
        execute_query_with_retry, transaction_context
    )
    USING_OPTIMIZED_POOL = True
except ImportError:
    from connection_pool_wrapper import (
        ConnectionPoolWrapper as ConnectionPool, 
        ConnectionManager
    )
    USING_OPTIMIZED_POOL = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f"logs/connection_pool_stress_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    ]
)

# Create directory if it doesn't exist
os.makedirs("logs", exist_ok=True)
os.makedirs("reports", exist_ok=True)

logger = logging.getLogger(__name__)

# Test parameters
DEFAULT_WORKERS = 20
DEFAULT_ITERATIONS = 50
DEFAULT_DURATION = 30  # seconds
DEFAULT_MAX_POOL_SIZE = 20
DEFAULT_TEST_TIMEOUT = 120  # seconds


class StressTestResults:
    """Class to collect and analyze stress test results."""
    
    def __init__(self):
        self.test_start_time = datetime.now()
        self.test_end_time = None
        self.total_queries = 0
        self.successful_queries = 0
        self.failed_queries = 0
        self.connection_errors = 0
        self.transaction_errors = 0
        self.query_times = []
        self.connection_times = []
        self.worker_results = []
        self.pool_metrics = {}
        self.pool_metrics_series = []
        self.error_counts = {
            'timeout': 0,
            'connection_reset': 0,
            'permission_denied': 0,
            'circuit_open': 0,
            'other': 0
        }
        
    def record_query_time(self, time_seconds: float):
        """Record a query execution time."""
        self.query_times.append(time_seconds)
        
    def record_connection_time(self, time_seconds: float):
        """Record a connection acquisition time."""
        self.connection_times.append(time_seconds)
        
    def add_worker_result(self, result: Dict[str, Any]):
        """Add results from a worker."""
        self.worker_results.append(result)
        self.total_queries += result.get('total_queries', 0)
        self.successful_queries += result.get('successful_queries', 0)
        self.failed_queries += result.get('failed_queries', 0)
        self.connection_errors += result.get('connection_errors', 0)
        self.transaction_errors += result.get('transaction_errors', 0)
        
        # Add query times from worker
        self.query_times.extend(result.get('query_times', []))
        self.connection_times.extend(result.get('connection_times', []))
        
        # Count error types
        for error_type, count in result.get('error_counts', {}).items():
            if error_type in self.error_counts:
                self.error_counts[error_type] += count
            else:
                self.error_counts[error_type] = count
                
    def record_pool_metrics(self, metrics: Dict[str, Any]):
        """Record pool metrics at a point in time."""
        # Store the metrics with a timestamp
        timestamped_metrics = {
            'timestamp': datetime.now().isoformat(),
            'metrics': metrics
        }
        self.pool_metrics_series.append(timestamped_metrics)
        
        # Store the most recent metrics
        self.pool_metrics = metrics
        
    def finalize(self):
        """Finalize results and calculate statistics."""
        self.test_end_time = datetime.now()
        
    def get_test_duration_seconds(self) -> float:
        """Get the test duration in seconds."""
        if self.test_end_time is None:
            return (datetime.now() - self.test_start_time).total_seconds()
        return (self.test_end_time - self.test_start_time).total_seconds()
        
    def get_queries_per_second(self) -> float:
        """Calculate queries per second."""
        duration = self.get_test_duration_seconds()
        if duration > 0:
            return self.successful_queries / duration
        return 0
        
    def get_success_rate(self) -> float:
        """Calculate success rate as percentage."""
        if self.total_queries > 0:
            return (self.successful_queries / self.total_queries) * 100
        return 0
        
    def get_query_time_stats(self) -> Dict[str, float]:
        """Get statistics about query execution times."""
        if not self.query_times:
            return {
                'min': 0,
                'max': 0,
                'avg': 0,
                'median': 0,
                'p95': 0,
                'p99': 0
            }
            
        try:
            sorted_times = sorted(self.query_times)
            p95_idx = int(len(sorted_times) * 0.95)
            p99_idx = int(len(sorted_times) * 0.99)
            
            return {
                'min': min(sorted_times),
                'max': max(sorted_times),
                'avg': statistics.mean(sorted_times),
                'median': statistics.median(sorted_times),
                'p95': sorted_times[p95_idx] if p95_idx < len(sorted_times) else sorted_times[-1],
                'p99': sorted_times[p99_idx] if p99_idx < len(sorted_times) else sorted_times[-1]
            }
        except Exception as e:
            logger.error(f"Error calculating query time stats: {e}")
            return {
                'min': 0,
                'max': 0,
                'avg': 0,
                'median': 0,
                'p95': 0,
                'p99': 0,
                'error': str(e)
            }
            
    def get_connection_time_stats(self) -> Dict[str, float]:
        """Get statistics about connection acquisition times."""
        if not self.connection_times:
            return {
                'min': 0,
                'max': 0,
                'avg': 0,
                'median': 0,
                'p95': 0,
                'p99': 0
            }
            
        try:
            sorted_times = sorted(self.connection_times)
            p95_idx = int(len(sorted_times) * 0.95)
            p99_idx = int(len(sorted_times) * 0.99)
            
            return {
                'min': min(sorted_times),
                'max': max(sorted_times),
                'avg': statistics.mean(sorted_times),
                'median': statistics.median(sorted_times),
                'p95': sorted_times[p95_idx] if p95_idx < len(sorted_times) else sorted_times[-1],
                'p99': sorted_times[p99_idx] if p99_idx < len(sorted_times) else sorted_times[-1]
            }
        except Exception as e:
            logger.error(f"Error calculating connection time stats: {e}")
            return {
                'min': 0,
                'max': 0,
                'avg': 0,
                'median': 0,
                'p95': 0,
                'p99': 0,
                'error': str(e)
            }
            
    def to_dict(self) -> Dict[str, Any]:
        """Convert results to a dictionary for serialization."""
        self.finalize()
        
        return {
            'test_start': self.test_start_time.isoformat(),
            'test_end': self.test_end_time.isoformat() if self.test_end_time else None,
            'test_duration_seconds': self.get_test_duration_seconds(),
            'using_optimized_pool': USING_OPTIMIZED_POOL,
            'total_queries': self.total_queries,
            'successful_queries': self.successful_queries,
            'failed_queries': self.failed_queries,
            'connection_errors': self.connection_errors,
            'transaction_errors': self.transaction_errors,
            'queries_per_second': self.get_queries_per_second(),
            'success_rate_percent': self.get_success_rate(),
            'query_time_stats': self.get_query_time_stats(),
            'connection_time_stats': self.get_connection_time_stats(),
            'error_counts': self.error_counts,
            'latest_pool_metrics': self.pool_metrics,
            'pool_metrics_series': self.pool_metrics_series[:10]  # Limit to avoid excessive size
        }
        
    def save_to_file(self, filename: str = None):
        """Save results to a JSON file."""
        if filename is None:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename = f"reports/connection_pool_stress_{timestamp}.json"
            
        with open(filename, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
            
        logger.info(f"Results saved to {filename}")
        return filename
        
    def generate_report(self, filename: str = None):
        """Generate a human-readable report."""
        if filename is None:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename = f"reports/connection_pool_stress_{timestamp}.md"
            
        results = self.to_dict()
        
        with open(filename, 'w') as f:
            f.write("# Connection Pool Stress Test Report\n\n")
            
            # Test Information
            f.write("## Test Information\n\n")
            f.write(f"- **Test Date**: {results['test_start']}\n")
            f.write(f"- **Test Duration**: {results['test_duration_seconds']:.2f} seconds\n")
            f.write(f"- **Pool Type**: {'Optimized' if results['using_optimized_pool'] else 'Standard'}\n\n")
            
            # Summary
            f.write("## Summary\n\n")
            f.write(f"- **Total Queries**: {results['total_queries']}\n")
            f.write(f"- **Successful Queries**: {results['successful_queries']} ({results['success_rate_percent']:.2f}%)\n")
            f.write(f"- **Failed Queries**: {results['failed_queries']}\n")
            f.write(f"- **Connection Errors**: {results['connection_errors']}\n")
            f.write(f"- **Transaction Errors**: {results['transaction_errors']}\n")
            f.write(f"- **Queries Per Second**: {results['queries_per_second']:.2f}\n\n")
            
            # Performance
            f.write("## Performance Metrics\n\n")
            
            # Query Times
            f.write("### Query Execution Times (seconds)\n\n")
            f.write("| Metric | Value |\n")
            f.write("|--------|-------|\n")
            for key, value in results['query_time_stats'].items():
                f.write(f"| {key.upper()} | {value:.6f} |\n")
            f.write("\n")
            
            # Connection Times
            f.write("### Connection Acquisition Times (seconds)\n\n")
            f.write("| Metric | Value |\n")
            f.write("|--------|-------|\n")
            for key, value in results['connection_time_stats'].items():
                f.write(f"| {key.upper()} | {value:.6f} |\n")
            f.write("\n")
            
            # Error Distribution
            f.write("### Error Distribution\n\n")
            f.write("| Error Type | Count |\n")
            f.write("|------------|-------|\n")
            for error_type, count in results['error_counts'].items():
                f.write(f"| {error_type} | {count} |\n")
            f.write("\n")
            
            # Pool Metrics
            f.write("## Pool Metrics\n\n")
            if results['latest_pool_metrics'] and 'pool' in results['latest_pool_metrics']:
                pool_metrics = results['latest_pool_metrics']['pool']
                f.write("| Metric | Value |\n")
                f.write("|--------|-------|\n")
                for key, value in pool_metrics.items():
                    if isinstance(value, (int, float)):
                        f.write(f"| {key} | {value} |\n")
                    else:
                        f.write(f"| {key} | {value} |\n")
                f.write("\n")
                
            # Circuit Breaker (if optimized pool)
            if USING_OPTIMIZED_POOL and 'circuit_breaker' in results['latest_pool_metrics']:
                f.write("### Circuit Breaker Status\n\n")
                circuit_metrics = results['latest_pool_metrics']['circuit_breaker']
                f.write("| Metric | Value |\n")
                f.write("|--------|-------|\n")
                for key, value in circuit_metrics.items():
                    f.write(f"| {key} | {value} |\n")
                f.write("\n")
                
            # Recommendations based on results
            f.write("## Recommendations\n\n")
            
            # Check if we had too many connection errors
            if results['connection_errors'] > 0.05 * results['total_queries']:
                f.write("- **High Connection Error Rate**: Consider increasing connection timeout or implementing more aggressive retry strategies.\n")
                
            # Check if p95 query time is much higher than median
            if results['query_time_stats']['p95'] > 5 * results['query_time_stats']['median']:
                f.write("- **Large Query Time Outliers**: Some queries are taking significantly longer than others. Consider investigating query patterns for optimization.\n")
                
            # Check pool utilization
            if 'latest_pool_metrics' in results and 'pool' in results['latest_pool_metrics'] and 'utilization' in results['latest_pool_metrics']['pool']:
                utilization = results['latest_pool_metrics']['pool']['utilization']
                if utilization > 0.8:
                    f.write("- **High Pool Utilization**: The connection pool was highly utilized. Consider increasing max_connections.\n")
                elif utilization < 0.3:
                    f.write("- **Low Pool Utilization**: The connection pool was underutilized. Consider reducing max_connections to free resources.\n")
                    
            f.write("\n")
            
            # Conclusion
            f.write("## Conclusion\n\n")
            if results['success_rate_percent'] > 98:
                f.write("The connection pool performed well under stress, with a high success rate and good performance metrics.\n")
            elif results['success_rate_percent'] > 90:
                f.write("The connection pool performed adequately, but there is room for optimization to improve reliability.\n")
            else:
                f.write("The connection pool struggled under stress. Consider implementing the recommendations above to improve reliability.\n")
                
        logger.info(f"Report generated at {filename}")
        return filename


def stress_test_worker(
    worker_id: int, 
    duration: int, 
    query_delay_range: Tuple[float, float] = (0.01, 0.3),
    simulate_errors: bool = False,
    use_transactions: bool = True
) -> Dict[str, Any]:
    """
    Worker function to stress test the connection pool.
    
    Args:
        worker_id: Unique identifier for the worker
        duration: How long to run the tests in seconds
        query_delay_range: Range of delays to simulate query processing (min, max)
        simulate_errors: Whether to simulate errors
        use_transactions: Whether to use transactions
        
    Returns:
        Dictionary of test results
    """
    results = {
        'worker_id': worker_id,
        'total_queries': 0,
        'successful_queries': 0,
        'failed_queries': 0,
        'connection_errors': 0,
        'transaction_errors': 0,
        'query_times': [],
        'connection_times': [],
        'error_counts': {
            'timeout': 0,
            'connection_reset': 0,
            'permission_denied': 0,
            'circuit_open': 0,
            'other': 0
        }
    }
    
    # Track start time
    start_time = time.time()
    end_time = start_time + duration
    
    iteration = 0
    
    while time.time() < end_time:
        iteration += 1
        results['total_queries'] += 1
        
        try:
            # Simulate query delay
            query_delay = random.uniform(query_delay_range[0], query_delay_range[1])
            
            # Simulate random errors if enabled
            if simulate_errors and random.random() < 0.05:  # 5% chance of error
                error_type = random.choice(['timeout', 'connection_reset', 'sql_error'])
                if error_type == 'timeout':
                    # Make query with extreme timeout
                    query_delay = 10  # 10 seconds, likely to timeout
                elif error_type == 'connection_reset':
                    # Force a connection reset error by using an intentionally malformed query
                    raise psycopg2.OperationalError("connection reset by peer")
                elif error_type == 'sql_error':
                    # Force a SQL error
                    execute_query_with_retry("SELECT invalid_column FROM nonexistent_table")
                    continue
            
            # Measure connection acquisition time
            conn_start = time.time()
            
            if use_transactions:
                # Use transaction context
                try:
                    with transaction_context() as conn:
                        # Record connection acquisition time
                        conn_time = time.time() - conn_start
                        results['connection_times'].append(conn_time)
                        
                        # Execute query
                        query_start = time.time()
                        with conn.cursor() as cursor:
                            # Use a simple sleep query to simulate work
                            cursor.execute(
                                "SELECT pg_sleep(%s)::TEXT AS sleep_result, %s AS worker_id, %s AS iteration",
                                (query_delay, worker_id, iteration)
                            )
                            
                            # Fetch results
                            result = cursor.fetchone()
                            
                        # Measure query execution time
                        query_time = time.time() - query_start
                        results['query_times'].append(query_time)
                        results['successful_queries'] += 1
                        
                except Exception as e:
                    error_str = str(e).lower()
                    results['transaction_errors'] += 1
                    
                    # Categorize errors
                    if 'timeout' in error_str or 'timed out' in error_str:
                        results['error_counts']['timeout'] += 1
                    elif 'reset' in error_str:
                        results['error_counts']['connection_reset'] += 1
                    elif 'permission' in error_str:
                        results['error_counts']['permission_denied'] += 1
                    elif 'circuit' in error_str and 'open' in error_str:
                        results['error_counts']['circuit_open'] += 1
                    else:
                        results['error_counts']['other'] += 1
                        
                    logger.debug(f"Worker {worker_id} transaction error: {e}")
            else:
                # Use simple query API
                try:
                    # Execute query with retry
                    query_start = time.time()
                    result = execute_query_with_retry(
                        "SELECT pg_sleep(%s)::TEXT AS sleep_result, %s AS worker_id, %s AS iteration",
                        {'sleep': query_delay, 'worker_id': worker_id, 'iteration': iteration}
                    )
                    
                    # Measure query execution time
                    query_time = time.time() - query_start
                    results['query_times'].append(query_time)
                    results['successful_queries'] += 1
                    
                except Exception as e:
                    error_str = str(e).lower()
                    results['failed_queries'] += 1
                    
                    # Categorize errors
                    if 'timeout' in error_str or 'timed out' in error_str:
                        results['error_counts']['timeout'] += 1
                    elif 'reset' in error_str:
                        results['error_counts']['connection_reset'] += 1
                    elif 'permission' in error_str:
                        results['error_counts']['permission_denied'] += 1
                    elif 'circuit' in error_str and 'open' in error_str:
                        results['error_counts']['circuit_open'] += 1
                    else:
                        results['error_counts']['other'] += 1
                        
                    logger.debug(f"Worker {worker_id} query error: {e}")
                    
        except Exception as e:
            # Connection error or other unexpected error
            error_str = str(e).lower()
            results['connection_errors'] += 1
            
            if 'timeout' in error_str or 'timed out' in error_str:
                results['error_counts']['timeout'] += 1
            elif 'reset' in error_str:
                results['error_counts']['connection_reset'] += 1
            elif 'permission' in error_str:
                results['error_counts']['permission_denied'] += 1
            elif 'circuit' in error_str and 'open' in error_str:
                results['error_counts']['circuit_open'] += 1
            else:
                results['error_counts']['other'] += 1
                
            logger.debug(f"Worker {worker_id} connection error: {e}")
        
        # Log progress
        if iteration % 10 == 0:
            elapsed = time.time() - start_time
            remaining = max(0, duration - elapsed)
            logger.debug(f"Worker {worker_id}: {iteration} queries, {remaining:.1f}s remaining")
            
    # Calculate actual duration
    actual_duration = time.time() - start_time
    results['duration'] = actual_duration
    
    logger.info(f"Worker {worker_id} completed {results['total_queries']} queries in {actual_duration:.2f}s "
               f"({results['successful_queries']} successful, {results['failed_queries']} failed, "
               f"{results['connection_errors']} connection errors)")
    
    return results


def stress_test_connection_pool(
    workers: int = DEFAULT_WORKERS,
    duration: int = DEFAULT_DURATION,
    optimized: bool = True,
    min_connections: int = 2,
    max_connections: int = DEFAULT_MAX_POOL_SIZE,
    simulate_errors: bool = False,
    use_transactions: bool = True,
    query_delay_range: Tuple[float, float] = (0.01, 0.5)
) -> Dict[str, Any]:
    """
    Run a stress test on the connection pool.
    
    Args:
        workers: Number of concurrent worker threads
        duration: How long to run the test in seconds
        optimized: Whether to use optimized pool (if available)
        min_connections: Minimum connections in pool
        max_connections: Maximum connections in pool
        simulate_errors: Whether to simulate errors
        use_transactions: Whether to use transactions
        query_delay_range: Range of delays to simulate query processing (min, max)
        
    Returns:
        Test results
    """
    # Check if we're using OptimizedConnectionPool vs ConnectionPoolWrapper
    global USING_OPTIMIZED_POOL
    
    # Get database config
    config = {
        'host': os.environ.get('SUPABASE_DB_HOST'),
        'port': int(os.environ.get('SUPABASE_DB_PORT', '5432')),
        'dbname': os.environ.get('SUPABASE_DB_NAME', 'postgres'),
        'user': os.environ.get('SUPABASE_DB_USER'),
        'password': os.environ.get('SUPABASE_DB_PASSWORD'),
        'min_connections': min_connections,
        'max_connections': max_connections,
        'connection_timeout': 30,  # 30 seconds
        'health_check_interval': 60  # 60 seconds
    }
    
    # Add optimized pool specific config
    if USING_OPTIMIZED_POOL and optimized:
        config.update({
            'circuit_window': 30,
            'circuit_threshold': 5,
            'circuit_timeout': 10,
            'max_retries': 3,
            'backoff_base': 2,
            'backoff_cap': 10,
            'jitter_factor': 0.1,
            'validation_query': 'SELECT 1'
        })
    
    # Initialize test results collector
    results = StressTestResults()
    
    # Initialize the appropriate pool based on availability and preference
    if USING_OPTIMIZED_POOL and optimized:
        # Use optimized pool
        pool = OptimizedConnectionPool.get_instance(config)
        logger.info(f"Using OptimizedConnectionPool with {min_connections}-{max_connections} connections")
    else:
        # Fall back to standard pool
        if USING_OPTIMIZED_POOL and not optimized:
            logger.info("Forced to use standard ConnectionPoolWrapper despite OptimizedConnectionPool availability")
            from connection_pool_wrapper import ConnectionPoolWrapper
            pool = ConnectionPoolWrapper.get_instance(config)
        else:
            # Standard pool is all we have
            pool = ConnectionPoolWrapper.get_instance(config)
        logger.info(f"Using ConnectionPoolWrapper with {min_connections}-{max_connections} connections")
        USING_OPTIMIZED_POOL = False
    
    logger.info(f"Starting stress test with {workers} workers for {duration} seconds")
    logger.info(f"Query delay range: {query_delay_range[0]:.2f}-{query_delay_range[1]:.2f} seconds")
    
    if simulate_errors:
        logger.info("Error simulation is enabled")
    
    if use_transactions:
        logger.info("Using transactions for queries")
    else:
        logger.info("Using simple queries (no transactions)")
    
    # Track start time
    start_time = time.time()
    
    # Start metrics collection thread
    stop_metrics = threading.Event()
    
    def metrics_collector():
        """Collect metrics periodically."""
        while not stop_metrics.is_set():
            try:
                # Get pool metrics if available
                if USING_OPTIMIZED_POOL:
                    metrics = pool.get_metrics()
                    results.record_pool_metrics(metrics)
                else:
                    metrics = {
                        'pool': {
                            'active_connections': getattr(pool, 'active_connections', 0),
                            'min_connections': getattr(pool, 'min_conn', min_connections),
                            'max_connections': getattr(pool, 'max_conn', max_connections)
                        }
                    }
                    results.record_pool_metrics(metrics)
            except Exception as e:
                logger.warning(f"Failed to collect metrics: {e}")
                
            # Sleep for a bit
            time.sleep(2)
    
    # Start metrics collector thread
    metrics_thread = threading.Thread(target=metrics_collector, daemon=True)
    metrics_thread.start()
    
    # Create and start worker threads
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        # Submit worker tasks
        futures = [
            executor.submit(
                stress_test_worker, 
                i, 
                duration,
                query_delay_range,
                simulate_errors,
                use_transactions
            ) 
            for i in range(workers)
        ]
        
        # Wait for all tasks to complete or timeout
        try:
            for future in concurrent.futures.as_completed(futures, timeout=DEFAULT_TEST_TIMEOUT):
                try:
                    worker_result = future.result()
                    # Add worker result to overall results
                    results.add_worker_result(worker_result)
                except Exception as e:
                    logger.error(f"Worker failed: {e}")
                    # Count as a worker error in results
                    results.failed_queries += 1
        except concurrent.futures.TimeoutError:
            logger.error(f"Test timeout after {DEFAULT_TEST_TIMEOUT} seconds!")
            for future in futures:
                if not future.done():
                    future.cancel()
    
    # Stop metrics collection
    stop_metrics.set()
    metrics_thread.join(timeout=5)
    
    # Finalize results
    results.finalize()
    
    # Calculate overall statistics
    total_time = time.time() - start_time
    logger.info(f"Stress test completed in {total_time:.2f} seconds")
    logger.info(f"Total queries: {results.total_queries}")
    logger.info(f"Successful queries: {results.successful_queries}")
    logger.info(f"Failed queries: {results.failed_queries}")
    logger.info(f"Connection errors: {results.connection_errors}")
    logger.info(f"Queries per second: {results.get_queries_per_second():.2f}")
    
    # Get query time statistics
    query_stats = results.get_query_time_stats()
    logger.info(f"Query time (avg): {query_stats['avg']:.6f}s, "
               f"(median): {query_stats['median']:.6f}s, "
               f"(p95): {query_stats['p95']:.6f}s")
    
    # Save results
    results_file = results.save_to_file()
    report_file = results.generate_report()
    
    logger.info(f"Results saved to {results_file}")
    logger.info(f"Report generated at {report_file}")
    
    return results.to_dict()


def run_full_test_suite():
    """
    Run a comprehensive test suite to evaluate connection pool performance.
    
    This runs multiple tests with different configurations to understand
    how the pool behaves under various conditions.
    """
    logger.info("Starting full connection pool test suite")
    
    # Define test scenarios
    scenarios = [
        {
            'name': 'Basic Load Test',
            'workers': 10,
            'duration': 15,
            'optimized': True,
            'min_connections': 2,
            'max_connections': 10,
            'simulate_errors': False,
            'use_transactions': True,
            'query_delay_range': (0.01, 0.2)
        },
        {
            'name': 'High Concurrency Test',
            'workers': 30,
            'duration': 15,
            'optimized': True,
            'min_connections': 5,
            'max_connections': 20,
            'simulate_errors': False,
            'use_transactions': True,
            'query_delay_range': (0.01, 0.2)
        },
        {
            'name': 'Error Resilience Test',
            'workers': 15,
            'duration': 15,
            'optimized': True,
            'min_connections': 3,
            'max_connections': 15,
            'simulate_errors': True,
            'use_transactions': True,
            'query_delay_range': (0.05, 0.5)
        },
        {
            'name': 'Comparison Test (Standard Pool)',
            'workers': 10,
            'duration': 15,
            'optimized': False,  # Force standard pool
            'min_connections': 2,
            'max_connections': 10,
            'simulate_errors': False,
            'use_transactions': True,
            'query_delay_range': (0.01, 0.2)
        }
    ]
    
    # Run each scenario
    results = {}
    for i, scenario in enumerate(scenarios):
        logger.info(f"Running test scenario {i+1}/{len(scenarios)}: {scenario['name']}")
        
        # Run the stress test
        test_result = stress_test_connection_pool(
            workers=scenario['workers'],
            duration=scenario['duration'],
            optimized=scenario['optimized'],
            min_connections=scenario['min_connections'],
            max_connections=scenario['max_connections'],
            simulate_errors=scenario['simulate_errors'],
            use_transactions=scenario['use_transactions'],
            query_delay_range=scenario['query_delay_range']
        )
        
        # Store results
        results[scenario['name']] = test_result
        
        # Wait a bit between tests
        if i < len(scenarios) - 1:
            logger.info("Waiting 5 seconds before next test scenario...")
            time.sleep(5)
    
    # Generate comparison report
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_file = f"reports/connection_pool_comparison_{timestamp}.md"
    
    with open(report_file, 'w') as f:
        f.write("# Connection Pool Test Suite Comparison Report\n\n")
        f.write(f"Test Date: {datetime.now().isoformat()}\n\n")
        
        # Comparative metrics table
        f.write("## Performance Comparison\n\n")
        f.write("| Scenario | QPS | Success Rate | Avg Query Time | p95 Query Time | Avg Conn Time | Errors |\n")
        f.write("|----------|-----|--------------|----------------|----------------|---------------|--------|\n")
        
        for name, result in results.items():
            qps = result['queries_per_second']
            success_rate = result['success_rate_percent']
            avg_query = result['query_time_stats']['avg']
            p95_query = result['query_time_stats']['p95']
            avg_conn = result['connection_time_stats']['avg']
            errors = result['failed_queries'] + result['connection_errors']
            
            f.write(f"| {name} | {qps:.2f} | {success_rate:.2f}% | {avg_query:.6f}s | "
                   f"{p95_query:.6f}s | {avg_conn:.6f}s | {errors} |\n")
                   
        f.write("\n")
        
        # Recommendations
        f.write("## Recommendations\n\n")
        
        # Compare optimized vs standard
        if 'Comparison Test (Standard Pool)' in results and 'Basic Load Test' in results:
            opt_result = results['Basic Load Test']
            std_result = results['Comparison Test (Standard Pool)']
            
            opt_qps = opt_result['queries_per_second']
            std_qps = std_result['queries_per_second']
            
            if opt_qps > std_qps * 1.2:
                f.write("- The optimized connection pool shows significantly better performance than the standard pool. ")
                f.write(f"It achieves {opt_qps:.2f} QPS vs {std_qps:.2f} QPS, a {((opt_qps/std_qps)-1)*100:.1f}% improvement.\n\n")
            elif opt_qps > std_qps:
                f.write("- The optimized connection pool performs slightly better than the standard pool. ")
                f.write(f"It achieves {opt_qps:.2f} QPS vs {std_qps:.2f} QPS, a {((opt_qps/std_qps)-1)*100:.1f}% improvement.\n\n")
            else:
                f.write("- The standard connection pool performs as well or better than the optimized pool in this test. ")
                f.write("Consider investigating why the optimized pool isn't providing the expected benefits.\n\n")
        
        # Pool sizing recommendation
        high_concurrency = results.get('High Concurrency Test', {})
        if high_concurrency:
            conn_time = high_concurrency['connection_time_stats']['avg']
            p95_conn_time = high_concurrency['connection_time_stats']['p95']
            pool_metrics = high_concurrency.get('latest_pool_metrics', {}).get('pool', {})
            utilization = pool_metrics.get('utilization', 0)
            
            if utilization > 0.85:
                f.write("- Under high concurrency, the connection pool shows high utilization (>85%). ")
                f.write(f"Consider increasing the maximum pool size above {pool_metrics.get('max_connections', 'current setting')}.\n\n")
            
            if p95_conn_time > conn_time * 3:
                f.write("- Connection acquisition time shows significant outliers under high concurrency. ")
                f.write(f"The p95 time ({p95_conn_time:.6f}s) is much higher than the average ({conn_time:.6f}s). ")
                f.write("This may indicate connection pool saturation or other issues.\n\n")
                
        # Error handling recommendation
        error_test = results.get('Error Resilience Test', {})
        if error_test:
            success_rate = error_test['success_rate_percent']
            if success_rate < 95:
                f.write("- The connection pool shows reduced resilience under error conditions, ")
                f.write(f"with a success rate of {success_rate:.2f}%. ")
                f.write("Consider enhancing error handling and retry mechanisms.\n\n")
            else:
                f.write("- The connection pool demonstrates good resilience under error conditions, ")
                f.write(f"maintaining a {success_rate:.2f}% success rate. ")
                f.write("The error handling and circuit breaker mechanisms are working effectively.\n\n")
                
        # Overall recommendation
        f.write("### Recommended Configuration\n\n")
        f.write("Based on the test results, the following connection pool configuration is recommended:\n\n")
        f.write("```python\n")
        f.write("connection_pool_config = {\n")
        
        # Determine recommended min_connections
        min_connections = 2
        if 'High Concurrency Test' in results:
            workers = scenarios[1]['workers']
            min_connections = max(min_connections, min(5, workers // 4))
            
        # Determine recommended max_connections
        max_connections = 10
        if 'High Concurrency Test' in results:
            workers = scenarios[1]['workers']
            high_concurrency = results['High Concurrency Test']
            pool_metrics = high_concurrency.get('latest_pool_metrics', {}).get('pool', {})
            utilization = pool_metrics.get('utilization', 0)
            
            if utilization > 0.7:
                # Need more connections
                max_connections = max(max_connections, min(50, workers))
            elif utilization < 0.3:
                # Could use fewer connections
                max_connections = max(min_connections * 2, workers // 2)
            else:
                # Current settings seem good
                max_connections = max(max_connections, min(30, workers))
                
        f.write(f"    'min_connections': {min_connections},\n")
        f.write(f"    'max_connections': {max_connections},\n")
        f.write("    'connection_timeout': 30,\n")
        f.write("    'connection_lifetime': 1800,  # 30 minutes\n")
        f.write("    'idle_timeout': 300,  # 5 minutes\n")
        f.write("    'health_check_interval': 60,  # 1 minute\n")
        
        # Add optimized pool settings if beneficial
        if 'Comparison Test (Standard Pool)' in results and 'Basic Load Test' in results:
            opt_result = results['Basic Load Test']
            std_result = results['Comparison Test (Standard Pool)']
            
            if opt_result['queries_per_second'] > std_result['queries_per_second']:
                f.write("    \n    # Optimized pool settings\n")
                f.write("    'retry_attempts': 3,\n")
                f.write("    'initial_retry_delay': 0.2,\n")
                f.write("    'max_retry_delay': 5,\n")
                f.write("    'retry_jitter_factor': 0.1,\n")
                f.write("    'validation_query': 'SELECT 1',\n")
                f.write("    'validation_timeout': 5,\n")
                f.write("    'circuit_breaker_threshold': 5,\n")
                f.write("    'circuit_breaker_timeout': 30,\n")
                f.write("    'circuit_breaker_reset': 2\n")
        
        f.write("}\n")
        f.write("```\n\n")
        
        # Next steps
        f.write("## Next Steps\n\n")
        f.write("1. Implement the recommended configuration changes\n")
        f.write("2. Run longer duration tests in production-like environment\n")
        f.write("3. Monitor pool metrics under real production load\n")
        f.write("4. Fine-tune parameters based on production insights\n")
    
    logger.info(f"Comparison report generated at {report_file}")
    return report_file
    

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Stress test connection pool performance")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS, 
                      help=f"Number of worker threads (default: {DEFAULT_WORKERS})")
    parser.add_argument("--duration", type=int, default=DEFAULT_DURATION,
                      help=f"Test duration in seconds (default: {DEFAULT_DURATION})")
    parser.add_argument("--standard", action="store_true",
                      help="Use standard connection pool instead of optimized (if available)")
    parser.add_argument("--min-connections", type=int, default=2,
                      help="Minimum connections in pool (default: 2)")
    parser.add_argument("--max-connections", type=int, default=DEFAULT_MAX_POOL_SIZE,
                      help=f"Maximum connections in pool (default: {DEFAULT_MAX_POOL_SIZE})")
    parser.add_argument("--simulate-errors", action="store_true",
                      help="Simulate random errors during testing")
    parser.add_argument("--no-transactions", action="store_true",
                      help="Don't use transactions for queries")
    parser.add_argument("--full-suite", action="store_true",
                      help="Run full test suite with multiple configurations")
    args = parser.parse_args()
    
    # Ensure logs directory exists
    os.makedirs('logs', exist_ok=True)
    os.makedirs('reports', exist_ok=True)
    
    # Environment variables check
    required_env = ['SUPABASE_DB_HOST', 'SUPABASE_DB_USER', 'SUPABASE_DB_PASSWORD']
    missing_env = [var for var in required_env if not os.environ.get(var)]
    
    if missing_env:
        logger.error(f"Missing required environment variables: {', '.join(missing_env)}")
        logger.error("Please set these environment variables or load them from .env file")
        
        # Try to load from .env file
        if os.path.exists('.env'):
            logger.info("Found .env file, attempting to load variables")
            with open('.env', 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    key, value = line.split('=', 1)
                    os.environ[key] = value
            
            # Check again
            missing_env = [var for var in required_env if not os.environ.get(var)]
            if missing_env:
                logger.error(f"Still missing required environment variables after loading .env: {', '.join(missing_env)}")
                sys.exit(1)
            else:
                logger.info("Successfully loaded environment variables from .env")
        else:
            sys.exit(1)
    
    # Run tests
    if args.full_suite:
        # Run full test suite
        report_file = run_full_test_suite()
        print(f"\nTest suite completed! Full report available at: {report_file}")
    else:
        # Run single test with specified parameters
        results = stress_test_connection_pool(
            workers=args.workers,
            duration=args.duration,
            optimized=not args.standard,
            min_connections=args.min_connections,
            max_connections=args.max_connections,
            simulate_errors=args.simulate_errors,
            use_transactions=not args.no_transactions,
            query_delay_range=(0.01, 0.3)
        )
        
        print("\nTest completed! Summary:")
        print(f"- Total queries: {results['total_queries']}")
        print(f"- Successful queries: {results['successful_queries']} ({results['success_rate_percent']:.2f}%)")
        print(f"- Queries per second: {results['queries_per_second']:.2f}")
        print(f"- Avg query time: {results['query_time_stats']['avg']:.6f}s")
        print(f"- P95 query time: {results['query_time_stats']['p95']:.6f}s")
        print(f"- Avg connection time: {results['connection_time_stats']['avg']:.6f}s")
        
        # Print report location
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        print(f"\nDetailed report available at: reports/connection_pool_stress_{timestamp}.md")


if __name__ == "__main__":
    main()