#!/usr/bin/env python3
"""
CryoProtect v2 - Test Connection Pooling

This script tests the connection pool under load conditions.
It simulates concurrent requests to verify the connection pool works correctly.
"""

import os
import sys
import time
import json
import logging
import argparse
import threading
import random
import statistics
from datetime import datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("test_connection_pooling.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("test_connection_pooling")

# Import connection pooling
try:
    from implement_connection_pooling import (
        initialize_connection_pool, execute_query,
        get_connection_pool_stats, shutdown_connection_pool
    )
except ImportError:
    logger.error("Error: implement_connection_pooling.py not found. Make sure it's in the same directory.")
    sys.exit(1)

# Import Supabase MCP tools
try:
    from supabase_mcp_tools import execute_sql_on_supabase
except ImportError:
    logger.error("Error: supabase_mcp_tools.py not found. Make sure it's in the same directory.")
    sys.exit(1)

def run_query(query, project_id=None, use_pool=True):
    """
    Run a query using either the connection pool or direct execution.
    
    Args:
        query: The SQL query to execute
        project_id: The Supabase project ID (required if use_pool=False)
        use_pool: Whether to use the connection pool
        
    Returns:
        Tuple of (result, execution_time)
    """
    start_time = time.time()
    
    try:
        if use_pool:
            result = execute_query(query)
        else:
            result = execute_sql_on_supabase(project_id, query)
        
        execution_time = time.time() - start_time
        return result, execution_time
    
    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"Error executing query: {str(e)}")
        return None, execution_time

def worker(worker_id, num_queries, query, project_id, use_pool):
    """
    Worker function to execute queries.
    
    Args:
        worker_id: ID of the worker
        num_queries: Number of queries to execute
        query: The SQL query to execute
        project_id: The Supabase project ID
        use_pool: Whether to use the connection pool
        
    Returns:
        List of execution times
    """
    execution_times = []
    errors = 0
    
    for i in range(num_queries):
        try:
            _, execution_time = run_query(query, project_id, use_pool)
            execution_times.append(execution_time)
            logger.debug(f"Worker {worker_id}, Query {i}: {execution_time:.4f}s")
        except Exception as e:
            errors += 1
            logger.error(f"Worker {worker_id}, Query {i} failed: {str(e)}")
    
    return {
        "worker_id": worker_id,
        "execution_times": execution_times,
        "errors": errors,
        "avg_time": statistics.mean(execution_times) if execution_times else 0,
        "min_time": min(execution_times) if execution_times else 0,
        "max_time": max(execution_times) if execution_times else 0
    }

def run_load_test(project_id, num_workers, queries_per_worker, use_pool=True):
    """
    Run a load test with multiple workers.
    
    Args:
        project_id: The Supabase project ID
        num_workers: Number of concurrent workers
        queries_per_worker: Number of queries per worker
        use_pool: Whether to use the connection pool
        
    Returns:
        Dictionary with test results
    """
    # Test query
    query = "SELECT COUNT(*) FROM public.molecules;"
    
    # Initialize connection pool if using it
    if use_pool:
        initialize_connection_pool(
            project_id=project_id,
            min_connections=2,
            max_connections=10,
            connection_timeout=30
        )
    
    # Run the load test
    start_time = time.time()
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i in range(num_workers):
            future = executor.submit(
                worker, i, queries_per_worker, query, project_id, use_pool
            )
            futures.append(future)
        
        # Wait for all workers to complete
        worker_results = [future.result() for future in futures]
    
    end_time = time.time()
    
    # Calculate statistics
    all_times = []
    total_errors = 0
    
    for result in worker_results:
        all_times.extend(result["execution_times"])
        total_errors += result["errors"]
    
    # Get pool stats if using the pool
    pool_stats = get_connection_pool_stats() if use_pool else None
    
    # Shutdown the pool if using it
    if use_pool:
        shutdown_connection_pool()
    
    # Calculate overall statistics
    total_queries = num_workers * queries_per_worker
    successful_queries = total_queries - total_errors
    
    stats = {
        "timestamp": datetime.now().isoformat(),
        "use_pool": use_pool,
        "num_workers": num_workers,
        "queries_per_worker": queries_per_worker,
        "total_queries": total_queries,
        "successful_queries": successful_queries,
        "total_time": end_time - start_time,
        "queries_per_second": successful_queries / (end_time - start_time),
        "success_rate": successful_queries / total_queries if total_queries > 0 else 0,
        "avg_time": statistics.mean(all_times) if all_times else 0,
        "min_time": min(all_times) if all_times else 0,
        "max_time": max(all_times) if all_times else 0,
        "median_time": statistics.median(all_times) if all_times else 0,
        "stdev_time": statistics.stdev(all_times) if len(all_times) > 1 else 0,
        "total_errors": total_errors,
        "worker_results": worker_results,
        "pool_stats": pool_stats
    }
    
    return stats

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Test connection pooling for CryoProtect database")
    parser.add_argument("--project-id", default="tsdlmynydfuypiugmkev",
                        help="Supabase project ID (default: tsdlmynydfuypiugmkev)")
    parser.add_argument("--workers", type=int, default=5,
                        help="Number of concurrent workers (default: 5)")
    parser.add_argument("--queries", type=int, default=20,
                        help="Number of queries per worker (default: 20)")
    parser.add_argument("--compare", action="store_true",
                        help="Compare performance with and without connection pooling")
    parser.add_argument("--output", default="connection_pool_test_results.json",
                        help="Output file for test results (default: connection_pool_test_results.json)")
    return parser.parse_args()

def main():
    """Main function to test connection pooling."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        logger.info(f"Testing connection pooling for project {args.project_id}")
        logger.info(f"Workers: {args.workers}, Queries per worker: {args.queries}")
        
        results = {}
        
        # Run test with connection pooling
        logger.info("Running test with connection pooling...")
        results["with_pool"] = run_load_test(
            project_id=args.project_id,
            num_workers=args.workers,
            queries_per_worker=args.queries,
            use_pool=True
        )
        
        logger.info(f"Test with pool: {results['with_pool']['queries_per_second']:.2f} queries/sec, " +
                   f"{results['with_pool']['success_rate'] * 100:.2f}% success rate")
        
        # Run test without connection pooling if requested
        if args.compare:
            logger.info("Running test without connection pooling...")
            results["without_pool"] = run_load_test(
                project_id=args.project_id,
                num_workers=args.workers,
                queries_per_worker=args.queries,
                use_pool=False
            )
            
            logger.info(f"Test without pool: {results['without_pool']['queries_per_second']:.2f} queries/sec, " +
                       f"{results['without_pool']['success_rate'] * 100:.2f}% success rate")
            
            # Calculate improvement
            with_pool_qps = results["with_pool"]["queries_per_second"]
            without_pool_qps = results["without_pool"]["queries_per_second"]
            
            if without_pool_qps > 0:
                improvement = (with_pool_qps - without_pool_qps) / without_pool_qps * 100
                logger.info(f"Performance improvement: {improvement:.2f}%")
                
                results["comparison"] = {
                    "improvement_percent": improvement,
                    "with_pool_qps": with_pool_qps,
                    "without_pool_qps": without_pool_qps
                }
        
        # Save results to file
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"Test results saved to {args.output}")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error testing connection pooling: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
