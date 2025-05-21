#!/usr/bin/env python3
"""
Test script for the optimized connection pool.

This script tests the OptimizedConnectionPool class with PostgreSQL connections.
"""

import os
import sys
import time
import random
import logging
import argparse
import threading
import psycopg2
from concurrent.futures import ThreadPoolExecutor
from psycopg2.extras import RealDictCursor
from typing import Dict, Any, List

# Import optimized connection pool
from optimized_connection_pool import (
    OptimizedConnectionPool, 
    retry_with_backoff,
    exponential_backoff
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class PostgresConnectionPool(OptimizedConnectionPool):
    """PostgreSQL implementation of the OptimizedConnectionPool."""
    
    def __init__(
        self,
        host: str,
        port: int,
        user: str,
        password: str,
        database: str,
        min_size: int = 2,
        max_size: int = 10,
        connection_timeout: int = 1800,
        validation_interval: int = 300,
        cleanup_interval: int = 600
    ):
        """
        Initialize PostgreSQL connection pool.
        
        Args:
            host: PostgreSQL host
            port: PostgreSQL port
            user: PostgreSQL user
            password: PostgreSQL password
            database: PostgreSQL database
            min_size: Minimum number of connections
            max_size: Maximum number of connections
            connection_timeout: Connection timeout in seconds
            validation_interval: Interval to validate idle connections
            cleanup_interval: Interval to cleanup expired connections
        """
        self.db_config = {
            'host': host,
            'port': port,
            'user': user,
            'password': password,
            'database': database
        }
        
        super().__init__(
            min_size=min_size,
            max_size=max_size,
            connection_timeout=connection_timeout,
            validation_interval=validation_interval,
            cleanup_interval=cleanup_interval,
            connection_factory=self._create_postgres_connection
        )
    
    def _create_postgres_connection(self) -> Dict[str, Any]:
        """Create a new PostgreSQL connection."""
        logger.debug(f"Creating new PostgreSQL connection to {self.db_config['host']}:{self.db_config['port']}")
        
        # Define a function to create the connection with retries
        def connect_with_retry():
            return psycopg2.connect(
                host=self.db_config['host'],
                port=self.db_config['port'],
                user=self.db_config['user'],
                password=self.db_config['password'],
                database=self.db_config['database'],
                cursor_factory=RealDictCursor
            )
        
        # Use retry_with_backoff for resilience
        conn = retry_with_backoff(
            connect_with_retry,
            max_retries=5,
            base_delay=0.5,
            max_delay=10.0,
            jitter=0.25,
            retry_on=(psycopg2.OperationalError, psycopg2.InterfaceError)
        )
        
        return conn
    
    def _validate_connection(self, connection: Dict[str, Any]) -> bool:
        """Validate PostgreSQL connection."""
        # First check the parent validation (timeout, errors)
        if not super()._validate_connection(connection):
            return False
        
        try:
            # Try to execute a simple query to check if the connection is alive
            conn = connection['data']
            cursor = conn.cursor()
            cursor.execute("SELECT 1")
            cursor.fetchone()
            cursor.close()
            return True
        except Exception as e:
            logger.warning(f"Connection validation failed: {str(e)}")
            return False

def run_test_query(pool, query, args=None):
    """Run a test query using the connection pool."""
    with pool.connection() as conn:
        cursor = conn.cursor()
        start_time = time.time()
        
        try:
            cursor.execute(query, args)
            result = cursor.fetchall()
            execution_time = time.time() - start_time
            
            return {
                'success': True,
                'rows': len(result),
                'execution_time': execution_time
            }
        except Exception as e:
            execution_time = time.time() - start_time
            
            return {
                'success': False,
                'error': str(e),
                'execution_time': execution_time
            }
        finally:
            cursor.close()

def worker(worker_id, pool, num_queries, queries):
    """Worker function for concurrent testing."""
    results = []
    
    for i in range(num_queries):
        query = random.choice(queries)
        logger.debug(f"Worker {worker_id} executing query {i+1}/{num_queries}: {query[:50]}...")
        
        try:
            result = run_test_query(pool, query)
            result['worker_id'] = worker_id
            result['query_index'] = i
            results.append(result)
            
            if result['success']:
                logger.debug(f"Worker {worker_id}, Query {i+1}: Success, {result['rows']} rows, {result['execution_time']:.4f}s")
            else:
                logger.warning(f"Worker {worker_id}, Query {i+1}: Failed, {result['error']}")
        except Exception as e:
            logger.error(f"Worker {worker_id}, Query {i+1}: Exception: {str(e)}")
            results.append({
                'worker_id': worker_id,
                'query_index': i,
                'success': False,
                'error': str(e),
                'execution_time': 0
            })
    
    return results

def run_load_test(pool, num_workers, queries_per_worker, test_queries):
    """Run load test with multiple workers."""
    logger.info(f"Starting load test with {num_workers} workers, {queries_per_worker} queries per worker")
    
    start_time = time.time()
    all_results = []
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        
        for i in range(num_workers):
            future = executor.submit(
                worker, i, pool, queries_per_worker, test_queries
            )
            futures.append(future)
        
        for future in futures:
            try:
                results = future.result()
                all_results.extend(results)
            except Exception as e:
                logger.error(f"Worker thread failed: {str(e)}")
    
    total_time = time.time() - start_time
    total_queries = len(all_results)
    successful_queries = sum(1 for r in all_results if r['success'])
    failed_queries = total_queries - successful_queries
    
    if total_queries > 0:
        success_rate = successful_queries / total_queries
        avg_execution_time = sum(r['execution_time'] for r in all_results) / total_queries
    else:
        success_rate = 0
        avg_execution_time = 0
    
    logger.info(f"Load test completed in {total_time:.2f}s")
    logger.info(f"Total queries: {total_queries}")
    logger.info(f"Successful queries: {successful_queries}")
    logger.info(f"Failed queries: {failed_queries}")
    logger.info(f"Success rate: {success_rate:.2%}")
    logger.info(f"Average execution time: {avg_execution_time:.4f}s")
    logger.info(f"Queries per second: {total_queries / total_time:.2f}")
    
    # Get and display pool stats
    stats = pool.get_stats()
    logger.info("Connection pool stats:")
    logger.info(f"  Created connections: {stats['created_connections']}")
    logger.info(f"  Reused connections: {stats['reused_connections']}")
    logger.info(f"  Discarded connections: {stats['discarded_connections']}")
    logger.info(f"  Idle connections: {stats['pool_state']['connections']['idle']}")
    logger.info(f"  In-use connections: {stats['pool_state']['connections']['in_use']}")
    logger.info(f"  Errors: {stats['errors']}")
    logger.info(f"  Average wait time: {stats['wait_time']['average']:.4f}s")
    logger.info(f"  Average query time: {stats['query_time']['average']:.4f}s")
    
    return {
        'total_time': total_time,
        'total_queries': total_queries,
        'successful_queries': successful_queries,
        'failed_queries': failed_queries,
        'success_rate': success_rate,
        'avg_execution_time': avg_execution_time,
        'queries_per_second': total_queries / total_time,
        'pool_stats': stats,
        'detailed_results': all_results
    }

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Test optimized connection pool')
    
    parser.add_argument('--host', default=os.environ.get('DB_HOST', 'localhost'),
                      help='PostgreSQL host (default: from DB_HOST env var or localhost)')
    parser.add_argument('--port', type=int, default=int(os.environ.get('DB_PORT', '5432')),
                      help='PostgreSQL port (default: from DB_PORT env var or 5432)')
    parser.add_argument('--user', default=os.environ.get('DB_USER', 'postgres'),
                      help='PostgreSQL user (default: from DB_USER env var or postgres)')
    parser.add_argument('--password', default=os.environ.get('DB_PASSWORD', ''),
                      help='PostgreSQL password (default: from DB_PASSWORD env var or empty)')
    parser.add_argument('--database', default=os.environ.get('DB_NAME', 'postgres'),
                      help='PostgreSQL database (default: from DB_NAME env var or postgres)')
    
    parser.add_argument('--min-connections', type=int, default=2,
                      help='Minimum number of connections in the pool (default: 2)')
    parser.add_argument('--max-connections', type=int, default=10,
                      help='Maximum number of connections in the pool (default: 10)')
    
    parser.add_argument('--workers', type=int, default=5,
                      help='Number of concurrent workers (default: 5)')
    parser.add_argument('--queries', type=int, default=10,
                      help='Number of queries per worker (default: 10)')
    
    parser.add_argument('--verbose', action='store_true',
                      help='Enable verbose logging')
    
    return parser.parse_args()

def main():
    """Main function."""
    args = parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Set up test queries
    test_queries = [
        "SELECT COUNT(*) FROM molecules",
        "SELECT COUNT(*) FROM mixtures",
        "SELECT COUNT(*) FROM molecular_properties",
        "SELECT * FROM molecules LIMIT 10",
        "SELECT * FROM mixtures LIMIT 10",
        "SELECT * FROM molecular_properties LIMIT 10",
        "SELECT m.id, m.name, COUNT(mp.id) FROM molecules m JOIN molecular_properties mp ON m.id = mp.molecule_id GROUP BY m.id, m.name LIMIT 10",
        "SELECT m.id, m.name, COUNT(mc.id) FROM mixtures m JOIN mixture_components mc ON m.id = mc.mixture_id GROUP BY m.id, m.name LIMIT 10"
    ]
    
    try:
        # Create connection pool
        pool = PostgresConnectionPool(
            host=args.host,
            port=args.port,
            user=args.user,
            password=args.password,
            database=args.database,
            min_size=args.min_connections,
            max_size=args.max_connections
        )
        
        logger.info(f"Created connection pool with min={args.min_connections}, max={args.max_connections}")
        
        # Run load test
        results = run_load_test(
            pool=pool,
            num_workers=args.workers,
            queries_per_worker=args.queries,
            test_queries=test_queries
        )
        
        # Shutdown pool
        logger.info("Shutting down connection pool")
        pool.shutdown()
        
        return 0
    
    except Exception as e:
        logger.error(f"Error in main: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())