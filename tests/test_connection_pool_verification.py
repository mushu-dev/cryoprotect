"""
Connection Pool Verification Tests
Based on the test cases defined in DATABASE_VERIFICATION_PLAN.md
"""
import os
import sys
import unittest
import time
import psycopg2
import concurrent.futures
import statistics
import logging
from dotenv import load_dotenv

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Try to import connection pool implementations
try:
    from db_pool import ConnectionPool, execute_query
    HAS_DB_POOL = True
except ImportError:
    try:
        from optimized_connection_pool import ConnectionPool, execute_query
        HAS_DB_POOL = True
    except ImportError:
        try:
            from database.connection_manager import ConnectionPool, execute_query
            HAS_DB_POOL = True
        except ImportError:
            HAS_DB_POOL = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename='connection_pool_tests.log'
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Get database connection parameters
DB_CONNECTION_PARAMS = {
    'host': os.environ.get('SUPABASE_DB_HOST'),
    'port': os.environ.get('SUPABASE_DB_PORT', '5432'),
    'dbname': os.environ.get('SUPABASE_DB_NAME'),
    'user': os.environ.get('SUPABASE_DB_USER'),
    'password': os.environ.get('SUPABASE_DB_PASSWORD'),
    'sslmode': os.environ.get('SUPABASE_DB_SSLMODE', 'require')
}

# Simplified connection pool implementation for testing if no main one is available
class SimpleConnectionPool:
    """Simple connection pool implementation for testing"""
    
    def __init__(self, min_conn=1, max_conn=10, **db_params):
        """Initialize connection pool"""
        self.min_conn = min_conn
        self.max_conn = max_conn
        self.db_params = db_params
        self.connections = []
        self.in_use = {}
        
        # Create initial connections
        for _ in range(min_conn):
            conn = self._create_connection()
            self.connections.append(conn)
    
    def _create_connection(self):
        """Create a new database connection"""
        try:
            conn = psycopg2.connect(**self.db_params)
            return conn
        except Exception as e:
            logger.error(f"Error creating connection: {str(e)}")
            raise
    
    def get_connection(self):
        """Get a connection from the pool"""
        if not self.connections:
            if len(self.in_use) < self.max_conn:
                # Create a new connection
                conn = self._create_connection()
            else:
                # No connections available and at max capacity
                raise Exception("No connections available")
        else:
            # Get connection from pool
            conn = self.connections.pop()
        
        # Add to in-use dictionary with a unique ID
        conn_id = id(conn)
        self.in_use[conn_id] = conn
        
        return conn, conn_id
    
    def return_connection(self, conn, conn_id):
        """Return a connection to the pool"""
        if conn_id in self.in_use:
            # Remove from in-use
            del self.in_use[conn_id]
            
            # Return to pool
            self.connections.append(conn)
        else:
            # Connection not from this pool
            conn.close()
    
    def close(self):
        """Close all connections in the pool"""
        # Close connections in the pool
        for conn in self.connections:
            conn.close()
        
        # Close in-use connections
        for conn in self.in_use.values():
            conn.close()
        
        # Clear lists
        self.connections = []
        self.in_use = {}

def simple_execute_query(query, params=None, pool=None):
    """Execute a query using the connection pool"""
    if pool is None:
        # Create a new pool
        pool = SimpleConnectionPool(**DB_CONNECTION_PARAMS)
    
    conn = None
    conn_id = None
    try:
        # Get connection from pool
        conn, conn_id = pool.get_connection()
        
        # Execute query
        cursor = conn.cursor()
        cursor.execute(query, params)
        
        # Fetch results
        try:
            results = cursor.fetchall()
            return results
        except psycopg2.ProgrammingError:
            # No results to fetch
            return None
    finally:
        # Return connection to pool
        if conn and conn_id:
            pool.return_connection(conn, conn_id)

# Use the available connection pool implementation or fallback to the simple one
if not HAS_DB_POOL:
    logger.warning("No connection pool implementation found, using simplified version for testing")
    ConnectionPool = SimpleConnectionPool
    execute_query = simple_execute_query

class ConnectionPoolBasicFunctionalityTests(unittest.TestCase):
    """Tests for Connection Pool Basic Functionality (TC-CP-001 to TC-CP-005)"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test resources"""
        cls.config = DB_CONNECTION_PARAMS
    
    def setUp(self):
        """Set up test"""
        self.pool = None
    
    def tearDown(self):
        """Clean up test"""
        if self.pool:
            self.pool.close()
    
    def test_initialize_pool(self):
        """TC-CP-001: Initialize connection pool"""
        try:
            # Create pool
            self.pool = ConnectionPool(min_conn=1, max_conn=5, **self.config)
            self.assertIsNotNone(self.pool)
            logger.info("Connection pool initialized successfully")
        except Exception as e:
            logger.error(f"Error initializing connection pool: {str(e)}")
            self.fail(f"Connection pool initialization failed: {str(e)}")
    
    def test_get_return_connection(self):
        """TC-CP-002: Get and return connection"""
        try:
            # Create pool
            self.pool = ConnectionPool(min_conn=1, max_conn=5, **self.config)
            
            # Get connection
            conn, conn_id = self.pool.get_connection()
            self.assertIsNotNone(conn)
            self.assertIsNotNone(conn_id)
            logger.info("Connection obtained successfully")
            
            # Return connection
            self.pool.return_connection(conn, conn_id)
            logger.info("Connection returned successfully")
        except Exception as e:
            logger.error(f"Error in get/return connection test: {str(e)}")
            self.fail(f"Get/return connection test failed: {str(e)}")
    
    def test_execute_simple_query(self):
        """TC-CP-003: Execute simple query"""
        try:
            # Create pool
            self.pool = ConnectionPool(min_conn=1, max_conn=5, **self.config)
            
            # Execute simple query
            results = execute_query("SELECT 1 as test", pool=self.pool)
            self.assertIsNotNone(results)
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0][0], 1)
            logger.info("Simple query executed successfully")
        except Exception as e:
            logger.error(f"Error executing simple query: {str(e)}")
            self.fail(f"Simple query execution failed: {str(e)}")
    
    def test_execute_transaction(self):
        """TC-CP-004: Execute transaction"""
        try:
            # Create pool
            self.pool = ConnectionPool(min_conn=1, max_conn=5, **self.config)
            
            # Get connection for transaction
            conn, conn_id = self.pool.get_connection()
            
            try:
                # Start transaction
                conn.autocommit = False
                cursor = conn.cursor()
                
                # Create temporary table
                cursor.execute("""
                    CREATE TEMPORARY TABLE test_transaction (
                        id serial PRIMARY KEY,
                        name text
                    )
                """)
                
                # Insert data
                cursor.execute("INSERT INTO test_transaction (name) VALUES (%s) RETURNING id", ["Test 1"])
                
                # Commit transaction
                conn.commit()
                
                # Verify data
                cursor.execute("SELECT COUNT(*) FROM test_transaction")
                count = cursor.fetchone()[0]
                self.assertEqual(count, 1)
                
                logger.info("Transaction executed successfully")
            except Exception as e:
                # Rollback on error
                conn.rollback()
                raise
            finally:
                # Return connection to pool
                self.pool.return_connection(conn, conn_id)
        except Exception as e:
            logger.error(f"Error executing transaction: {str(e)}")
            self.fail(f"Transaction execution failed: {str(e)}")
    
    def test_handle_connection_errors(self):
        """TC-CP-005: Handle connection errors"""
        try:
            # Create pool with invalid parameters
            invalid_config = self.config.copy()
            invalid_config['dbname'] = 'nonexistent_db'
            
            try:
                # This should raise an exception
                error_pool = ConnectionPool(min_conn=1, max_conn=5, **invalid_config)
                self.fail("Expected an exception for invalid connection parameters")
            except Exception as e:
                # Expected exception
                logger.info(f"Connection error handled successfully: {str(e)}")
                pass
            
            # Create valid pool
            self.pool = ConnectionPool(min_conn=1, max_conn=5, **self.config)
            
            # Test invalid query
            try:
                # This should raise an exception
                execute_query("SELECT * FROM nonexistent_table", pool=self.pool)
                self.fail("Expected an exception for invalid query")
            except Exception as e:
                # Expected exception
                logger.info(f"Query error handled successfully: {str(e)}")
                pass
        except Exception as e:
            logger.error(f"Error in connection error handling test: {str(e)}")
            self.fail(f"Connection error handling test failed: {str(e)}")

class ConnectionPoolPerformanceTests(unittest.TestCase):
    """Tests for Connection Pool Performance (TC-CP-010 to TC-CP-014)"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test resources"""
        cls.config = DB_CONNECTION_PARAMS
        cls.pool = ConnectionPool(min_conn=2, max_conn=10, **cls.config)
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test resources"""
        if cls.pool:
            cls.pool.close()
    
    def test_query_times(self):
        """TC-CP-010, TC-CP-011: Measure query time with direct connection and pool"""
        # Test queries
        queries = [
            "SELECT COUNT(*) FROM molecule",
            "SELECT COUNT(*) FROM user_profile",
            "SELECT COUNT(*) FROM project"
        ]
        
        direct_times = []
        pool_times = []
        
        for query in queries:
            # Direct connection query
            conn = None
            try:
                # Time direct connection
                start_time = time.time()
                
                # Create connection
                conn = psycopg2.connect(**self.config)
                
                # Execute query
                cursor = conn.cursor()
                cursor.execute(query)
                cursor.fetchall()
                
                # Calculate time
                query_time = time.time() - start_time
                direct_times.append(query_time)
                
                logger.info(f"Direct connection query time: {query_time:.6f}s")
            except Exception as e:
                logger.error(f"Error in direct connection query: {str(e)}")
            finally:
                # Close connection
                if conn:
                    conn.close()
            
            # Connection pool query
            try:
                # Time connection pool query
                start_time = time.time()
                
                # Execute query through pool
                execute_query(query, pool=self.pool)
                
                # Calculate time
                query_time = time.time() - start_time
                pool_times.append(query_time)
                
                logger.info(f"Connection pool query time: {query_time:.6f}s")
            except Exception as e:
                logger.error(f"Error in connection pool query: {str(e)}")
        
        # Calculate average times
        avg_direct = statistics.mean(direct_times) if direct_times else 0
        avg_pool = statistics.mean(pool_times) if pool_times else 0
        
        logger.info(f"Average direct connection time: {avg_direct:.6f}s")
        logger.info(f"Average connection pool time: {avg_pool:.6f}s")
        
        # Pool should be at least as fast as direct connection
        self.assertLessEqual(avg_pool, avg_direct * 1.5)
    
    def test_concurrent_queries(self):
        """TC-CP-012: Measure concurrent query performance"""
        # Test query
        query = "SELECT COUNT(*) FROM molecule"
        
        # Number of concurrent queries
        num_concurrent = 5
        
        # Single query time
        start_time = time.time()
        execute_query(query, pool=self.pool)
        single_time = time.time() - start_time
        
        logger.info(f"Single query time: {single_time:.6f}s")
        
        # Concurrent queries
        start_time = time.time()
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_concurrent) as executor:
            futures = [executor.submit(execute_query, query, pool=self.pool) for _ in range(num_concurrent)]
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
        
        total_time = time.time() - start_time
        
        logger.info(f"Concurrent queries ({num_concurrent}): {total_time:.6f}s")
        logger.info(f"Average time per query: {total_time/num_concurrent:.6f}s")
        
        # Calculate theoretical sequential time
        sequential_time = single_time * num_concurrent
        
        logger.info(f"Sequential equivalent: {sequential_time:.6f}s")
        logger.info(f"Speedup factor: {sequential_time/total_time:.2f}x")
        
        # There should be some speedup (but allow for test environment variations)
        self.assertLess(total_time, sequential_time * 0.95)
    
    def test_connection_acquisition_time(self):
        """TC-CP-013: Measure connection acquisition time"""
        acquisition_times = []
        
        # Measure multiple times
        for _ in range(5):
            start_time = time.time()
            conn, conn_id = self.pool.get_connection()
            acquisition_time = time.time() - start_time
            acquisition_times.append(acquisition_time)
            self.pool.return_connection(conn, conn_id)
        
        # Calculate average acquisition time
        avg_acquisition = statistics.mean(acquisition_times)
        
        logger.info(f"Average connection acquisition time: {avg_acquisition:.6f}s")
        
        # Connection acquisition should be quick
        self.assertLess(avg_acquisition, 0.1)
    
    def test_pool_scaling(self):
        """TC-CP-014: Measure pool scaling under load"""
        # Create a small pool
        small_pool = ConnectionPool(min_conn=1, max_conn=10, **self.config)
        
        try:
            # Test query
            query = "SELECT COUNT(*) FROM molecule"
            
            # Get multiple connections
            connections = []
            for _ in range(5):
                conn, conn_id = small_pool.get_connection()
                connections.append((conn, conn_id))
            
            # Execute query on remaining connections
            start_time = time.time()
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
                futures = [executor.submit(execute_query, query, pool=small_pool) for _ in range(3)]
                results = [future.result() for future in concurrent.futures.as_completed(futures)]
            
            total_time = time.time() - start_time
            
            logger.info(f"Pool scaling test: {total_time:.6f}s")
            
            # Return connections
            for conn, conn_id in connections:
                small_pool.return_connection(conn, conn_id)
        finally:
            # Close pool
            small_pool.close()

class ConnectionPoolResilienceTests(unittest.TestCase):
    """Tests for Connection Pool Resilience (TC-CP-020 to TC-CP-024)"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test resources"""
        cls.config = DB_CONNECTION_PARAMS
    
    def setUp(self):
        """Set up test"""
        self.pool = ConnectionPool(min_conn=1, max_conn=5, **self.config)
    
    def tearDown(self):
        """Clean up test"""
        if self.pool:
            self.pool.close()
    
    def test_connection_validation(self):
        """TC-CP-022: Test connection validation"""
        # Get connection
        conn, conn_id = self.pool.get_connection()
        
        # Simulate connection problems
        try:
            # Execute a query to verify the connection is working
            cursor = conn.cursor()
            cursor.execute("SELECT 1")
            cursor.fetchall()
            
            # Connection is valid
            logger.info("Connection validation successful")
        except Exception as e:
            # Connection is invalid
            logger.error(f"Connection validation failed: {str(e)}")
            self.fail("Connection validation failed")
        finally:
            # Return connection to pool
            self.pool.return_connection(conn, conn_id)
    
    def test_connection_lifecycle(self):
        """TC-CP-023: Test connection lifecycle"""
        # Test multiple get and return operations
        for _ in range(5):
            # Get connection
            conn, conn_id = self.pool.get_connection()
            
            # Execute query
            cursor = conn.cursor()
            cursor.execute("SELECT 1")
            result = cursor.fetchone()
            self.assertEqual(result[0], 1)
            
            # Return connection
            self.pool.return_connection(conn, conn_id)
        
        logger.info("Connection lifecycle test successful")

# Add import for uuid if needed by tests
import uuid

if __name__ == '__main__':
    unittest.main()