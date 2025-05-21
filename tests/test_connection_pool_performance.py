"""
Connection Pool Performance Tests
Focused specifically on measuring performance metrics and benchmarking the connection pool
"""
import os
import sys
import time
import statistics
import json
from datetime import datetime
import concurrent.futures
import logging
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
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
    filename='connection_pool_performance.log'
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

class ConnectionPoolPerformanceTester:
    """Class to test connection pool performance"""
    
    def __init__(self, db_params, min_conn=2, max_conn=10):
        """Initialize tester with database parameters and pool settings"""
        self.db_params = db_params
        self.min_conn = min_conn
        self.max_conn = max_conn
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'pool_settings': {
                'min_conn': min_conn,
                'max_conn': max_conn
            },
            'direct_connection': {
                'single_query': [],
                'multiple_queries': []
            },
            'connection_pool': {
                'single_query': [],
                'multiple_queries': [],
                'concurrent_queries': []
            },
            'acquisition_times': [],
            'pool_scaling': []
        }
    
    def test_direct_connection(self, query, iterations=5):
        """Test query performance with direct connection"""
        logger.info(f"Testing direct connection performance with {iterations} iterations")
        print(f"Testing direct connection performance...")
        
        # Test single query
        for i in range(iterations):
            conn = None
            try:
                start_time = time.time()
                
                # Create connection
                conn = psycopg2.connect(**self.db_params)
                
                # Execute query
                cursor = conn.cursor()
                cursor.execute(query)
                results = cursor.fetchall()
                
                query_time = time.time() - start_time
                self.results['direct_connection']['single_query'].append(query_time)
                
                logger.info(f"Direct connection query {i+1}: {query_time:.6f}s")
                print(f"  Query {i+1}: {query_time:.6f}s")
            except Exception as e:
                logger.error(f"Error in direct connection test: {str(e)}")
                print(f"  Error: {str(e)}")
            finally:
                if conn:
                    conn.close()
        
        # Calculate statistics
        self._add_statistics(self.results['direct_connection'], 'single_query')
        
        # Test multiple sequential queries
        num_queries = 10
        
        for i in range(iterations):
            conn = None
            try:
                start_time = time.time()
                
                # Create connection
                conn = psycopg2.connect(**self.db_params)
                cursor = conn.cursor()
                
                # Execute multiple queries sequentially
                for _ in range(num_queries):
                    cursor.execute(query)
                    results = cursor.fetchall()
                
                query_time = time.time() - start_time
                self.results['direct_connection']['multiple_queries'].append(query_time)
                
                logger.info(f"Direct connection multiple queries {i+1}: {query_time:.6f}s")
                print(f"  Multiple queries {i+1}: {query_time:.6f}s")
            except Exception as e:
                logger.error(f"Error in direct connection multiple queries test: {str(e)}")
                print(f"  Error: {str(e)}")
            finally:
                if conn:
                    conn.close()
        
        # Calculate statistics
        self._add_statistics(self.results['direct_connection'], 'multiple_queries')
    
    def test_connection_pool(self, query, iterations=5, concurrent_queries=10):
        """Test query performance with connection pool"""
        logger.info(f"Testing connection pool performance with {iterations} iterations")
        print(f"Testing connection pool performance...")
        
        # Create connection pool
        pool = ConnectionPool(min_conn=self.min_conn, max_conn=self.max_conn, **self.db_params)
        
        try:
            # Test single query
            for i in range(iterations):
                try:
                    start_time = time.time()
                    
                    # Execute query
                    execute_query(query, pool=pool)
                    
                    query_time = time.time() - start_time
                    self.results['connection_pool']['single_query'].append(query_time)
                    
                    logger.info(f"Connection pool query {i+1}: {query_time:.6f}s")
                    print(f"  Query {i+1}: {query_time:.6f}s")
                except Exception as e:
                    logger.error(f"Error in connection pool test: {str(e)}")
                    print(f"  Error: {str(e)}")
            
            # Calculate statistics
            self._add_statistics(self.results['connection_pool'], 'single_query')
            
            # Test multiple sequential queries
            num_queries = 10
            
            for i in range(iterations):
                try:
                    start_time = time.time()
                    
                    # Execute multiple queries sequentially
                    for _ in range(num_queries):
                        execute_query(query, pool=pool)
                    
                    query_time = time.time() - start_time
                    self.results['connection_pool']['multiple_queries'].append(query_time)
                    
                    logger.info(f"Connection pool multiple queries {i+1}: {query_time:.6f}s")
                    print(f"  Multiple queries {i+1}: {query_time:.6f}s")
                except Exception as e:
                    logger.error(f"Error in connection pool multiple queries test: {str(e)}")
                    print(f"  Error: {str(e)}")
            
            # Calculate statistics
            self._add_statistics(self.results['connection_pool'], 'multiple_queries')
            
            # Test concurrent queries
            for i in range(iterations):
                try:
                    start_time = time.time()
                    
                    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrent_queries) as executor:
                        futures = [executor.submit(execute_query, query, pool=pool) for _ in range(concurrent_queries)]
                        results = [future.result() for future in concurrent.futures.as_completed(futures)]
                    
                    query_time = time.time() - start_time
                    self.results['connection_pool']['concurrent_queries'].append(query_time)
                    
                    logger.info(f"Connection pool concurrent queries {i+1}: {query_time:.6f}s")
                    print(f"  Concurrent queries {i+1}: {query_time:.6f}s")
                except Exception as e:
                    logger.error(f"Error in connection pool concurrent queries test: {str(e)}")
                    print(f"  Error: {str(e)}")
            
            # Calculate statistics
            self._add_statistics(self.results['connection_pool'], 'concurrent_queries')
            
            # Calculate theoretical vs actual concurrent performance
            self._add_concurrent_performance_metrics()
            
            # Test connection acquisition time
            self._test_connection_acquisition(pool, iterations)
            
            # Test pool scaling
            self._test_pool_scaling(pool, iterations)
        
        finally:
            # Close pool
            pool.close()
    
    def _test_connection_acquisition(self, pool, iterations=5):
        """Test connection acquisition time"""
        logger.info(f"Testing connection acquisition time with {iterations} iterations")
        print(f"Testing connection acquisition time...")
        
        for i in range(iterations):
            try:
                start_time = time.time()
                conn, conn_id = pool.get_connection()
                acquisition_time = time.time() - start_time
                self.results['acquisition_times'].append(acquisition_time)
                
                logger.info(f"Connection acquisition {i+1}: {acquisition_time:.6f}s")
                print(f"  Acquisition {i+1}: {acquisition_time:.6f}s")
                
                # Return connection
                pool.return_connection(conn, conn_id)
            except Exception as e:
                logger.error(f"Error in connection acquisition test: {str(e)}")
                print(f"  Error: {str(e)}")
        
        # Calculate statistics
        self.results['acquisition_times_stats'] = {
            'mean': statistics.mean(self.results['acquisition_times']),
            'median': statistics.median(self.results['acquisition_times']),
            'min': min(self.results['acquisition_times']),
            'max': max(self.results['acquisition_times'])
        }
        
        logger.info(f"Connection acquisition stats: {self.results['acquisition_times_stats']}")
        print(f"  Average acquisition time: {self.results['acquisition_times_stats']['mean']:.6f}s")
    
    def _test_pool_scaling(self, pool, iterations=3):
        """Test pool scaling under load"""
        logger.info(f"Testing pool scaling with {iterations} iterations")
        print(f"Testing pool scaling under load...")
        
        for i in range(iterations):
            try:
                # Get multiple connections to force pool scaling
                connections = []
                start_time = time.time()
                
                # Get connections up to max_conn
                for _ in range(self.max_conn):
                    try:
                        conn, conn_id = pool.get_connection()
                        connections.append((conn, conn_id))
                    except Exception as e:
                        # Expected if pool reaches max_conn
                        logger.info(f"Pool reached max connections: {str(e)}")
                        break
                
                scaling_time = time.time() - start_time
                self.results['pool_scaling'].append({
                    'connections_obtained': len(connections),
                    'time': scaling_time
                })
                
                logger.info(f"Pool scaling {i+1}: {len(connections)} connections in {scaling_time:.6f}s")
                print(f"  Obtained {len(connections)} connections in {scaling_time:.6f}s")
                
                # Return connections
                for conn, conn_id in connections:
                    pool.return_connection(conn, conn_id)
            except Exception as e:
                logger.error(f"Error in pool scaling test: {str(e)}")
                print(f"  Error: {str(e)}")
    
    def _add_statistics(self, results_dict, key):
        """Add statistical analysis to results"""
        if len(results_dict[key]) > 0:
            results_dict[f"{key}_stats"] = {
                'mean': statistics.mean(results_dict[key]),
                'median': statistics.median(results_dict[key]),
                'min': min(results_dict[key]),
                'max': max(results_dict[key])
            }
    
    def _add_concurrent_performance_metrics(self):
        """Calculate and add concurrent performance metrics"""
        if (len(self.results['connection_pool']['single_query']) > 0 and 
            len(self.results['connection_pool']['concurrent_queries']) > 0):
            
            # Calculate average times
            avg_single = statistics.mean(self.results['connection_pool']['single_query'])
            avg_concurrent = statistics.mean(self.results['connection_pool']['concurrent_queries'])
            
            # Calculate theoretical sequential time for the same number of queries
            concurrent_count = 10  # default number of concurrent queries
            theoretical_sequential = avg_single * concurrent_count
            
            # Calculate speedup
            speedup = theoretical_sequential / avg_concurrent if avg_concurrent > 0 else 0
            
            # Add to results
            self.results['concurrent_performance'] = {
                'theoretical_sequential': theoretical_sequential,
                'actual_concurrent': avg_concurrent,
                'speedup': speedup
            }
            
            logger.info(f"Concurrent performance: {self.results['concurrent_performance']}")
            print(f"  Sequential equivalent: {theoretical_sequential:.6f}s")
            print(f"  Actual concurrent: {avg_concurrent:.6f}s")
            print(f"  Speedup: {speedup:.2f}x")
    
    def generate_report(self, output_file=None):
        """Generate a performance report"""
        if output_file is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = f"reports/connection_pool_performance_{timestamp}.json"
        
        # Ensure reports directory exists
        os.makedirs("reports", exist_ok=True)
        
        # Add summary to results
        self._add_summary()
        
        # Write results to file
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        logger.info(f"Performance report generated: {output_file}")
        print(f"\nPerformance report generated: {output_file}")
        
        # Also generate Markdown report
        markdown_file = output_file.replace('.json', '.md')
        self._generate_markdown_report(markdown_file)
        
        logger.info(f"Markdown report generated: {markdown_file}")
        print(f"Markdown report generated: {markdown_file}")
        
        return output_file
    
    def _add_summary(self):
        """Add summary section to results"""
        summary = {
            'direct_connection': {
                'single_query_avg': self.results['direct_connection'].get('single_query_stats', {}).get('mean', 0),
                'multiple_queries_avg': self.results['direct_connection'].get('multiple_queries_stats', {}).get('mean', 0),
            },
            'connection_pool': {
                'single_query_avg': self.results['connection_pool'].get('single_query_stats', {}).get('mean', 0),
                'multiple_queries_avg': self.results['connection_pool'].get('multiple_queries_stats', {}).get('mean', 0),
                'concurrent_queries_avg': self.results['connection_pool'].get('concurrent_queries_stats', {}).get('mean', 0),
            },
            'acquisition_time_avg': self.results.get('acquisition_times_stats', {}).get('mean', 0),
            'speedup': self.results.get('concurrent_performance', {}).get('speedup', 0),
        }
        
        # Calculate performance difference percentages
        if summary['direct_connection']['single_query_avg'] > 0:
            single_query_diff = (
                (summary['direct_connection']['single_query_avg'] - summary['connection_pool']['single_query_avg']) / 
                summary['direct_connection']['single_query_avg'] * 100
            )
            summary['single_query_improvement'] = single_query_diff
        
        if summary['direct_connection']['multiple_queries_avg'] > 0:
            multiple_queries_diff = (
                (summary['direct_connection']['multiple_queries_avg'] - summary['connection_pool']['multiple_queries_avg']) / 
                summary['direct_connection']['multiple_queries_avg'] * 100
            )
            summary['multiple_queries_improvement'] = multiple_queries_diff
        
        self.results['summary'] = summary
    
    def _generate_markdown_report(self, output_file):
        """Generate a markdown report"""
        with open(output_file, 'w') as f:
            f.write("# Connection Pool Performance Report\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Pool settings
            f.write("## Pool Settings\n\n")
            f.write(f"- Minimum Connections: {self.results['pool_settings']['min_conn']}\n")
            f.write(f"- Maximum Connections: {self.results['pool_settings']['max_conn']}\n\n")
            
            # Summary
            f.write("## Summary\n\n")
            
            summary = self.results.get('summary', {})
            
            f.write("| Metric | Direct Connection | Connection Pool | Improvement |\n")
            f.write("|--------|------------------|----------------|-------------|\n")
            
            single_query_direct = summary.get('direct_connection', {}).get('single_query_avg', 0)
            single_query_pool = summary.get('connection_pool', {}).get('single_query_avg', 0)
            single_query_imp = summary.get('single_query_improvement', 0)
            
            f.write(f"| Single Query | {single_query_direct:.6f}s | {single_query_pool:.6f}s | {single_query_imp:.2f}% |\n")
            
            multiple_queries_direct = summary.get('direct_connection', {}).get('multiple_queries_avg', 0)
            multiple_queries_pool = summary.get('connection_pool', {}).get('multiple_queries_avg', 0)
            multiple_queries_imp = summary.get('multiple_queries_improvement', 0)
            
            f.write(f"| Multiple Queries | {multiple_queries_direct:.6f}s | {multiple_queries_pool:.6f}s | {multiple_queries_imp:.2f}% |\n")
            
            concurrent_queries = summary.get('connection_pool', {}).get('concurrent_queries_avg', 0)
            speedup = summary.get('speedup', 0)
            
            f.write(f"| Concurrent Queries | N/A | {concurrent_queries:.6f}s | {speedup:.2f}x speedup |\n")
            
            acquisition_time = summary.get('acquisition_time_avg', 0)
            
            f.write(f"| Connection Acquisition | N/A | {acquisition_time:.6f}s | N/A |\n\n")
            
            # Detailed Results
            f.write("## Detailed Results\n\n")
            
            # Direct Connection
            f.write("### Direct Connection\n\n")
            
            # Single Query
            f.write("#### Single Query\n\n")
            f.write("| Iteration | Time (s) |\n")
            f.write("|-----------|----------|\n")
            
            for i, time_val in enumerate(self.results['direct_connection']['single_query']):
                f.write(f"| {i+1} | {time_val:.6f} |\n")
            
            f.write("\n")
            
            stats = self.results['direct_connection'].get('single_query_stats', {})
            f.write(f"- Mean: {stats.get('mean', 0):.6f}s\n")
            f.write(f"- Median: {stats.get('median', 0):.6f}s\n")
            f.write(f"- Min: {stats.get('min', 0):.6f}s\n")
            f.write(f"- Max: {stats.get('max', 0):.6f}s\n\n")
            
            # Multiple Queries
            f.write("#### Multiple Queries (10 queries)\n\n")
            f.write("| Iteration | Time (s) |\n")
            f.write("|-----------|----------|\n")
            
            for i, time_val in enumerate(self.results['direct_connection']['multiple_queries']):
                f.write(f"| {i+1} | {time_val:.6f} |\n")
            
            f.write("\n")
            
            stats = self.results['direct_connection'].get('multiple_queries_stats', {})
            f.write(f"- Mean: {stats.get('mean', 0):.6f}s\n")
            f.write(f"- Median: {stats.get('median', 0):.6f}s\n")
            f.write(f"- Min: {stats.get('min', 0):.6f}s\n")
            f.write(f"- Max: {stats.get('max', 0):.6f}s\n\n")
            
            # Connection Pool
            f.write("### Connection Pool\n\n")
            
            # Single Query
            f.write("#### Single Query\n\n")
            f.write("| Iteration | Time (s) |\n")
            f.write("|-----------|----------|\n")
            
            for i, time_val in enumerate(self.results['connection_pool']['single_query']):
                f.write(f"| {i+1} | {time_val:.6f} |\n")
            
            f.write("\n")
            
            stats = self.results['connection_pool'].get('single_query_stats', {})
            f.write(f"- Mean: {stats.get('mean', 0):.6f}s\n")
            f.write(f"- Median: {stats.get('median', 0):.6f}s\n")
            f.write(f"- Min: {stats.get('min', 0):.6f}s\n")
            f.write(f"- Max: {stats.get('max', 0):.6f}s\n\n")
            
            # Multiple Queries
            f.write("#### Multiple Queries (10 queries)\n\n")
            f.write("| Iteration | Time (s) |\n")
            f.write("|-----------|----------|\n")
            
            for i, time_val in enumerate(self.results['connection_pool']['multiple_queries']):
                f.write(f"| {i+1} | {time_val:.6f} |\n")
            
            f.write("\n")
            
            stats = self.results['connection_pool'].get('multiple_queries_stats', {})
            f.write(f"- Mean: {stats.get('mean', 0):.6f}s\n")
            f.write(f"- Median: {stats.get('median', 0):.6f}s\n")
            f.write(f"- Min: {stats.get('min', 0):.6f}s\n")
            f.write(f"- Max: {stats.get('max', 0):.6f}s\n\n")
            
            # Concurrent Queries
            f.write("#### Concurrent Queries (10 concurrent queries)\n\n")
            f.write("| Iteration | Time (s) |\n")
            f.write("|-----------|----------|\n")
            
            for i, time_val in enumerate(self.results['connection_pool']['concurrent_queries']):
                f.write(f"| {i+1} | {time_val:.6f} |\n")
            
            f.write("\n")
            
            stats = self.results['connection_pool'].get('concurrent_queries_stats', {})
            f.write(f"- Mean: {stats.get('mean', 0):.6f}s\n")
            f.write(f"- Median: {stats.get('median', 0):.6f}s\n")
            f.write(f"- Min: {stats.get('min', 0):.6f}s\n")
            f.write(f"- Max: {stats.get('max', 0):.6f}s\n\n")
            
            # Connection Acquisition
            f.write("### Connection Acquisition\n\n")
            f.write("| Iteration | Time (s) |\n")
            f.write("|-----------|----------|\n")
            
            for i, time_val in enumerate(self.results['acquisition_times']):
                f.write(f"| {i+1} | {time_val:.6f} |\n")
            
            f.write("\n")
            
            stats = self.results.get('acquisition_times_stats', {})
            f.write(f"- Mean: {stats.get('mean', 0):.6f}s\n")
            f.write(f"- Median: {stats.get('median', 0):.6f}s\n")
            f.write(f"- Min: {stats.get('min', 0):.6f}s\n")
            f.write(f"- Max: {stats.get('max', 0):.6f}s\n\n")
            
            # Pool Scaling
            f.write("### Pool Scaling\n\n")
            f.write("| Iteration | Connections | Time (s) |\n")
            f.write("|-----------|-------------|----------|\n")
            
            for i, scaling in enumerate(self.results['pool_scaling']):
                f.write(f"| {i+1} | {scaling['connections_obtained']} | {scaling['time']:.6f} |\n")
            
            f.write("\n")
            
            # Concurrent Performance
            f.write("### Concurrent Performance\n\n")
            
            concurrent_perf = self.results.get('concurrent_performance', {})
            theoretical = concurrent_perf.get('theoretical_sequential', 0)
            actual = concurrent_perf.get('actual_concurrent', 0)
            speedup = concurrent_perf.get('speedup', 0)
            
            f.write(f"- Theoretical Sequential Time (10 queries): {theoretical:.6f}s\n")
            f.write(f"- Actual Concurrent Time (10 queries): {actual:.6f}s\n")
            f.write(f"- Speedup Factor: {speedup:.2f}x\n\n")
            
            # Recommendation
            f.write("## Recommendations\n\n")
            
            if single_query_imp > 0 and multiple_queries_imp > 0 and speedup > 1:
                f.write("The connection pool shows significant performance improvements:\n\n")
                f.write(f"- {single_query_imp:.2f}% faster for single queries\n")
                f.write(f"- {multiple_queries_imp:.2f}% faster for multiple sequential queries\n")
                f.write(f"- {speedup:.2f}x speedup for concurrent queries\n\n")
                f.write("Recommendation: **Continue using the connection pool in production.**\n")
            elif speedup > 1:
                f.write("The connection pool shows mixed performance:\n\n")
                f.write(f"- {single_query_imp:.2f}% improvement for single queries\n")
                f.write(f"- {multiple_queries_imp:.2f}% improvement for multiple sequential queries\n")
                f.write(f"- {speedup:.2f}x speedup for concurrent queries\n\n")
                f.write("Recommendation: **Use the connection pool primarily for concurrent workloads.**\n")
            else:
                f.write("The connection pool shows limited performance benefits:\n\n")
                f.write(f"- {single_query_imp:.2f}% improvement for single queries\n")
                f.write(f"- {multiple_queries_imp:.2f}% improvement for multiple sequential queries\n")
                f.write(f"- {speedup:.2f}x speedup for concurrent queries\n\n")
                f.write("Recommendation: **Consider optimizing the connection pool or using direct connections for simple workloads.**\n")
            
            f.write("\nOptimal pool settings based on these tests:\n\n")
            f.write(f"- Minimum Connections: {self.results['pool_settings']['min_conn']}\n")
            f.write(f"- Maximum Connections: {self.results['pool_settings']['max_conn']}\n")

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Connection Pool Performance Tests')
    parser.add_argument('--query', type=str, default="SELECT COUNT(*) FROM molecule",
                        help='Query to use for tests')
    parser.add_argument('--iterations', type=int, default=5,
                        help='Number of iterations for each test')
    parser.add_argument('--concurrent', type=int, default=10,
                        help='Number of concurrent queries for concurrent test')
    parser.add_argument('--min-conn', type=int, default=2,
                        help='Minimum connections in pool')
    parser.add_argument('--max-conn', type=int, default=10,
                        help='Maximum connections in pool')
    parser.add_argument('--output', type=str, default=None,
                        help='Output file for results')
    return parser.parse_args()

def main():
    """Main function"""
    try:
        # Parse arguments
        args = parse_args()
        
        print("=================================================")
        print("Connection Pool Performance Tests")
        print("=================================================")
        print(f"Query: {args.query}")
        print(f"Iterations: {args.iterations}")
        print(f"Concurrent Queries: {args.concurrent}")
        print(f"Pool Settings: min_conn={args.min_conn}, max_conn={args.max_conn}")
        print("=================================================")
        
        # Create tester
        tester = ConnectionPoolPerformanceTester(
            DB_CONNECTION_PARAMS,
            min_conn=args.min_conn,
            max_conn=args.max_conn
        )
        
        # Run tests
        tester.test_direct_connection(args.query, args.iterations)
        tester.test_connection_pool(args.query, args.iterations, args.concurrent)
        
        # Generate report
        tester.generate_report(args.output)
        
        print("\nTests completed successfully!")
    
    except Exception as e:
        logger.error(f"Error running tests: {str(e)}", exc_info=True)
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()