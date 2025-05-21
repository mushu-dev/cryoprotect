# test_connection_pool.py
import logging
import sys
import threading
import time
from connection_pool_wrapper import ConnectionPoolWrapper, ConnectionManager
from config import get_db_config

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def run_test_query(thread_id):
    """Run a test query in a thread."""
    try:
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                cursor.execute("SELECT pg_sleep(0.5)")
                logger.info(f"Thread {thread_id} completed query")
    except Exception as e:
        logger.error(f"Thread {thread_id} query failed: {str(e)}")

def test_connection_pool():
    """Test the connection pool with concurrent connections."""
    # Create connection pool
    config = get_db_config()
    pool = ConnectionPoolWrapper.get_instance(config)
    
    # Test basic query
    with ConnectionManager() as conn:
        with conn.cursor() as cursor:
            cursor.execute("SELECT 1 as test")
            result = cursor.fetchone()
            logger.info(f"Basic query test: {result}")
    
    # Test multiple connections
    threads = []
    for i in range(5):
        thread = threading.Thread(target=run_test_query, args=(i,))
        threads.append(thread)
        thread.start()
    
    # Wait for threads to complete
    for thread in threads:
        thread.join()
    
    # Check active connections
    logger.info(f"Active connections: {pool.active_connections}")
    
    # Test health check
    pool._check_pool_health()
    
    return True

if __name__ == "__main__":
    logger.info("Testing connection pool...")
    test_connection_pool()
    logger.info("Connection pool test completed!")