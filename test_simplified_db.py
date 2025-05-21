#!/usr/bin/env python3
"""
Test script for simplified database module.
"""

import logging
import sys
import os
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from database import db, utils

def test_connection():
    """Test basic database connection."""
    logger.info("Testing database connection...")

    # Create config from Supabase settings
    config = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
        'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }

    logger.info(f"Connecting to: {config['host']}:{config['port']}/{config['database']} as {config['user']}")

    # Initialize database pool
    init_result = db.init_connection_pool(config=config)
    logger.info(f"Database initialization: {'Success' if init_result else 'Failed'}")

    if not init_result:
        logger.error("Failed to initialize database connection pool")
        return False

    # Test simple query
    success, message = db.test_connection()
    logger.info(f"Connection test: {message}")

    if not success:
        logger.error("Connection test failed")
        return False

    # Get list of tables
    tables = db.get_tables()
    if tables:
        logger.info(f"Found {len(tables)} tables in database")
        logger.info(f"Tables: {', '.join(tables[:10])}" + ("..." if len(tables) > 10 else ""))
    else:
        logger.warning("No tables found or could not retrieve table list")

    # Get more detailed stats
    try:
        stats = utils.test_database_connection()
        logger.info(f"Database size: {stats.get('database_size')}")
        logger.info(f"Table count: {stats.get('tables_count')}")
    except Exception as e:
        logger.error(f"Error getting database stats: {str(e)}")

    return success

def test_transaction():
    """Test transaction handling."""
    logger.info("Testing transaction handling...")
    
    try:
        # Execute a simple transaction
        with db.transaction() as cursor:
            cursor.execute("SELECT 1 as test")
            result = cursor.fetchone()
            
            if result and result['test'] == 1:
                logger.info("Transaction test successful")
                return True
            else:
                logger.error("Transaction test returned unexpected result")
                return False
    except Exception as e:
        logger.error(f"Transaction test failed: {str(e)}")
        return False

def test_batch_operations():
    """Test batch query execution."""
    logger.info("Testing batch operations...")
    
    try:
        queries = [
            "SELECT 1 as test1",
            "SELECT 2 as test2",
            "SELECT 3 as test3"
        ]
        
        results = db.execute_batch(queries)
        
        if (results and len(results) == 3 and 
                results[0][0]['test1'] == 1 and 
                results[1][0]['test2'] == 2 and 
                results[2][0]['test3'] == 3):
            logger.info("Batch operations test successful")
            return True
        else:
            logger.error("Batch operations test returned unexpected results")
            return False
    except Exception as e:
        logger.error(f"Batch operations test failed: {str(e)}")
        return False

def main():
    """Run all tests."""
    logger.info("Starting database module tests...")
    
    tests = [
        ("Connection Test", test_connection),
        ("Transaction Test", test_transaction),
        ("Batch Operations Test", test_batch_operations)
    ]
    
    success_count = 0
    
    for test_name, test_func in tests:
        logger.info(f"\n--- Running {test_name} ---")
        try:
            result = test_func()
            if result:
                logger.info(f"{test_name}: SUCCESS")
                success_count += 1
            else:
                logger.error(f"{test_name}: FAILED")
        except Exception as e:
            logger.error(f"{test_name}: ERROR - {str(e)}")
    
    logger.info(f"\nTest Summary: {success_count}/{len(tests)} tests passed")
    
    # Close connections
    db.close_all_connections()
    
    return success_count == len(tests)

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)