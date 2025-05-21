#!/usr/bin/env python3
"""
Simplified mock version of db_connection_utils.py for testing purposes.
"""

import logging
from contextlib import contextmanager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('db_connection_utils_mock')

@contextmanager
def get_db_connection():
    """
    Mock implementation of get_db_connection for testing.
    
    This function returns a mock connection object that can be used for testing.
    """
    logger.info("Getting mock database connection")
    
    class MockConnection:
        def cursor(self, *args, **kwargs):
            return MockCursor()
            
        def commit(self):
            logger.info("Mock commit")
            
        def rollback(self):
            logger.info("Mock rollback")
            
        def close(self):
            logger.info("Mock close")
    
    class MockCursor:
        def __enter__(self):
            return self
            
        def __exit__(self, exc_type, exc_val, exc_tb):
            pass
            
        def execute(self, query, params=None):
            logger.info(f"Mock execute: {query}")
            return []
            
        def fetchone(self):
            return {"id": "mock-id"}
            
        def fetchall(self):
            return [{"id": "mock-id"}]
    
    try:
        connection = MockConnection()
        yield connection
    finally:
        logger.info("Mock connection closed")

@contextmanager
def safe_transaction():
    """
    Mock implementation of safe_transaction for testing.
    """
    logger.info("Starting mock transaction")
    try:
        yield None
    except Exception as e:
        logger.error(f"Mock transaction error: {str(e)}")
        raise
    finally:
        logger.info("Ending mock transaction")