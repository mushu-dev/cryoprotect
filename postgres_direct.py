#!/usr/bin/env python3
"""
Mock implementation of postgres_direct.py for ChEMBL import script.
"""

import logging
from typing import Dict, List, Any, Optional, Tuple

logger = logging.getLogger(__name__)

class PostgresDirectConnection:
    """
    Mock implementation of PostgresDirectConnection for ChEMBL import script.
    """
    _instance = None
    
    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(PostgresDirectConnection, cls).__new__(cls)
            cls._instance.initialized = False
        return cls._instance
    
    def __init__(self, host=None, port=None, database=None, user=None, password=None):
        if self.initialized:
            return
            
        self.host = host or "localhost"
        self.port = port or 5432
        self.database = database or "postgres"
        self.user = user or "postgres"
        self.password = password or ""
        self.connection_pool = None
        self.initialized = True
        
        logger.info(f"Initialized PostgresDirectConnection to {self.host}:{self.port}/{self.database}")
        
    def get_connection(self):
        """Get a connection from the pool."""
        logger.info("Mock: Getting connection from pool")
        return MockConnection()
        
    def execute_query(self, query, params=None):
        """Execute a query and return the results."""
        logger.info(f"Mock: Executing query: {query}")
        return []
        
    def execute_batch(self, query, params_list):
        """Execute a batch of queries."""
        logger.info(f"Mock: Executing batch query: {query}")
        return []
        
    def transaction(self):
        """Return a transaction context manager."""
        return MockTransaction()
        
    def close(self):
        """Close all connections in the pool."""
        logger.info("Mock: Closing connection pool")
        
class MockConnection:
    """Mock database connection."""
    
    def cursor(self):
        """Return a cursor."""
        return MockCursor()
        
    def commit(self):
        """Commit the transaction."""
        logger.info("Mock: Committing transaction")
        
    def rollback(self):
        """Rollback the transaction."""
        logger.info("Mock: Rolling back transaction")
        
    def close(self):
        """Close the connection."""
        logger.info("Mock: Closing connection")
        
class MockCursor:
    """Mock database cursor."""
    
    def execute(self, query, params=None):
        """Execute a query."""
        logger.info(f"Mock: Executing query: {query}")
        
    def executemany(self, query, params_list):
        """Execute a batch of queries."""
        logger.info(f"Mock: Executing batch query: {query}")
        
    def fetchall(self):
        """Fetch all results."""
        return []
        
    def fetchone(self):
        """Fetch one result."""
        return None
        
    def close(self):
        """Close the cursor."""
        logger.info("Mock: Closing cursor")
        
class MockTransaction:
    """Mock transaction context manager."""
    
    def __enter__(self):
        logger.info("Mock: Entering transaction")
        return MockConnection()
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            logger.info("Mock: Committing transaction")
        else:
            logger.info(f"Mock: Rolling back transaction due to {exc_type.__name__}: {exc_val}")
        return False

# Create a singleton instance
db = PostgresDirectConnection()