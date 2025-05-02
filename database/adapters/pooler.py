#!/usr/bin/env python3
"""
Connection pooler adapter implementation for CryoProtect v2.

This module provides an enhanced connection pooling adapter for PostgreSQL databases.
"""

import os
import logging
import time
import threading
from datetime import datetime, timedelta
import psycopg2
from psycopg2.pool import ThreadedConnectionPool
from psycopg2.extras import RealDictCursor
from typing import Any, Dict, List, Optional, Union, Tuple, Set

from ..connection import ConnectionAdapter

logger = logging.getLogger(__name__)

class PoolerAdapter(ConnectionAdapter):
    """
    Enhanced connection pooling adapter implementation.
    
    This adapter provides advanced connection pooling features including:
    - Connection health monitoring
    - Automatic connection recycling
    - Connection timeout handling
    - Pool statistics tracking
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize adapter with configuration.
        
        Args:
            config: Dictionary containing connection parameters
        """
        self.config = config
        self.connection_pool = None
        self.pool_lock = threading.RLock()
        self.active_connections = set()
        self.last_connection_check = datetime.now()
        self.pool_stats = {
            'created': 0,
            'acquired': 0,
            'released': 0,
            'errors': 0,
            'recycled': 0
        }
        
        # Pool configuration
        self.pool_timeout = int(config.get('pool_timeout', 30))
        self.pool_recycle = int(config.get('pool_recycle', 1800))  # 30 minutes
        
    def connect(self) -> bool:
        """
        Establish connection pool to PostgreSQL.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        try:
            # Extract connection parameters
            host = self.config.get('host', 'localhost')
            port = self.config.get('port', 5432)
            dbname = self.config.get('database', 'postgres')
            user = self.config.get('user', 'postgres')
            password = self.config.get('password', '')
            min_conn = int(self.config.get('min_connections', 1))
            max_conn = int(self.config.get('max_connections', 10))
            
            # Create connection pool
            self.connection_pool = ThreadedConnectionPool(
                minconn=min_conn,
                maxconn=max_conn,
                host=host,
                port=port,
                dbname=dbname,
                user=user,
                password=password
            )
            
            # Start connection monitoring thread
            self._start_connection_monitor()
            
            logger.info(f"Connected to PostgreSQL with enhanced pooling at {host}:{port}/{dbname}")
            self.pool_stats['created'] = min_conn
            return True
        except Exception as e:
            logger.error(f"Failed to connect to PostgreSQL with enhanced pooling: {str(e)}")
            self.pool_stats['errors'] += 1
            return False
            
    def _start_connection_monitor(self):
        """Start a background thread to monitor connection health."""
        def monitor_connections():
            while self.connection_pool is not None:
                try:
                    self._check_connections()
                except Exception as e:
                    logger.error(f"Error in connection monitor: {str(e)}")
                time.sleep(60)  # Check every minute
                
        monitor_thread = threading.Thread(target=monitor_connections, daemon=True)
        monitor_thread.start()
        
    def _check_connections(self):
        """Check connection health and recycle old connections."""
        now = datetime.now()
        
        # Only check connections periodically
        if (now - self.last_connection_check).total_seconds() < self.pool_recycle / 2:
            return
            
        self.last_connection_check = now
        
        # Check if we need to recycle connections
        with self.pool_lock:
            if not self.connection_pool:
                return
                
            # Test a connection from the pool
            try:
                conn = self.connection_pool.getconn()
                with conn.cursor() as cursor:
                    cursor.execute("SELECT 1")
                    cursor.fetchone()
                self.connection_pool.putconn(conn)
            except Exception as e:
                logger.warning(f"Connection health check failed: {str(e)}")
                self.pool_stats['errors'] += 1
                
                # Try to recreate the pool
                try:
                    old_pool = self.connection_pool
                    self.connect()
                    old_pool.closeall()
                    self.pool_stats['recycled'] += 1
                    logger.info("Connection pool recycled after health check failure")
                except Exception as e2:
                    logger.error(f"Failed to recycle connection pool: {str(e2)}")
            
    def disconnect(self) -> bool:
        """
        Close all connections in the pool.
        
        Returns:
            bool: True if disconnection successful, False otherwise
        """
        try:
            with self.pool_lock:
                if self.connection_pool:
                    self.connection_pool.closeall()
                    self.connection_pool = None
                    self.active_connections.clear()
                    logger.info("Disconnected from PostgreSQL with enhanced pooling")
            return True
        except Exception as e:
            logger.error(f"Failed to disconnect from PostgreSQL with enhanced pooling: {str(e)}")
            self.pool_stats['errors'] += 1
            return False
            
    def reconnect(self) -> bool:
        """
        Re-establish a dropped connection to the database.
        
        This method will close any existing connection resources
        and then attempt to establish a new connection.
        
        Returns:
            bool: True if reconnection successful, False otherwise
        """
        logger.info("Attempting to reconnect to PostgreSQL with enhanced pooling")
        
        # First, ensure any existing connections are properly closed
        try:
            with self.pool_lock:
                if self.connection_pool:
                    self.connection_pool.closeall()
                    self.connection_pool = None
                    self.active_connections.clear()
                    logger.info("Closed existing connection pool")
        except Exception as e:
            logger.warning(f"Error closing existing connection pool: {str(e)}")
        
        # Attempt to establish a new connection
        return self.connect()
            
    def _get_connection(self):
        """
        Get a connection from the pool with timeout handling.
        
        Returns:
            Database connection
            
        Raises:
            Exception: If connection acquisition fails
        """
        with self.pool_lock:
            if not self.connection_pool:
                raise RuntimeError("Connection pool is not initialized")
                
            # Try to get a connection with timeout
            start_time = time.time()
            while True:
                try:
                    conn = self.connection_pool.getconn()
                    self.active_connections.add(conn)
                    self.pool_stats['acquired'] += 1
                    return conn
                except psycopg2.pool.PoolError:
                    # Pool is exhausted, wait and retry
                    if time.time() - start_time > self.pool_timeout:
                        self.pool_stats['errors'] += 1
                        raise TimeoutError("Timed out waiting for a database connection")
                    time.sleep(0.1)
                    
    def _release_connection(self, conn):
        """
        Release a connection back to the pool.
        
        Args:
            conn: Database connection
        """
        with self.pool_lock:
            if not self.connection_pool:
                return
                
            if conn in self.active_connections:
                self.active_connections.remove(conn)
                
            try:
                self.connection_pool.putconn(conn)
                self.pool_stats['released'] += 1
            except Exception as e:
                logger.error(f"Error releasing connection: {str(e)}")
                self.pool_stats['errors'] += 1
                
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """
        Execute SQL query and return results.
        
        Args:
            query: SQL query to execute
            params: Query parameters
            
        Returns:
            Query results
        """
        conn = None
        try:
            conn = self._get_connection()
            with conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute(query, params)
                
                if query.strip().upper().startswith('SELECT') or 'RETURNING' in query.upper():
                    result = cursor.fetchall()
                    return result
                else:
                    conn.commit()
                    return cursor.rowcount
        except Exception as e:
            if conn:
                conn.rollback()
            logger.error(f"Error executing query: {str(e)}")
            self.pool_stats['errors'] += 1
            raise
        finally:
            if conn:
                self._release_connection(conn)
                
    def execute_batch(self, queries: List[str]) -> List[Any]:
        """
        Execute multiple SQL queries and return results.
        
        Args:
            queries: List of SQL queries to execute
            
        Returns:
            List of query results
        """
        conn = None
        results = []
        
        try:
            conn = self._get_connection()
            with conn.cursor(cursor_factory=RealDictCursor) as cursor:
                for query in queries:
                    cursor.execute(query)
                    
                    if query.strip().upper().startswith('SELECT') or 'RETURNING' in query.upper():
                        results.append(cursor.fetchall())
                    else:
                        conn.commit()
                        results.append(cursor.rowcount)
                        
            return results
        except Exception as e:
            if conn:
                conn.rollback()
            logger.error(f"Error executing batch: {str(e)}")
            self.pool_stats['errors'] += 1
            raise
        finally:
            if conn:
                self._release_connection(conn)
                
    def begin_transaction(self) -> Any:
        """
        Begin a database transaction.
        
        Returns:
            Connection object representing the transaction
        """
        try:
            conn = self._get_connection()
            conn.autocommit = False
            return conn
        except Exception as e:
            logger.error(f"Error beginning transaction: {str(e)}")
            self.pool_stats['errors'] += 1
            raise
            
    def commit_transaction(self, transaction: Any) -> bool:
        """
        Commit a database transaction.
        
        Args:
            transaction: Connection object representing the transaction
            
        Returns:
            bool: True if commit successful, False otherwise
        """
        try:
            transaction.commit()
            self._release_connection(transaction)
            return True
        except Exception as e:
            logger.error(f"Error committing transaction: {str(e)}")
            self.pool_stats['errors'] += 1
            return False
            
    def rollback_transaction(self, transaction: Any) -> bool:
        """
        Rollback a database transaction.
        
        Args:
            transaction: Connection object representing the transaction
            
        Returns:
            bool: True if rollback successful, False otherwise
        """
        try:
            transaction.rollback()
            self._release_connection(transaction)
            return True
        except Exception as e:
            logger.error(f"Error rolling back transaction: {str(e)}")
            self.pool_stats['errors'] += 1
            return False
            
    def get_connection_info(self) -> Dict[str, Any]:
        """
        Get connection information.
        
        Returns:
            Dict with connection information
        """
        return {
            'type': 'pooler',
            'host': self.config.get('host', 'localhost'),
            'port': self.config.get('port', 5432),
            'database': self.config.get('database', 'postgres'),
            'user': self.config.get('user', 'postgres'),
            'pool_min_size': self.config.get('min_connections', 1),
            'pool_max_size': self.config.get('max_connections', 10),
            'pool_timeout': self.pool_timeout,
            'pool_recycle': self.pool_recycle,
            'active_connections': len(self.active_connections),
            'pool_stats': self.pool_stats,
            'connected': self.connection_pool is not None
        }
        
    def test_connection(self) -> Tuple[bool, str]:
        """
        Test database connection and return status with message.
        
        Returns:
            Tuple of (success: bool, message: str)
        """
        try:
            result = self.execute_query("SELECT 1 as test")
            if result and result[0]['test'] == 1:
                return True, "Connection successful"
            else:
                return False, "Connection test failed"
        except Exception as e:
            return False, f"Connection error: {str(e)}"
            
    def is_healthy(self) -> Tuple[bool, Dict[str, Any]]:
        """
        Perform a comprehensive health check on the database connection.
        
        This method checks various aspects of the connection health:
        - Basic connectivity
        - Connection latency
        - Transaction capability
        - Query execution capability
        - Connection pool status
        - Connection age and recycling status
        
        Returns:
            Tuple of (healthy: bool, health_metrics: Dict[str, Any])
        """
        metrics = {
            'basic_connectivity': False,
            'latency_ms': None,
            'transaction_capability': False,
            'query_capability': False,
            'pool_status': {
                'exists': self.connection_pool is not None,
                'active_connections': len(self.active_connections),
                'min_connections': self.config.get('min_connections', 1),
                'max_connections': self.config.get('max_connections', 10),
                'pool_timeout': self.pool_timeout,
                'pool_recycle': self.pool_recycle,
                'time_since_last_check': (datetime.now() - self.last_connection_check).total_seconds()
            },
            'stats': self.pool_stats.copy()
        }
        
        # Check basic connectivity with latency measurement
        try:
            start_time = time.time()
            success, message = self.test_connection()
            end_time = time.time()
            
            metrics['basic_connectivity'] = success
            metrics['latency_ms'] = round((end_time - start_time) * 1000, 2)
            metrics['connectivity_message'] = message
            
            if not success:
                return False, metrics
        except Exception as e:
            metrics['connectivity_error'] = str(e)
            return False, metrics
            
        # Check transaction capability
        try:
            # Begin a transaction
            transaction = self.begin_transaction()
            
            # Execute a simple query within the transaction
            with transaction.cursor() as cursor:
                cursor.execute("SELECT 1 as test")
                result = cursor.fetchone()
                
            # Commit the transaction
            self.commit_transaction(transaction)
            
            metrics['transaction_capability'] = True
        except Exception as e:
            metrics['transaction_error'] = str(e)
            try:
                # Attempt to rollback if transaction exists
                if 'transaction' in locals():
                    self.rollback_transaction(transaction)
            except:
                pass
                
        # Check query capability with a more complex query
        try:
            # Execute a query that checks database statistics
            result = self.execute_query(
                "SELECT count(*) as table_count FROM information_schema.tables WHERE table_schema = 'public'"
            )
            metrics['query_capability'] = True
            metrics['table_count'] = result[0]['table_count'] if result else 0
        except Exception as e:
            metrics['query_error'] = str(e)
            
        # Check connection pool status if it exists
        if self.connection_pool:
            try:
                with self.pool_lock:
                    metrics['pool_status']['used_connections'] = len(self.connection_pool._used)
                    metrics['pool_status']['unused_connections'] = len(self.connection_pool._unused)
                    metrics['pool_status']['total_connections'] = (
                        len(self.connection_pool._used) + len(self.connection_pool._unused)
                    )
            except Exception as e:
                metrics['pool_status']['error'] = str(e)
                
        # Check if pool needs recycling based on age
        pool_age = (datetime.now() - self.last_connection_check).total_seconds()
        metrics['pool_status']['pool_age_seconds'] = pool_age
        metrics['pool_status']['recycle_needed'] = pool_age > self.pool_recycle
        
        # Determine overall health
        is_healthy = (
            metrics['basic_connectivity'] and
            metrics['transaction_capability'] and
            metrics['query_capability'] and
            metrics['pool_status']['exists'] and
            not metrics['pool_status']['recycle_needed']
        )
        
        return is_healthy, metrics