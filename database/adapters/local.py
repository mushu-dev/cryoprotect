#!/usr/bin/env python3
"""
Local PostgreSQL adapter implementation for CryoProtect v2.

This module provides a connection adapter for local PostgreSQL databases.
"""

import os
import logging
import time
import json
import psycopg2
from psycopg2.pool import ThreadedConnectionPool
from psycopg2.extras import RealDictCursor
from typing import Any, Dict, List, Optional, Union, Tuple

from ..connection import ConnectionAdapter, generate_connection_request_id

logger = logging.getLogger(__name__)

class LocalAdapter(ConnectionAdapter):
    """
    Local PostgreSQL database adapter implementation.
    
    This adapter provides connection to a local PostgreSQL database
    with connection pooling and comprehensive error handling.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize adapter with configuration.
        
        Args:
            config: Dictionary containing connection parameters
        """
        self.config = config
        self.connection_pool = None
        
    def connect(self) -> bool:
        """
        Establish connection pool to local PostgreSQL.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        request_id = generate_connection_request_id()
        logger.info("[ConnReq:%s] Attempting to connect to local PostgreSQL", request_id)
        
        start_time = time.time()
        
        try:
            # Extract connection parameters
            host = self.config.get('host', 'localhost')
            port = self.config.get('port', 5432)
            dbname = self.config.get('database', 'postgres')
            user = self.config.get('user', 'postgres')
            password = self.config.get('password', '')
            min_conn = int(self.config.get('min_connections', 1))
            max_conn = int(self.config.get('max_connections', 10))
            
            # Log connection attempt (without password)
            logger.debug("[ConnReq:%s] Connection parameters: host=%s, port=%d, dbname=%s, user=%s, min_conn=%d, max_conn=%d",
                        request_id, host, port, dbname, user, min_conn, max_conn)
            
            # Create connection pool
            pool_start_time = time.time()
            self.connection_pool = ThreadedConnectionPool(
                minconn=min_conn,
                maxconn=max_conn,
                host=host,
                port=port,
                dbname=dbname,
                user=user,
                password=password
            )
            pool_creation_time = time.time() - pool_start_time
            
            logger.info("[ConnReq:%s] Created connection pool to local PostgreSQL at %s:%d/%s (duration=%.3fs)",
                       request_id, host, port, dbname, pool_creation_time)
            
            # Verify connection by executing a simple query
            test_start_time = time.time()
            success, message = self.test_connection()
            test_duration = time.time() - test_start_time
            
            if success:
                logger.info("[ConnReq:%s] Connection test successful (duration=%.3fs): %s",
                           request_id, test_duration, message)
            else:
                logger.warning("[ConnReq:%s] Connection test failed (duration=%.3fs): %s",
                              request_id, test_duration, message)
                self.connection_pool = None
                return False
            
            total_duration = time.time() - start_time
            logger.info("[ConnReq:%s] Successfully connected to local PostgreSQL (total_duration=%.3fs)",
                       request_id, total_duration)
            
            return True
        except psycopg2.OperationalError as e:
            total_duration = time.time() - start_time
            logger.error("[ConnReq:%s] Operational error connecting to local PostgreSQL (duration=%.3fs): %s",
                        request_id, total_duration, str(e))
            self.connection_pool = None
            return False
        except Exception as e:
            total_duration = time.time() - start_time
            logger.error("[ConnReq:%s] Failed to connect to local PostgreSQL (duration=%.3fs): %s",
                        request_id, total_duration, str(e))
            self.connection_pool = None
            return False
            
    def disconnect(self) -> bool:
        """
        Close all connections in the pool.
        
        Returns:
            bool: True if disconnection successful, False otherwise
        """
        request_id = generate_connection_request_id()
        logger.info("[ConnReq:%s] Disconnecting from local PostgreSQL", request_id)
        
        start_time = time.time()
        
        if not self.connection_pool:
            logger.debug("[ConnReq:%s] No active connection pool to disconnect", request_id)
            return True
            
        try:
            # Get connection pool stats before closing
            try:
                used_connections = len(self.connection_pool._used)
                unused_connections = len(self.connection_pool._unused)
                total_connections = used_connections + unused_connections
                
                logger.debug("[ConnReq:%s] Connection pool stats before disconnect: used=%d, unused=%d, total=%d",
                            request_id, used_connections, unused_connections, total_connections)
            except Exception as stats_error:
                logger.debug("[ConnReq:%s] Could not get connection pool stats: %s",
                            request_id, str(stats_error))
            
            # Close all connections
            self.connection_pool.closeall()
            
            duration = time.time() - start_time
            logger.info("[ConnReq:%s] Successfully disconnected from local PostgreSQL (duration=%.3fs)",
                       request_id, duration)
            
            self.connection_pool = None
            return True
        except Exception as e:
            duration = time.time() - start_time
            logger.error("[ConnReq:%s] Failed to disconnect from local PostgreSQL (duration=%.3fs): %s",
                        request_id, duration, str(e))
            return False
            
    def reconnect(self) -> bool:
        """
        Re-establish a dropped connection to the database.
        
        This method will close any existing connection resources
        and then attempt to establish a new connection.
        
        Returns:
            bool: True if reconnection successful, False otherwise
        """
        request_id = generate_connection_request_id()
        logger.info("[ConnReq:%s] Attempting to reconnect to local PostgreSQL", request_id)
        
        start_time = time.time()
        
        # First, ensure any existing connections are properly closed
        try:
            if self.connection_pool:
                # Get connection pool stats before closing
                try:
                    used_connections = len(self.connection_pool._used)
                    unused_connections = len(self.connection_pool._unused)
                    total_connections = used_connections + unused_connections
                    
                    logger.debug("[ConnReq:%s] Connection pool stats before reconnect: used=%d, unused=%d, total=%d",
                                request_id, used_connections, unused_connections, total_connections)
                except Exception as stats_error:
                    logger.debug("[ConnReq:%s] Could not get connection pool stats: %s",
                                request_id, str(stats_error))
                
                close_start = time.time()
                self.connection_pool.closeall()
                close_duration = time.time() - close_start
                
                logger.info("[ConnReq:%s] Closed existing connection pool (duration=%.3fs)",
                           request_id, close_duration)
        except Exception as e:
            logger.warning("[ConnReq:%s] Error closing existing connection pool: %s",
                          request_id, str(e))
        
        # Set connection pool to None to ensure we create a new one
        self.connection_pool = None
        
        # Attempt to establish a new connection
        connect_start = time.time()
        result = self.connect()
        connect_duration = time.time() - connect_start
        
        total_duration = time.time() - start_time
        
        if result:
            logger.info("[ConnReq:%s] Successfully reconnected to local PostgreSQL (total_duration=%.3fs)",
                       request_id, total_duration)
        else:
            logger.error("[ConnReq:%s] Failed to reconnect to local PostgreSQL (total_duration=%.3fs)",
                        request_id, total_duration)
        
        return result
            
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
            conn = self.connection_pool.getconn()
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
            raise
        finally:
            if conn:
                self.connection_pool.putconn(conn)
                
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
            conn = self.connection_pool.getconn()
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
            raise
        finally:
            if conn:
                self.connection_pool.putconn(conn)
                
    def begin_transaction(self) -> Any:
        """
        Begin a database transaction.
        
        Returns:
            Connection object representing the transaction
        """
        try:
            conn = self.connection_pool.getconn()
            conn.autocommit = False
            return conn
        except Exception as e:
            logger.error(f"Error beginning transaction: {str(e)}")
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
            self.connection_pool.putconn(transaction)
            return True
        except Exception as e:
            logger.error(f"Error committing transaction: {str(e)}")
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
            self.connection_pool.putconn(transaction)
            return True
        except Exception as e:
            logger.error(f"Error rolling back transaction: {str(e)}")
            return False
            
    def get_connection_info(self) -> Dict[str, Any]:
        """
        Get connection information.
        
        Returns:
            Dict with connection information
        """
        return {
            'type': 'local_postgresql',
            'host': self.config.get('host', 'localhost'),
            'port': self.config.get('port', 5432),
            'database': self.config.get('database', 'postgres'),
            'user': self.config.get('user', 'postgres'),
            'pool_min_size': self.config.get('min_connections', 1),
            'pool_max_size': self.config.get('max_connections', 10),
            'connected': self.connection_pool is not None
        }
        
    def test_connection(self) -> Tuple[bool, str]:
        """
        Test database connection and return status with message.
        
        Returns:
            Tuple of (success: bool, message: str)
        """
        request_id = generate_connection_request_id()
        logger.debug("[ConnReq:%s] Testing local PostgreSQL connection", request_id)
        
        start_time = time.time()
        
        if not self.connection_pool:
            duration = time.time() - start_time
            logger.warning("[ConnReq:%s] Connection test failed: No connection pool exists (duration=%.3fs)",
                          request_id, duration)
            return False, "No connection pool exists"
        
        try:
            query_start = time.time()
            result = self.execute_query("SELECT 1 as test")
            query_duration = time.time() - query_start
            
            if result and result[0]['test'] == 1:
                duration = time.time() - start_time
                logger.debug("[ConnReq:%s] Connection test successful (duration=%.3fs, query_time=%.3fs)",
                            request_id, duration, query_duration)
                return True, "Connection successful"
            else:
                duration = time.time() - start_time
                logger.warning("[ConnReq:%s] Connection test failed: Unexpected result (duration=%.3fs)",
                              request_id, duration)
                return False, "Connection test failed: Unexpected result"
        except psycopg2.OperationalError as e:
            duration = time.time() - start_time
            logger.error("[ConnReq:%s] Connection test failed: Operational error (duration=%.3fs): %s",
                        request_id, duration, str(e))
            return False, f"Connection operational error: {str(e)}"
        except Exception as e:
            duration = time.time() - start_time
            logger.error("[ConnReq:%s] Connection test failed: Unexpected error (duration=%.3fs): %s",
                        request_id, duration, str(e))
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
        
        Returns:
            Tuple of (healthy: bool, health_metrics: Dict[str, Any])
        """
        request_id = generate_connection_request_id()
        logger.info("[ConnReq:%s] Performing comprehensive health check on local PostgreSQL connection",
                   request_id)
        
        health_check_start = time.time()
        
        metrics = {
            'basic_connectivity': False,
            'latency_ms': None,
            'transaction_capability': False,
            'query_capability': False,
            'pool_status': {
                'exists': self.connection_pool is not None,
                'open_connections': 0,
                'min_connections': self.config.get('min_connections', 1),
                'max_connections': self.config.get('max_connections', 10)
            },
            'timestamp': time.time(),
            'check_id': request_id
        }
        
        # Check basic connectivity with latency measurement
        logger.debug("[ConnReq:%s] Checking basic connectivity", request_id)
        try:
            conn_start_time = time.time()
            success, message = self.test_connection()
            conn_end_time = time.time()
            conn_duration = conn_end_time - conn_start_time
            
            metrics['basic_connectivity'] = success
            metrics['latency_ms'] = round(conn_duration * 1000, 2)
            metrics['connectivity_message'] = message
            
            if success:
                logger.debug("[ConnReq:%s] Basic connectivity check passed (latency=%.2fms)",
                            request_id, metrics['latency_ms'])
            else:
                logger.warning("[ConnReq:%s] Basic connectivity check failed: %s (latency=%.2fms)",
                              request_id, message, metrics['latency_ms'])
                return False, metrics
        except Exception as e:
            logger.error("[ConnReq:%s] Error during basic connectivity check: %s",
                        request_id, str(e))
            metrics['connectivity_error'] = str(e)
            return False, metrics
            
        # Check transaction capability
        logger.debug("[ConnReq:%s] Checking transaction capability", request_id)
        try:
            # Begin a transaction
            tx_start_time = time.time()
            transaction = self.begin_transaction()
            
            # Execute a simple query within the transaction
            with transaction.cursor() as cursor:
                cursor.execute("SELECT 1 as test")
                result = cursor.fetchone()
                
            # Commit the transaction
            self.commit_transaction(transaction)
            tx_duration = time.time() - tx_start_time
            
            metrics['transaction_capability'] = True
            metrics['transaction_duration_ms'] = round(tx_duration * 1000, 2)
            logger.debug("[ConnReq:%s] Transaction capability check passed (duration=%.2fms)",
                        request_id, metrics['transaction_duration_ms'])
        except Exception as e:
            logger.warning("[ConnReq:%s] Transaction capability check failed: %s",
                          request_id, str(e))
            metrics['transaction_error'] = str(e)
            try:
                # Attempt to rollback if transaction exists
                if 'transaction' in locals():
                    self.rollback_transaction(transaction)
                    logger.debug("[ConnReq:%s] Transaction rolled back after error", request_id)
            except Exception as rollback_error:
                logger.error("[ConnReq:%s] Error rolling back transaction: %s",
                            request_id, str(rollback_error))
                
        # Check query capability with a more complex query
        logger.debug("[ConnReq:%s] Checking query capability", request_id)
        try:
            # Execute a query that checks database statistics
            query_start_time = time.time()
            result = self.execute_query(
                "SELECT count(*) as table_count FROM information_schema.tables WHERE table_schema = 'public'"
            )
            query_duration = time.time() - query_start_time
            
            metrics['query_capability'] = True
            metrics['query_duration_ms'] = round(query_duration * 1000, 2)
            metrics['table_count'] = result[0]['table_count'] if result else 0
            
            logger.debug("[ConnReq:%s] Query capability check passed (duration=%.2fms, table_count=%d)",
                        request_id, metrics['query_duration_ms'], metrics['table_count'])
        except Exception as e:
            logger.warning("[ConnReq:%s] Query capability check failed: %s",
                          request_id, str(e))
            metrics['query_error'] = str(e)
            
        # Check connection pool status if it exists
        logger.debug("[ConnReq:%s] Checking connection pool status", request_id)
        if self.connection_pool:
            try:
                metrics['pool_status']['open_connections'] = len(self.connection_pool._used)
                metrics['pool_status']['available_connections'] = len(self.connection_pool._unused)
                metrics['pool_status']['total_connections'] = (
                    len(self.connection_pool._used) + len(self.connection_pool._unused)
                )
                
                logger.debug("[ConnReq:%s] Connection pool status: open=%d, available=%d, total=%d",
                            request_id,
                            metrics['pool_status']['open_connections'],
                            metrics['pool_status']['available_connections'],
                            metrics['pool_status']['total_connections'])
            except Exception as e:
                logger.warning("[ConnReq:%s] Error checking connection pool status: %s",
                              request_id, str(e))
                metrics['pool_status']['error'] = str(e)
        else:
            logger.warning("[ConnReq:%s] Connection pool does not exist", request_id)
                
        # Determine overall health
        is_healthy = (
            metrics['basic_connectivity'] and
            metrics['transaction_capability'] and
            metrics['query_capability'] and
            metrics['pool_status']['exists']
        )
        
        # Calculate total health check duration
        total_duration = time.time() - health_check_start
        metrics['total_check_duration_ms'] = round(total_duration * 1000, 2)
        
        if is_healthy:
            logger.info("[ConnReq:%s] Health check passed (duration=%.2fms)",
                       request_id, metrics['total_check_duration_ms'])
        else:
            logger.warning("[ConnReq:%s] Health check failed (duration=%.2fms): %s",
                          request_id, metrics['total_check_duration_ms'],
                          json.dumps({k: v for k, v in metrics.items()
                                     if k not in ['timestamp', 'check_id']}))
        
        return is_healthy, metrics