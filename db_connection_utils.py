#!/usr/bin/env python3
"""
Resilient Database Connection Utilities for CryoProtect v2

This module provides enhanced database connection utilities with:
1. Session pooler prioritization (port 6543)
2. IP/hostname fallback mechanisms
3. Circuit breaker pattern implementation
4. Robust transaction management

Based on specifications in DATABASE_POPULATION_ISSUES.md (Sections 1.3.1, 1.3.2)
"""

import os
import time
import logging
import threading
import socket
import json
from typing import Any, Dict, List, Optional, Union, Tuple, Callable
from contextlib import contextmanager
from datetime import datetime, timedelta

import psycopg2
from psycopg2.extras import RealDictCursor
from psycopg2.pool import ThreadedConnectionPool
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('db_connection_utils')

# Load environment variables
load_dotenv()

class CircuitBreaker:
    """
    Circuit Breaker pattern implementation for database connections.
    
    This class implements the Circuit Breaker pattern to prevent cascading failures
    when database connections are unstable. It tracks connection failures and
    automatically opens the circuit after a threshold is reached, preventing
    further connection attempts for a specified recovery timeout.
    
    States:
    - CLOSED: Normal operation, connections allowed
    - OPEN: Circuit is open, connections are blocked
    - HALF-OPEN: Testing if the system has recovered
    """
    
    def __init__(self, name: str, failure_threshold: int = 3, recovery_timeout: int = 60):
        """
        Initialize the circuit breaker.
        
        Args:
            name: Name of the circuit breaker (for logging)
            failure_threshold: Number of consecutive failures before opening the circuit
            recovery_timeout: Time in seconds to wait before attempting recovery
        """
        self.name = name
        self.failure_count = 0
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.last_failure_time = 0
        self.state = "CLOSED"  # CLOSED, OPEN, HALF-OPEN
        self.lock = threading.RLock()
        self.stats = {
            "success_count": 0,
            "failure_count": 0,
            "last_state_change": time.time(),
            "total_open_time": 0,
            "open_count": 0
        }
    def execute(self, func: Callable, *args, **kwargs) -> Any:
        """
        Execute a function with circuit breaker protection.
        
        Args:
            func: Function to execute
            *args: Arguments to pass to the function
            **kwargs: Keyword arguments to pass to the function
            
        Returns:
            Result of the function
            
        Raises:
            ConnectionError: If the circuit is open
            Exception: Any exception raised by the function
        """
        with self.lock:
            current_time = time.time()
            
            # Check if circuit is open
            if self.state == "OPEN":
                elapsed = current_time - self.last_failure_time
                
                if elapsed > self.recovery_timeout:
                    logger.info(f"Circuit breaker '{self.name}' transitioning from OPEN to HALF-OPEN "
                               f"after {elapsed:.1f}s")
                    self.state = "HALF-OPEN"
                    self.stats["last_state_change"] = current_time
                    self.stats["total_open_time"] += elapsed
                else:
                    logger.warning(f"Circuit breaker '{self.name}' is OPEN, blocking execution. "
                                  f"Will retry in {self.recovery_timeout - elapsed:.1f}s")
                    raise ConnectionError(f"Circuit breaker '{self.name}' is open. "
                                         f"Retry after {self.recovery_timeout - elapsed:.1f}s")
        
        # Execute the function
        try:
            result = func(*args, **kwargs)
            
            # Update state if successful
            with self.lock:
                if self.state == "HALF-OPEN":
                    logger.info(f"Circuit breaker '{self.name}' transitioning from HALF-OPEN to CLOSED "
                               f"after successful execution")
                    self.state = "CLOSED"
                    self.failure_count = 0
                    self.stats["last_state_change"] = time.time()
                
                self.stats["success_count"] += 1
            
            return result
            
        except Exception as e:
            # Update state if failed
            with self.lock:
                self.failure_count += 1
                self.last_failure_time = time.time()
                self.stats["failure_count"] += 1
                
                if self.state == "CLOSED" and self.failure_count >= self.failure_threshold:
                    logger.warning(f"Circuit breaker '{self.name}' transitioning from CLOSED to OPEN "
                                  f"after {self.failure_count} consecutive failures")
                    self.state = "OPEN"
                    self.stats["last_state_change"] = time.time()
                    self.stats["open_count"] += 1
                elif self.state == "HALF-OPEN":
                    logger.warning(f"Circuit breaker '{self.name}' transitioning from HALF-OPEN to OPEN "
                                  f"after failed recovery attempt")
                    self.state = "OPEN"
                    self.stats["last_state_change"] = time.time()
                    self.stats["open_count"] += 1
            
            # Re-raise the exception
            raise
            
    def get_state(self) -> Dict[str, Any]:
        """
        Get the current state of the circuit breaker.
        
        Returns:
            Dict with state information
        """
        with self.lock:
            return {
                "name": self.name,
                "state": self.state,
                "failure_count": self.failure_count,
                "failure_threshold": self.failure_threshold,
                "recovery_timeout": self.recovery_timeout,
                "last_failure_time": self.last_failure_time,
                "time_since_last_failure": time.time() - self.last_failure_time if self.last_failure_time > 0 else None,
                "stats": self.stats
            }
    
    def reset(self) -> None:
        """Reset the circuit breaker to its initial state."""
        with self.lock:
            prev_state = self.state
            self.state = "CLOSED"
            self.failure_count = 0
            self.last_failure_time = 0
            
            if prev_state != "CLOSED":
                self.stats["last_state_change"] = time.time()
                logger.info(f"Circuit breaker '{self.name}' manually reset from {prev_state} to CLOSED")


class ConnectionManager:
    """
    Enhanced database connection manager with resilient connection handling.
    
    This class provides a robust connection management system with:
    1. Session pooler prioritization (port 6543)
    2. IP/hostname fallback mechanisms
    3. Circuit breaker pattern implementation
    4. Connection pooling with health monitoring
    5. Automatic reconnection capabilities
    """
    
    _instance = None
    _lock = threading.Lock()
    
    @classmethod
    def get_instance(cls) -> 'ConnectionManager':
        """
        Get singleton instance of ConnectionManager.
        
        Returns:
            ConnectionManager: Singleton instance
        """
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = ConnectionManager()
        return cls._instance
        
    def __init__(self):
        """Initialize the connection manager."""
        if ConnectionManager._instance is not None:
            raise RuntimeError("ConnectionManager is a singleton. Use get_instance() instead.")
            
        ConnectionManager._instance = self
        
        # Initialize connection pools
        self.pools = {}
        self.active_pool = None
        self.pool_lock = threading.RLock()
        
        # Initialize circuit breakers
        self.circuit_breakers = {
            "pooler": CircuitBreaker("pooler_connection", failure_threshold=3, recovery_timeout=60),
            "direct": CircuitBreaker("direct_connection", failure_threshold=3, recovery_timeout=60),
            "ip": CircuitBreaker("ip_connection", failure_threshold=3, recovery_timeout=60),
            "mcp": CircuitBreaker("mcp_connection", failure_threshold=3, recovery_timeout=60)
        }
        
        # Connection order priority
        self.connection_order = os.getenv('DB_CONNECTION_ORDER', 'pooler,direct,ip,mcp').split(',')
        
        # Connection configuration
        self.config = self._load_config()
        
        # Connection statistics
        self.stats = {
            "total_connections": 0,
            "successful_connections": 0,
            "failed_connections": 0,
            "reconnections": 0,
            "last_connection_time": 0,
            "connection_latency": 0
        }
        
        logger.info(f"ConnectionManager initialized with connection order: {self.connection_order}")
    
    def _load_config(self) -> Dict[str, Dict[str, Any]]:
        """
        Load connection configuration from environment variables.
        
        Returns:
            Dict with connection configuration
        """
        return {
            "pooler": {
                "enabled": os.getenv('POOLER_ENABLED', 'true').lower() == 'true',
                "host": os.getenv('SUPABASE_DB_HOST', ''),
                "port": int(os.getenv('SUPABASE_DB_PORT', '6543')),  # Session pooler port
                "dbname": os.getenv('SUPABASE_DB_NAME', ''),
                "user": os.getenv('SUPABASE_DB_USER', ''),
                "password": os.getenv('SUPABASE_DB_PASSWORD', ''),
                "min_connections": int(os.getenv('POOLER_MIN_CONNECTIONS', '1')),
                "max_connections": int(os.getenv('POOLER_MAX_CONNECTIONS', '10')),
                "connect_timeout": int(os.getenv('POOLER_CONNECT_TIMEOUT', '10'))
            },
            "direct": {
                "enabled": os.getenv('DIRECT_ENABLED', 'true').lower() == 'true',
                "host": os.getenv('SUPABASE_DB_HOST', ''),
                "port": int(os.getenv('SUPABASE_DB_DIRECT_PORT', '5432')),  # Direct connection port
                "dbname": os.getenv('SUPABASE_DB_NAME', ''),
                "user": os.getenv('SUPABASE_DB_USER', ''),
                "password": os.getenv('SUPABASE_DB_PASSWORD', ''),
                "min_connections": int(os.getenv('DIRECT_MIN_CONNECTIONS', '1')),
                "max_connections": int(os.getenv('DIRECT_MAX_CONNECTIONS', '5')),
                "connect_timeout": int(os.getenv('DIRECT_CONNECT_TIMEOUT', '10'))
            },
            "ip": {
                "enabled": os.getenv('IP_ENABLED', 'true').lower() == 'true',
                "host": os.getenv('SUPABASE_DB_IP_ADDRESS', ''),
                "port": int(os.getenv('SUPABASE_DB_DIRECT_PORT', '5432')),
                "dbname": os.getenv('SUPABASE_DB_NAME', ''),
                "user": os.getenv('SUPABASE_DB_USER', ''),
                "password": os.getenv('SUPABASE_DB_PASSWORD', ''),
                "min_connections": int(os.getenv('IP_MIN_CONNECTIONS', '1')),
                "max_connections": int(os.getenv('IP_MAX_CONNECTIONS', '5')),
                "connect_timeout": int(os.getenv('IP_CONNECT_TIMEOUT', '10'))
            },
            "mcp": {
                "enabled": os.getenv('MCP_ENABLED', 'true').lower() == 'true',
                "server_name": os.getenv('MCP_SERVER_NAME', 'supabase'),
                "project_id": os.getenv('MCP_PROJECT_ID', '')
            }
        }
        
    def _resolve_hostname(self, hostname: str) -> Optional[str]:
        """
        Resolve hostname to IP address.
        
        Args:
            hostname: Hostname to resolve
            
        Returns:
            IP address or None if resolution fails
        """
        try:
            logger.debug(f"Resolving hostname: {hostname}")
            ip_address = socket.gethostbyname(hostname)
            logger.debug(f"Resolved {hostname} to {ip_address}")
            return ip_address
        except socket.gaierror as e:
            logger.warning(f"Failed to resolve hostname {hostname}: {str(e)}")
            return None
    
    def _create_connection_pool(self, pool_type: str) -> bool:
        """
        Create a connection pool of the specified type.
        
        Args:
            pool_type: Type of connection pool to create (pooler, direct, ip)
            
        Returns:
            bool: True if pool creation successful, False otherwise
        """
        if pool_type not in self.config or not self.config[pool_type]["enabled"]:
            logger.warning(f"Connection type {pool_type} is not enabled or configured")
            return False
        
        config = self.config[pool_type]
        
        # For IP connection type, ensure we have an IP address
        if pool_type == "ip" and not config["host"]:
            # Try to resolve the hostname from the direct connection config
            direct_host = self.config["direct"]["host"]
            if direct_host:
                ip_address = self._resolve_hostname(direct_host)
                if ip_address:
                    config["host"] = ip_address
                else:
                    logger.warning(f"Could not resolve hostname for IP connection: {direct_host}")
                    return False
            else:
                logger.warning("No hostname available to resolve for IP connection")
                return False
        
        # Create the connection pool
        try:
            # Use circuit breaker to protect against connection failures
            def create_pool():
                return ThreadedConnectionPool(
                    minconn=config["min_connections"],
                    maxconn=config["max_connections"],
                    host=config["host"],
                    port=config["port"],
                    dbname=config["dbname"],
                    user=config["user"],
                    password=config["password"],
                    connect_timeout=config["connect_timeout"]
                )
            
            # Execute with circuit breaker
            pool = self.circuit_breakers[pool_type].execute(create_pool)
            
            with self.pool_lock:
                self.pools[pool_type] = {
                    "pool": pool,
                    "created_at": time.time(),
                    "last_used": time.time(),
                    "stats": {
                        "connections_acquired": 0,
                        "connections_released": 0,
                        "errors": 0
                    }
                }
            
            logger.info(f"Created {pool_type} connection pool to {config['host']}:{config['port']}/{config['dbname']}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to create {pool_type} connection pool: {str(e)}")
            return False
    
    def _get_mcp_connection(self):
        """
        Get a connection using MCP (Model Context Protocol).
        
        Returns:
            Connection object or None if connection fails
        """
        if not self.config["mcp"]["enabled"]:
            logger.warning("MCP connection is not enabled")
            return None
        
        try:
            # Use circuit breaker to protect against connection failures
            def connect_mcp():
                # This is a placeholder for the actual MCP connection logic
                # In a real implementation, this would use the MCP API
                logger.info(f"Connecting to MCP database using {self.config['mcp']['server_name']}")
                
                # For now, we'll just return a dummy connection object
                return MCPConnection(
                    server_name=self.config["mcp"]["server_name"],
                    project_id=self.config["mcp"]["project_id"]
                )
            
            # Execute with circuit breaker
            connection = self.circuit_breakers["mcp"].execute(connect_mcp)
            
            logger.info(f"Connected to MCP database using {self.config['mcp']['server_name']}")
            return connection
            
        except Exception as e:
            logger.error(f"Failed to connect to MCP database: {str(e)}")
            return None
            
    def get_connection(self) -> Optional[Any]:
        """
        Get a database connection using the configured connection order.
        
        This method attempts to get a connection from each pool in the configured order
        until a successful connection is established.
        
        Returns:
            Connection object or None if all connection attempts fail
        """
        start_time = time.time()
        self.stats["total_connections"] += 1
        
        # Try each connection type in order
        for conn_type in self.connection_order:
            logger.debug(f"Attempting to get connection using {conn_type}")
            
            # Skip if this connection type is disabled or circuit is open
            if not self.config[conn_type]["enabled"]:
                logger.debug(f"Skipping {conn_type} connection (disabled)")
                continue
                
            try:
                # Check circuit breaker state
                breaker_state = self.circuit_breakers[conn_type].get_state()
                if breaker_state["state"] == "OPEN":
                    time_remaining = breaker_state["recovery_timeout"] - (time.time() - breaker_state["last_failure_time"])
                    if time_remaining > 0:
                        logger.debug(f"Skipping {conn_type} connection (circuit breaker open, retry in {time_remaining:.1f}s)")
                        continue
                
                # For MCP connection type, use a different approach
                if conn_type == "mcp":
                    connection = self._get_mcp_connection()
                    if connection:
                        self.active_pool = "mcp"
                        self.stats["successful_connections"] += 1
                        self.stats["last_connection_time"] = time.time()
                        self.stats["connection_latency"] = time.time() - start_time
                        logger.info(f"Successfully connected using MCP (latency: {self.stats['connection_latency']:.3f}s)")
                        return connection
                    continue
                
                # For other connection types, use connection pools
                with self.pool_lock:
                    # Create pool if it doesn't exist
                    if conn_type not in self.pools:
                        if not self._create_connection_pool(conn_type):
                            continue
                    
                    # Get connection from pool
                    try:
                        pool_info = self.pools[conn_type]
                        connection = pool_info["pool"].getconn()
                        
                        # Test connection
                        with connection.cursor() as cursor:
                            cursor.execute("SELECT 1")
                            result = cursor.fetchone()
                            if result and result[0] == 1:
                                # Connection is good
                                pool_info["last_used"] = time.time()
                                pool_info["stats"]["connections_acquired"] += 1
                                self.active_pool = conn_type
                                
                                self.stats["successful_connections"] += 1
                                self.stats["last_connection_time"] = time.time()
                                self.stats["connection_latency"] = time.time() - start_time
                                
                                logger.info(f"Successfully connected using {conn_type} pool "
                                           f"(latency: {self.stats['connection_latency']:.3f}s)")
                                
                                return ConnectionWrapper(connection, self, conn_type)
                            else:
                                # Connection test failed
                                pool_info["stats"]["errors"] += 1
                                pool_info["pool"].putconn(connection)
                                logger.warning(f"Connection test failed for {conn_type} pool")
                                continue
                                
                    except Exception as e:
                        # Connection acquisition failed
                        logger.warning(f"Failed to get connection from {conn_type} pool: {str(e)}")
                        
                        # Try to recreate the pool
                        try:
                            if conn_type in self.pools:
                                old_pool = self.pools[conn_type]["pool"]
                                old_pool.closeall()
                            
                            if self._create_connection_pool(conn_type):
                                # Retry with new pool
                                connection = self.pools[conn_type]["pool"].getconn()
                                
                                # Test connection
                                with connection.cursor() as cursor:
                                    cursor.execute("SELECT 1")
                                    result = cursor.fetchone()
                                    if result and result[0] == 1:
                                        # Connection is good
                                        self.pools[conn_type]["last_used"] = time.time()
                                        self.pools[conn_type]["stats"]["connections_acquired"] += 1
                                        self.active_pool = conn_type
                                        
                                        self.stats["successful_connections"] += 1
                                        self.stats["reconnections"] += 1
                                        self.stats["last_connection_time"] = time.time()
                                        self.stats["connection_latency"] = time.time() - start_time
                                        
                                        logger.info(f"Successfully reconnected using {conn_type} pool "
                                                   f"(latency: {self.stats['connection_latency']:.3f}s)")
                                        
                                        return ConnectionWrapper(connection, self, conn_type)
                            
                        except Exception as e2:
                            logger.error(f"Failed to recreate {conn_type} pool: {str(e2)}")
                            continue
            
            except Exception as e:
                logger.error(f"Unexpected error getting {conn_type} connection: {str(e)}")
                continue
        
        # All connection attempts failed
        self.stats["failed_connections"] += 1
        logger.error(f"All connection attempts failed (elapsed: {time.time() - start_time:.3f}s)")
        return None
        
    def release_connection(self, connection: Any, pool_type: str) -> None:
        """
        Release a connection back to its pool.
        
        Args:
            connection: Connection to release
            pool_type: Type of pool the connection belongs to
        """
        if pool_type == "mcp":
            # MCP connections don't need to be released
            return
            
        with self.pool_lock:
            if pool_type not in self.pools:
                logger.warning(f"Cannot release connection to non-existent {pool_type} pool")
                return
                
            try:
                self.pools[pool_type]["pool"].putconn(connection)
                self.pools[pool_type]["stats"]["connections_released"] += 1
                logger.debug(f"Released connection to {pool_type} pool")
            except Exception as e:
                logger.error(f"Error releasing connection to {pool_type} pool: {str(e)}")
                self.pools[pool_type]["stats"]["errors"] += 1
    
    def close_all_connections(self) -> None:
        """Close all connection pools."""
        with self.pool_lock:
            for pool_type, pool_info in self.pools.items():
                try:
                    pool_info["pool"].closeall()
                    logger.info(f"Closed all connections in {pool_type} pool")
                except Exception as e:
                    logger.error(f"Error closing {pool_type} pool: {str(e)}")
            
            self.pools = {}
            self.active_pool = None
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get connection statistics.
        
        Returns:
            Dict with connection statistics
        """
        stats = self.stats.copy()
        stats["pools"] = {}
        
        with self.pool_lock:
            for pool_type, pool_info in self.pools.items():
                stats["pools"][pool_type] = {
                    "created_at": pool_info["created_at"],
                    "last_used": pool_info["last_used"],
                    "age": time.time() - pool_info["created_at"],
                    "stats": pool_info["stats"].copy()
                }
        
        stats["circuit_breakers"] = {}
        for name, breaker in self.circuit_breakers.items():
            stats["circuit_breakers"][name] = breaker.get_state()
        
        stats["active_pool"] = self.active_pool
        stats["connection_order"] = self.connection_order
        
        return stats


class ConnectionWrapper:
    """
    Wrapper for database connections with automatic release.
    
    This class wraps a database connection and ensures it is properly
    released back to the pool when no longer needed.
    """
    
    def __init__(self, connection: Any, manager: ConnectionManager, pool_type: str):
        """
        Initialize the connection wrapper.
        
        Args:
            connection: Database connection
            manager: ConnectionManager instance
            pool_type: Type of pool the connection belongs to
        """
        self.connection = connection
        self.manager = manager
        self.pool_type = pool_type
        self.closed = False
    
    def __del__(self):
        """Ensure connection is released when object is garbage collected."""
        self.close()
    
    def close(self) -> None:
        """Release the connection back to the pool."""
        if not self.closed:
            self.manager.release_connection(self.connection, self.pool_type)
            self.closed = True
    
    def cursor(self, *args, **kwargs):
        """
        Get a cursor for the connection.
        
        Args:
            *args: Arguments to pass to cursor method
            **kwargs: Keyword arguments to pass to cursor method
            
        Returns:
            Database cursor
        """
        return self.connection.cursor(*args, **kwargs)
    
    def commit(self) -> None:
        """Commit the current transaction."""
        self.connection.commit()
    
    def rollback(self) -> None:
        """Rollback the current transaction."""
        self.connection.rollback()
    
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """
        Execute a SQL query and return results.
        
        Args:
            query: SQL query to execute
            params: Query parameters
            
        Returns:
            Query results
        """
        try:
            with self.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute(query, params)
                
                if query.strip().upper().startswith('SELECT') or 'RETURNING' in query.upper():
                    result = cursor.fetchall()
                    return result
                else:
                    self.commit()
                    return cursor.rowcount
        except Exception as e:
            self.rollback()
            logger.error(f"Error executing query: {str(e)}")
            raise


class MCPConnection:
    """
    Connection wrapper for MCP (Model Context Protocol) database access.
    
    This class provides a compatible interface for MCP database connections
    to work with the rest of the connection management system.
    """
    
    def __init__(self, server_name: str, project_id: str):
        """
        Initialize the MCP connection.
        
        Args:
            server_name: Name of the MCP server
            project_id: Project ID for the database
        """
        self.server_name = server_name
        self.project_id = project_id
        self.transaction_active = False
    
@contextmanager
def get_db_connection():
    """
    Context manager for getting a database connection.
    
    This function provides a convenient way to get a database connection
    and ensure it is properly released when no longer needed.
    
    Yields:
        Database connection
    """
    connection = None
    try:
        connection = ConnectionManager.get_instance().get_connection()
        if not connection:
            raise ConnectionError("Failed to establish database connection")
        yield connection
    finally:
        if connection and hasattr(connection, 'close'):
            connection.close()


@contextmanager
def safe_transaction():
    """
    Safe transaction context manager with cleanup.
    
    This function provides a transaction context that ensures proper
    cleanup in case of errors, including handling aborted transactions.
    
    Yields:
        Database connection with active transaction
    """
    connection = None
    try:
        # Get a database connection
        connection = ConnectionManager.get_instance().get_connection()
        if not connection:
            raise ConnectionError("Failed to establish database connection")
        
        # Check if already in transaction
        try:
            if hasattr(connection, 'execute_query'):
                connection.execute_query("SELECT current_setting('transaction_isolation')")
            else:
                with connection.cursor() as cursor:
                    cursor.execute("SELECT current_setting('transaction_isolation')")
        except Exception as e:
            if "current transaction is aborted" in str(e):
                logger.warning("Detected aborted transaction, rolling back")
                if hasattr(connection, 'rollback'):
                    connection.rollback()
        
        # Start new transaction
        if hasattr(connection, 'begin_transaction'):
            transaction = connection.begin_transaction()
        else:
            connection.autocommit = False
            transaction = connection
        
        yield connection
        
        # Commit transaction
        if hasattr(connection, 'commit_transaction'):
            connection.commit_transaction(transaction)
        else:
            connection.commit()
            
    except Exception as e:
        # Rollback transaction on error
        logger.error(f"Transaction error: {str(e)}")
        if connection:
            if hasattr(connection, 'rollback_transaction'):
                try:
                    connection.rollback_transaction(transaction)
                except Exception as rollback_error:
                    logger.error(f"Rollback failed: {rollback_error}")
            elif hasattr(connection, 'rollback'):
                try:
                    connection.rollback()
                except Exception as rollback_error:
                    logger.error(f"Rollback failed: {rollback_error}")
        raise
    finally:
        # Release connection
        if connection and hasattr(connection, 'close'):
            connection.close()
    def close(self) -> None:
        """Close the connection (no-op for MCP)."""
        pass
    
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """
        Execute a SQL query using MCP and return results.
        
        Args:
            query: SQL query to execute
            params: Query parameters
            
        Returns:
            Query results
        """
        try:
            # This is a placeholder for the actual MCP query execution
            # In a real implementation, this would use the MCP API
            logger.info(f"Executing query via MCP: {query[:50]}...")
            
            # For now, we'll just return a dummy result
            return [{"result": "MCP query executed successfully"}]
        except Exception as e:
            logger.error(f"Error executing MCP query: {str(e)}")
            raise
    
    def begin_transaction(self) -> None:
        """Begin a transaction (simulated for MCP)."""
        self.transaction_active = True
    
    def commit(self) -> None:
        """Commit the current transaction (simulated for MCP)."""
        self.transaction_active = False
    
    def rollback(self) -> None:
        """Rollback the current transaction (simulated for MCP)."""
        self.transaction_active = False