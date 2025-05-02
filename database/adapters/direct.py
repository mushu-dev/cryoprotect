#!/usr/bin/env python3
"""
Direct PostgreSQL adapter implementation for CryoProtect v2.

This module provides a direct connection adapter for PostgreSQL databases
with DNS/IP resolution fallback.
"""

import os
import logging
import socket
import psycopg2
from psycopg2.pool import ThreadedConnectionPool
from psycopg2.extras import RealDictCursor
from typing import Any, Dict, List, Optional, Union, Tuple
import subprocess

from ..connection import ConnectionAdapter

logger = logging.getLogger(__name__)

class DirectAdapter(ConnectionAdapter):
    """
    Direct PostgreSQL database adapter implementation.
    
    This adapter provides direct connection to a PostgreSQL database
    with DNS/IP resolution fallback for improved reliability.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize adapter with configuration.
        
        Args:
            config: Dictionary containing connection parameters
        """
        self.config = config
        self.connection_pool = None
        self.ip_address = None
        
    def _resolve_hostname(self, hostname: str) -> Optional[str]:
        """
        Resolve hostname to IP address with multiple methods.
        
        Args:
            hostname: Hostname to resolve
            
        Returns:
            IP address or None if resolution fails
        """
        methods = [
            self._resolve_with_socket,
            self._resolve_with_nslookup,
            self._resolve_with_dig,
            self._resolve_with_alternative_dns
        ]
        
        for method in methods:
            try:
                ip = method(hostname)
                if ip:
                    logger.info(f"Resolved {hostname} to {ip} using {method.__name__}")
                    return ip
            except Exception as e:
                logger.debug(f"Failed to resolve {hostname} using {method.__name__}: {str(e)}")
                
        return None
        
    def _resolve_with_socket(self, hostname: str) -> Optional[str]:
        """Resolve hostname using socket."""
        try:
            return socket.gethostbyname(hostname)
        except:
            return None
            
    def _resolve_with_nslookup(self, hostname: str) -> Optional[str]:
        """Resolve hostname using nslookup."""
        try:
            output = subprocess.check_output(["nslookup", hostname], universal_newlines=True)
            for line in output.splitlines():
                if "Address:" in line and not "localhost" in line:
                    return line.split("Address:")[1].strip()
            return None
        except:
            return None
            
    def _resolve_with_dig(self, hostname: str) -> Optional[str]:
        """Resolve hostname using dig."""
        try:
            output = subprocess.check_output(["dig", "+short", hostname], universal_newlines=True)
            if output:
                return output.strip()
            return None
        except:
            return None
            
    def _resolve_with_alternative_dns(self, hostname: str) -> Optional[str]:
        """Resolve hostname using alternative DNS servers."""
        dns_servers = ["8.8.8.8", "1.1.1.1", "9.9.9.9"]
        
        for dns in dns_servers:
            try:
                output = subprocess.check_output(
                    ["nslookup", hostname, dns],
                    universal_newlines=True
                )
                for line in output.splitlines():
                    if "Address:" in line and not dns in line and not "localhost" in line:
                        return line.split("Address:")[1].strip()
            except:
                continue
                
        return None
        
    def connect(self) -> bool:
        """
        Establish connection pool to PostgreSQL.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        try:
            # Extract connection parameters
            host = self.config.get('host')
            port = self.config.get('port', 5432)
            dbname = self.config.get('database', 'postgres')
            user = self.config.get('user')
            password = self.config.get('password')
            min_conn = int(self.config.get('min_connections', 1))
            max_conn = int(self.config.get('max_connections', 10))
            
            if not host or not user or not password:
                logger.error("Missing required connection parameters")
                return False
                
            # Try direct hostname connection
            try:
                self.connection_pool = ThreadedConnectionPool(
                    minconn=min_conn,
                    maxconn=max_conn,
                    host=host,
                    port=port,
                    dbname=dbname,
                    user=user,
                    password=password
                )
                logger.info(f"Connected to PostgreSQL at {host}:{port}/{dbname}")
                return True
            except (psycopg2.OperationalError, socket.gaierror) as e:
                logger.warning(f"Failed to connect with hostname: {str(e)}")
                
                # Resolve hostname to IP
                self.ip_address = self._resolve_hostname(host)
                if not self.ip_address:
                    logger.error(f"Failed to resolve hostname {host} to IP address")
                    return False
                    
                # Try connection with IP address
                try:
                    self.connection_pool = ThreadedConnectionPool(
                        minconn=min_conn,
                        maxconn=max_conn,
                        host=self.ip_address,
                        port=port,
                        dbname=dbname,
                        user=user,
                        password=password
                    )
                    logger.info(f"Connected to PostgreSQL at {self.ip_address}:{port}/{dbname} (IP fallback)")
                    return True
                except Exception as e2:
                    logger.error(f"Failed to connect with IP fallback: {str(e2)}")
                    return False
        except Exception as e:
            logger.error(f"Failed to connect to PostgreSQL: {str(e)}")
            return False
            
    def disconnect(self) -> bool:
        """
        Close all connections in the pool.
        
        Returns:
            bool: True if disconnection successful, False otherwise
        """
        try:
            if self.connection_pool:
                self.connection_pool.closeall()
                logger.info("Disconnected from PostgreSQL")
            return True
        except Exception as e:
            logger.error(f"Failed to disconnect from PostgreSQL: {str(e)}")
            return False
            
    def reconnect(self) -> bool:
        """
        Re-establish a dropped connection to the database.
        
        This method will close any existing connection resources
        and then attempt to establish a new connection with DNS/IP fallback.
        
        Returns:
            bool: True if reconnection successful, False otherwise
        """
        logger.info("Attempting to reconnect to PostgreSQL with direct connection")
        
        # First, ensure any existing connections are properly closed
        try:
            if self.connection_pool:
                self.connection_pool.closeall()
                self.connection_pool = None
                logger.info("Closed existing connection pool")
        except Exception as e:
            logger.warning(f"Error closing existing connection pool: {str(e)}")
        
        # Reset IP address to force DNS resolution again
        # This helps if the IP address has changed or if DNS issues were temporary
        self.ip_address = None
        
        # Attempt to establish a new connection
        return self.connect()
            
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
            'type': 'direct',
            'host': self.config.get('host'),
            'ip_address': self.ip_address,
            'port': self.config.get('port', 5432),
            'database': self.config.get('database', 'postgres'),
            'user': self.config.get('user'),
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
        - DNS resolution status
        
        Returns:
            Tuple of (healthy: bool, health_metrics: Dict[str, Any])
        """
        import time
        
        metrics = {
            'basic_connectivity': False,
            'latency_ms': None,
            'transaction_capability': False,
            'query_capability': False,
            'dns_resolution': {
                'hostname': self.config.get('host'),
                'resolved_ip': self.ip_address,
                'resolution_method': None
            },
            'pool_status': {
                'exists': self.connection_pool is not None,
                'min_connections': self.config.get('min_connections', 1),
                'max_connections': self.config.get('max_connections', 10)
            }
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
                metrics['pool_status']['used_connections'] = len(self.connection_pool._used)
                metrics['pool_status']['unused_connections'] = len(self.connection_pool._unused)
                metrics['pool_status']['total_connections'] = (
                    len(self.connection_pool._used) + len(self.connection_pool._unused)
                )
            except Exception as e:
                metrics['pool_status']['error'] = str(e)
                
        # Check DNS resolution - verify the hostname still resolves to the same IP
        try:
            hostname = self.config.get('host')
            if hostname:
                current_ip = self._resolve_with_socket(hostname)
                
                if current_ip:
                    metrics['dns_resolution']['current_ip'] = current_ip
                    metrics['dns_resolution']['ip_match'] = (current_ip == self.ip_address) if self.ip_address else False
                    metrics['dns_resolution']['resolution_method'] = '_resolve_with_socket'
                else:
                    # Try alternative resolution methods
                    for method_name in ['_resolve_with_nslookup', '_resolve_with_dig', '_resolve_with_alternative_dns']:
                        method = getattr(self, method_name)
                        current_ip = method(hostname)
                        if current_ip:
                            metrics['dns_resolution']['current_ip'] = current_ip
                            metrics['dns_resolution']['ip_match'] = (current_ip == self.ip_address) if self.ip_address else False
                            metrics['dns_resolution']['resolution_method'] = method_name
                            break
                            
                # If we couldn't resolve the hostname at all, that's a problem
                if 'current_ip' not in metrics['dns_resolution']:
                    metrics['dns_resolution']['resolution_failed'] = True
        except Exception as e:
            metrics['dns_resolution']['error'] = str(e)
            
        # Determine overall health
        is_healthy = (
            metrics['basic_connectivity'] and
            metrics['transaction_capability'] and
            metrics['query_capability'] and
            metrics['pool_status']['exists'] and
            ('resolution_failed' not in metrics['dns_resolution'] or not metrics['dns_resolution']['resolution_failed'])
        )
        
        return is_healthy, metrics