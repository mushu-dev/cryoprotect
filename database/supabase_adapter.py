import os
import logging
import socket
import psycopg2
from psycopg2.pool import ThreadedConnectionPool
from psycopg2.extras import RealDictCursor
from typing import Any, Dict, List, Optional, Union, Tuple
import subprocess

from .adapter import DatabaseAdapter

logger = logging.getLogger(__name__)

class SupabaseDirectAdapter(DatabaseAdapter):
    """Supabase direct connection adapter implementation."""
    
    def __init__(self, config: Dict[str, Any]):
        """Initialize adapter with configuration."""
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
        Establish connection pool to Supabase PostgreSQL.
        
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
                logger.error("Missing required Supabase connection parameters")
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
                logger.info(f"Connected to Supabase PostgreSQL at {host}:{port}/{dbname}")
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
                    logger.info(f"Connected to Supabase PostgreSQL at {self.ip_address}:{port}/{dbname} (IP fallback)")
                    return True
                except Exception as e2:
                    logger.error(f"Failed to connect with IP fallback: {str(e2)}")
                    return False
        except Exception as e:
            logger.error(f"Failed to connect to Supabase PostgreSQL: {str(e)}")
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
                logger.info("Disconnected from Supabase PostgreSQL")
            return True
        except Exception as e:
            logger.error(f"Failed to disconnect from Supabase PostgreSQL: {str(e)}")
            return False
            
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
            'type': 'supabase_direct',
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