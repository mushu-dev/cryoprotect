# Database Connection Implementation Plan

## Current Issue Analysis

The project is currently blocked due to persistent Supabase connection issues:

1. **DNS Resolution Failures**: Unable to resolve `db.tsdlmynydfuypiugmkev.supabase.co`
2. **Authentication Failures**: Local PostgreSQL authentication with `postgres` user failing
3. **IP Fallback Issues**: Connection attempts to fallback IP address (142.250.72.78) timing out
4. **MCP Inconsistency**: MCP direct SQL execution works while script connections fail

These issues are preventing database population verification and blocking all subsequent tasks.

## Implementation Strategy

We will implement a flexible database adapter that can seamlessly switch between:
- Local PostgreSQL database (for development)
- Supabase direct connection (for production)
- Supabase MCP connection (as final fallback)

### 1. Database Adapter Pattern

```
├── database/
│   ├── __init__.py
│   ├── adapter.py            # Abstract adapter interface
│   ├── local_adapter.py      # Local PostgreSQL implementation  
│   ├── supabase_adapter.py   # Supabase direct implementation
│   ├── mcp_adapter.py        # MCP fallback implementation
│   ├── connection_manager.py # Manages connection lifecycle
│   └── config.py             # Connection configuration
```

### 2. Connection Manager Implementation

The Connection Manager will:
- Determine which adapter to use based on environment
- Manage connection pool lifecycle
- Handle failover between adapters
- Provide unified query interface
- Implement transaction support

## Implementation Tasks

### Task 1: Create Database Adapter Interface

**File**: `database/adapter.py`

```python
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Union, Tuple

class DatabaseAdapter(ABC):
    """Abstract database adapter interface."""
    
    @abstractmethod
    def connect(self) -> bool:
        """Establish connection to the database."""
        pass
        
    @abstractmethod
    def disconnect(self) -> bool:
        """Close database connection."""
        pass
        
    @abstractmethod
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """Execute SQL query and return results."""
        pass
        
    @abstractmethod
    def execute_batch(self, queries: List[str]) -> List[Any]:
        """Execute multiple SQL queries and return results."""
        pass
        
    @abstractmethod
    def begin_transaction(self) -> Any:
        """Begin a database transaction."""
        pass
        
    @abstractmethod
    def commit_transaction(self, transaction: Any) -> bool:
        """Commit a database transaction."""
        pass
        
    @abstractmethod
    def rollback_transaction(self, transaction: Any) -> bool:
        """Rollback a database transaction."""
        pass
        
    @abstractmethod
    def get_connection_info(self) -> Dict[str, Any]:
        """Get connection information."""
        pass
        
    @abstractmethod
    def test_connection(self) -> Tuple[bool, str]:
        """Test database connection and return status with message."""
        pass
```

### Task 2: Implement Local PostgreSQL Adapter

**File**: `database/local_adapter.py`

```python
import os
import logging
import psycopg2
from psycopg2.pool import ThreadedConnectionPool
from psycopg2.extras import RealDictCursor
from typing import Any, Dict, List, Optional, Union, Tuple

from .adapter import DatabaseAdapter

logger = logging.getLogger(__name__)

class LocalPostgreSQLAdapter(DatabaseAdapter):
    """Local PostgreSQL database adapter implementation."""
    
    def __init__(self, config: Dict[str, Any]):
        """Initialize adapter with configuration."""
        self.config = config
        self.connection_pool = None
        
    def connect(self) -> bool:
        """
        Establish connection pool to local PostgreSQL.
        
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
            logger.info(f"Connected to local PostgreSQL at {host}:{port}/{dbname}")
            return True
        except Exception as e:
            logger.error(f"Failed to connect to local PostgreSQL: {str(e)}")
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
                logger.info("Disconnected from local PostgreSQL")
            return True
        except Exception as e:
            logger.error(f"Failed to disconnect from local PostgreSQL: {str(e)}")
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
        try:
            result = self.execute_query("SELECT 1 as test")
            if result and result[0]['test'] == 1:
                return True, "Connection successful"
            else:
                return False, "Connection test failed"
        except Exception as e:
            return False, f"Connection error: {str(e)}"
```

### Task 3: Implement Supabase Direct Adapter

**File**: `database/supabase_adapter.py`

```python
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
```

### Task 4: Implement MCP Adapter

**File**: `database/mcp_adapter.py`

```python
import os
import logging
import json
from typing import Any, Dict, List, Optional, Union, Tuple

from .adapter import DatabaseAdapter

logger = logging.getLogger(__name__)

class MCPAdapter(DatabaseAdapter):
    """MCP-based database adapter implementation."""
    
    def __init__(self, config: Dict[str, Any]):
        """Initialize adapter with configuration."""
        self.config = config
        self.project_id = config.get('project_id')
        self.connected = False
        
        # Import MCP functions
        try:
            from use_mcp_tool import execute_sql, get_project_id
            self.execute_sql_through_mcp = execute_sql
            self.get_project_id_for_mcp = get_project_id
            
            if not self.project_id:
                self.project_id = self.get_project_id_for_mcp()
        except ImportError:
            logger.error("Failed to import MCP tools. MCP adapter will not work.")
            self.execute_sql_through_mcp = None
            self.get_project_id_for_mcp = None
        
    def connect(self) -> bool:
        """
        Establish connection to MCP.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        try:
            if not self.execute_sql_through_mcp or not self.project_id:
                logger.error("MCP tools not available or project ID not set")
                return False
                
            # Test connection
            result = self.execute_sql_through_mcp("SELECT 1 as test", self.project_id)
            if isinstance(result, list) and len(result) > 0 and result[0].get('test') == 1:
                self.connected = True
                logger.info(f"Connected to Supabase via MCP (project_id: {self.project_id})")
                return True
            else:
                logger.error(f"Failed to connect to Supabase via MCP: {result}")
                return False
        except Exception as e:
            logger.error(f"Failed to connect to Supabase via MCP: {str(e)}")
            return False
            
    def disconnect(self) -> bool:
        """
        Close MCP connection (no-op since MCP is stateless).
        
        Returns:
            bool: Always True
        """
        self.connected = False
        return True
            
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """
        Execute SQL query through MCP and return results.
        
        Args:
            query: SQL query to execute
            params: Query parameters (not used, must be embedded in query)
            
        Returns:
            Query results
        """
        if not self.execute_sql_through_mcp or not self.project_id:
            raise ValueError("MCP tools not available or project ID not set")
            
        try:
            # Handle params by formatting query (not ideal but MCP doesn't support params)
            if params:
                if isinstance(params, dict):
                    # Convert dict params to SQL params
                    for key, value in params.items():
                        placeholder = f"%({key})s"
                        if isinstance(value, str):
                            formatted_value = f"'{value.replace("'", "''")}'"
                        elif value is None:
                            formatted_value = "NULL"
                        elif isinstance(value, bool):
                            formatted_value = "TRUE" if value else "FALSE"
                        else:
                            formatted_value = str(value)
                        query = query.replace(placeholder, formatted_value)
                else:
                    # Convert tuple params to SQL params
                    for i, value in enumerate(params):
                        placeholder = f"%s"
                        if isinstance(value, str):
                            formatted_value = f"'{value.replace("'", "''")}'"
                        elif value is None:
                            formatted_value = "NULL"
                        elif isinstance(value, bool):
                            formatted_value = "TRUE" if value else "FALSE"
                        else:
                            formatted_value = str(value)
                        query = query.replace(placeholder, formatted_value, 1)
            
            result = self.execute_sql_through_mcp(query, self.project_id)
            return result
        except Exception as e:
            logger.error(f"Error executing query via MCP: {str(e)}")
            raise
                
    def execute_batch(self, queries: List[str]) -> List[Any]:
        """
        Execute multiple SQL queries through MCP and return results.
        
        Args:
            queries: List of SQL queries to execute
            
        Returns:
            List of query results
        """
        if not self.execute_sql_through_mcp or not self.project_id:
            raise ValueError("MCP tools not available or project ID not set")
            
        try:
            # Combine queries into a transaction
            transaction_queries = ["BEGIN;"] + queries + ["COMMIT;"]
            combined_query = "\n".join(transaction_queries)
            
            result = self.execute_sql_through_mcp(combined_query, self.project_id)
            
            # Parse results (this is approximate since MCP returns a combined result)
            if isinstance(result, list) and len(result) > 0:
                return [result]  # Return full result as single batch result
            return [None]
        except Exception as e:
            logger.error(f"Error executing batch via MCP: {str(e)}")
            raise
                
    def begin_transaction(self) -> Any:
        """
        Begin a database transaction.
        
        Returns:
            List to collect transaction queries
        """
        return []
            
    def commit_transaction(self, transaction: List[str]) -> bool:
        """
        Commit a database transaction.
        
        Args:
            transaction: List of queries to execute in transaction
            
        Returns:
            bool: True if commit successful, False otherwise
        """
        try:
            # Combine queries into a transaction
            transaction_queries = ["BEGIN;"] + transaction + ["COMMIT;"]
            combined_query = "\n".join(transaction_queries)
            
            self.execute_sql_through_mcp(combined_query, self.project_id)
            return True
        except Exception as e:
            logger.error(f"Error committing transaction via MCP: {str(e)}")
            return False
            
    def rollback_transaction(self, transaction: List[str]) -> bool:
        """
        Rollback a database transaction (no-op for MCP).
        
        Args:
            transaction: Transaction object (ignored)
            
        Returns:
            bool: Always True
        """
        # No need to do anything, just clear the transaction list
        transaction.clear()
        return True
            
    def get_connection_info(self) -> Dict[str, Any]:
        """
        Get connection information.
        
        Returns:
            Dict with connection information
        """
        return {
            'type': 'mcp',
            'project_id': self.project_id,
            'connected': self.connected
        }
        
    def test_connection(self) -> Tuple[bool, str]:
        """
        Test database connection and return status with message.
        
        Returns:
            Tuple of (success: bool, message: str)
        """
        try:
            if not self.execute_sql_through_mcp or not self.project_id:
                return False, "MCP tools not available or project ID not set"
                
            result = self.execute_sql_through_mcp("SELECT 1 as test", self.project_id)
            if isinstance(result, list) and len(result) > 0 and result[0].get('test') == 1:
                return True, "Connection successful"
            else:
                return False, f"Connection test failed: {result}"
        except Exception as e:
            return False, f"Connection error: {str(e)}"
```

### Task 5: Implement Connection Manager

**File**: `database/connection_manager.py`

```python
import os
import logging
from typing import Any, Dict, List, Optional, Union, Tuple
from dotenv import load_dotenv

from .adapter import DatabaseAdapter
from .local_adapter import LocalPostgreSQLAdapter
from .supabase_adapter import SupabaseDirectAdapter
from .mcp_adapter import MCPAdapter

logger = logging.getLogger(__name__)

class ConnectionManager:
    """Database connection manager."""
    
    _instance = None
    
    @classmethod
    def get_instance(cls) -> 'ConnectionManager':
        """Get singleton instance of ConnectionManager."""
        if cls._instance is None:
            cls._instance = ConnectionManager()
        return cls._instance
    
    def __init__(self):
        """Initialize connection manager."""
        # Load environment variables
        load_dotenv()
        
        self.adapters = {}
        self.primary_adapter = None
        self.fallback_adapter = None
        self.active_adapter = None
        
        # Initialize configuration
        self.init_configuration()
        
    def init_configuration(self):
        """Initialize configuration from environment variables."""
        # Load environment variables
        self.db_config = self._load_db_config()
        
        # Determine connection mode
        self.connection_mode = os.getenv('DB_CONNECTION_MODE', 'auto').lower()
        
        # Create adapters based on configuration
        if self.connection_mode == 'local' or self.connection_mode == 'auto':
            self.adapters['local'] = LocalPostgreSQLAdapter(self.db_config.get('local', {}))
            
        if self.connection_mode == 'supabase' or self.connection_mode == 'auto':
            self.adapters['supabase'] = SupabaseDirectAdapter(self.db_config.get('supabase', {}))
            
        if self.connection_mode == 'mcp' or self.connection_mode == 'auto':
            self.adapters['mcp'] = MCPAdapter(self.db_config.get('mcp', {}))
            
        # Set primary and fallback adapters
        if self.connection_mode == 'auto':
            self.primary_adapter = 'supabase'
            self.fallback_adapter = 'local'
            self.last_fallback = 'mcp'
        else:
            self.primary_adapter = self.connection_mode
            self.fallback_adapter = None
            self.last_fallback = None
            
    def _load_db_config(self) -> Dict[str, Dict[str, Any]]:
        """
        Load database configuration from environment variables.
        
        Returns:
            Dict with configuration for all adapter types
        """
        # Supabase configuration
        supabase_config = {
            'host': os.getenv('SUPABASE_DB_HOST') or os.getenv('DB_HOST'),
            'port': os.getenv('SUPABASE_DB_PORT') or os.getenv('DB_PORT', '5432'),
            'database': os.getenv('SUPABASE_DB_NAME') or os.getenv('DB_NAME', 'postgres'),
            'user': os.getenv('SUPABASE_DB_USER') or os.getenv('DB_USER'),
            'password': os.getenv('SUPABASE_DB_PASSWORD') or os.getenv('DB_PASSWORD'),
            'min_connections': os.getenv('SUPABASE_DB_MIN_CONNECTIONS', '1'),
            'max_connections': os.getenv('SUPABASE_DB_MAX_CONNECTIONS', '10')
        }
        
        # Local configuration
        local_config = {
            'host': os.getenv('LOCAL_DB_HOST', 'localhost'),
            'port': os.getenv('LOCAL_DB_PORT', '5432'),
            'database': os.getenv('LOCAL_DB_NAME', 'postgres'),
            'user': os.getenv('LOCAL_DB_USER', 'postgres'),
            'password': os.getenv('LOCAL_DB_PASSWORD', ''),
            'min_connections': os.getenv('LOCAL_DB_MIN_CONNECTIONS', '1'),
            'max_connections': os.getenv('LOCAL_DB_MAX_CONNECTIONS', '5')
        }
        
        # MCP configuration
        mcp_config = {
            'project_id': os.getenv('SUPABASE_PROJECT_ID'),
        }
        
        return {
            'supabase': supabase_config,
            'local': local_config,
            'mcp': mcp_config
        }
        
    def connect(self) -> bool:
        """
        Connect to database using primary adapter or fallbacks.
        
        Returns:
            bool: True if connected successfully, False otherwise
        """
        # Try primary adapter first
        adapter_order = [self.primary_adapter]
        
        # Add fallbacks if in auto mode
        if self.connection_mode == 'auto' and self.fallback_adapter:
            adapter_order.append(self.fallback_adapter)
            
        if self.connection_mode == 'auto' and self.last_fallback:
            adapter_order.append(self.last_fallback)
            
        # Try adapters in order
        for adapter_name in adapter_order:
            if adapter_name not in self.adapters:
                continue
                
            adapter = self.adapters[adapter_name]
            try:
                if adapter.connect():
                    self.active_adapter = adapter_name
                    logger.info(f"Connected using {adapter_name} adapter")
                    return True
            except Exception as e:
                logger.warning(f"Failed to connect with {adapter_name} adapter: {str(e)}")
                
        logger.error("Failed to connect to database with any adapter")
        return False
        
    def disconnect(self) -> bool:
        """
        Disconnect from database.
        
        Returns:
            bool: True if disconnection successful, False otherwise
        """
        success = True
        
        # Disconnect all adapters
        for name, adapter in self.adapters.items():
            try:
                adapter.disconnect()
            except Exception as e:
                logger.error(f"Error disconnecting {name} adapter: {str(e)}")
                success = False
                
        self.active_adapter = None
        return success
        
    def get_active_adapter(self) -> Optional[DatabaseAdapter]:
        """
        Get the currently active adapter.
        
        Returns:
            Active DatabaseAdapter or None if not connected
        """
        if not self.active_adapter:
            return None
        return self.adapters.get(self.active_adapter)
        
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """
        Execute SQL query using the active adapter.
        
        Args:
            query: SQL query to execute
            params: Query parameters
            
        Returns:
            Query results
        """
        adapter = self.get_active_adapter()
        if not adapter:
            if not self.connect():
                raise ConnectionError("Not connected to database")
            adapter = self.get_active_adapter()
            
        return adapter.execute_query(query, params)
        
    def execute_batch(self, queries: List[str]) -> List[Any]:
        """
        Execute multiple SQL queries using the active adapter.
        
        Args:
            queries: List of SQL queries to execute
            
        Returns:
            List of query results
        """
        adapter = self.get_active_adapter()
        if not adapter:
            if not self.connect():
                raise ConnectionError("Not connected to database")
            adapter = self.get_active_adapter()
            
        return adapter.execute_batch(queries)
        
    def begin_transaction(self) -> Any:
        """
        Begin a database transaction.
        
        Returns:
            Transaction object
        """
        adapter = self.get_active_adapter()
        if not adapter:
            if not self.connect():
                raise ConnectionError("Not connected to database")
            adapter = self.get_active_adapter()
            
        return adapter.begin_transaction()
        
    def commit_transaction(self, transaction: Any) -> bool:
        """
        Commit a database transaction.
        
        Args:
            transaction: Transaction object
            
        Returns:
            bool: True if commit successful, False otherwise
        """
        adapter = self.get_active_adapter()
        if not adapter:
            raise ConnectionError("Not connected to database")
            
        return adapter.commit_transaction(transaction)
        
    def rollback_transaction(self, transaction: Any) -> bool:
        """
        Rollback a database transaction.
        
        Args:
            transaction: Transaction object
            
        Returns:
            bool: True if rollback successful, False otherwise
        """
        adapter = self.get_active_adapter()
        if not adapter:
            raise ConnectionError("Not connected to database")
            
        return adapter.rollback_transaction(transaction)
        
    def test_all_connections(self) -> Dict[str, Tuple[bool, str]]:
        """
        Test all configured database connections.
        
        Returns:
            Dict mapping adapter names to (success, message) tuples
        """
        results = {}
        
        for name, adapter in self.adapters.items():
            try:
                logger.info(f"Testing {name} adapter connection...")
                
                # Try to connect if not already connected
                if name != self.active_adapter:
                    try:
                        adapter.connect()
                    except:
                        pass
                        
                success, message = adapter.test_connection()
                results[name] = (success, message)
                
                logger.info(f"Connection test for {name} adapter: {'SUCCESS' if success else 'FAILED'} - {message}")
                
                # Disconnect if we just connected for testing
                if name != self.active_adapter:
                    try:
                        adapter.disconnect()
                    except:
                        pass
            except Exception as e:
                results[name] = (False, f"Error testing connection: {str(e)}")
                logger.error(f"Error testing {name} adapter: {str(e)}")
                
        return results
        
    def get_connection_info(self) -> Dict[str, Any]:
        """
        Get connection information for all adapters.
        
        Returns:
            Dict with connection information
        """
        info = {
            'connection_mode': self.connection_mode,
            'primary_adapter': self.primary_adapter,
            'fallback_adapter': self.fallback_adapter,
            'active_adapter': self.active_adapter,
            'adapters': {}
        }
        
        for name, adapter in self.adapters.items():
            try:
                info['adapters'][name] = adapter.get_connection_info()
            except Exception as e:
                info['adapters'][name] = {'error': str(e)}
                
        return info
```

### Task 6: Create Database Initialization Script

**File**: `database/init_local_db.py`

```python
#!/usr/bin/env python3
"""
Initialize local PostgreSQL database for development.

This script:
1. Creates required database and schema
2. Applies all migrations
3. Creates test data for development
"""

import os
import sys
import logging
import argparse
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from dotenv import load_dotenv
import glob
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Initialize local PostgreSQL database')
    parser.add_argument('--reset', action='store_true', help='Reset database if it exists')
    parser.add_argument('--skip-migrations', action='store_true', help='Skip applying migrations')
    parser.add_argument('--skip-test-data', action='store_true', help='Skip creating test data')
    return parser.parse_args()

def get_connection_params():
    """Get database connection parameters from environment."""
    load_dotenv()
    
    # Get connection parameters
    params = {
        'host': os.getenv('LOCAL_DB_HOST', 'localhost'),
        'port': os.getenv('LOCAL_DB_PORT', '5432'),
        'user': os.getenv('LOCAL_DB_USER', 'postgres'),
        'password': os.getenv('LOCAL_DB_PASSWORD', '')
    }
    
    # Database name
    db_name = os.getenv('LOCAL_DB_NAME', 'cryoprotect')
    
    return params, db_name

def create_database(params, db_name, reset=False):
    """
    Create database if it doesn't exist.
    
    Args:
        params: Connection parameters
        db_name: Database name
        reset: Whether to reset database if it exists
        
    Returns:
        bool: True if database created, False if it already exists
    """
    conn = None
    try:
        # Connect to PostgreSQL server
        conn = psycopg2.connect(**params)
        conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        
        with conn.cursor() as cursor:
            # Check if database exists
            cursor.execute("SELECT 1 FROM pg_database WHERE datname = %s", (db_name,))
            exists = cursor.fetchone() is not None
            
            if exists:
                if reset:
                    logger.info(f"Resetting database '{db_name}'...")
                    # Terminate all connections
                    cursor.execute(f"SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE datname = '{db_name}'")
                    # Drop database
                    cursor.execute(f"DROP DATABASE {db_name}")
                    logger.info(f"Database '{db_name}' dropped")
                else:
                    logger.info(f"Database '{db_name}' already exists")
                    return False
            
            # Create database
            cursor.execute(f"CREATE DATABASE {db_name}")
            logger.info(f"Database '{db_name}' created")
            
            return True
    except Exception as e:
        logger.error(f"Error creating database: {str(e)}")
        return False
    finally:
        if conn:
            conn.close()

def apply_migrations(params, db_name):
    """
    Apply all migrations in order.
    
    Args:
        params: Connection parameters
        db_name: Database name
        
    Returns:
        bool: True if migrations applied successfully, False otherwise
    """
    conn = None
    try:
        # Connect to database
        conn_params = params.copy()
        conn_params['dbname'] = db_name
        conn = psycopg2.connect(**conn_params)
        
        # Get migration files
        migrations_dir = os.path.join(os.getcwd(), 'migrations')
        migration_files = glob.glob(os.path.join(migrations_dir, '*.sql'))
        
        # Sort migration files by numeric prefix
        def get_migration_number(file_path):
            match = re.search(r'^(\d+)', os.path.basename(file_path))
            return int(match.group(1)) if match else 0
            
        migration_files.sort(key=get_migration_number)
        
        if not migration_files:
            logger.warning("No migration files found")
            return True
            
        # Apply migrations
        for migration_file in migration_files:
            logger.info(f"Applying migration: {os.path.basename(migration_file)}")
            
            with conn.cursor() as cursor:
                # Create migrations table if it doesn't exist
                cursor.execute("""
                    CREATE TABLE IF NOT EXISTS migrations (
                        id SERIAL PRIMARY KEY,
                        filename VARCHAR(255) NOT NULL,
                        applied_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
                    )
                """)
                
                # Check if migration already applied
                cursor.execute("SELECT 1 FROM migrations WHERE filename = %s", (os.path.basename(migration_file),))
                if cursor.fetchone():
                    logger.info(f"Migration {os.path.basename(migration_file)} already applied")
                    continue
                    
                # Read and apply migration
                with open(migration_file, 'r') as f:
                    sql = f.read()
                    cursor.execute(sql)
                    
                # Record migration
                cursor.execute("INSERT INTO migrations (filename) VALUES (%s)", (os.path.basename(migration_file),))
                
                # Commit transaction
                conn.commit()
                
        logger.info("All migrations applied successfully")
        return True
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error applying migrations: {str(e)}")
        return False
    finally:
        if conn:
            conn.close()

def create_test_data(params, db_name):
    """
    Create test data for development.
    
    Args:
        params: Connection parameters
        db_name: Database name
        
    Returns:
        bool: True if test data created successfully, False otherwise
    """
    conn = None
    try:
        # Connect to database
        conn_params = params.copy()
        conn_params['dbname'] = db_name
        conn = psycopg2.connect(**conn_params)
        
        logger.info("Creating test data...")
        
        with conn.cursor() as cursor:
            # Create test user
            cursor.execute("""
                INSERT INTO auth.users (id, email, encrypted_password, email_confirmed_at, created_at, updated_at)
                VALUES (
                    '00000000-0000-0000-0000-000000000000',
                    'test@example.com',
                    '$2a$10$abcdefghijklmnopqrstuvwxyz',
                    NOW(),
                    NOW(),
                    NOW()
                )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test user profile
            cursor.execute("""
                INSERT INTO user_profile (id, auth_user_id, display_name, email, affiliation, created_at, updated_at)
                VALUES (
                    '00000000-0000-0000-0000-000000000001',
                    '00000000-0000-0000-0000-000000000000',
                    'Test User',
                    'test@example.com',
                    'Test Organization',
                    NOW(),
                    NOW()
                )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test molecules
            cursor.execute("""
                INSERT INTO molecules (
                    id, name, formula, molecular_weight, smiles, inchi, inchi_key,
                    chembl_id, pubchem_cid, data_source, created_at, updated_at
                )
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000010',
                        'Dimethyl sulfoxide',
                        'C2H6OS',
                        78.13,
                        'CS(=O)C',
                        'InChI=1S/C2H6OS/c1-4(2)3/h1-2H3',
                        'IAZDPXIOMUYVGZ-UHFFFAOYSA-N',
                        'CHEMBL422',
                        '679',
                        'test_data',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000011',
                        'Glycerol',
                        'C3H8O3',
                        92.09,
                        'C(C(CO)O)O',
                        'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2',
                        'PEDCQBHIVMGVHV-UHFFFAOYSA-N',
                        'CHEMBL692',
                        '753',
                        'test_data',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000012',
                        'Trehalose',
                        'C12H22O11',
                        342.30,
                        'C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O',
                        'InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12-10(19)6(15)5(14)3(2-20)22-12/h3-19H,1-2H2/t3-,4+,5-,6-,7+,8+,9-,10-,11-,12+/m1/s1',
                        'HBVRQNNGHBPKND-DZGCQCFKSA-N',
                        'CHEMBL15132',
                        '7427',
                        'test_data',
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test property types
            cursor.execute("""
                INSERT INTO property_types (id, name, description, data_type, unit, created_at, updated_at)
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000020',
                        'logP',
                        'Octanol-water partition coefficient',
                        'numeric',
                        '',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000021',
                        'molecular_weight',
                        'Molecular weight',
                        'numeric',
                        'g/mol',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000022',
                        'h_bond_donors',
                        'Number of hydrogen bond donors',
                        'numeric',
                        '',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000023',
                        'h_bond_acceptors',
                        'Number of hydrogen bond acceptors',
                        'numeric',
                        '',
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test molecular properties
            cursor.execute("""
                INSERT INTO molecular_properties (
                    id, molecule_id, property_type_id, value, source, confidence, created_at, updated_at
                )
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000030',
                        '00000000-0000-0000-0000-000000000010',
                        '00000000-0000-0000-0000-000000000020',
                        -1.35,
                        'test_data',
                        0.95,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000031',
                        '00000000-0000-0000-0000-000000000010',
                        '00000000-0000-0000-0000-000000000021',
                        78.13,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000032',
                        '00000000-0000-0000-0000-000000000010',
                        '00000000-0000-0000-0000-000000000022',
                        0,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000033',
                        '00000000-0000-0000-0000-000000000010',
                        '00000000-0000-0000-0000-000000000023',
                        1,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000034',
                        '00000000-0000-0000-0000-000000000011',
                        '00000000-0000-0000-0000-000000000020',
                        -1.76,
                        'test_data',
                        0.9,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000035',
                        '00000000-0000-0000-0000-000000000011',
                        '00000000-0000-0000-0000-000000000021',
                        92.09,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000036',
                        '00000000-0000-0000-0000-000000000012',
                        '00000000-0000-0000-0000-000000000020',
                        -3.77,
                        'test_data',
                        0.85,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000037',
                        '00000000-0000-0000-0000-000000000012',
                        '00000000-0000-0000-0000-000000000021',
                        342.30,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test mixtures
            cursor.execute("""
                INSERT INTO mixtures (id, name, description, creator_id, created_at, updated_at)
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000040',
                        'DMSO-Glycerol Mix',
                        'A test mixture of DMSO and glycerol',
                        '00000000-0000-0000-0000-000000000001',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000041',
                        'Trehalose Solution',
                        'A test mixture with trehalose',
                        '00000000-0000-0000-0000-000000000001',
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test mixture components
            cursor.execute("""
                INSERT INTO mixture_components (
                    id, mixture_id, molecule_id, concentration, concentration_unit, created_at, updated_at
                )
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000050',
                        '00000000-0000-0000-0000-000000000040',
                        '00000000-0000-0000-0000-000000000010',
                        10.0,
                        'percent_v_v',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000051',
                        '00000000-0000-0000-0000-000000000040',
                        '00000000-0000-0000-0000-000000000011',
                        5.0,
                        'percent_v_v',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000052',
                        '00000000-0000-0000-0000-000000000041',
                        '00000000-0000-0000-0000-000000000012',
                        0.3,
                        'molar',
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Commit transaction
            conn.commit()
            
        logger.info("Test data created successfully")
        return True
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error creating test data: {str(e)}")
        return False
    finally:
        if conn:
            conn.close()

def main():
    """Main function."""
    args = parse_args()
    
    # Get connection parameters
    params, db_name = get_connection_params()
    
    # Create database
    created = create_database(params, db_name, args.reset)
    
    # Apply migrations
    if not args.skip_migrations:
        if not apply_migrations(params, db_name):
            logger.error("Failed to apply migrations")
            return 1
            
    # Create test data
    if not args.skip_test_data:
        if not create_test_data(params, db_name):
            logger.error("Failed to create test data")
            return 1
            
    logger.info("Database initialization completed successfully")
    return 0

if __name__ == '__main__':
    sys.exit(main())
```

### Task 7: Create Database Helper Functions

**File**: `database/utils.py`

```python
#!/usr/bin/env python3
"""
Database utility functions for CryoProtect v2.
"""

import os
import logging
from typing import Any, Dict, List, Optional, Union, Tuple, Callable
import time
import functools

from ..database.connection_manager import ConnectionManager

logger = logging.getLogger(__name__)

def get_db():
    """
    Get database connection manager.
    
    Returns:
        ConnectionManager: Database connection manager
    """
    return ConnectionManager.get_instance()

def with_connection(f: Callable) -> Callable:
    """
    Decorator to ensure database connection is established before function execution.
    
    Args:
        f: Function to decorate
        
    Returns:
        Decorated function
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        db = get_db()
        
        # Ensure connection
        if not db.get_active_adapter():
            if not db.connect():
                raise ConnectionError("Failed to connect to database")
                
        return f(*args, **kwargs)
    return wrapper

def with_retry(max_retries: int = 3, backoff: float = 1.5) -> Callable:
    """
    Decorator to retry database operations with exponential backoff.
    
    Args:
        max_retries: Maximum number of retry attempts
        backoff: Backoff factor for retry delays
        
    Returns:
        Decorator function
    """
    def decorator(f: Callable) -> Callable:
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            last_exception = None
            retry_count = 0
            
            while retry_count <= max_retries:
                try:
                    return f(*args, **kwargs)
                except Exception as e:
                    retry_count += 1
                    last_exception = e
                    
                    # Check if we should retry
                    if retry_count > max_retries:
                        break
                        
                    # Calculate backoff time
                    delay = backoff ** (retry_count - 1)
                    
                    logger.warning(
                        f"Retrying database operation after error: {str(e)}. "
                        f"Retry {retry_count}/{max_retries} in {delay:.2f}s..."
                    )
                    
                    time.sleep(delay)
                    
            # Re-raise the last exception
            raise last_exception
        return wrapper
    return decorator

def with_transaction(f: Callable) -> Callable:
    """
    Decorator to execute function in a database transaction.
    
    Args:
        f: Function to decorate
        
    Returns:
        Decorated function
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        db = get_db()
        
        # Ensure connection
        if not db.get_active_adapter():
            if not db.connect():
                raise ConnectionError("Failed to connect to database")
                
        # Begin transaction
        transaction = db.begin_transaction()
        
        try:
            # Execute function
            result = f(*args, transaction=transaction, **kwargs)
            
            # Commit transaction
            db.commit_transaction(transaction)
            
            return result
        except Exception as e:
            # Rollback transaction
            db.rollback_transaction(transaction)
            raise
    return wrapper

@with_connection
def execute_query(query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
    """
    Execute SQL query.
    
    Args:
        query: SQL query to execute
        params: Query parameters
        
    Returns:
        Query results
    """
    db = get_db()
    return db.execute_query(query, params)

@with_connection
def execute_batch(queries: List[str]) -> List[Any]:
    """
    Execute multiple SQL queries.
    
    Args:
        queries: List of SQL queries to execute
        
    Returns:
        List of query results
    """
    db = get_db()
    return db.execute_batch(queries)

@with_connection
def get_molecule_by_id(molecule_id: str) -> Optional[Dict[str, Any]]:
    """
    Get molecule by ID.
    
    Args:
        molecule_id: Molecule ID
        
    Returns:
        Molecule data or None if not found
    """
    db = get_db()
    results = db.execute_query(
        "SELECT * FROM molecules WHERE id = %s",
        (molecule_id,)
    )
    
    if results and len(results) > 0:
        return results[0]
    return None

@with_connection
def get_molecule_properties(molecule_id: str) -> List[Dict[str, Any]]:
    """
    Get properties for a molecule.
    
    Args:
        molecule_id: Molecule ID
        
    Returns:
        List of property data
    """
    db = get_db()
    return db.execute_query(
        """
        SELECT mp.*, pt.name as property_name, pt.unit 
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE mp.molecule_id = %s
        """,
        (molecule_id,)
    )

@with_connection
def get_molecules_by_inchikey(inchi_key: str) -> List[Dict[str, Any]]:
    """
    Get molecules by InChI Key.
    
    Args:
        inchi_key: InChI Key
        
    Returns:
        List of molecule data
    """
    db = get_db()
    return db.execute_query(
        "SELECT * FROM molecules WHERE inchi_key = %s",
        (inchi_key,)
    )

@with_connection
def insert_molecule(
    name: str,
    formula: Optional[str] = None,
    molecular_weight: Optional[float] = None,
    smiles: Optional[str] = None,
    inchi: Optional[str] = None,
    inchi_key: Optional[str] = None,
    chembl_id: Optional[str] = None,
    pubchem_cid: Optional[str] = None,
    data_source: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Insert new molecule.
    
    Args:
        name: Molecule name
        formula: Molecular formula
        molecular_weight: Molecular weight
        smiles: SMILES string
        inchi: InChI string
        inchi_key: InChI Key
        chembl_id: ChEMBL ID
        pubchem_cid: PubChem CID
        data_source: Data source
        
    Returns:
        Inserted molecule data
    """
    db = get_db()
    result = db.execute_query(
        """
        INSERT INTO molecules (
            name, formula, molecular_weight, smiles, inchi, inchi_key,
            chembl_id, pubchem_cid, data_source, created_at, updated_at
        )
        VALUES (
            %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW()
        )
        RETURNING *
        """,
        (
            name, formula, molecular_weight, smiles, inchi, inchi_key,
            chembl_id, pubchem_cid, data_source
        )
    )
    
    if result and len(result) > 0:
        return result[0]
    return None

@with_connection
def update_molecule(
    molecule_id: str,
    name: Optional[str] = None,
    formula: Optional[str] = None,
    molecular_weight: Optional[float] = None,
    smiles: Optional[str] = None,
    inchi: Optional[str] = None,
    inchi_key: Optional[str] = None,
    chembl_id: Optional[str] = None,
    pubchem_cid: Optional[str] = None,
    data_source: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Update molecule.
    
    Args:
        molecule_id: Molecule ID
        name: Molecule name
        formula: Molecular formula
        molecular_weight: Molecular weight
        smiles: SMILES string
        inchi: InChI string
        inchi_key: InChI Key
        chembl_id: ChEMBL ID
        pubchem_cid: PubChem CID
        data_source: Data source
        
    Returns:
        Updated molecule data
    """
    # Build update fields
    update_fields = []
    params = []
    
    if name is not None:
        update_fields.append("name = %s")
        params.append(name)
        
    if formula is not None:
        update_fields.append("formula = %s")
        params.append(formula)
        
    if molecular_weight is not None:
        update_fields.append("molecular_weight = %s")
        params.append(molecular_weight)
        
    if smiles is not None:
        update_fields.append("smiles = %s")
        params.append(smiles)
        
    if inchi is not None:
        update_fields.append("inchi = %s")
        params.append(inchi)
        
    if inchi_key is not None:
        update_fields.append("inchi_key = %s")
        params.append(inchi_key)
        
    if chembl_id is not None:
        update_fields.append("chembl_id = %s")
        params.append(chembl_id)
        
    if pubchem_cid is not None:
        update_fields.append("pubchem_cid = %s")
        params.append(pubchem_cid)
        
    if data_source is not None:
        update_fields.append("data_source = %s")
        params.append(data_source)
        
    # Add updated_at field
    update_fields.append("updated_at = NOW()")
    
    # If no fields to update, return current molecule
    if not update_fields:
        return get_molecule_by_id(molecule_id)
        
    # Build query
    query = f"""
        UPDATE molecules
        SET {", ".join(update_fields)}
        WHERE id = %s
        RETURNING *
    """
    
    # Add molecule_id to params
    params.append(molecule_id)
    
    # Execute query
    db = get_db()
    result = db.execute_query(query, tuple(params))
    
    if result and len(result) > 0:
        return result[0]
    return None

@with_connection
def set_molecule_property(
    molecule_id: str,
    property_type_id: str,
    value: Union[float, str, bool],
    source: Optional[str] = None,
    confidence: Optional[float] = None
) -> Optional[Dict[str, Any]]:
    """
    Set molecule property.
    
    Args:
        molecule_id: Molecule ID
        property_type_id: Property type ID
        value: Property value
        source: Data source
        confidence: Confidence value
        
    Returns:
        Property data
    """
    db = get_db()
    
    # Check if property exists
    existing = db.execute_query(
        """
        SELECT * FROM molecular_properties
        WHERE molecule_id = %s AND property_type_id = %s
        """,
        (molecule_id, property_type_id)
    )
    
    if existing and len(existing) > 0:
        # Update existing property
        result = db.execute_query(
            """
            UPDATE molecular_properties
            SET value = %s, source = %s, confidence = %s, updated_at = NOW()
            WHERE molecule_id = %s AND property_type_id = %s
            RETURNING *
            """,
            (value, source, confidence, molecule_id, property_type_id)
        )
    else:
        # Insert new property
        result = db.execute_query(
            """
            INSERT INTO molecular_properties
            (molecule_id, property_type_id, value, source, confidence, created_at, updated_at)
            VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
            RETURNING *
            """,
            (molecule_id, property_type_id, value, source, confidence)
        )
    
    if result and len(result) > 0:
        return result[0]
    return None

@with_connection
def get_or_create_property_type(
    name: str,
    description: Optional[str] = None,
    data_type: str = 'numeric',
    unit: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Get or create property type.
    
    Args:
        name: Property type name
        description: Property type description
        data_type: Data type (numeric, text, boolean)
        unit: Unit of measurement
        
    Returns:
        Property type data
    """
    db = get_db()
    
    # Check if property type exists
    existing = db.execute_query(
        "SELECT * FROM property_types WHERE name = %s",
        (name,)
    )
    
    if existing and len(existing) > 0:
        return existing[0]
        
    # Create new property type
    result = db.execute_query(
        """
        INSERT INTO property_types
        (name, description, data_type, unit, created_at, updated_at)
        VALUES (%s, %s, %s, %s, NOW(), NOW())
        RETURNING *
        """,
        (name, description, data_type, unit)
    )
    
    if result and len(result) > 0:
        return result[0]
    return None

@with_connection
def test_database_connection() -> Dict[str, Any]:
    """
    Test database connection and return detailed status.
    
    Returns:
        Dict with connection status and test results
    """
    db = get_db()
    adapter_info = db.get_connection_info()
    test_results = db.test_all_connections()
    
    # Try simple query
    query_result = None
    query_error = None
    
    try:
        result = db.execute_query("SELECT 1 as test")
        if result and len(result) > 0:
            query_result = "Success"
    except Exception as e:
        query_error = str(e)
    
    return {
        'connection_info': adapter_info,
        'test_results': test_results,
        'query_test': {
            'result': query_result,
            'error': query_error
        }
    }
```

## Implementation Plan

### Week 1: Core Implementation

Day 1-2: Implement adapter pattern
- Create abstract DatabaseAdapter interface
- Implement LocalPostgreSQLAdapter
- Implement SupabaseDirectAdapter
- Implement MCPAdapter

Day 3-4: Implement connection management
- Create ConnectionManager
- Implement adapter selection logic
- Add fallback mechanisms
- Test all connection methods

Day 5: Implement database utilities
- Create helper functions for common operations
- Add transaction support
- Implement retry mechanisms
- Create diagnostic tools

### Week 2: Integration & Testing

Day 1-2: Update database population scripts
- Modify ChEMBL import script to use new adapter
- Update PubChem import script to use new adapter
- Enhance verification script with better diagnostics
- Test with local PostgreSQL database

Day 3-4: Create database initialization tools
- Implement local database initialization script
- Create migration application tool
- Add test data generation
- Document setup process

Day 5: Comprehensive testing
- Test connection fallback mechanisms
- Verify data population with all adapters
- Measure performance improvements
- Create validation report

## Success Criteria

1. **Connection Reliability**
   - Zero connection failures during normal operation
   - Automatic failover between connection methods
   - Clear error reporting for diagnostic purposes
   - Self-healing connection management

2. **Local Development Support**
   - Fully functional local PostgreSQL environment
   - Identical schema and behavior to production
   - Easy setup with initialization script
   - Comprehensive test data

3. **Performance Improvements**
   - 50%+ faster database operations compared to MCP
   - Efficient connection pooling
   - Transaction support for bulk operations
   - Optimized query execution

4. **Code Quality**
   - Clear abstraction via adapter pattern
   - Comprehensive error handling
   - Well-documented interfaces
   - High test coverage