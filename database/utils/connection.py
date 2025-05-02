"""
Database connection utilities.

This module provides connection management functions for database operations.
"""

import os
import logging
import contextlib
from typing import Dict, Optional, Any, List
from functools import wraps

# Import Supabase client
try:
    from supabase import create_client, Client
except ImportError:
    raise ImportError("supabase package not installed. Install with 'pip install supabase'")

logger = logging.getLogger(__name__)

def get_connection_config() -> Dict[str, str]:
    """
    Get database connection configuration from environment variables.
    
    Returns:
        Dictionary with connection parameters
    """
    return {
        'url': os.environ.get('SUPABASE_URL'),
        'key': os.environ.get('SUPABASE_KEY'),
        'service_role': os.environ.get('SUPABASE_SERVICE_ROLE_KEY'),
        'database_url': os.environ.get('DATABASE_URL')
    }

class SupabaseClientWrapper:
    """
    Wrapper for Supabase client that adds fallback methods for missing functionality.
    """
    
    def __init__(self, client: Client):
        """
        Initialize the wrapper with a Supabase client.
        
        Args:
            client: Supabase client instance
        """
        self._client = client
        
    def __getattr__(self, name):
        """
        Forward attribute access to the wrapped client.
        
        Args:
            name: Attribute name
            
        Returns:
            The requested attribute from the wrapped client
        """
        return getattr(self._client, name)
    
    def sql(self, query: str):
        """
        Execute a SQL query.
        
        This method is not natively available in the Supabase client,
        but is required by the migration management module.
        
        Args:
            query: SQL query to execute
            
        Returns:
            Query builder with execute method
        """
        logger.info("Using sql method wrapper")
        
        # Create a wrapper for the execute method
        class ExecuteWrapper:
            def __init__(self, client, query):
                self.client = client
                self.query = query
                
            def execute(self):
                try:
                    # Execute the query using RPC with exec_sql function
                    return self.client.rpc('exec_sql', {'query': self.query}).execute()
                except Exception as e:
                    logger.error(f"Error executing SQL query: {str(e)}")
                    raise
        
        return ExecuteWrapper(self._client, query)
    
    def rpc(self, function_name: str, params: Dict = None):
        """
        Call a remote procedure.
        
        Provides a fallback for the 'has_table' function if it's not available.
        
        Args:
            function_name: Name of the RPC function to call
            params: Parameters to pass to the function
            
        Returns:
            RPC query builder
        """
        if function_name == 'has_table' and params and 'table_name' in params:
            # Fallback implementation for has_table
            return self._has_table_fallback(params['table_name'])
        
        # Use the standard RPC method
        return self._client.rpc(function_name, params)
    
    def _has_table_fallback(self, table_name: str):
        """
        Fallback implementation for the has_table RPC function.
        
        Args:
            table_name: Name of the table to check
            
        Returns:
            Query builder with execute method that returns a result with data
        """
        logger.info(f"Using fallback implementation for has_table({table_name})")
        
        # Create a query to check if the table exists
        query = f"""
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.tables
            WHERE table_schema = 'public'
            AND table_name = '{table_name}'
        ) AS table_exists;
        """
        
        # Create a wrapper for the result
        class HasTableResult:
            def __init__(self, exists: bool):
                self.data = [exists]
                
        # Create a wrapper for the execute method
        class ExecuteWrapper:
            def __init__(self, client, query):
                self.client = client
                self.query = query
                
            def execute(self):
                try:
                    # Execute the query using RPC with exec_sql function
                    response = self.client.rpc('exec_sql', {'query': self.query}).execute()
                    
                    # Extract the result
                    if hasattr(response, 'data') and response.data:
                        # The response format might vary, try different approaches
                        
                        # Log the response data for debugging
                        logger.info(f"Response data: {response.data}")
                        
                        # If it's a list of dictionaries
                        if isinstance(response.data, list) and len(response.data) > 0:
                            if 'table_exists' in response.data[0]:
                                exists = response.data[0]['table_exists']
                                return HasTableResult(exists)
                            elif 'exists' in response.data[0]:
                                exists = response.data[0]['exists']
                                return HasTableResult(exists)
                        
                        # If it's a dictionary with a 'success' key
                        if isinstance(response.data, dict) and 'success' in response.data:
                            # If success is True, assume the table exists
                            return HasTableResult(response.data['success'])
                        
                        # Try to interpret the response as a boolean
                        try:
                            # If the first element is a boolean or can be converted to one
                            if isinstance(response.data, list) and len(response.data) > 0:
                                if isinstance(response.data[0], bool):
                                    return HasTableResult(response.data[0])
                                elif str(response.data[0]).lower() in ('true', 't', 'yes', 'y', '1'):
                                    return HasTableResult(True)
                                elif str(response.data[0]).lower() in ('false', 'f', 'no', 'n', '0'):
                                    return HasTableResult(False)
                        except (IndexError, ValueError, TypeError):
                            pass
                        
                        # If all else fails, try a direct query to the table
                        try:
                            # Try to query the table directly
                            table_response = self.client.table(self.table_name).select('count').limit(1).execute()
                            # If we get here without an error, the table exists
                            return HasTableResult(True)
                        except Exception:
                            # If we get an error, the table probably doesn't exist
                            return HasTableResult(False)
                    
                    return HasTableResult(False)
                except Exception as e:
                    logger.error(f"Error in has_table fallback: {str(e)}")
                    return HasTableResult(False)
        
        return ExecuteWrapper(self._client, query)
        
        return ExecuteWrapper(self._client, query)

def create_connection(config: Optional[Dict] = None) -> Any:
    """
    Create a database connection using provided or environment config.
    
    This function returns a Supabase client that supports the methods
    required by the migration management module (.sql(), .table(), .rpc()).
    
    Args:
        config: Optional configuration dictionary
        
    Returns:
        Database connection object (Supabase client)
        
    Raises:
        ValueError: If required configuration is missing
        ConnectionError: If connection fails
    """
    if config is None:
        config = get_connection_config()
    
    url = config.get('url')
    key = config.get('key') or config.get('service_role')
    
    if not url or not key:
        logger.error("Missing required Supabase configuration (url and key)")
        raise ValueError("Missing required Supabase configuration (url and key)")
    
    try:
        logger.info("Creating database connection to Supabase")
        client = create_client(url, key)
        
        # Wrap the client to provide fallbacks for missing functionality
        wrapped_client = SupabaseClientWrapper(client)
        
        return wrapped_client
    except Exception as e:
        logger.error(f"Failed to create Supabase connection: {str(e)}")
        raise ConnectionError(f"Failed to create Supabase connection: {str(e)}")

@contextlib.contextmanager
def supabase_connection():
    """
    Context manager for a Supabase database connection.
    Yields a connection object for use in a with-statement.
    
    Example:
        with supabase_connection() as conn:
            result = conn.table('migrations').select('*').execute()
    
    Yields:
        A Supabase client connection
    """
    conn = None
    try:
        conn = create_connection()
        yield conn
    finally:
        # No explicit cleanup needed for Supabase client
        pass

# Connection pool support
_connection_pool = None
_pool_initialized = False

def initialize_connection_pool(min_connections: int = 2, max_connections: int = 10):
    """
    Initialize the database connection pool for Supabase.
    
    Args:
        min_connections: Minimum number of connections to maintain
        max_connections: Maximum number of connections allowed
    """
    global _connection_pool, _pool_initialized
    
    if _pool_initialized:
        logger.warning("Connection pool already initialized")
        return
    
    try:
        # Try to import the connection pool wrapper
        from connection_pool_wrapper import initialize_supabase_pool
        
        # Get connection config
        config = get_connection_config()
        
        # Initialize the pool
        _connection_pool = initialize_supabase_pool(
            supabase_url=config.get('url'),
            supabase_key=config.get('key') or config.get('service_role'),
            min_connections=min_connections,
            max_connections=max_connections
        )
        
        _pool_initialized = True
        logger.info(f"Initialized Supabase connection pool (min={min_connections}, max={max_connections})")
    except ImportError:
        logger.warning("connection_pool_wrapper not found. Connection pooling not available.")
        _pool_initialized = False
    except Exception as e:
        logger.error(f"Error initializing connection pool: {str(e)}")
        _pool_initialized = False

def close_connections():
    """
    Close all active database connections in the pool.
    """
    global _connection_pool, _pool_initialized
    
    if not _pool_initialized or _connection_pool is None:
        logger.warning("Connection pool not initialized")
        return
    
    try:
        # Try to import the connection pool wrapper
        from connection_pool_wrapper import shutdown_supabase_pool
        
        # Shutdown the pool
        shutdown_supabase_pool()
        
        _connection_pool = None
        _pool_initialized = False
        logger.info("Closed all connections in the Supabase connection pool")
    except ImportError:
        logger.warning("connection_pool_wrapper not found. Connection pooling not available.")
    except Exception as e:
        logger.error(f"Error closing connections: {str(e)}")

def get_supabase_service_client():
    """
    Get a Supabase service client for advanced operations.
    
    Returns:
        An instance of the Supabase service client.
    """
    config = get_connection_config()
    
    # Ensure we use the service role key
    service_role = config.get('service_role')
    if not service_role:
        logger.warning("SUPABASE_SERVICE_ROLE_KEY not set in environment")
    
    # Create a new config with the service role key
    service_config = {
        'url': config.get('url'),
        'key': service_role or config.get('key')
    }
    
    return create_connection(service_config)