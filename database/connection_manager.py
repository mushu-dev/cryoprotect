#!/usr/bin/env python3
"""
Database connection manager for CryoProtect v2.

This module provides a singleton connection manager that handles
database connections using various adapters.
"""

import os
import logging
import importlib
from typing import Any, Dict, List, Optional, Union, Tuple, Type
from dotenv import load_dotenv

from .adapter import DatabaseAdapter
from .local_adapter import LocalPostgreSQLAdapter
from .supabase_adapter import SupabaseDirectAdapter
from .mcp_adapter import MCPAdapter
from .enhanced_convex_adapter import ConvexAdapter

logger = logging.getLogger(__name__)

class ConnectionManager:
    """
    Singleton database connection manager.
    
    This class manages database connections using various adapters
    and provides a unified interface for database operations.
    """
    
    _instance = None
    
    @classmethod
    def get_instance(cls):
        """
        Get singleton instance of ConnectionManager.
        
        Returns:
            ConnectionManager: Singleton instance
        """
        if cls._instance is None:
            cls._instance = ConnectionManager()
        return cls._instance
    
    def __init__(self):
        """Initialize connection manager."""
        if ConnectionManager._instance is not None:
            raise RuntimeError("ConnectionManager is a singleton. Use get_instance() instead.")
            
        ConnectionManager._instance = self
        
        self.adapters = {}
        self.active_adapter = None
        self.config = self._load_config()
        
    def _load_config(self) -> Dict[str, Any]:
        """
        Load configuration from environment variables.
        
        Returns:
            Dict with configuration values
        """
        # Load environment variables
        load_dotenv()
        
        # Determine environment
        environment = os.getenv('ENVIRONMENT', 'development')
        
        # Load adapter configurations
        config = {
            'environment': environment,
            'adapters': {
                'local': {
                    'enabled': os.getenv('LOCAL_DB_ENABLED', 'true').lower() == 'true',
                    'host': os.getenv('LOCAL_DB_HOST', 'localhost'),
                    'port': int(os.getenv('LOCAL_DB_PORT', '5432')),
                    'database': os.getenv('LOCAL_DB_NAME', 'postgres'),
                    'user': os.getenv('LOCAL_DB_USER', 'postgres'),
                    'password': os.getenv('LOCAL_DB_PASSWORD', ''),
                    'min_connections': int(os.getenv('LOCAL_DB_MIN_CONNECTIONS', '1')),
                    'max_connections': int(os.getenv('LOCAL_DB_MAX_CONNECTIONS', '10'))
                },
                'supabase': {
                    'enabled': os.getenv('SUPABASE_DB_ENABLED', 'true').lower() == 'true',
                    'host': os.getenv('SUPABASE_DB_HOST', ''),
                    'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
                    'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
                    'user': os.getenv('SUPABASE_DB_USER', ''),
                    'password': os.getenv('SUPABASE_DB_PASSWORD', ''),
                    'min_connections': int(os.getenv('SUPABASE_DB_MIN_CONNECTIONS', '1')),
                    'max_connections': int(os.getenv('SUPABASE_DB_MAX_CONNECTIONS', '10')),
                    'ip_address': os.getenv('SUPABASE_DB_IP_ADDRESS', '')
                },
                'convex': {
                    'enabled': os.getenv('CONVEX_DB_ENABLED', 'false').lower() == 'true',
                    'url': os.getenv('CONVEX_URL', ''),
                    'key': os.getenv('CONVEX_DEPLOYMENT_KEY', ''),
                    'timeout': int(os.getenv('CONVEX_TIMEOUT', '30')),
                    'retry_count': int(os.getenv('CONVEX_RETRY_COUNT', '3')),
                    'circuit_breaker_threshold': int(os.getenv('CONVEX_CIRCUIT_BREAKER_THRESHOLD', '5')),
                    'circuit_breaker_timeout': int(os.getenv('CONVEX_CIRCUIT_BREAKER_TIMEOUT', '60'))
                },
                'mcp': {
                    'enabled': os.getenv('MCP_DB_ENABLED', 'true').lower() == 'true',
                    'server_name': os.getenv('MCP_SERVER_NAME', 'supabase'),
                    'timeout': int(os.getenv('MCP_TIMEOUT', '30'))
                }
            },
            'adapter_order': os.getenv('DB_ADAPTER_ORDER', 'local,supabase,convex,mcp').split(',')
        }
        
        return config
    
    def connect(self) -> bool:
        """
        Connect to database using available adapters.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        # Try adapters in configured order
        for adapter_name in self.config['adapter_order']:
            adapter_config = self.config['adapters'].get(adapter_name)
            
            if not adapter_config or not adapter_config.get('enabled', False):
                logger.debug(f"Skipping disabled adapter: {adapter_name}")
                continue
                
            logger.info(f"Trying to connect using {adapter_name} adapter")
            
            # Initialize adapter if not already initialized
            if adapter_name not in self.adapters:
                self._initialize_adapter(adapter_name, adapter_config)
                
            adapter = self.adapters.get(adapter_name)
            if not adapter:
                logger.warning(f"Failed to initialize {adapter_name} adapter")
                continue
                
            # Try to connect
            if adapter.connect():
                self.active_adapter = adapter_name
                logger.info(f"Successfully connected using {adapter_name} adapter")
                return True
                
            logger.warning(f"Failed to connect using {adapter_name} adapter")
            
        logger.error("Failed to connect using any adapter")
        return False
        
    def _initialize_adapter(self, adapter_name: str, adapter_config: Dict[str, Any]) -> None:
        """
        Initialize database adapter.
        
        Args:
            adapter_name: Name of the adapter
            adapter_config: Adapter configuration
        """
        try:
            if adapter_name == 'local':
                self.adapters[adapter_name] = LocalPostgreSQLAdapter(adapter_config)
            elif adapter_name == 'supabase':
                self.adapters[adapter_name] = SupabaseDirectAdapter(adapter_config)
            elif adapter_name == 'convex':
                self.adapters[adapter_name] = ConvexAdapter(adapter_config)
            elif adapter_name == 'mcp':
                self.adapters[adapter_name] = MCPAdapter(adapter_config)
            else:
                logger.warning(f"Unknown adapter type: {adapter_name}")
        except Exception as e:
            logger.error(f"Error initializing {adapter_name} adapter: {str(e)}")
            
    def disconnect(self) -> bool:
        """
        Disconnect from all database adapters.
        
        Returns:
            bool: True if disconnection successful, False otherwise
        """
        success = True
        
        for adapter_name, adapter in self.adapters.items():
            try:
                if adapter.disconnect():
                    logger.info(f"Disconnected from {adapter_name} adapter")
                else:
                    logger.warning(f"Failed to disconnect from {adapter_name} adapter")
                    success = False
            except Exception as e:
                logger.error(f"Error disconnecting from {adapter_name} adapter: {str(e)}")
                success = False
                
        self.active_adapter = None
        return success
        
    def get_active_adapter(self) -> Optional[DatabaseAdapter]:
        """
        Get the currently active database adapter.
        
        Returns:
            DatabaseAdapter or None if no active adapter
        """
        if not self.active_adapter:
            return None
            
        return self.adapters.get(self.active_adapter)
        
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """
        Execute SQL query using active adapter.
        
        Args:
            query: SQL query to execute
            params: Query parameters
            
        Returns:
            Query results
            
        Raises:
            ConnectionError: If no active adapter
        """
        adapter = self.get_active_adapter()
        if not adapter:
            if not self.connect():
                raise ConnectionError("No active database adapter")
            adapter = self.get_active_adapter()
            
        return adapter.execute_query(query, params)
        
    def execute_batch(self, queries: List[str]) -> List[Any]:
        """
        Execute multiple SQL queries using active adapter.
        
        Args:
            queries: List of SQL queries to execute
            
        Returns:
            List of query results
            
        Raises:
            ConnectionError: If no active adapter
        """
        adapter = self.get_active_adapter()
        if not adapter:
            if not self.connect():
                raise ConnectionError("No active database adapter")
            adapter = self.get_active_adapter()
            
        return adapter.execute_batch(queries)
        
    def begin_transaction(self) -> Any:
        """
        Begin a database transaction using active adapter.
        
        Returns:
            Transaction object
            
        Raises:
            ConnectionError: If no active adapter
        """
        adapter = self.get_active_adapter()
        if not adapter:
            if not self.connect():
                raise ConnectionError("No active database adapter")
            adapter = self.get_active_adapter()
            
        return adapter.begin_transaction()
        
    def commit_transaction(self, transaction: Any) -> bool:
        """
        Commit a database transaction using active adapter.
        
        Args:
            transaction: Transaction object
            
        Returns:
            bool: True if commit successful, False otherwise
            
        Raises:
            ConnectionError: If no active adapter
        """
        adapter = self.get_active_adapter()
        if not adapter:
            raise ConnectionError("No active database adapter")
            
        return adapter.commit_transaction(transaction)
        
    def rollback_transaction(self, transaction: Any) -> bool:
        """
        Rollback a database transaction using active adapter.
        
        Args:
            transaction: Transaction object
            
        Returns:
            bool: True if rollback successful, False otherwise
            
        Raises:
            ConnectionError: If no active adapter
        """
        adapter = self.get_active_adapter()
        if not adapter:
            raise ConnectionError("No active database adapter")
            
        return adapter.rollback_transaction(transaction)
        
    def get_connection_info(self) -> Dict[str, Any]:
        """
        Get connection information for all adapters.
        
        Returns:
            Dict with connection information
        """
        info = {
            'environment': self.config['environment'],
            'active_adapter': self.active_adapter,
            'adapters': {}
        }
        
        for adapter_name, adapter in self.adapters.items():
            try:
                adapter_info = adapter.get_connection_info()
                info['adapters'][adapter_name] = adapter_info
            except Exception as e:
                info['adapters'][adapter_name] = {
                    'error': str(e)
                }
                
        return info
        
    def test_all_connections(self) -> Dict[str, Tuple[bool, str]]:
        """
        Test all database connections.
        
        Returns:
            Dict with adapter names as keys and (success, message) tuples as values
        """
        results = {}
        
        for adapter_name, adapter in self.adapters.items():
            try:
                success, message = adapter.test_connection()
                results[adapter_name] = (success, message)
            except Exception as e:
                results[adapter_name] = (False, str(e))
                
        return results