#!/usr/bin/env python3
"""
Revised Database Connection System for CryoProtect v2.

This module provides a connection factory and adapter interface for
database connections with automatic fallback logic.
"""

import os
import logging
import importlib
import time
import json
import uuid
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Union, Tuple, Type
from dotenv import load_dotenv

# Configure logger
logger = logging.getLogger(__name__)

# Create a connection request ID generator for tracking connection lifecycle
def generate_connection_request_id():
    """Generate a unique ID for tracking connection requests through the lifecycle."""
    return str(uuid.uuid4())[:8]

class ConnectionAdapter(ABC):
    """
    Abstract connection adapter interface.
    
    This class defines the interface that all database connection adapters
    must implement. It provides methods for connection management, query
    execution, transaction handling, and connection information.
    """
    
    @abstractmethod
    def connect(self) -> bool:
        """
        Establish connection to the database.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        pass
        
    @abstractmethod
    def disconnect(self) -> bool:
        """
        Close database connection.
        
        Returns:
            bool: True if disconnection successful, False otherwise
        """
        pass
        
    @abstractmethod
    def reconnect(self) -> bool:
        """
        Re-establish a dropped connection to the database.
        
        This method should handle cleaning up any existing connection resources
        before attempting to establish a new connection.
        
        Returns:
            bool: True if reconnection successful, False otherwise
        """
        pass
        
    @abstractmethod
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """
        Execute SQL query and return results.
        
        Args:
            query: SQL query to execute
            params: Query parameters
            
        Returns:
            Query results
        """
        pass
        
    @abstractmethod
    def execute_batch(self, queries: List[str]) -> List[Any]:
        """
        Execute multiple SQL queries and return results.
        
        Args:
            queries: List of SQL queries to execute
            
        Returns:
            List of query results
        """
        pass
        
    @abstractmethod
    def begin_transaction(self) -> Any:
        """
        Begin a database transaction.
        
        Returns:
            Transaction object
        """
        pass
        
    @abstractmethod
    def commit_transaction(self, transaction: Any) -> bool:
        """
        Commit a database transaction.
        
        Args:
            transaction: Transaction object
            
        Returns:
            bool: True if commit successful, False otherwise
        """
        pass
        
    @abstractmethod
    def rollback_transaction(self, transaction: Any) -> bool:
        """
        Rollback a database transaction.
        
        Args:
            transaction: Transaction object
            
        Returns:
            bool: True if rollback successful, False otherwise
        """
        pass
        
    @abstractmethod
    def get_connection_info(self) -> Dict[str, Any]:
        """
        Get connection information.
        
        Returns:
            Dict with connection information
        """
        pass
        
    @abstractmethod
    def test_connection(self) -> Tuple[bool, str]:
        """
        Test database connection and return status with message.
        
        This is a basic connectivity test that should be fast and lightweight.
        
        Returns:
            Tuple of (success: bool, message: str)
        """
        pass
        
    @abstractmethod
    def is_healthy(self) -> Tuple[bool, Dict[str, Any]]:
        """
        Perform a comprehensive health check on the database connection.
        
        This method should check various aspects of the connection health:
        - Basic connectivity
        - Connection latency
        - Transaction capability
        - Query execution capability
        - Connection pool status (if applicable)
        
        Returns:
            Tuple of (healthy: bool, health_metrics: Dict[str, Any])
            where health_metrics contains detailed information about the connection health
        """
        pass


class ConnectionFactory:
    """
    Factory class for creating and managing database connection adapters.
    
    This class provides a unified interface for obtaining database connections
    with automatic fallback logic. It manages multiple adapter types and
    attempts to connect using each adapter in a configured order until a
    successful connection is established.
    """
    
    _instance = None
    
    @classmethod
    def get_instance(cls):
        """
        Get singleton instance of ConnectionFactory.
        
        Returns:
            ConnectionFactory: Singleton instance
        """
        if cls._instance is None:
            cls._instance = ConnectionFactory()
        return cls._instance
    
    def __init__(self):
        """Initialize connection factory."""
        if ConnectionFactory._instance is not None:
            raise RuntimeError("ConnectionFactory is a singleton. Use get_instance() instead.")
            
        ConnectionFactory._instance = self
        
        self.adapters = {}
        self.active_adapter = None
        self.config = self._load_config()
        self.connection_attempts = {}  # Track connection attempts for backoff
        self.health_check_results = {}  # Store recent health check results
        self.connection_stats = {
            "total_requests": 0,
            "successful_connections": 0,
            "failed_connections": 0,
            "reconnections": 0,
            "fallbacks": 0,
            "adapter_stats": {}
        }
        
        logger.info("ConnectionFactory initialized with adapter order: %s",
                   self.config.get('adapter_order', []))
        
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
                'pooler': {
                    'enabled': os.getenv('POOLER_DB_ENABLED', 'true').lower() == 'true',
                    'host': os.getenv('POOLER_DB_HOST', 'localhost'),
                    'port': int(os.getenv('POOLER_DB_PORT', '5432')),
                    'database': os.getenv('POOLER_DB_NAME', 'postgres'),
                    'user': os.getenv('POOLER_DB_USER', 'postgres'),
                    'password': os.getenv('POOLER_DB_PASSWORD', ''),
                    'min_connections': int(os.getenv('POOLER_DB_MIN_CONNECTIONS', '1')),
                    'max_connections': int(os.getenv('POOLER_DB_MAX_CONNECTIONS', '10')),
                    'pool_timeout': int(os.getenv('POOLER_DB_POOL_TIMEOUT', '30')),
                    'pool_recycle': int(os.getenv('POOLER_DB_POOL_RECYCLE', '1800'))
                },
                'direct': {
                    'enabled': os.getenv('DIRECT_DB_ENABLED', 'true').lower() == 'true',
                    'host': os.getenv('DIRECT_DB_HOST', ''),
                    'port': int(os.getenv('DIRECT_DB_PORT', '5432')),
                    'database': os.getenv('DIRECT_DB_NAME', 'postgres'),
                    'user': os.getenv('DIRECT_DB_USER', ''),
                    'password': os.getenv('DIRECT_DB_PASSWORD', ''),
                    'min_connections': int(os.getenv('DIRECT_DB_MIN_CONNECTIONS', '1')),
                    'max_connections': int(os.getenv('DIRECT_DB_MAX_CONNECTIONS', '10')),
                    'ip_address': os.getenv('DIRECT_DB_IP_ADDRESS', '')
                },
                'mcp': {
                    'enabled': os.getenv('MCP_DB_ENABLED', 'true').lower() == 'true',
                    'server_name': os.getenv('MCP_SERVER_NAME', 'supabase'),
                    'project_id': os.getenv('MCP_PROJECT_ID', ''),
                    'timeout': int(os.getenv('MCP_TIMEOUT', '30'))
                }
            },
            'adapter_order': os.getenv('DB_ADAPTER_ORDER', 'local,pooler,direct,mcp').split(',')
        }
        
        return config
    
    def get_connection(self) -> Optional[ConnectionAdapter]:
        """
        Get a database connection using available adapters with automatic fallback.
        
        This method attempts to connect using each adapter in the configured order
        until a successful connection is established. If a connection is already
        active, it returns the existing connection after validating its health.
        
        Returns:
            ConnectionAdapter or None if no connection could be established
        """
        request_id = generate_connection_request_id()
        start_time = time.time()
        self.connection_stats["total_requests"] += 1
        
        logger.info("[ConnReq:%s] Connection request started", request_id)
        
        # Return existing connection if available and healthy
        if self.active_adapter and self.active_adapter in self.adapters:
            adapter = self.adapters[self.active_adapter]
            adapter_name = self.active_adapter
            
            logger.debug("[ConnReq:%s] Testing existing %s connection", request_id, adapter_name)
            
            # Test if connection is still valid
            test_start = time.time()
            success, message = adapter.test_connection()
            test_duration = time.time() - test_start
            
            if success:
                logger.info("[ConnReq:%s] Using existing %s connection (test_duration=%.3fs)",
                           request_id, adapter_name, test_duration)
                
                # Update adapter stats
                if adapter_name not in self.connection_stats["adapter_stats"]:
                    self.connection_stats["adapter_stats"][adapter_name] = {
                        "connection_attempts": 0,
                        "connection_successes": 0,
                        "connection_failures": 0,
                        "reconnection_attempts": 0,
                        "reconnection_successes": 0,
                        "last_connection_time": time.time()
                    }
                
                return adapter
            else:
                logger.warning("[ConnReq:%s] Existing %s connection is no longer valid: %s",
                              request_id, adapter_name, message)
                
                # Attempt to reconnect the current adapter
                logger.info("[ConnReq:%s] Attempting to reconnect %s adapter",
                           request_id, adapter_name)
                
                # Update reconnection stats
                self.connection_stats["reconnections"] += 1
                if adapter_name in self.connection_stats["adapter_stats"]:
                    self.connection_stats["adapter_stats"][adapter_name]["reconnection_attempts"] += 1
                
                if self._attempt_reconnect(self.active_adapter, request_id):
                    logger.info("[ConnReq:%s] Successfully reconnected %s adapter",
                               request_id, adapter_name)
                    
                    # Update reconnection success stats
                    if adapter_name in self.connection_stats["adapter_stats"]:
                        self.connection_stats["adapter_stats"][adapter_name]["reconnection_successes"] += 1
                        self.connection_stats["adapter_stats"][adapter_name]["last_connection_time"] = time.time()
                    
                    return self.adapters[self.active_adapter]
                
                # If reconnection failed, clear active adapter and try fallback
                logger.warning("[ConnReq:%s] Reconnection to %s adapter failed, trying fallback",
                              request_id, adapter_name)
                self.active_adapter = None
                self.connection_stats["fallbacks"] += 1
        
        # Try adapters in configured order
        logger.info("[ConnReq:%s] Trying adapters in configured order: %s",
                   request_id, self.config['adapter_order'])
        
        for adapter_name in self.config['adapter_order']:
            adapter_config = self.config['adapters'].get(adapter_name)
            
            if not adapter_config or not adapter_config.get('enabled', False):
                logger.debug("[ConnReq:%s] Skipping disabled adapter: %s",
                            request_id, adapter_name)
                continue
                
            # Check if we should apply backoff for this adapter
            if not self._should_attempt_connection(adapter_name):
                logger.debug("[ConnReq:%s] Skipping %s adapter due to backoff policy",
                            request_id, adapter_name)
                continue
                
            logger.info("[ConnReq:%s] Trying to connect using %s adapter",
                       request_id, adapter_name)
            
            # Initialize adapter if not already initialized
            if adapter_name not in self.adapters:
                logger.debug("[ConnReq:%s] Initializing %s adapter", request_id, adapter_name)
                self._initialize_adapter(adapter_name, adapter_config)
                
            adapter = self.adapters.get(adapter_name)
            if not adapter:
                logger.warning("[ConnReq:%s] Failed to initialize %s adapter",
                              request_id, adapter_name)
                self._record_connection_attempt(adapter_name, False)
                
                # Update adapter stats
                if adapter_name not in self.connection_stats["adapter_stats"]:
                    self.connection_stats["adapter_stats"][adapter_name] = {
                        "connection_attempts": 1,
                        "connection_successes": 0,
                        "connection_failures": 1,
                        "reconnection_attempts": 0,
                        "reconnection_successes": 0
                    }
                else:
                    self.connection_stats["adapter_stats"][adapter_name]["connection_attempts"] += 1
                    self.connection_stats["adapter_stats"][adapter_name]["connection_failures"] += 1
                
                continue
                
            # Update adapter stats before connection attempt
            if adapter_name not in self.connection_stats["adapter_stats"]:
                self.connection_stats["adapter_stats"][adapter_name] = {
                    "connection_attempts": 1,
                    "connection_successes": 0,
                    "connection_failures": 0,
                    "reconnection_attempts": 0,
                    "reconnection_successes": 0
                }
            else:
                self.connection_stats["adapter_stats"][adapter_name]["connection_attempts"] += 1
            
            # Try to connect
            connect_start = time.time()
            connection_success = adapter.connect()
            connect_duration = time.time() - connect_start
            
            if connection_success:
                self.active_adapter = adapter_name
                self._record_connection_attempt(adapter_name, True)
                
                # Update connection success stats
                self.connection_stats["successful_connections"] += 1
                self.connection_stats["adapter_stats"][adapter_name]["connection_successes"] += 1
                self.connection_stats["adapter_stats"][adapter_name]["last_connection_time"] = time.time()
                
                logger.info("[ConnReq:%s] Successfully connected using %s adapter (duration=%.3fs)",
                           request_id, adapter_name, connect_duration)
                
                # Log connection details
                conn_info = adapter.get_connection_info()
                logger.info("[ConnReq:%s] Connection details: %s",
                           request_id, json.dumps(conn_info))
                
                total_duration = time.time() - start_time
                logger.info("[ConnReq:%s] Connection request completed successfully in %.3fs",
                           request_id, total_duration)
                
                return adapter
                
            # Update connection failure stats
            self.connection_stats["failed_connections"] += 1
            self.connection_stats["adapter_stats"][adapter_name]["connection_failures"] += 1
            
            logger.warning("[ConnReq:%s] Failed to connect using %s adapter (duration=%.3fs)",
                          request_id, adapter_name, connect_duration)
            self._record_connection_attempt(adapter_name, False)
            
        total_duration = time.time() - start_time
        logger.error("[ConnReq:%s] Failed to connect using any adapter after %.3fs",
                    request_id, total_duration)
        
        # Log connection statistics
        logger.info("[ConnReq:%s] Connection statistics: %s",
                   request_id, json.dumps({
                       "total_requests": self.connection_stats["total_requests"],
                       "successful_connections": self.connection_stats["successful_connections"],
                       "failed_connections": self.connection_stats["failed_connections"],
                       "reconnections": self.connection_stats["reconnections"],
                       "fallbacks": self.connection_stats["fallbacks"]
                   }))
        
        return None
        
    def _should_attempt_connection(self, adapter_name: str) -> bool:
        """
        Determine if a connection attempt should be made based on backoff policy.
        
        Args:
            adapter_name: Name of the adapter
            
        Returns:
            bool: True if connection attempt should be made, False otherwise
        """
        import time
        
        if adapter_name not in self.connection_attempts:
            return True
            
        attempts = self.connection_attempts[adapter_name]
        if attempts['success']:
            return True
            
        # Calculate backoff time based on number of consecutive failures
        # Using exponential backoff with a maximum of 5 minutes
        backoff_seconds = min(2 ** attempts['consecutive_failures'], 300)
        
        # Check if enough time has passed since the last attempt
        time_since_last_attempt = time.time() - attempts['last_attempt_time']
        return time_since_last_attempt >= backoff_seconds
        
    def _record_connection_attempt(self, adapter_name: str, success: bool) -> None:
        """
        Record a connection attempt for backoff calculation.
        
        Args:
            adapter_name: Name of the adapter
            success: Whether the connection attempt was successful
        """
        import time
        
        if adapter_name not in self.connection_attempts:
            self.connection_attempts[adapter_name] = {
                'success': False,
                'consecutive_failures': 0,
                'last_attempt_time': 0
            }
            
        attempts = self.connection_attempts[adapter_name]
        attempts['last_attempt_time'] = time.time()
        
        if success:
            attempts['success'] = True
            attempts['consecutive_failures'] = 0
        else:
            attempts['success'] = False
            attempts['consecutive_failures'] += 1
            
    def _attempt_reconnect(self, adapter_name: str, request_id: str = None) -> bool:
        """
        Attempt to reconnect a specific adapter.
        
        Args:
            adapter_name: Name of the adapter to reconnect
            request_id: Optional connection request ID for logging
            
        Returns:
            bool: True if reconnection successful, False otherwise
        """
        if not request_id:
            request_id = generate_connection_request_id()
            
        if adapter_name not in self.adapters:
            logger.warning("[ConnReq:%s] Cannot reconnect non-existent adapter: %s",
                          request_id, adapter_name)
            return False
            
        adapter = self.adapters[adapter_name]
        
        # Record the reconnection attempt
        self._record_connection_attempt(adapter_name, False)
        
        # Attempt reconnection
        try:
            logger.info("[ConnReq:%s] Executing reconnect() on %s adapter",
                       request_id, adapter_name)
            
            reconnect_start = time.time()
            reconnect_success = adapter.reconnect()
            reconnect_duration = time.time() - reconnect_start
            
            if reconnect_success:
                self._record_connection_attempt(adapter_name, True)
                logger.info("[ConnReq:%s] Successfully reconnected %s adapter (duration=%.3fs)",
                           request_id, adapter_name, reconnect_duration)
                
                # Log connection details after successful reconnection
                try:
                    conn_info = adapter.get_connection_info()
                    logger.info("[ConnReq:%s] Reconnection details: %s",
                               request_id, json.dumps(conn_info))
                except Exception as e:
                    logger.warning("[ConnReq:%s] Failed to get connection info after reconnection: %s",
                                  request_id, str(e))
                
                return True
            else:
                logger.warning("[ConnReq:%s] Reconnect() on %s adapter returned False (duration=%.3fs)",
                              request_id, adapter_name, reconnect_duration)
        except Exception as e:
            logger.error("[ConnReq:%s] Error reconnecting %s adapter: %s",
                        request_id, adapter_name, str(e))
            
        return False
        
    def _initialize_adapter(self, adapter_name: str, adapter_config: Dict[str, Any], request_id: str = None) -> None:
        """
        Initialize database adapter.
        
        Args:
            adapter_name: Name of the adapter
            adapter_config: Adapter configuration
            request_id: Optional connection request ID for logging
        """
        if not request_id:
            request_id = generate_connection_request_id()
            
        logger.info("[ConnReq:%s] Initializing %s adapter", request_id, adapter_name)
        
        # Log sanitized configuration (without sensitive data)
        safe_config = adapter_config.copy()
        if 'password' in safe_config:
            safe_config['password'] = '******'
        logger.debug("[ConnReq:%s] %s adapter configuration: %s",
                    request_id, adapter_name, json.dumps(safe_config))
        
        init_start = time.time()
        try:
            # Import the appropriate adapter class
            if adapter_name == 'local':
                from .adapters.local import LocalAdapter
                logger.debug("[ConnReq:%s] Importing LocalAdapter class", request_id)
                self.adapters[adapter_name] = LocalAdapter(adapter_config)
            elif adapter_name == 'pooler':
                from .adapters.pooler import PoolerAdapter
                logger.debug("[ConnReq:%s] Importing PoolerAdapter class", request_id)
                self.adapters[adapter_name] = PoolerAdapter(adapter_config)
            elif adapter_name == 'direct':
                from .adapters.direct import DirectAdapter
                logger.debug("[ConnReq:%s] Importing DirectAdapter class", request_id)
                self.adapters[adapter_name] = DirectAdapter(adapter_config)
            elif adapter_name == 'mcp':
                from .adapters.mcp import MCPAdapter
                logger.debug("[ConnReq:%s] Importing MCPAdapter class", request_id)
                self.adapters[adapter_name] = MCPAdapter(adapter_config)
            else:
                logger.warning("[ConnReq:%s] Unknown adapter type: %s", request_id, adapter_name)
                return
                
            init_duration = time.time() - init_start
            logger.info("[ConnReq:%s] Successfully initialized %s adapter (duration=%.3fs)",
                       request_id, adapter_name, init_duration)
            
        except ImportError as e:
            init_duration = time.time() - init_start
            logger.error("[ConnReq:%s] Import error initializing %s adapter: %s (duration=%.3fs)",
                        request_id, adapter_name, str(e), init_duration)
            # Track adapter initialization failures in stats
            if adapter_name not in self.connection_stats["adapter_stats"]:
                self.connection_stats["adapter_stats"][adapter_name] = {
                    "connection_attempts": 0,
                    "connection_successes": 0,
                    "connection_failures": 0,
                    "reconnection_attempts": 0,
                    "reconnection_successes": 0,
                    "initialization_failures": 1
                }
            else:
                if "initialization_failures" not in self.connection_stats["adapter_stats"][adapter_name]:
                    self.connection_stats["adapter_stats"][adapter_name]["initialization_failures"] = 1
                else:
                    self.connection_stats["adapter_stats"][adapter_name]["initialization_failures"] += 1
                    
        except Exception as e:
            init_duration = time.time() - init_start
            logger.error("[ConnReq:%s] Error initializing %s adapter: %s (duration=%.3fs)",
                        request_id, adapter_name, str(e), init_duration)
            # Track adapter initialization failures in stats
            if adapter_name not in self.connection_stats["adapter_stats"]:
                self.connection_stats["adapter_stats"][adapter_name] = {
                    "connection_attempts": 0,
                    "connection_successes": 0,
                    "connection_failures": 0,
                    "reconnection_attempts": 0,
                    "reconnection_successes": 0,
                    "initialization_failures": 1
                }
            else:
                if "initialization_failures" not in self.connection_stats["adapter_stats"][adapter_name]:
                    self.connection_stats["adapter_stats"][adapter_name]["initialization_failures"] = 1
                else:
                    self.connection_stats["adapter_stats"][adapter_name]["initialization_failures"] += 1
            
    def close_all_connections(self) -> bool:
        """
        Close all database connections.
        
        Returns:
            bool: True if all disconnections successful, False otherwise
        """
        request_id = generate_connection_request_id()
        success = True
        
        logger.info("[ConnReq:%s] Closing all database connections", request_id)
        start_time = time.time()
        
        adapter_results = {}
        
        for adapter_name, adapter in self.adapters.items():
            try:
                logger.debug("[ConnReq:%s] Attempting to disconnect %s adapter",
                            request_id, adapter_name)
                
                disconnect_start = time.time()
                disconnect_success = adapter.disconnect()
                disconnect_duration = time.time() - disconnect_start
                
                adapter_results[adapter_name] = {
                    "success": disconnect_success,
                    "duration": disconnect_duration
                }
                
                if disconnect_success:
                    logger.info("[ConnReq:%s] Disconnected from %s adapter (duration=%.3fs)",
                               request_id, adapter_name, disconnect_duration)
                else:
                    logger.warning("[ConnReq:%s] Failed to disconnect from %s adapter (duration=%.3fs)",
                                  request_id, adapter_name, disconnect_duration)
                    success = False
            except Exception as e:
                disconnect_duration = time.time() - disconnect_start if 'disconnect_start' in locals() else 0
                logger.error("[ConnReq:%s] Error disconnecting from %s adapter: %s (duration=%.3fs)",
                            request_id, adapter_name, str(e), disconnect_duration)
                
                adapter_results[adapter_name] = {
                    "success": False,
                    "error": str(e),
                    "duration": disconnect_duration
                }
                
                success = False
                
        total_duration = time.time() - start_time
        
        if success:
            logger.info("[ConnReq:%s] Successfully closed all database connections (duration=%.3fs)",
                       request_id, total_duration)
        else:
            logger.warning("[ConnReq:%s] Failed to close some database connections (duration=%.3fs): %s",
                          request_id, total_duration, json.dumps(adapter_results))
            
        self.active_adapter = None
        return success
        
    def get_active_adapter_name(self) -> Optional[str]:
        """
        Get the name of the currently active adapter.
        
        Returns:
            str or None if no active adapter
        """
        return self.active_adapter
        
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
        request_id = generate_connection_request_id()
        results = {}
        
        logger.info("[ConnReq:%s] Testing all database connections", request_id)
        start_time = time.time()
        
        for adapter_name, adapter in self.adapters.items():
            try:
                logger.debug("[ConnReq:%s] Testing %s adapter connection",
                            request_id, adapter_name)
                
                test_start = time.time()
                success, message = adapter.test_connection()
                test_duration = time.time() - test_start
                
                results[adapter_name] = (success, message)
                
                if success:
                    logger.info("[ConnReq:%s] %s adapter connection test successful (duration=%.3fs): %s",
                               request_id, adapter_name, test_duration, message)
                else:
                    logger.warning("[ConnReq:%s] %s adapter connection test failed (duration=%.3fs): %s",
                                  request_id, adapter_name, test_duration, message)
            except Exception as e:
                test_duration = time.time() - test_start if 'test_start' in locals() else 0
                logger.error("[ConnReq:%s] Error testing %s adapter connection: %s (duration=%.3fs)",
                            request_id, adapter_name, str(e), test_duration)
                results[adapter_name] = (False, str(e))
        
        total_duration = time.time() - start_time
        success_count = sum(1 for result in results.values() if result[0])
        
        logger.info("[ConnReq:%s] Connection test results: %d/%d successful (duration=%.3fs)",
                   request_id, success_count, len(results), total_duration)
                
        return results
        
    def check_all_connections_health(self) -> Dict[str, Dict[str, Any]]:
        """
        Check health of all database connections.
        
        This method performs a comprehensive health check on all adapters
        and stores the results for future reference.
        
        Returns:
            Dict with adapter names as keys and health check results as values
        """
        request_id = generate_connection_request_id()
        results = {}
        
        logger.info("[ConnReq:%s] Performing comprehensive health check on all database connections",
                   request_id)
        check_start_time = time.time()
        
        healthy_adapters = 0
        unhealthy_adapters = 0
        
        for adapter_name, adapter in self.adapters.items():
            try:
                logger.debug("[ConnReq:%s] Checking health of %s adapter",
                            request_id, adapter_name)
                
                start_time = time.time()
                healthy, metrics = adapter.is_healthy()
                elapsed_time = time.time() - start_time
                
                result = {
                    'healthy': healthy,
                    'metrics': metrics,
                    'check_time': time.time(),
                    'check_duration': elapsed_time
                }
                
                results[adapter_name] = result
                self.health_check_results[adapter_name] = result
                
                if healthy:
                    healthy_adapters += 1
                    logger.info("[ConnReq:%s] %s adapter is healthy (duration=%.3fs)",
                               request_id, adapter_name, elapsed_time)
                    logger.debug("[ConnReq:%s] %s adapter health metrics: %s",
                                request_id, adapter_name, json.dumps(metrics))
                else:
                    unhealthy_adapters += 1
                    logger.warning("[ConnReq:%s] %s adapter is unhealthy (duration=%.3fs): %s",
                                  request_id, adapter_name, elapsed_time, json.dumps(metrics))
                    
                    if adapter_name == self.active_adapter:
                        logger.error("[ConnReq:%s] Active adapter %s is unhealthy, may need fallback",
                                    request_id, adapter_name)
                    
            except Exception as e:
                unhealthy_adapters += 1
                logger.error("[ConnReq:%s] Error checking health of %s adapter: %s",
                            request_id, adapter_name, str(e))
                
                results[adapter_name] = {
                    'healthy': False,
                    'error': str(e),
                    'check_time': time.time()
                }
                self.health_check_results[adapter_name] = results[adapter_name]
        
        total_duration = time.time() - check_start_time
        logger.info("[ConnReq:%s] Health check complete: %d healthy, %d unhealthy adapters (duration=%.3fs)",
                   request_id, healthy_adapters, unhealthy_adapters, total_duration)
                
        return results
        
    def validate_connection_before_use(self, adapter: ConnectionAdapter, adapter_name: str) -> bool:
        """
        Validate connection health before use and attempt reconnection if needed.
        
        Args:
            adapter: The adapter to validate
            adapter_name: Name of the adapter
            
        Returns:
            bool: True if connection is valid or was successfully reconnected, False otherwise
        """
        request_id = generate_connection_request_id()
        
        logger.info("[ConnReq:%s] Validating %s adapter connection before use",
                   request_id, adapter_name)
        start_time = time.time()
        
        try:
            # First try a simple connection test
            logger.debug("[ConnReq:%s] Performing basic connection test on %s adapter",
                        request_id, adapter_name)
            
            test_start = time.time()
            success, message = adapter.test_connection()
            test_duration = time.time() - test_start
            
            if success:
                logger.info("[ConnReq:%s] %s adapter connection is valid (duration=%.3fs)",
                           request_id, adapter_name, test_duration)
                return True
                
            # If simple test fails, attempt reconnection
            logger.warning("[ConnReq:%s] Connection validation failed for %s adapter: %s (duration=%.3fs)",
                          request_id, adapter_name, message, test_duration)
            
            logger.info("[ConnReq:%s] Attempting reconnection for %s adapter",
                       request_id, adapter_name)
            
            if self._attempt_reconnect(adapter_name, request_id):
                logger.info("[ConnReq:%s] Successfully reconnected %s adapter",
                           request_id, adapter_name)
                return True
                
            # If reconnection fails, try fallback
            logger.error("[ConnReq:%s] Failed to reconnect %s adapter, attempting fallback",
                        request_id, adapter_name)
            
            self.active_adapter = None
            self.connection_stats["fallbacks"] += 1
            
            fallback_start = time.time()
            new_adapter = self.get_connection()
            fallback_duration = time.time() - fallback_start
            
            if new_adapter is not None:
                logger.info("[ConnReq:%s] Successfully established fallback connection (duration=%.3fs)",
                           request_id, fallback_duration)
            else:
                logger.error("[ConnReq:%s] Failed to establish fallback connection (duration=%.3fs)",
                            request_id, fallback_duration)
                
            total_duration = time.time() - start_time
            logger.info("[ConnReq:%s] Connection validation completed in %.3fs with result: %s",
                       request_id, total_duration, new_adapter is not None)
                
            return new_adapter is not None
            
        except Exception as e:
            total_duration = time.time() - start_time
            logger.error("[ConnReq:%s] Error validating connection for %s adapter: %s (duration=%.3fs)",
                        request_id, adapter_name, str(e), total_duration)
            return False


# Convenience function to get a database connection
def get_db_connection() -> Optional[ConnectionAdapter]:
    """
    Get a database connection using the ConnectionFactory.
    
    This is a convenience function that uses the ConnectionFactory
    to obtain a database connection with automatic fallback logic.
    
    Returns:
        ConnectionAdapter or None if no connection could be established
    """
    request_id = generate_connection_request_id()
    logger.debug("[ConnReq:%s] get_db_connection() called", request_id)
    
    start_time = time.time()
    factory = ConnectionFactory.get_instance()
    adapter = factory.get_connection()
    
    duration = time.time() - start_time
    if adapter:
        adapter_name = factory.get_active_adapter_name()
        logger.debug("[ConnReq:%s] get_db_connection() returned %s adapter (duration=%.3fs)",
                    request_id, adapter_name, duration)
    else:
        logger.warning("[ConnReq:%s] get_db_connection() failed to return an adapter (duration=%.3fs)",
                      request_id, duration)
    
    return adapter


# Convenience function to close all database connections
def close_all_db_connections() -> bool:
    """
    Close all database connections.
    
    This is a convenience function that uses the ConnectionFactory
    to close all database connections.
    
    Returns:
        bool: True if all disconnections successful, False otherwise
    """
    request_id = generate_connection_request_id()
    logger.debug("[ConnReq:%s] close_all_db_connections() called", request_id)
    
    start_time = time.time()
    factory = ConnectionFactory.get_instance()
    result = factory.close_all_connections()
    
    duration = time.time() - start_time
    if result:
        logger.debug("[ConnReq:%s] close_all_db_connections() succeeded (duration=%.3fs)",
                    request_id, duration)
    else:
        logger.warning("[ConnReq:%s] close_all_db_connections() failed to close some connections (duration=%.3fs)",
                      request_id, duration)
    
    return result


# Convenience function to get connection information
def get_db_connection_info() -> Dict[str, Any]:
    """
    Get connection information for all adapters.
    
    This is a convenience function that uses the ConnectionFactory
    to get connection information for all adapters.
    
    Returns:
        Dict with connection information
    """
    request_id = generate_connection_request_id()
    logger.debug("[ConnReq:%s] get_db_connection_info() called", request_id)
    
    start_time = time.time()
    factory = ConnectionFactory.get_instance()
    info = factory.get_connection_info()
    
    duration = time.time() - start_time
    logger.debug("[ConnReq:%s] get_db_connection_info() completed (duration=%.3fs)",
                request_id, duration)
    
    # Log active adapter if available
    if info.get('active_adapter'):
        logger.debug("[ConnReq:%s] Active adapter: %s",
                    request_id, info['active_adapter'])
    
    return info


# Convenience function to test all connections
def test_all_db_connections() -> Dict[str, Tuple[bool, str]]:
    """
    Test all database connections.
    
    This is a convenience function that uses the ConnectionFactory
    to test all database connections.
    
    Returns:
        Dict with adapter names as keys and (success, message) tuples as values
    """
    request_id = generate_connection_request_id()
    logger.debug("[ConnReq:%s] test_all_db_connections() called", request_id)
    
    start_time = time.time()
    factory = ConnectionFactory.get_instance()
    results = factory.test_all_connections()
    
    duration = time.time() - start_time
    success_count = sum(1 for result in results.values() if result[0])
    
    logger.debug("[ConnReq:%s] test_all_db_connections() completed with %d/%d successful connections (duration=%.3fs)",
                request_id, success_count, len(results), duration)
    
    return results


# Convenience function to check health of all connections
def check_all_db_connections_health() -> Dict[str, Dict[str, Any]]:
    """
    Check health of all database connections.
    
    This is a convenience function that uses the ConnectionFactory
    to check the health of all database connections.
    
    Returns:
        Dict with adapter names as keys and health check results as values
    """
    request_id = generate_connection_request_id()
    logger.debug("[ConnReq:%s] check_all_db_connections_health() called", request_id)
    
    start_time = time.time()
    factory = ConnectionFactory.get_instance()
    results = factory.check_all_connections_health()
    
    duration = time.time() - start_time
    healthy_count = sum(1 for result in results.values() if result.get('healthy', False))
    
    logger.debug("[ConnReq:%s] check_all_db_connections_health() completed with %d/%d healthy connections (duration=%.3fs)",
                request_id, healthy_count, len(results), duration)
    
    # Log any unhealthy adapters
    if healthy_count < len(results):
        unhealthy_adapters = [name for name, result in results.items() if not result.get('healthy', False)]
        logger.warning("[ConnReq:%s] Unhealthy adapters detected: %s",
                      request_id, ', '.join(unhealthy_adapters))
    
    return results


# Convenience function to validate connection before use
def validate_db_connection(adapter: ConnectionAdapter, adapter_name: str) -> bool:
    """
    Validate connection health before use and attempt reconnection if needed.
    
    This is a convenience function that uses the ConnectionFactory
    to validate a connection before use.
    
    Args:
        adapter: The adapter to validate
        adapter_name: Name of the adapter
        
    Returns:
        bool: True if connection is valid or was successfully reconnected, False otherwise
    """
    request_id = generate_connection_request_id()
    logger.debug("[ConnReq:%s] validate_db_connection() called for %s adapter",
                request_id, adapter_name)
    
    start_time = time.time()
    factory = ConnectionFactory.get_instance()
    result = factory.validate_connection_before_use(adapter, adapter_name)
    
    duration = time.time() - start_time
    if result:
        logger.debug("[ConnReq:%s] %s adapter connection validated successfully (duration=%.3fs)",
                    request_id, adapter_name, duration)
    else:
        logger.warning("[ConnReq:%s] %s adapter connection validation failed (duration=%.3fs)",
                      request_id, adapter_name, duration)
    
    return result