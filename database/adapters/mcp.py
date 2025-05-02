#!/usr/bin/env python3
"""
MCP adapter implementation for CryoProtect v2.

This module provides a connection adapter for PostgreSQL databases
through the Model Context Protocol (MCP) server.
"""

import os
import logging
import json
from typing import Any, Dict, List, Optional, Union, Tuple

from ..connection import ConnectionAdapter

logger = logging.getLogger(__name__)

class MCPAdapter(ConnectionAdapter):
    """
    MCP-based database adapter implementation.
    
    This adapter provides connection to a PostgreSQL database through
    the Model Context Protocol (MCP) server, which allows for executing
    SQL queries without direct database access.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize adapter with configuration.
        
        Args:
            config: Dictionary containing connection parameters
        """
        self.config = config
        self.server_name = config.get('server_name', 'supabase')
        self.project_id = config.get('project_id')
        self.timeout = int(config.get('timeout', 30))
        self.connected = False
        
        # Import MCP functions
        try:
            # First try to import from the global namespace (if available)
            try:
                from use_mcp_tool import use_mcp_tool
                self.use_mcp_tool = use_mcp_tool
            except ImportError:
                # Fall back to a local implementation
                logger.warning("Could not import use_mcp_tool from global namespace, using local implementation")
                self.use_mcp_tool = self._local_mcp_tool
                
            # Get project ID if not provided
            if not self.project_id:
                self.project_id = self._get_project_id()
        except Exception as e:
            logger.error(f"Failed to initialize MCP tools: {str(e)}")
            self.use_mcp_tool = None
        
    def _local_mcp_tool(self, server_name: str, tool_name: str, arguments: Dict[str, Any]) -> Any:
        """
        Local implementation of use_mcp_tool for testing.
        
        Args:
            server_name: Name of the MCP server
            tool_name: Name of the tool to execute
            arguments: Tool arguments
            
        Returns:
            Tool result
        """
        logger.warning("Using local MCP tool implementation (for testing only)")
        if tool_name == "execute_sql":
            # Mock implementation for testing
            query = arguments.get("query", "")
            if "SELECT 1" in query:
                return [{"test": 1}]
            return []
        elif tool_name == "get_project_id":
            return "test-project-id"
        else:
            raise ValueError(f"Unknown tool: {tool_name}")
            
    def _get_project_id(self) -> str:
        """
        Get Supabase project ID using MCP.
        
        Returns:
            Project ID
        """
        try:
            result = self.use_mcp_tool(
                self.server_name,
                "get_project_id",
                {}
            )
            return result
        except Exception as e:
            logger.error(f"Failed to get project ID: {str(e)}")
            return None
            
    def _execute_sql_through_mcp(self, query: str) -> Any:
        """
        Execute SQL query through MCP.
        
        Args:
            query: SQL query to execute
            
        Returns:
            Query results
        """
        try:
            result = self.use_mcp_tool(
                self.server_name,
                "execute_sql",
                {
                    "project_id": self.project_id,
                    "query": query
                }
            )
            return result
        except Exception as e:
            logger.error(f"Failed to execute SQL through MCP: {str(e)}")
            raise
        
    def connect(self) -> bool:
        """
        Establish connection to MCP.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        try:
            if not self.use_mcp_tool or not self.project_id:
                logger.error("MCP tools not available or project ID not set")
                return False
                
            # Test connection
            result = self._execute_sql_through_mcp("SELECT 1 as test")
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
        
    def reconnect(self) -> bool:
        """
        Re-establish a connection to the MCP server.
        
        Since MCP connections are stateless, this simply tests the connection
        and updates the connection status.
        
        Returns:
            bool: True if reconnection successful, False otherwise
        """
        logger.info("Attempting to reconnect to Supabase via MCP")
        
        # Reset connection status
        self.connected = False
        
        # Attempt to get project ID if not set
        if not self.project_id:
            self.project_id = self._get_project_id()
            if not self.project_id:
                logger.error("Failed to get project ID during reconnection")
                return False
        
        # Test connection
        return self.connect()
            
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """
        Execute SQL query through MCP and return results.
        
        Args:
            query: SQL query to execute
            params: Query parameters (not used, must be embedded in query)
            
        Returns:
            Query results
        """
        if not self.use_mcp_tool or not self.project_id:
            raise ValueError("MCP tools not available or project ID not set")
            
        try:
            # Handle params by formatting query (not ideal but MCP doesn't support params)
            if params:
                if isinstance(params, dict):
                    # Convert dict params to SQL params
                    for key, value in params.items():
                        placeholder = f"%({key})s"
                        if isinstance(value, str):
                            # Escape single quotes for SQL
                            escaped_value = value.replace("'", "''")
                            formatted_value = f"'{escaped_value}'"
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
                        placeholder = "%s"
                        if isinstance(value, str):
                            # Escape single quotes for SQL
                            escaped_value = value.replace("'", "''")
                            formatted_value = f"'{escaped_value}'"
                        elif value is None:
                            formatted_value = "NULL"
                        elif isinstance(value, bool):
                            formatted_value = "TRUE" if value else "FALSE"
                        else:
                            formatted_value = str(value)
                        query = query.replace(placeholder, formatted_value, 1)
            
            result = self._execute_sql_through_mcp(query)
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
        if not self.use_mcp_tool or not self.project_id:
            raise ValueError("MCP tools not available or project ID not set")
            
        try:
            # Combine queries into a transaction
            transaction_queries = ["BEGIN;"] + queries + ["COMMIT;"]
            combined_query = "\n".join(transaction_queries)
            
            result = self._execute_sql_through_mcp(combined_query)
            
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
            
            self._execute_sql_through_mcp(combined_query)
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
            'server_name': self.server_name,
            'project_id': self.project_id,
            'timeout': self.timeout,
            'connected': self.connected
        }
        
    def test_connection(self) -> Tuple[bool, str]:
        """
        Test database connection and return status with message.
        
        Returns:
            Tuple of (success: bool, message: str)
        """
        try:
            if not self.use_mcp_tool or not self.project_id:
                return False, "MCP tools not available or project ID not set"
                
            result = self._execute_sql_through_mcp("SELECT 1 as test")
            if isinstance(result, list) and len(result) > 0 and result[0].get('test') == 1:
                return True, "Connection successful"
            else:
                return False, f"Connection test failed: {result}"
        except Exception as e:
            error_message = str(e)
            return False, f"Connection error: {error_message}"
            
    def is_healthy(self) -> Tuple[bool, Dict[str, Any]]:
        """
        Perform a comprehensive health check on the MCP connection.
        
        This method checks various aspects of the connection health:
        - Basic connectivity to MCP server
        - Connection latency
        - Query execution capability
        - Project ID validity
        
        Returns:
            Tuple of (healthy: bool, health_metrics: Dict[str, Any])
        """
        import time
        
        metrics = {
            'basic_connectivity': False,
            'latency_ms': None,
            'query_capability': False,
            'mcp_status': {
                'server_name': self.server_name,
                'project_id': self.project_id,
                'timeout': self.timeout,
                'connected': self.connected
            }
        }
        
        # Check if MCP tools are available
        if not self.use_mcp_tool:
            metrics['mcp_status']['error'] = "MCP tools not available"
            return False, metrics
            
        # Check if project ID is set
        if not self.project_id:
            metrics['mcp_status']['error'] = "Project ID not set"
            return False, metrics
            
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
            
        # Check query capability with a more complex query
        try:
            # Execute a query that checks database statistics
            result = self.execute_query(
                "SELECT count(*) as table_count FROM information_schema.tables WHERE table_schema = 'public'"
            )
            metrics['query_capability'] = True
            metrics['table_count'] = result[0]['table_count'] if result and len(result) > 0 else 0
        except Exception as e:
            metrics['query_error'] = str(e)
            
        # Check transaction capability
        try:
            # Create a transaction
            transaction = self.begin_transaction()
            
            # Add a query to the transaction
            transaction.append("SELECT 1 as test")
            
            # Commit the transaction
            commit_success = self.commit_transaction(transaction)
            
            metrics['transaction_capability'] = commit_success
        except Exception as e:
            metrics['transaction_error'] = str(e)
            
        # Check project validity by getting project info
        try:
            project_info = self.use_mcp_tool(
                self.server_name,
                "get_project",
                {"id": self.project_id}
            )
            metrics['project_valid'] = project_info is not None
            metrics['project_info'] = {
                'exists': project_info is not None
            }
            if project_info:
                # Extract relevant project information without sensitive data
                if isinstance(project_info, dict):
                    metrics['project_info']['name'] = project_info.get('name')
                    metrics['project_info']['region'] = project_info.get('region')
                    metrics['project_info']['status'] = project_info.get('status')
        except Exception as e:
            metrics['project_error'] = str(e)
            
        # Determine overall health
        is_healthy = (
            metrics['basic_connectivity'] and
            metrics['query_capability'] and
            metrics.get('transaction_capability', False) and
            metrics.get('project_valid', False)
        )
        
        return is_healthy, metrics