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
            error_message = str(e)
            return False, f"Connection error: {error_message}"