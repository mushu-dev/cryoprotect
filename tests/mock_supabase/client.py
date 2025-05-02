"""
CryoProtect Analyzer - Mock Supabase Client

This module provides a mock implementation of the Supabase client.
"""

from unittest.mock import MagicMock
from .query import MockQueryBuilder
from .rpc import MockRPC
from .auth import MockAuth

class MockSupabaseClient:
    """Mock implementation of the Supabase client."""
    
    def __init__(self, url, key, options=None):
        self.url = url
        self.key = key
        self.options = options or {}
        self.auth = MockAuth()
        # Add rpcs attribute for test compatibility
        self.rpcs = {
            'calculate_cryoprotectant_score': lambda *args, **kwargs: 85.5,
            'import_molecule_from_pubchem': lambda *args, **kwargs: {},
            'compare_prediction_with_experiment': lambda *args, **kwargs: {},
            # Add more RPCs as needed for test coverage
        }
        
        # Create MagicMock objects for table and rpc methods
        self._table_mock = MagicMock()
        self._rpc_mock = MagicMock()
    
    def table(self, table_name):
        """
        Get a query builder for a table.
        
        Args:
            table_name: Name of the table
            
        Returns:
            MagicMock that can be used for chaining method calls
        """
        # In real tests, we'll use the mock's return_value
        # But for backward compatibility, we'll also set up the actual implementation
        actual_builder = MockQueryBuilder(table_name)
        
        # Configure the mock to return the actual builder when needed
        self._table_mock.actual_builder = actual_builder
        
        # Return the mock for chaining
        return self._table_mock
    
    def rpc(self, function_name, params=None):
        """
        Call a stored procedure.
        
        Args:
            function_name: Name of the function to call
            params: Parameters to pass to the function
            
        Returns:
            MagicMock that can be used for chaining method calls
        """
        # In real tests, we'll use the mock's return_value
        # But for backward compatibility, we'll also set up the actual implementation
        actual_rpc = MockRPC(function_name, params)
        
        # Configure the mock to return the actual RPC when needed
        self._rpc_mock.actual_rpc = actual_rpc
        self._rpc_mock.function_name = function_name
        self._rpc_mock.params = params
        
        # Return the mock for chaining
        return self._rpc_mock
    
    def from_(self, table_name):
        """
        Alias for table().
        
        Args:
            table_name: Name of the table
            
        Returns:
            MagicMock for the table
        """
        return self.table(table_name)


def create_mock_client(url="https://mock.supabase.co", key="mock-key", options=None):
    """
    Create a mock Supabase client.
    
    Args:
        url: Supabase URL (default: "https://mock.supabase.co")
        key: Supabase key (default: "mock-key")
        options: Additional options
        
    Returns:
        MockSupabaseClient
    """
    return MockSupabaseClient(url, key, options)