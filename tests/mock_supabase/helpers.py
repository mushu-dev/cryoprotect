"""
CryoProtect Analyzer - Mock Supabase Helpers

This module provides helper functions for using the mock Supabase client in tests.
"""

import unittest
from unittest.mock import patch, MagicMock
from .client import create_mock_client
from .data import reset_mock_data, load_test_data

def patch_supabase(load_data=True):
    """
    Decorator to patch the Supabase client in tests.
    
    Args:
        load_data: Whether to load test data (default: True)
        
    Returns:
        Decorator function
    """
    def decorator(test_method):
        @patch('api.utils.create_client')
        def wrapper(self, mock_create_client, *args, **kwargs):
            # Create a MagicMock instead of the actual client
            mock_client = MagicMock()
            
            # Set up common mock methods
            mock_table = MagicMock()
            mock_client.table.return_value = mock_table
            
            mock_rpc = MagicMock()
            mock_client.rpc.return_value = mock_rpc
            
            # Set the mock client as the return value
            mock_create_client.return_value = mock_client
            
            # Reset and optionally load test data
            reset_mock_data()
            if load_data:
                load_test_data()
            
            # Call the test method
            return test_method(self, mock_client, *args, **kwargs)
        
        return wrapper
    
    return decorator


class MockSupabaseTestCase(unittest.TestCase):
    """Base test case with Supabase mocking."""
    
    def setUp(self):
        """Set up the test case."""
        # Create patcher
        self.supabase_patcher = patch('api.utils.create_client')
        
        # Start patcher and get mock
        self.mock_create_client = self.supabase_patcher.start()
        
        # Create a MagicMock instead of the actual client
        self.mock_client = MagicMock()
        
        # Set up common mock methods
        mock_table = MagicMock()
        self.mock_client.table.return_value = mock_table
        
        mock_rpc = MagicMock()
        self.mock_client.rpc.return_value = mock_rpc
        
        # Set the mock client as the return value
        self.mock_create_client.return_value = self.mock_client
        
        # Reset and load test data
        reset_mock_data()
        load_test_data()
    
    def tearDown(self):
        """Tear down the test case."""
        # Stop patcher
        self.supabase_patcher.stop()


def mock_rpc_function(name, result=None, error=None):
    """
    Create a mock RPC function.
    
    Args:
        name: Name of the RPC function
        result: Result to return
        error: Error to return
        
    Returns:
        Mock RPC function
    """
    from .rpc import register_rpc_function
    
    def mock_func(params):
        if error:
            raise Exception(error)
        return result
    
    register_rpc_function(name, mock_func)
    return mock_func


def configure_test_scenario(scenario='success'):
    """
    Configure a test scenario.
    
    Args:
        scenario: Scenario to configure ('success', 'failure', 'timeout')
    """
    from .data import _mock_rpcs
    
    if scenario == 'failure':
        # Configure all RPCs to fail
        for name in list(_mock_rpcs.keys()):
            mock_rpc_function(name, error="Simulated failure")
    
    elif scenario == 'timeout':
        # Configure all RPCs to timeout
        for name in list(_mock_rpcs.keys()):
            def timeout_func(params):
                import time
                time.sleep(10)  # Long sleep to simulate timeout
                return None
            
            _mock_rpcs[name] = timeout_func
    
    # For 'success', we don't need to do anything as the default behavior is success