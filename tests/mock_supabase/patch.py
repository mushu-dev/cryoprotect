"""
CryoProtect Analyzer - Mock Supabase Patch

This module provides functions to patch the Supabase client in the application.
"""

import unittest
from unittest.mock import patch
from flask import g
from .client import create_mock_client
from .data import reset_mock_data, load_test_data

def patch_supabase_in_app(app, load_data=True):
    """
    Patch the Supabase client in a Flask application.
    
    This function patches the get_supabase_client function in api.utils
    to return a mock Supabase client. It also adds a teardown function
    to reset the mock data after each request.
    
    Args:
        app: Flask application
        load_data: Whether to load test data (default: True)
    """
    # Import here to avoid circular imports
    from api.utils import get_supabase_client as original_get_supabase_client
    
    # Create mock client
    mock_client = create_mock_client()
    
    # Reset and optionally load test data
    reset_mock_data()
    if load_data:
        load_test_data()
    
    # Define patched function
    def patched_get_supabase_client():
        if not hasattr(g, 'supabase'):
            g.supabase = mock_client
        return g.supabase
    
    # Patch the function
    app.config['TESTING'] = True
    app.get_supabase_client = patched_get_supabase_client
    
    # Add teardown function
    @app.teardown_request
    def teardown_request(exception=None):
        if hasattr(g, 'supabase'):
            del g.supabase


def patch_supabase_globally(load_data=True):
    """
    Patch the Supabase client globally.
    
    This function patches the create_client function in the supabase module
    to return a mock Supabase client. This is useful for testing code that
    directly imports and uses the Supabase client.
    
    Args:
        load_data: Whether to load test data (default: True)
        
    Returns:
        Tuple of (patcher, mock_client)
    """
    # Create mock client
    mock_client = create_mock_client()
    
    # Reset and optionally load test data
    reset_mock_data()
    if load_data:
        load_test_data()
    
    # Create patcher
    patcher = patch('supabase.create_client', return_value=mock_client)
    
    # Start patcher
    patcher.start()
    
    return patcher, mock_client


def patch_supabase_in_utils(load_data=True):
    """
    Patch the Supabase client in the utils module.
    
    This function patches the get_supabase_client function in api.utils
    to return a mock Supabase client. This is useful for testing code that
    uses the get_supabase_client function.
    
    Args:
        load_data: Whether to load test data (default: True)
        
    Returns:
        Tuple of (patcher, mock_client)
    """
    # Create mock client
    mock_client = create_mock_client()
    
    # Reset and optionally load test data
    reset_mock_data()
    if load_data:
        load_test_data()
    
    # Create patcher for create_client
    create_client_patcher = patch('supabase.create_client', return_value=mock_client)
    
    # Start patcher
    create_client_patcher.start()
    
    return create_client_patcher, mock_client


class SupabasePatchMixin:
    """Mixin to patch Supabase in test cases."""
    
    def setUp_supabase_patch(self, load_data=True):
        """
        Set up Supabase patching.
        
        Args:
            load_data: Whether to load test data (default: True)
        """
        # Create patcher
        self.supabase_patcher, self.mock_client = patch_supabase_globally(load_data)
    
    def tearDown_supabase_patch(self):
        """Tear down Supabase patching."""
        # Stop patcher
        self.supabase_patcher.stop()