"""
CryoProtect Analyzer - Mock Supabase

This module provides mock implementations of the Supabase client for offline testing.
It simulates database operations without requiring an actual Supabase connection.
"""

from .client import MockSupabaseClient, create_mock_client
from .data import reset_mock_data, load_test_data, DEFAULT_USER_ID
from .query import MockQueryBuilder, MockResponse
from .rpc import MockRPC, register_rpc_function
from .auth import MockAuth, MockMFA
from .helpers import (
    patch_supabase, 
    MockSupabaseTestCase, 
    mock_rpc_function, 
    configure_test_scenario
)
from .patch import (
    patch_supabase_in_app,
    patch_supabase_globally,
    patch_supabase_in_utils,
    SupabasePatchMixin
)

__all__ = [
    # Client
    'MockSupabaseClient',
    'create_mock_client',
    
    # Data
    'reset_mock_data',
    'load_test_data',
    'DEFAULT_USER_ID',
    
    # Query
    'MockQueryBuilder',
    'MockResponse',
    
    # RPC
    'MockRPC',
    'register_rpc_function',
    
    # Auth
    'MockAuth',
    'MockMFA',
    
    # Helpers
    'patch_supabase',
    'MockSupabaseTestCase',
    'mock_rpc_function',
    'configure_test_scenario',
    
    # Patch
    'patch_supabase_in_app',
    'patch_supabase_globally',
    'patch_supabase_in_utils',
    'SupabasePatchMixin'
]