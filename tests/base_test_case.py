"""
CryoProtect Analyzer - Base Test Case

This module provides a base test case class that sets up the Flask app context
and provides utilities for mocking Supabase.
"""

import unittest
from flask import g, current_app
from app import create_app
from api import init_app
from tests.mock_supabase.helpers import patch_supabase, MockSupabaseTestCase

class BaseTestCase(unittest.TestCase):
    """Base test case with Flask app context."""
    
    @classmethod
    def setUpClass(cls):
        """Set up Flask app and app context once for all tests."""
        # Create the Flask app with testing config
        cls.app = create_app(testing=True)
        
        # Push the app context for the duration of all tests
        cls.app_context = cls.app.app_context()
        cls.app_context.push()
        
        # Create the test client
        cls.client = cls.app.test_client()
    
    @classmethod
    def tearDownClass(cls):
        """Pop the Flask app context after all tests."""
        cls.app_context.pop()

class MockSupabaseBaseTestCase(MockSupabaseTestCase):
    """Base test case with Flask app context and Supabase mocking."""
    
    @classmethod
    def setUpClass(cls):
        """Set up Flask app and app context once for all tests."""
        # Create the Flask app with testing config
        cls.app = create_app(testing=True)
        
        # Push the app context for the duration of all tests
        cls.app_context = cls.app.app_context()
        cls.app_context.push()
        
        # Create the test client
        cls.client = cls.app.test_client()
    
    def setUp(self):
        """Set up test case."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Ensure we have a valid app context
        if not current_app:
            self.app_context = self.app.app_context()
            self.app_context.push()
    
    def tearDown(self):
        """Tear down test case."""
        # Call the parent tearDown
        super().tearDown()
    
    @classmethod
    def tearDownClass(cls):
        """Pop the Flask app context after all tests."""
        cls.app_context.pop()