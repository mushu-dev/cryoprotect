"""
Unit tests for the pagination utilities module.

These tests verify that the pagination utilities correctly:
- Parse pagination parameters from requests
- Calculate pagination metadata
- Format paginated responses
- Handle edge cases and validation
"""

import unittest
from unittest.mock import patch, MagicMock
import json
from flask import Flask, request
import pytest

from api.pagination_utils import (
    parse_pagination_params,
    get_pagination_metadata,
    paginate_query_response,
    get_total_count
)


class TestPaginationUtils(unittest.TestCase):
    """Test suite for pagination utilities."""
    
    def setUp(self):
        """Set up test environment."""
        self.app = Flask(__name__)
        self.app_context = self.app.app_context()
        self.app_context.push()
        self.request_context = None
    
    def tearDown(self):
        """Clean up test environment."""
        if self.request_context:
            self.request_context.pop()
        self.app_context.pop()
    
    def test_parse_pagination_params_default(self):
        """Test parsing pagination parameters with default values."""
        with self.app.test_request_context('/?'):
            limit, offset = parse_pagination_params()
            self.assertEqual(limit, 100)
            self.assertEqual(offset, 0)
    
    def test_parse_pagination_params_custom(self):
        """Test parsing pagination parameters with custom values."""
        with self.app.test_request_context('/?limit=50&offset=150'):
            limit, offset = parse_pagination_params()
            self.assertEqual(limit, 50)
            self.assertEqual(offset, 150)
    
    def test_parse_pagination_params_max_limit(self):
        """Test that limit is capped at max_limit."""
        with self.app.test_request_context('/?limit=1000'):
            limit, offset = parse_pagination_params(max_limit=500)
            self.assertEqual(limit, 500)
            self.assertEqual(offset, 0)
    
    def test_parse_pagination_params_page_per_page(self):
        """Test parsing pagination using page and per_page parameters."""
        with self.app.test_request_context('/?page=3&per_page=25'):
            limit, offset = parse_pagination_params()
            self.assertEqual(limit, 25)
            self.assertEqual(offset, 50)  # (3-1) * 25 = 50
    
    def test_parse_pagination_params_invalid_limit(self):
        """Test handling invalid limit parameter."""
        with self.app.test_request_context('/?limit=-10'):
            with self.assertRaises(ValueError):
                parse_pagination_params()
    
    def test_parse_pagination_params_invalid_offset(self):
        """Test handling invalid offset parameter."""
        with self.app.test_request_context('/?offset=-5'):
            with self.assertRaises(ValueError):
                parse_pagination_params()
    
    def test_parse_pagination_params_invalid_page(self):
        """Test handling invalid page parameter."""
        with self.app.test_request_context('/?page=-1&per_page=10'):
            limit, offset = parse_pagination_params()
            self.assertEqual(limit, 10)
            self.assertEqual(offset, 0)  # Defaults to 0 for invalid page
    
    def test_get_pagination_metadata_first_page(self):
        """Test generating pagination metadata for the first page."""
        metadata = get_pagination_metadata(
            total_count=100,
            limit=10,
            offset=0,
            resource_endpoint="/api/molecules"
        )
        
        self.assertEqual(metadata["total_items"], 100)
        self.assertEqual(metadata["total_pages"], 10)
        self.assertEqual(metadata["current_page"], 1)
        self.assertEqual(metadata["per_page"], 10)
        self.assertEqual(metadata["offset"], 0)
        
        # Check links
        self.assertIn("first", metadata["links"])
        self.assertNotIn("prev", metadata["links"])
        self.assertIn("next", metadata["links"])
        self.assertIn("last", metadata["links"])
        
        # Verify link URLs
        self.assertIn("/api/molecules?limit=10&offset=0", metadata["links"]["first"])
        self.assertIn("/api/molecules?limit=10&offset=10", metadata["links"]["next"])
        self.assertIn("/api/molecules?limit=10&offset=90", metadata["links"]["last"])
    
    def test_get_pagination_metadata_middle_page(self):
        """Test generating pagination metadata for a middle page."""
        metadata = get_pagination_metadata(
            total_count=100,
            limit=10,
            offset=30,
            resource_endpoint="/api/molecules"
        )
        
        self.assertEqual(metadata["total_items"], 100)
        self.assertEqual(metadata["total_pages"], 10)
        self.assertEqual(metadata["current_page"], 4)
        self.assertEqual(metadata["per_page"], 10)
        self.assertEqual(metadata["offset"], 30)
        
        # Check links
        self.assertIn("first", metadata["links"])
        self.assertIn("prev", metadata["links"])
        self.assertIn("next", metadata["links"])
        self.assertIn("last", metadata["links"])
        
        # Verify link URLs
        self.assertIn("/api/molecules?limit=10&offset=0", metadata["links"]["first"])
        self.assertIn("/api/molecules?limit=10&offset=20", metadata["links"]["prev"])
        self.assertIn("/api/molecules?limit=10&offset=40", metadata["links"]["next"])
        self.assertIn("/api/molecules?limit=10&offset=90", metadata["links"]["last"])
    
    def test_get_pagination_metadata_last_page(self):
        """Test generating pagination metadata for the last page."""
        metadata = get_pagination_metadata(
            total_count=100,
            limit=10,
            offset=90,
            resource_endpoint="/api/molecules"
        )
        
        self.assertEqual(metadata["total_items"], 100)
        self.assertEqual(metadata["total_pages"], 10)
        self.assertEqual(metadata["current_page"], 10)
        self.assertEqual(metadata["per_page"], 10)
        self.assertEqual(metadata["offset"], 90)
        
        # Check links
        self.assertIn("first", metadata["links"])
        self.assertIn("prev", metadata["links"])
        self.assertNotIn("next", metadata["links"])
        self.assertIn("last", metadata["links"])
        
        # Verify link URLs
        self.assertIn("/api/molecules?limit=10&offset=0", metadata["links"]["first"])
        self.assertIn("/api/molecules?limit=10&offset=80", metadata["links"]["prev"])
        self.assertIn("/api/molecules?limit=10&offset=90", metadata["links"]["last"])
    
    def test_get_pagination_metadata_with_query_params(self):
        """Test generating pagination metadata with additional query parameters."""
        metadata = get_pagination_metadata(
            total_count=100,
            limit=10,
            offset=0,
            resource_endpoint="/api/molecules",
            query_params={"name": "aspirin", "sort": "name"}
        )
        
        # Verify query parameters are included in links
        self.assertIn("name=aspirin", metadata["links"]["first"])
        self.assertIn("sort=name", metadata["links"]["first"])
        self.assertIn("name=aspirin", metadata["links"]["next"])
        self.assertIn("sort=name", metadata["links"]["next"])
    
    def test_paginate_query_response(self):
        """Test formatting a paginated response."""
        data = [{"id": 1, "name": "Item 1"}, {"id": 2, "name": "Item 2"}]
        
        response = paginate_query_response(
            query_result=data,
            total_count=100,
            limit=10,
            offset=0,
            resource_endpoint="/api/items"
        )
        
        # Check response structure
        self.assertIn("data", response)
        self.assertIn("pagination", response)
        
        # Check data
        self.assertEqual(response["data"], data)
        
        # Check pagination metadata
        self.assertEqual(response["pagination"]["total_items"], 100)
        self.assertEqual(response["pagination"]["total_pages"], 10)
        self.assertEqual(response["pagination"]["current_page"], 1)
    
    @patch('api.pagination_utils.supabase')
    def test_get_total_count(self, mock_supabase):
        """Test getting total count from a table."""
        # Mock the Supabase client and response
        mock_query = MagicMock()
        mock_query.select.return_value = mock_query
        mock_query.execute.return_value = MagicMock(count=42)
        
        mock_supabase.table.return_value = mock_query
        
        # Call the function
        count = get_total_count(mock_supabase, "test_table")
        
        # Verify the result
        self.assertEqual(count, 42)
        
        # Verify the correct methods were called
        mock_supabase.table.assert_called_once_with("test_table")
        mock_query.select.assert_called_once_with("*", count="exact")
    
    @patch('api.pagination_utils.supabase')
    def test_get_total_count_with_filters(self, mock_supabase):
        """Test getting total count with filters."""
        # Mock the Supabase client and response
        mock_query = MagicMock()
        mock_query.select.return_value = mock_query
        mock_query.eq.return_value = mock_query
        mock_query.execute.return_value = MagicMock(count=10)
        
        mock_supabase.table.return_value = mock_query
        
        # Call the function with filters
        filters = {"status": "active"}
        count = get_total_count(mock_supabase, "test_table", filters)
        
        # Verify the result
        self.assertEqual(count, 10)
        
        # Verify the correct methods were called
        mock_supabase.table.assert_called_once_with("test_table")
        mock_query.select.assert_called_once_with("*", count="exact")
        mock_query.eq.assert_called_once_with("status", "active")


if __name__ == '__main__':
    unittest.main()