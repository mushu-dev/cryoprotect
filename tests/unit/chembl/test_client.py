"""
Unit tests for the enhanced ChEMBL client.

Tests the weekday-specific rate limiting, resilient caching, and circuit breaker pattern.
"""

import unittest
import time
import json
import os
import shutil
import tempfile
from datetime import datetime
from unittest import mock
from unittest.mock import patch, MagicMock, PropertyMock

import requests
import pytest

from chembl.client import ResilientChEMBLClient
from chembl.utils import CircuitBreakerError


class TestResilientChEMBLClient(unittest.TestCase):
    """Test the ResilientChEMBLClient class."""

    def setUp(self):
        """Set up the test environment."""
        # Create a temporary directory for cache
        self.temp_dir = tempfile.mkdtemp()
        
        # Initialize client with test settings
        self.client = ResilientChEMBLClient(
            cache_dir=self.temp_dir,
            weekday_requests_per_second=10.0,  # High value for faster tests
            weekend_requests_per_second=15.0,
            monday_requests_per_second=5.0,
            max_retries=2,
            failure_threshold=2,
            recovery_timeout=1,  # Short timeout for faster tests
            cache_ttl=3600,  # 1 hour for faster tests
            error_cache_ttl=60,  # 1 minute for faster tests
            memory_cache_size=10
        )

    def tearDown(self):
        """Clean up the test environment."""
        # Remove the temporary directory
        shutil.rmtree(self.temp_dir)

    @patch('chembl.rate_limiter.datetime')
    def test_monday_rate_limiting(self, mock_datetime):
        """Test that stricter rate limits are applied on Mondays."""
        # Mock Monday (weekday 0)
        mock_datetime.now.return_value = MagicMock(weekday=lambda: 0)
        
        # Get the rate limiter stats before any requests
        stats_before = self.client.get_rate_limiter_stats()
        
        # Make a request
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = {"molecules": []}
            mock_get.return_value = mock_response
            
            start_time = time.time()
            self.client._make_request("molecule", {"limit": 1})
            
            # Get the rate limiter stats after the request
            stats_after = self.client.get_rate_limiter_stats()
            
            # Verify that Monday rate limit was applied
            self.assertEqual(stats_after['day_type'], 'monday')
            self.assertEqual(stats_after['requests_per_second']['monday'], 5.0)
            
            # Verify the current delay is using the Monday setting
            self.assertAlmostEqual(stats_after['current_delay'], 1.0/5.0, places=2)

    @patch('chembl.rate_limiter.datetime')
    def test_weekday_rate_limiting(self, mock_datetime):
        """Test that regular weekday rate limits are applied Tuesday-Friday."""
        # Mock Tuesday (weekday 1)
        mock_datetime.now.return_value = MagicMock(weekday=lambda: 1)
        
        # Make a request
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = {"molecules": []}
            mock_get.return_value = mock_response
            
            self.client._make_request("molecule", {"limit": 1})
            
            # Get the rate limiter stats after the request
            stats = self.client.get_rate_limiter_stats()
            
            # Verify that weekday rate limit was applied
            self.assertEqual(stats['day_type'], 'weekday')
            self.assertEqual(stats['requests_per_second']['weekday'], 10.0)
            
            # Verify the current delay is using the weekday setting
            self.assertAlmostEqual(stats['current_delay'], 1.0/10.0, places=2)

    @patch('chembl.rate_limiter.datetime')
    def test_weekend_rate_limiting(self, mock_datetime):
        """Test that weekend rate limits are applied on Saturday-Sunday."""
        # Mock Saturday (weekday 5)
        mock_datetime.now.return_value = MagicMock(weekday=lambda: 5)
        
        # Make a request
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = {"molecules": []}
            mock_get.return_value = mock_response
            
            self.client._make_request("molecule", {"limit": 1})
            
            # Get the rate limiter stats after the request
            stats = self.client.get_rate_limiter_stats()
            
            # Verify that weekend rate limit was applied
            self.assertEqual(stats['day_type'], 'weekend')
            self.assertEqual(stats['requests_per_second']['weekend'], 15.0)
            
            # Verify the current delay is using the weekend setting
            self.assertAlmostEqual(stats['current_delay'], 1.0/15.0, places=2)

    def test_resilient_caching_success(self):
        """Test that successful responses are cached correctly."""
        # Mock a successful response
        test_data = {"molecule_chembl_id": "CHEMBL25", "pref_name": "Aspirin"}
        
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = test_data
            mock_get.return_value = mock_response
            
            # First request should hit the API
            result1 = self.client.get_molecule_by_chembl_id("CHEMBL25")
            
            # Second request should hit the cache
            result2 = self.client.get_molecule_by_chembl_id("CHEMBL25")
            
            # Verify that the API was only called once
            self.assertEqual(mock_get.call_count, 1)
            
            # Verify that both results are the same
            self.assertEqual(result1, result2)
            
            # Verify cache stats
            cache_stats = self.client.get_cache_stats()
            self.assertEqual(cache_stats["memory_hits"], 1)

    def test_resilient_caching_error(self):
        """Test that error responses are cached with a shorter TTL."""
        # First request will fail
        with patch('requests.get') as mock_get:
            mock_get.side_effect = requests.RequestException("API Error")
            
            # This should fail and cache the error
            result1 = self.client.get_molecule_by_chembl_id("CHEMBL25")
            
            # Verify that it's an error response
            self.assertIn("Error", result1)
            
            # Second request should hit the error cache
            result2 = self.client.get_molecule_by_chembl_id("CHEMBL25")
            
            # Verify that the API was only called once
            self.assertEqual(mock_get.call_count, 1)
            
            # Verify that both results are the same
            self.assertEqual(result1["Error"], result2["Error"])
            
            # Verify cache stats
            cache_stats = self.client.get_cache_stats()
            self.assertGreaterEqual(cache_stats["errors"], 1)

    def test_circuit_breaker_opens(self):
        """Test that the circuit breaker opens after multiple failures."""
        # Mock multiple failures
        with patch('requests.get') as mock_get:
            mock_get.side_effect = requests.RequestException("API Error")
            
            # Make enough requests to open the circuit breaker
            for _ in range(3):  # More than failure_threshold
                try:
                    self.client._make_request("molecule", {"limit": 1})
                except Exception:
                    pass
            
            # Verify that the circuit is open
            circuit_stats = self.client.get_circuit_breaker_stats()
            self.assertEqual(circuit_stats["state"], "open")

    def test_circuit_breaker_fallback_to_cache(self):
        """Test that when circuit is open, client falls back to cache."""
        # First, cache a successful response
        test_data = {"molecule_chembl_id": "CHEMBL25", "pref_name": "Aspirin"}
        
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = test_data
            mock_get.return_value = mock_response
            
            # Cache the data
            self.client.get_molecule_by_chembl_id("CHEMBL25")
        
        # Now force the circuit to open
        self.client.circuit_breaker.state = "open"
        
        # Try to get the data again
        with patch('requests.get') as mock_get:
            # This should not be called because circuit is open
            result = self.client.get_molecule_by_chembl_id("CHEMBL25")
            
            # Verify that the API was not called
            mock_get.assert_not_called()
            
            # Verify that we got the cached data
            self.assertEqual(result["ChEMBL ID"], "CHEMBL25")

    def test_circuit_breaker_half_open(self):
        """Test that the circuit breaker transitions to half-open after timeout."""
        # Set up the circuit breaker
        self.client.circuit_breaker.state = "open"
        self.client.circuit_breaker.last_failure_time = time.time() - 2  # More than recovery_timeout
        
        # Mock a successful response
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = {"molecule_chembl_id": "CHEMBL25"}
            mock_get.return_value = mock_response
            
            # This should transition to half-open and then closed
            self.client._make_request("molecule", {"limit": 1})
            
            # Verify that the circuit is closed
            circuit_stats = self.client.get_circuit_breaker_stats()
            self.assertEqual(circuit_stats["state"], "closed")

    def test_reset_circuit_breaker(self):
        """Test manually resetting the circuit breaker."""
        # Set the circuit breaker to open
        self.client.circuit_breaker.state = "open"
        self.client.circuit_breaker.failures = 5
        
        # Reset the circuit breaker
        self.client.reset_circuit_breaker()
        
        # Verify that the circuit is closed
        circuit_stats = self.client.get_circuit_breaker_stats()
        self.assertEqual(circuit_stats["state"], "closed")
        self.assertEqual(circuit_stats["failures"], 0)


if __name__ == '__main__':
    unittest.main()