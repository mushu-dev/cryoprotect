"""
Unit tests for the PubChem utility functions.
"""

import time
import unittest
from unittest import mock

from .utils import (
    CircuitBreaker,
    CircuitBreakerError,
    CircuitState,
    retry_with_backoff,
    create_request_key
)

class TestCircuitBreaker(unittest.TestCase):
    """Test cases for CircuitBreaker."""
    
    def setUp(self):
        """Set up test environment."""
        self.circuit_breaker = CircuitBreaker(
            failure_threshold=2,
            recovery_timeout=1,
            expected_exceptions=(ValueError,),
            name="test_circuit"
        )
    
    def test_initial_state(self):
        """Test initial circuit breaker state."""
        self.assertEqual(self.circuit_breaker.state, CircuitState.CLOSED)
        self.assertEqual(self.circuit_breaker.failure_count, 0)
    
    def test_successful_call(self):
        """Test successful function call."""
        # Define a test function
        def test_func():
            return "success"
        
        # Wrap with circuit breaker
        wrapped_func = self.circuit_breaker(test_func)
        
        # Call the function
        result = wrapped_func()
        
        # Check result
        self.assertEqual(result, "success")
        
        # Check circuit breaker state
        self.assertEqual(self.circuit_breaker.state, CircuitState.CLOSED)
        self.assertEqual(self.circuit_breaker.failure_count, 0)
        
        # Check stats
        stats = self.circuit_breaker.get_stats()
        self.assertEqual(stats["stats"]["success_count"], 1)
        self.assertEqual(stats["stats"]["failure_count"], 0)
    
    def test_failed_call(self):
        """Test failed function call."""
        # Define a test function that raises an exception
        def test_func():
            raise ValueError("Test error")
        
        # Wrap with circuit breaker
        wrapped_func = self.circuit_breaker(test_func)
        
        # Call the function (should raise the exception)
        with self.assertRaises(ValueError):
            wrapped_func()
        
        # Check circuit breaker state
        self.assertEqual(self.circuit_breaker.state, CircuitState.CLOSED)
        self.assertEqual(self.circuit_breaker.failure_count, 1)
        
        # Check stats
        stats = self.circuit_breaker.get_stats()
        self.assertEqual(stats["stats"]["success_count"], 0)
        self.assertEqual(stats["stats"]["failure_count"], 1)
    
    def test_circuit_open(self):
        """Test circuit opening after threshold failures."""
        # Define a test function that raises an exception
        def test_func():
            raise ValueError("Test error")
        
        # Wrap with circuit breaker
        wrapped_func = self.circuit_breaker(test_func)
        
        # Call the function multiple times to exceed threshold
        for _ in range(2):
            with self.assertRaises(ValueError):
                wrapped_func()
        
        # Check circuit breaker state
        self.assertEqual(self.circuit_breaker.state, CircuitState.OPEN)
        
        # Call the function again (should raise CircuitBreakerError)
        with self.assertRaises(CircuitBreakerError):
            wrapped_func()
        
        # Check stats
        stats = self.circuit_breaker.get_stats()
        self.assertEqual(stats["stats"]["rejected_count"], 1)
    
    def test_circuit_recovery(self):
        """Test circuit recovery after timeout."""
        # Define a test function that raises an exception
        def test_func():
            raise ValueError("Test error")
        
        # Wrap with circuit breaker
        wrapped_func = self.circuit_breaker(test_func)
        
        # Call the function multiple times to exceed threshold
        for _ in range(2):
            with self.assertRaises(ValueError):
                wrapped_func()
        
        # Check circuit breaker state
        self.assertEqual(self.circuit_breaker.state, CircuitState.OPEN)
        
        # Wait for recovery timeout
        time.sleep(1.1)
        
        # Define a successful function
        def success_func():
            return "success"
        
        # Replace the wrapped function
        self.circuit_breaker.call = lambda func, *args, **kwargs: self.circuit_breaker.call(success_func)
        
        # Call the function again (should enter half-open state and succeed)
        result = wrapped_func()
        
        # Check result
        self.assertEqual(result, "success")
        
        # Check circuit breaker state
        self.assertEqual(self.circuit_breaker.state, CircuitState.CLOSED)
    
    def test_reset(self):
        """Test circuit breaker reset."""
        # Define a test function that raises an exception
        def test_func():
            raise ValueError("Test error")
        
        # Wrap with circuit breaker
        wrapped_func = self.circuit_breaker(test_func)
        
        # Call the function multiple times to exceed threshold
        for _ in range(2):
            with self.assertRaises(ValueError):
                wrapped_func()
        
        # Check circuit breaker state
        self.assertEqual(self.circuit_breaker.state, CircuitState.OPEN)
        
        # Reset the circuit breaker
        self.circuit_breaker.reset()
        
        # Check circuit breaker state
        self.assertEqual(self.circuit_breaker.state, CircuitState.CLOSED)
        self.assertEqual(self.circuit_breaker.failure_count, 0)


class TestRetryWithBackoff(unittest.TestCase):
    """Test cases for retry_with_backoff decorator."""
    
    def test_successful_call(self):
        """Test successful function call."""
        # Define a test function
        @retry_with_backoff(max_retries=3, base_delay=0.1)
        def test_func():
            return "success"
        
        # Call the function
        result = test_func()
        
        # Check result
        self.assertEqual(result, "success")
    
    def test_retry_and_succeed(self):
        """Test function that fails initially but succeeds after retries."""
        # Mock time.sleep to avoid actual delays
        with mock.patch('time.sleep'):
            # Define a test function that fails twice then succeeds
            attempt = [0]
            
            @retry_with_backoff(max_retries=3, base_delay=0.1)
            def test_func():
                attempt[0] += 1
                if attempt[0] <= 2:
                    raise ValueError(f"Attempt {attempt[0]} failed")
                return "success"
            
            # Call the function
            result = test_func()
            
            # Check result
            self.assertEqual(result, "success")
            self.assertEqual(attempt[0], 3)
    
    def test_max_retries_exceeded(self):
        """Test function that always fails and exceeds max retries."""
        # Mock time.sleep to avoid actual delays
        with mock.patch('time.sleep'):
            # Define a test function that always fails
            @retry_with_backoff(max_retries=3, base_delay=0.1)
            def test_func():
                raise ValueError("Always fails")
            
            # Call the function (should raise the exception after max retries)
            with self.assertRaises(ValueError):
                test_func()
    
    def test_backoff_timing(self):
        """Test that backoff timing is calculated correctly."""
        # Define a test function that always fails
        @retry_with_backoff(
            max_retries=3,
            base_delay=1.0,
            max_delay=10.0,
            backoff_factor=2.0,
            jitter=False
        )
        def test_func():
            raise ValueError("Always fails")
        
        # Mock time.sleep to track calls
        with mock.patch('time.sleep') as mock_sleep:
            # Call the function (should raise the exception after max retries)
            with self.assertRaises(ValueError):
                test_func()
            
            # Check that sleep was called with increasing delays
            expected_calls = [mock.call(1.0), mock.call(2.0), mock.call(4.0)]
            self.assertEqual(mock_sleep.call_args_list, expected_calls)


class TestCreateRequestKey(unittest.TestCase):
    """Test cases for create_request_key function."""
    
    def test_url_only(self):
        """Test key creation with URL only."""
        url = "https://example.com/api"
        key = create_request_key(url)
        self.assertEqual(key, url)
    
    def test_url_with_params(self):
        """Test key creation with URL and parameters."""
        url = "https://example.com/api"
        params = {"param1": "value1", "param2": "value2"}
        key = create_request_key(url, params)
        self.assertEqual(key, "https://example.com/api?param1=value1&param2=value2")
    
    def test_params_order(self):
        """Test that parameter order doesn't affect the key."""
        url = "https://example.com/api"
        params1 = {"param1": "value1", "param2": "value2"}
        params2 = {"param2": "value2", "param1": "value1"}
        key1 = create_request_key(url, params1)
        key2 = create_request_key(url, params2)
        self.assertEqual(key1, key2)


if __name__ == "__main__":
    unittest.main()