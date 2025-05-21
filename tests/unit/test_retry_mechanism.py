#!/usr/bin/env python3
"""
Unit tests for the retry mechanism.

These tests verify that the retry mechanism correctly retries operations
that fail with the specified exceptions, and gives up after the maximum
number of retries.
"""

import unittest
import time
from unittest.mock import patch, MagicMock
from api.resiliency.retry import retry_with_backoff

class TestRetryMechanism(unittest.TestCase):
    """Tests for the retry mechanism."""
    
    def test_successful_operation(self):
        """Test that a successful operation is not retried."""
        # Create a mock function that succeeds
        mock_func = MagicMock(return_value="success")
        
        # Apply the retry decorator
        decorated_func = retry_with_backoff()(mock_func)
        
        # Call the decorated function
        result = decorated_func()
        
        # Verify that the function was called exactly once
        mock_func.assert_called_once()
        
        # Verify that the result is correct
        self.assertEqual(result, "success")
    
    def test_failed_operation_retries(self):
        """Test that a failed operation is retried."""
        # Create a mock function that fails twice then succeeds
        mock_func = MagicMock(side_effect=[
            ConnectionError("First failure"),
            ConnectionError("Second failure"),
            "success"
        ])
        
        # Apply the retry decorator with quick retry times for testing
        decorated_func = retry_with_backoff(
            max_retries=3,
            base_delay=0.01,
            max_delay=0.1,
            backoff_factor=1
        )(mock_func)
        
        # Call the decorated function
        result = decorated_func()
        
        # Verify that the function was called exactly three times
        self.assertEqual(mock_func.call_count, 3)
        
        # Verify that the result is correct
        self.assertEqual(result, "success")
    
    def test_max_retries_reached(self):
        """Test that the function gives up after max_retries."""
        # Create a mock function that always fails
        mock_func = MagicMock(side_effect=ConnectionError("Always fails"))
        
        # Apply the retry decorator with quick retry times for testing
        decorated_func = retry_with_backoff(
            max_retries=3,
            base_delay=0.01,
            max_delay=0.1,
            backoff_factor=1
        )(mock_func)
        
        # Call the decorated function and expect an exception
        with self.assertRaises(ConnectionError):
            decorated_func()
        
        # Verify that the function was called exactly four times (initial + 3 retries)
        self.assertEqual(mock_func.call_count, 4)
    
    def test_specific_exceptions(self):
        """Test that only specified exceptions trigger retries."""
        # Create a mock function that raises different exceptions
        mock_func = MagicMock(side_effect=[
            ConnectionError("Retryable error"),
            ValueError("Non-retryable error"),
            "success"
        ])
        
        # Apply the retry decorator with specific exceptions
        decorated_func = retry_with_backoff(
            max_retries=3,
            exceptions=(ConnectionError,),
            base_delay=0.01,
            max_delay=0.1,
            backoff_factor=1
        )(mock_func)
        
        # Call the decorated function and expect ValueError
        with self.assertRaises(ValueError):
            decorated_func()
        
        # Verify that the function was called exactly twice
        # (once for ConnectionError which is retried, and once for ValueError which is not)
        self.assertEqual(mock_func.call_count, 2)
    
    def test_exponential_backoff(self):
        """Test that retries use exponential backoff."""
        # Create a mock sleep function to capture delay times
        mock_sleep = MagicMock()
        
        # Create a mock function that always fails
        mock_func = MagicMock(side_effect=[
            ConnectionError("First failure"),
            ConnectionError("Second failure"),
            ConnectionError("Third failure"),
            ConnectionError("Fourth failure")
        ])
        
        # Apply the retry decorator
        decorated_func = retry_with_backoff(
            max_retries=3,
            base_delay=0.1,
            max_delay=10.0,
            backoff_factor=2,
            jitter=0  # Disable jitter for predictable testing
        )(mock_func)
        
        # Patch time.sleep to record delays
        with patch('time.sleep', mock_sleep):
            # Call the decorated function and expect an exception
            with self.assertRaises(ConnectionError):
                decorated_func()
        
        # Verify that the function was called exactly four times
        self.assertEqual(mock_func.call_count, 4)
        
        # Verify that sleep was called with exponentially increasing delays
        mock_sleep.assert_any_call(0.1)  # First retry: base_delay
        mock_sleep.assert_any_call(0.2)  # Second retry: base_delay * backoff_factor
        mock_sleep.assert_any_call(0.4)  # Third retry: base_delay * backoff_factor^2
    
    def test_giveup_after(self):
        """Test that the function gives up after the specified time limit."""
        # Create a mock function that always fails
        mock_func = MagicMock(side_effect=ConnectionError("Always fails"))
        
        # Apply the retry decorator with a time limit
        decorated_func = retry_with_backoff(
            max_retries=10,  # High value to ensure giveup_after triggers first
            base_delay=0.1,
            giveup_after=0.5  # Give up after 0.5 seconds
        )(mock_func)
        
        # Mock time.time to return increasing values
        with patch('time.time', side_effect=[0, 0.2, 0.4, 0.6, 0.8]):
            # Call the decorated function and expect an exception
            with self.assertRaises(ConnectionError):
                decorated_func()
        
        # Verify that the function was called at most 3 times
        # (0.6 > 0.5, so should give up after the third call)
        self.assertLessEqual(mock_func.call_count, 3)
    
    def test_on_retry_callback(self):
        """Test that the on_retry callback is called before each retry."""
        # Create a mock function that fails twice then succeeds
        mock_func = MagicMock(side_effect=[
            ConnectionError("First failure"),
            ConnectionError("Second failure"),
            "success"
        ])
        
        # Create a mock callback
        mock_callback = MagicMock()
        
        # Apply the retry decorator with the callback
        decorated_func = retry_with_backoff(
            max_retries=3,
            base_delay=0.01,
            on_retry=mock_callback
        )(mock_func)
        
        # Call the decorated function
        result = decorated_func()
        
        # Verify that the function was called exactly three times
        self.assertEqual(mock_func.call_count, 3)
        
        # Verify that the callback was called exactly twice (once before each retry)
        self.assertEqual(mock_callback.call_count, 2)
        
        # Verify that the result is correct
        self.assertEqual(result, "success")
    
    def test_retry_condition(self):
        """Test that the retry_condition function controls when to retry."""
        # Create a mock function that raises different ConnectionErrors
        mock_func = MagicMock(side_effect=[
            ConnectionError("Retry this"),
            ConnectionError("Don't retry this"),
            ConnectionError("Retry this"),
            "success"
        ])
        
        # Create a retry condition that only retries on specific error messages
        def retry_condition(e):
            return "Retry this" in str(e)
        
        # Apply the retry decorator with the condition
        decorated_func = retry_with_backoff(
            max_retries=3,
            base_delay=0.01,
            retry_condition=retry_condition
        )(mock_func)
        
        # Call the decorated function and expect the non-retryable error
        with self.assertRaises(ConnectionError) as cm:
            decorated_func()
        
        # Verify that the function was called exactly twice
        self.assertEqual(mock_func.call_count, 2)
        
        # Verify that the correct error was raised
        self.assertEqual(str(cm.exception), "Don't retry this")

if __name__ == '__main__':
    unittest.main()