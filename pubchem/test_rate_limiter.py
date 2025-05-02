"""
Unit tests for the PubChem rate limiter module.
"""

import time
import unittest
import datetime
from unittest import mock

from .rate_limiter import AdaptiveRateLimiter

class TestAdaptiveRateLimiter(unittest.TestCase):
    """Test cases for AdaptiveRateLimiter."""
    
    def setUp(self):
        """Set up test environment."""
        self.weekday_rps = 2.0
        self.weekend_rps = 5.0
        self.rate_limiter = AdaptiveRateLimiter(
            weekday_requests_per_second=self.weekday_rps,
            weekend_requests_per_second=self.weekend_rps,
            burst_size=5
        )
    
    def test_is_weekend(self):
        """Test weekend detection."""
        # Mock datetime.now to return a weekday (Wednesday)
        wednesday = datetime.datetime(2025, 4, 23)  # A Wednesday
        with mock.patch('datetime.datetime') as mock_datetime:
            mock_datetime.now.return_value = wednesday
            self.assertFalse(self.rate_limiter._is_weekend())
        
        # Mock datetime.now to return a weekend (Saturday)
        saturday = datetime.datetime(2025, 4, 26)  # A Saturday
        with mock.patch('datetime.datetime') as mock_datetime:
            mock_datetime.now.return_value = saturday
            self.assertTrue(self.rate_limiter._is_weekend())
    
    def test_get_current_delay_weekday(self):
        """Test delay calculation on weekdays."""
        # Mock _is_weekend to return False (weekday)
        with mock.patch.object(self.rate_limiter, '_is_weekend', return_value=False):
            # Mock _is_weekday_hours to return True
            with mock.patch.object(self.rate_limiter, '_is_weekday_hours', return_value=True):
                delay = self.rate_limiter._get_current_delay()
                self.assertEqual(delay, 1.0 / self.weekday_rps)
    
    def test_get_current_delay_weekend(self):
        """Test delay calculation on weekends."""
        # Mock _is_weekend to return True (weekend)
        with mock.patch.object(self.rate_limiter, '_is_weekend', return_value=True):
            delay = self.rate_limiter._get_current_delay()
            self.assertEqual(delay, 1.0 / self.weekend_rps)
    
    def test_wait_no_delay_with_tokens(self):
        """Test wait with available tokens."""
        # Set up rate limiter with tokens
        self.rate_limiter.tokens = 5
        
        # Mock time.sleep to track calls
        with mock.patch('time.sleep') as mock_sleep:
            # Call wait
            start_time = time.time()
            self.rate_limiter.wait()
            elapsed_time = time.time() - start_time
            
            # Check that sleep was not called
            mock_sleep.assert_not_called()
            
            # Check that elapsed time is minimal
            self.assertLess(elapsed_time, 0.1)
            
            # Check that a token was consumed
            self.assertEqual(self.rate_limiter.tokens, 4)
    
    def test_wait_with_delay(self):
        """Test wait with no tokens available."""
        # Set up rate limiter with no tokens
        self.rate_limiter.tokens = 0
        
        # Mock _get_current_delay to return a fixed delay
        test_delay = 0.1
        with mock.patch.object(self.rate_limiter, '_get_current_delay', return_value=test_delay):
            # Mock time.sleep to track calls
            with mock.patch('time.sleep') as mock_sleep:
                # Call wait
                self.rate_limiter.wait()
                
                # Check that sleep was called with the correct delay
                mock_sleep.assert_called_once_with(test_delay)
    
    def test_burst_handling(self):
        """Test handling of request bursts."""
        # Set up rate limiter with burst size
        burst_size = 5
        rate_limiter = AdaptiveRateLimiter(
            weekday_requests_per_second=10.0,  # High RPS to minimize natural delays
            weekend_requests_per_second=10.0,
            burst_size=burst_size
        )
        
        # Mock time.sleep to avoid actual delays
        with mock.patch('time.sleep'):
            # Make burst_size requests (should use tokens without delay)
            for _ in range(burst_size):
                rate_limiter.wait()
            
            # Check that tokens are depleted
            self.assertEqual(rate_limiter.tokens, 0)
            
            # Make one more request (should require delay)
            with mock.patch('time.sleep') as mock_sleep:
                rate_limiter.wait()
                mock_sleep.assert_called_once()
    
    def test_token_refill(self):
        """Test token refill based on time elapsed."""
        # Set up rate limiter with no tokens
        self.rate_limiter.tokens = 0
        self.rate_limiter.last_refill_time = time.time() - 1.0  # 1 second ago
        
        # Mock _get_current_delay to return a fixed delay
        with mock.patch.object(self.rate_limiter, '_get_current_delay', return_value=0.5):  # 2 tokens per second
            # Mock time.sleep to avoid actual delays
            with mock.patch('time.sleep'):
                # Call wait
                self.rate_limiter.wait()
                
                # Check that tokens were refilled (1 second * 2 tokens/second = 2 tokens)
                # But 1 token was used, so 1 should remain
                self.assertEqual(self.rate_limiter.tokens, 1)
    
    def test_stats_tracking(self):
        """Test statistics tracking."""
        # Reset stats
        self.rate_limiter.stats = {
            "total_requests": 0,
            "weekday_requests": 0,
            "weekend_requests": 0,
            "total_wait_time": 0.0
        }
        
        # Mock _is_weekend to alternate between weekday and weekend
        side_effects = [False, True]
        with mock.patch.object(self.rate_limiter, '_is_weekend', side_effect=side_effects):
            # Mock time.sleep to avoid actual delays
            with mock.patch('time.sleep'):
                # Make requests
                self.rate_limiter.wait()  # Weekday
                self.rate_limiter.wait()  # Weekend
                
                # Check stats
                stats = self.rate_limiter.get_stats()
                self.assertEqual(stats["total_requests"], 2)
                self.assertEqual(stats["weekday_requests"], 1)
                self.assertEqual(stats["weekend_requests"], 1)


if __name__ == "__main__":
    unittest.main()