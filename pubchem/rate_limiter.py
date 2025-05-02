"""
Adaptive rate limiter for PubChem API.

This module provides rate limiting functionality that adapts based on
the day of the week (weekday vs weekend) to respect PubChem API limits.
"""

import time
import logging
import datetime
from typing import Optional, Dict, Any
import threading

logger = logging.getLogger(__name__)

class AdaptiveRateLimiter:
    """
    Adaptive rate limiter for PubChem API.
    
    Features:
    - Different rate limits for weekdays and weekends
    - Thread-safe implementation
    - Configurable rate limits
    """
    
    def __init__(
        self,
        weekday_requests_per_second: float = 2.0,  # 2 requests per second on weekdays
        weekend_requests_per_second: float = 5.0,  # 5 requests per second on weekends
        burst_size: int = 10,
        weekday_hours: tuple = (0, 24)  # Default: all day is considered "weekday hours"
    ):
        """
        Initialize the rate limiter.
        
        Args:
            weekday_requests_per_second: Maximum requests per second on weekdays
            weekend_requests_per_second: Maximum requests per second on weekends
            burst_size: Maximum number of requests that can be made in a burst
            weekday_hours: Tuple of (start_hour, end_hour) for weekday rate limiting
        """
        self.weekday_delay = 1.0 / weekday_requests_per_second if weekday_requests_per_second > 0 else 0
        self.weekend_delay = 1.0 / weekend_requests_per_second if weekend_requests_per_second > 0 else 0
        self.burst_size = burst_size
        self.weekday_hours = weekday_hours
        
        self.tokens = burst_size
        self.last_refill_time = time.time()
        self.lock = threading.RLock()
        
        # Stats
        self.stats = {
            "total_requests": 0,
            "weekday_requests": 0,
            "weekend_requests": 0,
            "total_wait_time": 0.0
        }
    
    def _is_weekend(self) -> bool:
        """
        Check if current time is considered a weekend.
        
        Returns:
            True if current time is a weekend, False otherwise
        """
        now = datetime.datetime.now()
        # 5 = Saturday, 6 = Sunday
        return now.weekday() >= 5
    
    def _is_weekday_hours(self) -> bool:
        """
        Check if current time is within weekday hours.
        
        Returns:
            True if current time is within weekday hours, False otherwise
        """
        now = datetime.datetime.now()
        current_hour = now.hour
        start_hour, end_hour = self.weekday_hours
        
        return start_hour <= current_hour < end_hour
    
    def _get_current_delay(self) -> float:
        """
        Get the current delay based on day of week and time.
        
        Returns:
            Delay in seconds between requests
        """
        if self._is_weekend():
            return self.weekend_delay
        else:
            if self._is_weekday_hours():
                return self.weekday_delay
            else:
                # Outside of weekday hours, use weekend rate
                return self.weekend_delay
    
    def wait(self) -> None:
        """
        Wait for the appropriate amount of time to respect rate limits.
        
        This method implements a token bucket algorithm for rate limiting.
        """
        with self.lock:
            # Update stats
            self.stats["total_requests"] += 1
            if self._is_weekend():
                self.stats["weekend_requests"] += 1
            else:
                self.stats["weekday_requests"] += 1
            
            # Calculate tokens to add based on time elapsed
            now = time.time()
            time_elapsed = now - self.last_refill_time
            self.last_refill_time = now
            
            # Add tokens based on time elapsed and current rate
            current_delay = self._get_current_delay()
            tokens_to_add = time_elapsed / current_delay if current_delay > 0 else self.burst_size
            self.tokens = min(self.burst_size, self.tokens + tokens_to_add)
            
            # If we have at least one token, use it and return immediately
            if self.tokens >= 1:
                self.tokens -= 1
                return
            
            # Otherwise, calculate wait time
            wait_time = (1 - self.tokens) * current_delay
            self.tokens = 0
            
            # Update stats
            self.stats["total_wait_time"] += wait_time
            
            # Log if wait time is significant
            if wait_time > 1.0:
                logger.info(f"Rate limiting: waiting for {wait_time:.2f} seconds")
            
            # Wait
            time.sleep(wait_time)
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get rate limiter statistics.
        
        Returns:
            Dictionary with rate limiter statistics
        """
        with self.lock:
            return self.stats.copy()