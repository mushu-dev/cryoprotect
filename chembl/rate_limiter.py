"""
Adaptive rate limiter for ChEMBL API.

This module provides a rate limiter that adapts to the ChEMBL API's rate limits,
slowing down during weekdays and speeding up during weekends.
"""

import time
import logging
from datetime import datetime
from typing import Dict, Any

logger = logging.getLogger(__name__)

class AdaptiveRateLimiter:
    """
    Adaptive rate limiter for ChEMBL API.
    
    Features:
    - Different rate limits for weekdays and weekends
    - Adaptive delay based on API response times
    - Automatic backoff on rate limit errors
    """
    
    def __init__(
        self,
        weekday_requests_per_second: float = 5.0,
        weekend_requests_per_second: float = 5.0,
        monday_requests_per_second: float = 3.0,  # Stricter limit for Mondays
        min_delay: float = 0.1,
        max_delay: float = 5.0,
        backoff_factor: float = 1.5,
        memory_threshold_percent: float = 80.0
    ):
        """
        Initialize the rate limiter.
        
        Args:
            weekday_requests_per_second: Maximum requests per second on weekdays (Tue-Fri)
            weekend_requests_per_second: Maximum requests per second on weekends
            monday_requests_per_second: Maximum requests per second on Mondays (stricter)
            min_delay: Minimum delay between requests in seconds
            max_delay: Maximum delay between requests in seconds
            backoff_factor: Factor to increase delay by on rate limit errors
            memory_threshold_percent: Memory usage threshold (as percentage) to trigger additional delay
        """
        self.weekday_delay = 1.0 / weekday_requests_per_second
        self.weekend_delay = 1.0 / weekend_requests_per_second
        self.monday_delay = 1.0 / monday_requests_per_second
        self.min_delay = min_delay
        self.max_delay = max_delay
        self.backoff_factor = backoff_factor
        self.memory_threshold_percent = memory_threshold_percent
        
        self.current_delay = self.weekday_delay
        self.last_request_time = 0.0
        self.request_count = 0
        self.rate_limit_errors = 0
        self.memory_slowdowns = 0
        
        logger.info(f"AdaptiveRateLimiter initialized (weekday: {self.weekday_delay:.2f}s, weekend: {self.weekend_delay:.2f}s, memory threshold: {self.memory_threshold_percent}%)")
    
    def wait(self) -> None:
        """
        Wait for the appropriate amount of time before making a request.
        Includes memory usage monitoring to slow down if memory usage is high.
        """
        # Update current delay based on day of week
        now = datetime.now()
        weekday = now.weekday()
        
        # Apply different rate limits based on day of week
        if weekday == 0:  # Monday
            self.current_delay = self.monday_delay
            logger.debug("Using Monday rate limit")
        elif weekday >= 5:  # Weekend (5 = Saturday, 6 = Sunday)
            self.current_delay = self.weekend_delay
            logger.debug("Using weekend rate limit")
        else:  # Regular weekday (Tuesday-Friday)
            self.current_delay = self.weekday_delay
            logger.debug("Using weekday rate limit")
        
        # Apply backoff if there have been rate limit errors
        if self.rate_limit_errors > 0:
            self.current_delay *= (self.backoff_factor ** self.rate_limit_errors)
            self.current_delay = min(self.current_delay, self.max_delay)
        
        # Check memory usage and slow down if necessary
        memory_usage_percent = self._get_memory_usage_percent()
        if memory_usage_percent > self.memory_threshold_percent:
            memory_factor = 1.0 + ((memory_usage_percent - self.memory_threshold_percent) / 20.0)
            self.current_delay *= memory_factor
            self.memory_slowdowns += 1
            logger.warning(f"Memory usage at {memory_usage_percent:.1f}% (above {self.memory_threshold_percent}% threshold). "
                          f"Increasing delay by factor of {memory_factor:.2f} to {self.current_delay:.2f}s")
        
        # Ensure minimum delay
        self.current_delay = max(self.current_delay, self.min_delay)
        # Ensure maximum delay
        self.current_delay = min(self.current_delay, self.max_delay)
        
        # Calculate time to wait
        elapsed = time.time() - self.last_request_time
        wait_time = max(0, self.current_delay - elapsed)
        
        if wait_time > 0:
            time.sleep(wait_time)
        
        self.last_request_time = time.time()
        self.request_count += 1
    
    def report_rate_limit_error(self) -> None:
        """
        Report a rate limit error to increase backoff.
        """
        self.rate_limit_errors += 1
        logger.warning(f"Rate limit error reported. Increasing delay to {self.current_delay * self.backoff_factor:.2f}s")
    
    def report_success(self) -> None:
        """
        Report a successful request to gradually reduce backoff.
        """
        if self.rate_limit_errors > 0:
            self.rate_limit_errors -= 0.5
            logger.info(f"Successful request. Reducing rate limit error count to {self.rate_limit_errors}")
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get rate limiter statistics.
        
        Returns:
            Dictionary with rate limiter statistics
        """
        memory_usage = self._get_memory_usage_percent()
        now = datetime.now()
        weekday = now.weekday()
        day_type = "monday" if weekday == 0 else "weekend" if weekday >= 5 else "weekday"
        
        return {
            "request_count": self.request_count,
            "current_delay": self.current_delay,
            "rate_limit_errors": self.rate_limit_errors,
            "day_type": day_type,
            "day_of_week": weekday,
            "memory_usage_percent": memory_usage,
            "memory_slowdowns": self.memory_slowdowns,
            "requests_per_second": {
                "monday": 1.0 / self.monday_delay,
                "weekday": 1.0 / self.weekday_delay,
                "weekend": 1.0 / self.weekend_delay
            }
        }
        
    def _get_memory_usage_percent(self) -> float:
        """
        Get current memory usage as a percentage.
        
        Returns:
            Memory usage percentage (0-100)
        """
        try:
            import psutil
            process = psutil.Process()
            memory_info = process.memory_info()
            memory_percent = process.memory_percent()
            logger.debug(f"Memory usage: {memory_percent:.1f}% ({memory_info.rss / (1024 * 1024):.1f} MB)")
            return memory_percent
        except ImportError:
            logger.warning("psutil not installed, cannot monitor memory usage")
            return 0.0
        except Exception as e:
            logger.warning(f"Error getting memory usage: {e}")
            return 0.0