"""
Simple rate limiter for PubChem API.

This module provides a basic rate limiting functionality for the PubChem API
to avoid exceeding the API rate limits.
"""

import time
import logging
from typing import Dict, Any

logger = logging.getLogger(__name__)

class RateLimiter:
    """
    Simple rate limiter for PubChem API.
    
    This is a simplified version of AdaptiveRateLimiter for use with
    the enhanced PubChem importer.
    """
    
    def __init__(self, requests_per_minute: float = 5.0):
        """
        Initialize the rate limiter.
        
        Args:
            requests_per_minute: Maximum requests per minute
        """
        self.delay = 60.0 / requests_per_minute if requests_per_minute > 0 else 0
        self.last_request_time = 0
        
        # Stats
        self.stats = {
            "total_requests": 0,
            "total_wait_time": 0.0
        }
    
    def wait(self) -> None:
        """
        Wait for the appropriate amount of time to respect rate limits.
        """
        # Update stats
        self.stats["total_requests"] += 1
        
        # Calculate time since last request
        now = time.time()
        time_since_last_request = now - self.last_request_time
        
        # If we need to wait
        if time_since_last_request < self.delay:
            wait_time = self.delay - time_since_last_request
            
            # Update stats
            self.stats["total_wait_time"] += wait_time
            
            # Log if wait time is significant
            if wait_time > 1.0:
                logger.debug(f"Rate limiting: waiting for {wait_time:.2f} seconds")
            
            # Wait
            time.sleep(wait_time)
            
            # Update last request time
            self.last_request_time = time.time()
        else:
            # No need to wait
            self.last_request_time = now
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get rate limiter statistics.
        
        Returns:
            Dictionary with rate limiter statistics
        """
        return self.stats.copy()