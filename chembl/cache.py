"""
Caching system for ChEMBL API responses.

This module provides a two-level caching system (memory and disk) for ChEMBL API responses
to reduce API calls and improve performance.
"""

import os
import json
import time
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List, Union

logger = logging.getLogger(__name__)

class ChEMBLCache:
    """
    Two-level cache for ChEMBL API responses.
    
    Features:
    - In-memory LRU cache for fast access to frequently used data
    - Disk-based persistent cache for long-term storage
    - TTL (time-to-live) for cache entries
    - Error handling with separate TTL for error responses
    - Resilient to disk I/O failures
    """
    
    def __init__(
        self,
        cache_dir: str = "cache/chembl",
        memory_size: int = 1000,
        ttl: int = 86400 * 30,  # 30 days
        error_ttl: int = 3600,  # 1 hour for error responses
        max_retries: int = 3    # Maximum retries for disk operations
    ):
        """
        Initialize the cache.
        
        Args:
            cache_dir: Directory to store cache files
            memory_size: Maximum number of items to keep in memory cache
            ttl: Time-to-live for cache entries in seconds
            error_ttl: Time-to-live for error responses in seconds (shorter)
            max_retries: Maximum retries for disk operations
        """
        self.cache_dir = Path(cache_dir)
        self.memory_size = memory_size
        self.ttl = ttl
        self.error_ttl = error_ttl
        self.max_retries = max_retries
        self.memory_cache: Dict[str, Dict[str, Any]] = {}
        self.access_times: Dict[str, float] = {}
        self.stats = {
            "memory_hits": 0,
            "disk_hits": 0,
            "misses": 0,
            "sets": 0,
            "errors": 0,
            "disk_failures": 0
        }
        
        # Create cache directory if it doesn't exist
        os.makedirs(self.cache_dir, exist_ok=True)
        
        logger.info(f"ChEMBL cache initialized at {self.cache_dir}")
    
    def _get_cache_path(self, key: str) -> Path:
        """Get the file path for a cache key."""
        # Use first 2 chars of key as subdirectory to avoid too many files in one directory
        subdir = key[:2] if len(key) >= 2 else "xx"
        subdir_path = self.cache_dir / subdir
        os.makedirs(subdir_path, exist_ok=True)
        return subdir_path / f"{key}.json"
    
    def get(self, key: str) -> Optional[Dict[str, Any]]:
        """
        Get a value from the cache.
        
        Args:
            key: Cache key
            
        Returns:
            Cached value or None if not found or expired
        """
        # Check memory cache first
        if key in self.memory_cache:
            entry = self.memory_cache[key]
            
            # Determine TTL based on whether this is an error response
            current_ttl = self.error_ttl if self._is_error_response(entry.get("data", {})) else self.ttl
            
            # Check if entry is expired
            if "timestamp" in entry and time.time() - entry["timestamp"] <= current_ttl:
                # Update access time
                self.access_times[key] = time.time()
                self.stats["memory_hits"] += 1
                return entry["data"]
            else:
                # Remove expired entry
                del self.memory_cache[key]
                if key in self.access_times:
                    del self.access_times[key]
        
        # Check disk cache with retries
        for attempt in range(self.max_retries):
            try:
                cache_path = self._get_cache_path(key)
                if cache_path.exists():
                    with open(cache_path, "r") as f:
                        entry = json.load(f)
                    
                    # Determine TTL based on whether this is an error response
                    current_ttl = self.error_ttl if self._is_error_response(entry.get("data", {})) else self.ttl
                    
                    # Check if entry is expired
                    if "timestamp" in entry and time.time() - entry["timestamp"] <= current_ttl:
                        # Add to memory cache
                        self._add_to_memory_cache(key, entry)
                        self.stats["disk_hits"] += 1
                        return entry["data"]
                    else:
                        # Remove expired entry
                        try:
                            cache_path.unlink(missing_ok=True)
                        except Exception:
                            # Ignore errors when removing expired entries
                            pass
                # If we get here without exceptions, break the retry loop
                break
            except Exception as e:
                if attempt < self.max_retries - 1:
                    logger.warning(f"Error reading cache file {key} (attempt {attempt+1}/{self.max_retries}): {str(e)}")
                    time.sleep(0.1 * (2 ** attempt))  # Exponential backoff
                else:
                    logger.error(f"Failed to read cache file after {self.max_retries} attempts: {str(e)}")
                    self.stats["disk_failures"] += 1
        
        self.stats["misses"] += 1
        return None
    
    def set(self, key: str, value: Dict[str, Any]) -> None:
        """
        Set a value in the cache.
        
        Args:
            key: Cache key
            value: Value to cache
        """
        # Track error responses
        is_error = self._is_error_response(value)
        if is_error:
            self.stats["errors"] += 1
        
        entry = {
            "timestamp": time.time(),
            "data": value,
            "is_error": is_error
        }
        
        # Add to memory cache
        self._add_to_memory_cache(key, entry)
        
        # Add to disk cache with retries
        for attempt in range(self.max_retries):
            try:
                cache_path = self._get_cache_path(key)
                with open(cache_path, "w") as f:
                    json.dump(entry, f)
                # If successful, break the retry loop
                break
            except Exception as e:
                if attempt < self.max_retries - 1:
                    logger.warning(f"Error writing cache file {key} (attempt {attempt+1}/{self.max_retries}): {str(e)}")
                    time.sleep(0.1 * (2 ** attempt))  # Exponential backoff
                else:
                    logger.error(f"Failed to write cache file after {self.max_retries} attempts: {str(e)}")
                    self.stats["disk_failures"] += 1
        
        self.stats["sets"] += 1
    
    def _add_to_memory_cache(self, key: str, entry: Dict[str, Any]) -> None:
        """Add an entry to the memory cache, evicting LRU items if necessary."""
        # Evict least recently used items if cache is full
        if len(self.memory_cache) >= self.memory_size and key not in self.memory_cache:
            # Find least recently used key
            if self.access_times:
                lru_key = min(self.access_times, key=self.access_times.get)
                del self.memory_cache[lru_key]
                del self.access_times[lru_key]
        
        # Add to memory cache
        self.memory_cache[key] = entry
        self.access_times[key] = time.time()
    
    def clear(self) -> None:
        """Clear the cache (both memory and disk)."""
        # Clear memory cache
        self.memory_cache.clear()
        self.access_times.clear()
        
        # Clear disk cache
        try:
            for subdir in self.cache_dir.iterdir():
                if subdir.is_dir():
                    for cache_file in subdir.iterdir():
                        if cache_file.is_file() and cache_file.suffix == ".json":
                            cache_file.unlink(missing_ok=True)
        except Exception as e:
            logger.warning(f"Error clearing disk cache: {str(e)}")
        
        logger.info("Cache cleared")
    
    def get_stats(self) -> Dict[str, int]:
        """
        Get cache statistics.
        
        Returns:
            Dictionary with cache statistics
        """
        return {
            "memory_hits": self.stats["memory_hits"],
            "disk_hits": self.stats["disk_hits"],
            "misses": self.stats["misses"],
            "sets": self.stats["sets"],
            "errors": self.stats["errors"],
            "disk_failures": self.stats["disk_failures"],
            "memory_size": len(self.memory_cache),
            "memory_capacity": self.memory_size,
            "hit_rate": self._calculate_hit_rate()
        }
    
    def _calculate_hit_rate(self) -> float:
        """Calculate the cache hit rate as a percentage."""
        total_requests = self.stats["memory_hits"] + self.stats["disk_hits"] + self.stats["misses"]
        if total_requests == 0:
            return 0.0
        return 100.0 * (self.stats["memory_hits"] + self.stats["disk_hits"]) / total_requests
    
    def _is_error_response(self, data: Dict[str, Any]) -> bool:
        """
        Check if a response contains an error.
        
        Args:
            data: Response data
            
        Returns:
            True if the response contains an error, False otherwise
        """
        # Check for common error indicators in the response
        if not isinstance(data, dict):
            return False
        
        # Check for explicit error field
        if "Error" in data:
            return True
        
        # Check for HTTP error status
        if "status" in data and isinstance(data["status"], int) and data["status"] >= 400:
            return True
            
        # Check for error message
        if "error" in data or "errors" in data:
            return True
            
        return False