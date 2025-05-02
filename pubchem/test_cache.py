"""
Unit tests for the PubChem cache module.
"""

import os
import time
import json
import shutil
import unittest
from pathlib import Path
from unittest import mock

from .cache import PubChemCache

class TestPubChemCache(unittest.TestCase):
    """Test cases for PubChemCache."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_cache_dir = "test_cache"
        self.cache = PubChemCache(cache_dir=self.test_cache_dir, memory_size=10, ttl=1)
    
    def tearDown(self):
        """Clean up test environment."""
        if os.path.exists(self.test_cache_dir):
            shutil.rmtree(self.test_cache_dir)
    
    def test_memory_cache(self):
        """Test in-memory caching."""
        # Set a value
        self.cache.set("test_key", "test_value")
        
        # Get the value
        value = self.cache.get("test_key")
        
        # Check that the value was retrieved from memory
        self.assertEqual(value, "test_value")
        self.assertEqual(self.cache.stats["memory_hits"], 1)
        self.assertEqual(self.cache.stats["disk_hits"], 0)
        self.assertEqual(self.cache.stats["misses"], 0)
    
    def test_disk_cache(self):
        """Test disk-based caching."""
        # Set a value
        self.cache.set("test_key", "test_value")
        
        # Create a new cache instance (to clear memory cache)
        new_cache = PubChemCache(cache_dir=self.test_cache_dir, memory_size=10, ttl=1)
        
        # Get the value
        value = new_cache.get("test_key")
        
        # Check that the value was retrieved from disk
        self.assertEqual(value, "test_value")
        self.assertEqual(new_cache.stats["memory_hits"], 0)
        self.assertEqual(new_cache.stats["disk_hits"], 1)
        self.assertEqual(new_cache.stats["misses"], 0)
    
    def test_cache_expiration(self):
        """Test cache expiration."""
        # Set a value
        self.cache.set("test_key", "test_value")
        
        # Wait for TTL to expire
        time.sleep(1.1)
        
        # Get the value
        value = self.cache.get("test_key")
        
        # Check that the value was not found
        self.assertIsNone(value)
        self.assertEqual(self.cache.stats["memory_hits"], 0)
        self.assertEqual(self.cache.stats["disk_hits"], 0)
        self.assertEqual(self.cache.stats["misses"], 1)
    
    def test_cache_clear(self):
        """Test cache clearing."""
        # Set some values
        self.cache.set("test_key1", "test_value1")
        self.cache.set("test_key2", "test_value2")
        
        # Clear the cache
        self.cache.clear()
        
        # Get the values
        value1 = self.cache.get("test_key1")
        value2 = self.cache.get("test_key2")
        
        # Check that the values were not found
        self.assertIsNone(value1)
        self.assertIsNone(value2)
        self.assertEqual(self.cache.stats["memory_hits"], 0)
        self.assertEqual(self.cache.stats["disk_hits"], 0)
        self.assertEqual(self.cache.stats["misses"], 2)
    
    def test_cache_stats(self):
        """Test cache statistics."""
        # Set a value
        self.cache.set("test_key", "test_value")
        
        # Get the value multiple times
        self.cache.get("test_key")
        self.cache.get("test_key")
        self.cache.get("nonexistent_key")
        
        # Check stats
        stats = self.cache.get_stats()
        self.assertEqual(stats["memory_hits"], 2)
        self.assertEqual(stats["disk_hits"], 0)
        self.assertEqual(stats["misses"], 1)
        self.assertEqual(stats["writes"], 1)
    
    def test_disk_cache_error_handling(self):
        """Test disk cache error handling."""
        # Mock open to raise an exception
        with mock.patch("builtins.open", side_effect=IOError("Test error")):
            # Set a value (should fail silently)
            self.cache.set("test_key", "test_value")
            
            # Get the value (should return None)
            value = self.cache.get("test_key")
            
            # Check that the value was not found
            self.assertIsNone(value)
            self.assertEqual(self.cache.stats["memory_hits"], 0)
            self.assertEqual(self.cache.stats["disk_hits"], 0)
            self.assertEqual(self.cache.stats["misses"], 1)
    
    def test_memory_cache_size_limit(self):
        """Test memory cache size limit."""
        # Set more values than the memory cache size
        for i in range(20):
            self.cache.set(f"test_key{i}", f"test_value{i}")
        
        # Get a value that should have been evicted from memory
        value = self.cache.get("test_key0")
        
        # Check that the value was retrieved from disk, not memory
        self.assertEqual(value, "test_value0")
        self.assertEqual(self.cache.stats["memory_hits"], 0)
        self.assertEqual(self.cache.stats["disk_hits"], 1)
        self.assertEqual(self.cache.stats["misses"], 0)


if __name__ == "__main__":
    unittest.main()