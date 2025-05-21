"""
Tests for the caching system.

This module tests the various cache backends and the cache manager,
ensuring proper functionality of memory and disk caching.
"""

import os
import time
import json
import pickle
import tempfile
import unittest
import asyncio
import logging
import shutil
from typing import Dict, Any, Optional

# Import module to test
from unified_importer.core.caching import (
    CacheEntry,
    MemoryCache,
    DiskCache,
    CacheManager,
    AsyncCacheManager,
    create_cache_key
)


class TestCacheEntry(unittest.TestCase):
    """Tests for the CacheEntry class."""
    
    def test_init(self):
        """Test initialization of cache entry."""
        entry = CacheEntry(key="test", value={"data": 123})
        self.assertEqual(entry.key, "test")
        self.assertEqual(entry.value, {"data": 123})
        self.assertIsNone(entry.expires_at)
        self.assertEqual(entry.size_bytes, 0)
        self.assertEqual(entry.hit_count, 0)
        
        # Check that created_at and last_accessed are set
        self.assertIsInstance(entry.created_at, float)
        self.assertIsInstance(entry.last_accessed, float)
    
    def test_is_expired(self):
        """Test expiration checking."""
        # Entry with no expiration
        entry1 = CacheEntry(key="test1", value="value1")
        self.assertFalse(entry1.is_expired())
        
        # Entry with future expiration
        entry2 = CacheEntry(key="test2", value="value2", expires_at=time.time() + 60)
        self.assertFalse(entry2.is_expired())
        
        # Entry with past expiration
        entry3 = CacheEntry(key="test3", value="value3", expires_at=time.time() - 60)
        self.assertTrue(entry3.is_expired())
    
    def test_time_to_expiry(self):
        """Test time to expiry calculation."""
        # Entry with no expiration
        entry1 = CacheEntry(key="test1", value="value1")
        self.assertIsNone(entry1.time_to_expiry())
        
        # Entry with future expiration
        entry2 = CacheEntry(key="test2", value="value2", expires_at=time.time() + 60)
        self.assertGreater(entry2.time_to_expiry(), 0)
        self.assertLessEqual(entry2.time_to_expiry(), 60)
        
        # Entry with past expiration
        entry3 = CacheEntry(key="test3", value="value3", expires_at=time.time() - 60)
        self.assertEqual(entry3.time_to_expiry(), 0)
    
    def test_age(self):
        """Test age calculation."""
        entry = CacheEntry(key="test", value="value")
        time.sleep(0.01)  # Ensure some time passes
        self.assertGreater(entry.age(), 0)
    
    def test_record_hit(self):
        """Test hit recording."""
        entry = CacheEntry(key="test", value="value")
        self.assertEqual(entry.hit_count, 0)
        
        # Record hits
        entry.record_hit()
        self.assertEqual(entry.hit_count, 1)
        
        entry.record_hit()
        self.assertEqual(entry.hit_count, 2)
    
    def test_to_dict(self):
        """Test conversion to dictionary."""
        entry = CacheEntry(key="test", value="value")
        data_dict = entry.to_dict()
        
        self.assertEqual(data_dict["key"], "test")
        self.assertEqual(data_dict["value"], "value")
        self.assertIsInstance(data_dict["created_at"], float)
    
    def test_from_dict(self):
        """Test creation from dictionary."""
        data_dict = {
            "key": "test",
            "value": "value",
            "created_at": time.time(),
            "expires_at": None,
            "size_bytes": 100,
            "source": "test_source",
            "hit_count": 5,
            "last_accessed": time.time()
        }
        
        entry = CacheEntry.from_dict(data_dict)
        self.assertEqual(entry.key, "test")
        self.assertEqual(entry.value, "value")
        self.assertEqual(entry.size_bytes, 100)
        self.assertEqual(entry.source, "test_source")
        self.assertEqual(entry.hit_count, 5)


class TestMemoryCache(unittest.TestCase):
    """Tests for the MemoryCache class."""
    
    def setUp(self):
        """Set up test case."""
        self.cache = MemoryCache(
            max_size_bytes=1024 * 1024,  # 1MB
            max_items=10,
            default_ttl=60.0,
            cleanup_interval=1.0
        )
    
    def tearDown(self):
        """Clean up after test case."""
        self.cache.close()
    
    def test_put_get(self):
        """Test putting and getting items."""
        # Put an item
        entry = CacheEntry(key="test", value="value")
        success = self.cache.put(entry)
        self.assertTrue(success)
        
        # Get the item
        result = self.cache.get("test")
        self.assertIsNotNone(result)
        self.assertEqual(result.key, "test")
        self.assertEqual(result.value, "value")
        self.assertEqual(result.hit_count, 1)
    
    def test_get_nonexistent(self):
        """Test getting a nonexistent item."""
        result = self.cache.get("nonexistent")
        self.assertIsNone(result)
    
    def test_delete(self):
        """Test deleting an item."""
        # Put an item
        entry = CacheEntry(key="test", value="value")
        self.cache.put(entry)
        
        # Delete the item
        success = self.cache.delete("test")
        self.assertTrue(success)
        
        # Try to get the deleted item
        result = self.cache.get("test")
        self.assertIsNone(result)
        
        # Try to delete a nonexistent item
        success = self.cache.delete("nonexistent")
        self.assertFalse(success)
    
    def test_clear(self):
        """Test clearing the cache."""
        # Put some items
        for i in range(5):
            entry = CacheEntry(key=f"test{i}", value=f"value{i}")
            self.cache.put(entry)
        
        # Clear the cache
        success = self.cache.clear()
        self.assertTrue(success)
        
        # Try to get the items
        for i in range(5):
            result = self.cache.get(f"test{i}")
            self.assertIsNone(result)
    
    def test_expiration(self):
        """Test item expiration."""
        # Put an item with a short TTL
        entry = CacheEntry(key="test", value="value", expires_at=time.time() + 0.1)
        self.cache.put(entry)
        
        # Get the item before expiration
        result = self.cache.get("test")
        self.assertIsNotNone(result)
        
        # Wait for the item to expire
        time.sleep(0.2)
        
        # Try to get the expired item
        result = self.cache.get("test")
        self.assertIsNone(result)
    
    def test_cleanup(self):
        """Test automatic cleanup of expired items."""
        # Put items with short TTLs
        for i in range(5):
            entry = CacheEntry(key=f"test{i}", value=f"value{i}", expires_at=time.time() + 0.1)
            self.cache.put(entry)
        
        # Wait for the items to expire and cleanup to run
        time.sleep(1.5)
        
        # Check stats
        stats = self.cache.get_stats()
        self.assertGreaterEqual(stats["expirations"], 5)
        
        # Try to get the expired items
        for i in range(5):
            result = self.cache.get(f"test{i}")
            self.assertIsNone(result)
    
    def test_eviction(self):
        """Test LRU eviction when cache is full."""
        # Fill the cache to capacity
        for i in range(10):  # max_items = 10
            entry = CacheEntry(key=f"test{i}", value=f"value{i}", size_bytes=100)
            self.cache.put(entry)
        
        # All items should be in the cache
        for i in range(10):
            result = self.cache.get(f"test{i}")
            self.assertIsNotNone(result)
        
        # Put a new item, which should evict the oldest
        entry = CacheEntry(key="new", value="new_value", size_bytes=100)
        self.cache.put(entry)
        
        # The new item should be in the cache
        result = self.cache.get("new")
        self.assertIsNotNone(result)
        
        # At least one old item should be evicted
        evicted_count = 0
        for i in range(10):
            if self.cache.get(f"test{i}") is None:
                evicted_count += 1
        
        self.assertGreater(evicted_count, 0)
    
    def test_size_limit(self):
        """Test size-based eviction."""
        # Put a large item
        large_data = "x" * 900 * 1024  # 900KB
        entry = CacheEntry(key="large", value=large_data, size_bytes=900 * 1024)
        self.cache.put(entry)
        
        # Put more items until we exceed the size limit
        for i in range(5):
            entry = CacheEntry(key=f"test{i}", value=f"value{i}", size_bytes=50 * 1024)
            self.cache.put(entry)
        
        # Check that some items were evicted
        self.assertLess(len(self.cache.cache), 6)
    
    def test_get_stats(self):
        """Test statistics collection."""
        # Put some items
        for i in range(5):
            entry = CacheEntry(key=f"test{i}", value=f"value{i}")
            self.cache.put(entry)
        
        # Get some items multiple times
        for _ in range(3):
            self.cache.get("test0")
        
        self.cache.get("test1")
        
        # Get a nonexistent item
        self.cache.get("nonexistent")
        
        # Check stats
        stats = self.cache.get_stats()
        self.assertEqual(stats["insertions"], 5)
        self.assertEqual(stats["hits"], 4)
        self.assertEqual(stats["misses"], 1)
        self.assertEqual(stats["items"], 5)


class TestDiskCache(unittest.TestCase):
    """Tests for the DiskCache class."""
    
    def setUp(self):
        """Set up test case."""
        # Create a temporary directory for the cache
        self.temp_dir = tempfile.mkdtemp()
        
        self.cache = DiskCache(
            cache_dir=self.temp_dir,
            max_size_bytes=1024 * 1024,  # 1MB
            max_items=10,
            default_ttl=60.0,
            cleanup_interval=1.0
        )
    
    def tearDown(self):
        """Clean up after test case."""
        self.cache.close()
        shutil.rmtree(self.temp_dir)
    
    def test_put_get(self):
        """Test putting and getting items."""
        # Put an item
        entry = CacheEntry(key="test", value="value")
        success = self.cache.put(entry)
        self.assertTrue(success)
        
        # Get the item
        result = self.cache.get("test")
        self.assertIsNotNone(result)
        self.assertEqual(result.key, "test")
        self.assertEqual(result.value, "value")
        self.assertEqual(result.hit_count, 1)
    
    def test_get_nonexistent(self):
        """Test getting a nonexistent item."""
        result = self.cache.get("nonexistent")
        self.assertIsNone(result)
    
    def test_delete(self):
        """Test deleting an item."""
        # Put an item
        entry = CacheEntry(key="test", value="value")
        self.cache.put(entry)
        
        # Delete the item
        success = self.cache.delete("test")
        self.assertTrue(success)
        
        # Try to get the deleted item
        result = self.cache.get("test")
        self.assertIsNone(result)
    
    def test_clear(self):
        """Test clearing the cache."""
        # Put some items
        for i in range(5):
            entry = CacheEntry(key=f"test{i}", value=f"value{i}")
            self.cache.put(entry)
        
        # Clear the cache
        success = self.cache.clear()
        self.assertTrue(success)
        
        # Try to get the items
        for i in range(5):
            result = self.cache.get(f"test{i}")
            self.assertIsNone(result)
    
    def test_expiration(self):
        """Test item expiration."""
        # Put an item with a short TTL
        entry = CacheEntry(key="test", value="value", expires_at=time.time() + 0.1)
        self.cache.put(entry)
        
        # Get the item before expiration
        result = self.cache.get("test")
        self.assertIsNotNone(result)
        
        # Wait for the item to expire
        time.sleep(0.2)
        
        # Try to get the expired item
        result = self.cache.get("test")
        self.assertIsNone(result)
    
    def test_complex_data(self):
        """Test caching complex data structures."""
        # Create a complex nested structure
        complex_data = {
            "string": "value",
            "number": 123,
            "list": [1, 2, 3],
            "dict": {"a": 1, "b": 2},
            "nested": {
                "list": [{"x": 1}, {"y": 2}],
                "tuple": (1, 2, 3)
            }
        }
        
        # Put the complex data
        entry = CacheEntry(key="complex", value=complex_data)
        self.cache.put(entry)
        
        # Get the complex data
        result = self.cache.get("complex")
        self.assertIsNotNone(result)
        
        # Check that the complex data is intact
        self.assertEqual(result.value["string"], "value")
        self.assertEqual(result.value["number"], 123)
        self.assertEqual(result.value["list"], [1, 2, 3])
        self.assertEqual(result.value["dict"], {"a": 1, "b": 2})
        self.assertEqual(result.value["nested"]["list"][0]["x"], 1)
    
    def test_get_stats(self):
        """Test statistics collection."""
        # Put some items
        for i in range(5):
            entry = CacheEntry(key=f"test{i}", value=f"value{i}")
            self.cache.put(entry)
        
        # Get some items multiple times
        for _ in range(3):
            self.cache.get("test0")
        
        self.cache.get("test1")
        
        # Get a nonexistent item
        self.cache.get("nonexistent")
        
        # Check stats
        stats = self.cache.get_stats()
        self.assertEqual(stats["insertions"], 5)
        self.assertEqual(stats["hits"], 4)
        self.assertEqual(stats["misses"], 1)
        self.assertEqual(stats["items"], 5)
    
    def test_persistence(self):
        """Test cache persistence between instances."""
        # Put an item
        entry = CacheEntry(key="test", value="value")
        self.cache.put(entry)
        
        # Close the cache
        self.cache.close()
        
        # Create a new cache with the same directory
        new_cache = DiskCache(
            cache_dir=self.temp_dir,
            max_size_bytes=1024 * 1024,
            max_items=10,
            default_ttl=60.0,
            cleanup_interval=1.0
        )
        
        try:
            # Try to get the item from the new cache
            result = new_cache.get("test")
            self.assertIsNotNone(result)
            self.assertEqual(result.value, "value")
        finally:
            new_cache.close()


class TestCacheManager(unittest.TestCase):
    """Tests for the CacheManager class."""
    
    def setUp(self):
        """Set up test case."""
        # Create a temporary directory for the disk cache
        self.temp_dir = tempfile.mkdtemp()
        
        # Create cache backends
        self.memory_cache = MemoryCache(
            max_size_bytes=1024 * 1024,
            max_items=10,
            default_ttl=60.0,
            cleanup_interval=1.0
        )
        
        self.disk_cache = DiskCache(
            cache_dir=self.temp_dir,
            max_size_bytes=1024 * 1024,
            max_items=10,
            default_ttl=60.0,
            cleanup_interval=1.0
        )
        
        # Create cache manager
        self.manager = CacheManager(
            memory_cache=self.memory_cache,
            disk_cache=self.disk_cache
        )
    
    def tearDown(self):
        """Clean up after test case."""
        self.manager.close()
        shutil.rmtree(self.temp_dir)
    
    def test_get_backend(self):
        """Test getting the appropriate backend for a policy."""
        # Default policy should use memory cache
        backend = self.manager.get_backend("default")
        self.assertIs(backend, self.memory_cache)
        
        # PubChem policy should use disk cache
        backend = self.manager.get_backend("pubchem")
        self.assertIs(backend, self.disk_cache)
        
        # ChEMBL policy should use disk cache
        backend = self.manager.get_backend("chembl")
        self.assertIs(backend, self.disk_cache)
        
        # Transient policy should use memory cache
        backend = self.manager.get_backend("transient")
        self.assertIs(backend, self.memory_cache)
        
        # Unknown policy should use default (memory cache)
        backend = self.manager.get_backend("unknown")
        self.assertIs(backend, self.memory_cache)
    
    def test_get_ttl(self):
        """Test getting the TTL for a policy."""
        # Default policy TTL
        ttl = self.manager.get_ttl("default")
        self.assertEqual(ttl, 3600.0)  # 1 hour
        
        # PubChem policy TTL
        ttl = self.manager.get_ttl("pubchem")
        self.assertEqual(ttl, 86400.0 * 7)  # 7 days
        
        # Unknown policy should use default TTL
        ttl = self.manager.get_ttl("unknown")
        self.assertEqual(ttl, 3600.0)
    
    def test_put_get(self):
        """Test putting and getting items with different policies."""
        # Put items with different policies
        self.manager.put("key1", "value1", "default")
        self.manager.put("key2", "value2", "pubchem")
        
        # Get items
        value1 = self.manager.get("key1", "default")
        value2 = self.manager.get("key2", "pubchem")
        
        self.assertEqual(value1, "value1")
        self.assertEqual(value2, "value2")
        
        # Items should be in the correct backends
        self.assertIsNotNone(self.memory_cache.get("key1"))
        self.assertIsNone(self.memory_cache.get("key2"))
        
        self.assertIsNotNone(self.disk_cache.get("key2"))
        self.assertIsNone(self.disk_cache.get("key1"))
    
    def test_custom_ttl(self):
        """Test setting a custom TTL."""
        # Put an item with a custom TTL
        self.manager.put("key", "value", "default", ttl=0.1)
        
        # Get the item before expiration
        value = self.manager.get("key", "default")
        self.assertEqual(value, "value")
        
        # Wait for the item to expire
        time.sleep(0.2)
        
        # Try to get the expired item
        value = self.manager.get("key", "default")
        self.assertIsNone(value)
    
    def test_delete(self):
        """Test deleting items with different policies."""
        # Put items with different policies
        self.manager.put("key1", "value1", "default")
        self.manager.put("key2", "value2", "pubchem")
        
        # Delete items
        success1 = self.manager.delete("key1", "default")
        success2 = self.manager.delete("key2", "pubchem")
        
        self.assertTrue(success1)
        self.assertTrue(success2)
        
        # Try to get the deleted items
        value1 = self.manager.get("key1", "default")
        value2 = self.manager.get("key2", "pubchem")
        
        self.assertIsNone(value1)
        self.assertIsNone(value2)
    
    def test_clear(self):
        """Test clearing caches."""
        # Put items with different policies
        self.manager.put("key1", "value1", "default")
        self.manager.put("key2", "value2", "pubchem")
        
        # Clear only the memory cache
        self.manager.clear("default")
        
        # Item in memory cache should be gone
        value1 = self.manager.get("key1", "default")
        self.assertIsNone(value1)
        
        # Item in disk cache should still be there
        value2 = self.manager.get("key2", "pubchem")
        self.assertEqual(value2, "value2")
        
        # Put a new item in memory cache
        self.manager.put("key3", "value3", "default")
        
        # Clear all caches
        self.manager.clear()
        
        # All items should be gone
        value2 = self.manager.get("key2", "pubchem")
        value3 = self.manager.get("key3", "default")
        
        self.assertIsNone(value2)
        self.assertIsNone(value3)
    
    def test_get_stats(self):
        """Test getting statistics."""
        # Put some items
        self.manager.put("key1", "value1", "default")
        self.manager.put("key2", "value2", "pubchem")
        
        # Get items
        self.manager.get("key1", "default")
        self.manager.get("key2", "pubchem")
        
        # Get stats
        stats = self.manager.get_stats()
        
        # Should have stats for both backends
        self.assertIn("memory", stats)
        self.assertIn("disk", stats)
        
        # Memory cache should have 1 item
        self.assertEqual(stats["memory"]["items"], 1)
        
        # Disk cache should have 1 item
        self.assertEqual(stats["disk"]["items"], 1)


class TestAsyncCacheManager(unittest.IsolatedAsyncioTestCase):
    """Tests for the AsyncCacheManager class."""
    
    async def asyncSetUp(self):
        """Set up test case."""
        # Create a temporary directory for the disk cache
        self.temp_dir = tempfile.mkdtemp()
        
        # Create cache backends
        self.memory_cache = MemoryCache(
            max_size_bytes=1024 * 1024,
            max_items=10,
            default_ttl=60.0,
            cleanup_interval=1.0
        )
        
        self.disk_cache = DiskCache(
            cache_dir=self.temp_dir,
            max_size_bytes=1024 * 1024,
            max_items=10,
            default_ttl=60.0,
            cleanup_interval=1.0
        )
        
        # Create cache manager
        self.cache_manager = CacheManager(
            memory_cache=self.memory_cache,
            disk_cache=self.disk_cache
        )
        
        # Create async cache manager
        self.async_manager = AsyncCacheManager(
            cache_manager=self.cache_manager
        )
    
    async def asyncTearDown(self):
        """Clean up after test case."""
        self.async_manager.close()
        shutil.rmtree(self.temp_dir)
    
    async def test_put_get(self):
        """Test putting and getting items asynchronously."""
        # Put items with different policies
        await self.async_manager.put("key1", "value1", "default")
        await self.async_manager.put("key2", "value2", "pubchem")
        
        # Get items
        value1 = await self.async_manager.get("key1", "default")
        value2 = await self.async_manager.get("key2", "pubchem")
        
        self.assertEqual(value1, "value1")
        self.assertEqual(value2, "value2")
    
    async def test_delete(self):
        """Test deleting items asynchronously."""
        # Put items
        await self.async_manager.put("key1", "value1", "default")
        
        # Delete items
        success = await self.async_manager.delete("key1", "default")
        
        self.assertTrue(success)
        
        # Try to get the deleted item
        value = await self.async_manager.get("key1", "default")
        
        self.assertIsNone(value)
    
    async def test_clear(self):
        """Test clearing caches asynchronously."""
        # Put items
        await self.async_manager.put("key1", "value1", "default")
        
        # Clear all caches
        success = await self.async_manager.clear()
        
        self.assertTrue(success)
        
        # Try to get the item
        value = await self.async_manager.get("key1", "default")
        
        self.assertIsNone(value)
    
    async def test_get_stats(self):
        """Test getting statistics asynchronously."""
        # Put some items
        await self.async_manager.put("key1", "value1", "default")
        await self.async_manager.put("key2", "value2", "pubchem")
        
        # Get items
        await self.async_manager.get("key1", "default")
        await self.async_manager.get("key2", "pubchem")
        
        # Get stats
        stats = await self.async_manager.get_stats()
        
        # Should have stats for both backends
        self.assertIn("memory", stats)
        self.assertIn("disk", stats)


class TestCacheKeyCreation(unittest.TestCase):
    """Tests for the cache key creation utility."""
    
    def test_create_cache_key(self):
        """Test creating cache keys from parts."""
        # Simple key
        key = create_cache_key("test", ["a", "b", "c"])
        self.assertEqual(key, "test:a:b:c")
        
        # Key with different types
        key = create_cache_key("test", ["a", 123, True, None])
        self.assertEqual(key, "test:a:123:True:None")
        
        # Empty parts
        key = create_cache_key("test", [])
        self.assertEqual(key, "test:")
        
        # Long key that gets hashed
        long_part = "x" * 1000
        key = create_cache_key("test", [long_part], max_length=50)
        self.assertEqual(len(key), 50)
        self.assertTrue(key.startswith("test:"))


if __name__ == "__main__":
    unittest.main()