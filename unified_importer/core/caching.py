"""
Caching system for the unified molecular importer.

This module provides caching capabilities to improve the performance of
the importer by reducing redundant API calls and database operations.
It implements various caching strategies including memory caching,
disk-based caching, and persistent caching between runs.
"""

import os
import time
import json
import pickle
import hashlib
import logging
import threading
import asyncio
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Set, Tuple, Callable, TypeVar, Generic
from dataclasses import dataclass, field
from datetime import datetime, timedelta
import sqlite3

# Type variables for generic functions
T = TypeVar('T')  # Cache key type (usually str)
V = TypeVar('V')  # Cache value type (Any)


@dataclass
class CacheEntry(Generic[T, V]):
    """
    Data structure for a cache entry with metadata.
    
    Attributes:
        key: Cache key
        value: Cached value
        created_at: Timestamp when the entry was created
        expires_at: Timestamp when the entry expires (None if no expiration)
        size_bytes: Approximate size of the entry in bytes
        source: Source of the data (e.g., "pubchem", "chembl")
        hit_count: Number of times this entry has been accessed
        last_accessed: Timestamp when the entry was last accessed
    """
    key: T
    value: V
    created_at: float = field(default_factory=time.time)
    expires_at: Optional[float] = None
    size_bytes: int = 0
    source: str = ""
    hit_count: int = 0
    last_accessed: float = field(default_factory=time.time)
    
    def is_expired(self) -> bool:
        """Check if the cache entry is expired."""
        if self.expires_at is None:
            return False
        return time.time() > self.expires_at
    
    def time_to_expiry(self) -> Optional[float]:
        """Get the time remaining until expiry in seconds."""
        if self.expires_at is None:
            return None
        return max(0.0, self.expires_at - time.time())
    
    def age(self) -> float:
        """Get the age of the cache entry in seconds."""
        return time.time() - self.created_at
    
    def record_hit(self) -> None:
        """Record a cache hit."""
        self.hit_count += 1
        self.last_accessed = time.time()
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the cache entry to a dictionary."""
        return {
            "key": self.key,
            "value": self.value,
            "created_at": self.created_at,
            "expires_at": self.expires_at,
            "size_bytes": self.size_bytes,
            "source": self.source,
            "hit_count": self.hit_count,
            "last_accessed": self.last_accessed
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'CacheEntry':
        """Create a cache entry from a dictionary."""
        return cls(
            key=data["key"],
            value=data["value"],
            created_at=data["created_at"],
            expires_at=data["expires_at"],
            size_bytes=data["size_bytes"],
            source=data["source"],
            hit_count=data["hit_count"],
            last_accessed=data["last_accessed"]
        )


class CacheBackend(Generic[T, V]):
    """Abstract base class for cache backends."""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the cache backend.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
    
    def get(self, key: T) -> Optional[CacheEntry[T, V]]:
        """
        Get a value from the cache.
        
        Args:
            key: Cache key
            
        Returns:
            Cache entry or None if not found
        """
        raise NotImplementedError("Subclasses must implement get")
    
    def put(self, entry: CacheEntry[T, V]) -> bool:
        """
        Store a value in the cache.
        
        Args:
            entry: Cache entry to store
            
        Returns:
            True if the value was stored successfully
        """
        raise NotImplementedError("Subclasses must implement put")
    
    def delete(self, key: T) -> bool:
        """
        Delete a value from the cache.
        
        Args:
            key: Cache key
            
        Returns:
            True if the value was deleted successfully
        """
        raise NotImplementedError("Subclasses must implement delete")
    
    def clear(self) -> bool:
        """
        Clear all values from the cache.
        
        Returns:
            True if the cache was cleared successfully
        """
        raise NotImplementedError("Subclasses must implement clear")
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about the cache.
        
        Returns:
            Dictionary with cache statistics
        """
        raise NotImplementedError("Subclasses must implement get_stats")
    
    def close(self) -> None:
        """Close the cache and release resources."""
        pass


class MemoryCache(CacheBackend[T, V]):
    """
    In-memory cache backend.
    
    This class provides a thread-safe in-memory cache with
    size limits, automatic expiration, and LRU eviction.
    """
    
    def __init__(
        self,
        max_size_bytes: int = 100 * 1024 * 1024,  # 100MB
        max_items: int = 10000,
        default_ttl: Optional[float] = 3600.0,  # 1 hour
        cleanup_interval: float = 60.0,  # 1 minute
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the memory cache.
        
        Args:
            max_size_bytes: Maximum cache size in bytes
            max_items: Maximum number of items in the cache
            default_ttl: Default time-to-live in seconds
            cleanup_interval: Interval for cleanup in seconds
            logger: Logger instance
        """
        super().__init__(logger)
        self.max_size_bytes = max_size_bytes
        self.max_items = max_items
        self.default_ttl = default_ttl
        self.cleanup_interval = cleanup_interval
        
        # Cache state
        self.cache: Dict[T, CacheEntry[T, V]] = {}
        self.size_bytes = 0
        self.lock = threading.RLock()
        
        # Cache statistics
        self.stats = {
            "hits": 0,
            "misses": 0,
            "insertions": 0,
            "evictions": 0,
            "expirations": 0,
            "size_bytes": 0,
            "items": 0,
            "last_cleanup": 0.0
        }
        
        # Start cleanup thread
        self._stop_event = threading.Event()
        self._cleanup_thread = threading.Thread(
            target=self._cleanup_loop,
            daemon=True,
            name="memory-cache-cleanup"
        )
        self._cleanup_thread.start()
    
    def get(self, key: T) -> Optional[CacheEntry[T, V]]:
        """
        Get a value from the cache.
        
        Args:
            key: Cache key
            
        Returns:
            Cache entry or None if not found
        """
        with self.lock:
            entry = self.cache.get(key)
            
            if entry is None:
                self.stats["misses"] += 1
                return None
            
            # Check if expired
            if entry.is_expired():
                self.cache.pop(key, None)
                self.size_bytes -= entry.size_bytes
                self.stats["expirations"] += 1
                self.stats["misses"] += 1
                self.stats["size_bytes"] = self.size_bytes
                self.stats["items"] = len(self.cache)
                return None
            
            # Record hit
            entry.record_hit()
            self.stats["hits"] += 1
            return entry
    
    def put(self, entry: CacheEntry[T, V]) -> bool:
        """
        Store a value in the cache.
        
        Args:
            entry: Cache entry to store
            
        Returns:
            True if the value was stored successfully
        """
        # Calculate size if not provided
        if entry.size_bytes <= 0:
            try:
                # Try to estimate size using pickle
                serialized = pickle.dumps(entry.value)
                entry.size_bytes = len(serialized)
            except Exception:
                # Fallback to a conservative estimate
                entry.size_bytes = 1024  # 1KB
        
        # Set expiry if not provided
        if entry.expires_at is None and self.default_ttl is not None:
            entry.expires_at = time.time() + self.default_ttl
        
        with self.lock:
            # Check if we need to evict items
            while (len(self.cache) >= self.max_items or 
                   self.size_bytes + entry.size_bytes > self.max_size_bytes):
                if not self._evict_one():
                    # Can't evict any more items
                    self.logger.warning(
                        f"Cache full, can't store item of size {entry.size_bytes} bytes"
                    )
                    return False
            
            # Store the entry
            old_entry = self.cache.get(entry.key)
            if old_entry:
                self.size_bytes -= old_entry.size_bytes
            
            self.cache[entry.key] = entry
            self.size_bytes += entry.size_bytes
            
            # Update stats
            self.stats["insertions"] += 1
            self.stats["size_bytes"] = self.size_bytes
            self.stats["items"] = len(self.cache)
            
            return True
    
    def delete(self, key: T) -> bool:
        """
        Delete a value from the cache.
        
        Args:
            key: Cache key
            
        Returns:
            True if the value was deleted successfully
        """
        with self.lock:
            entry = self.cache.pop(key, None)
            if entry is None:
                return False
            
            self.size_bytes -= entry.size_bytes
            self.stats["size_bytes"] = self.size_bytes
            self.stats["items"] = len(self.cache)
            return True
    
    def clear(self) -> bool:
        """
        Clear all values from the cache.
        
        Returns:
            True if the cache was cleared successfully
        """
        with self.lock:
            self.cache.clear()
            self.size_bytes = 0
            self.stats["size_bytes"] = 0
            self.stats["items"] = 0
            return True
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about the cache.
        
        Returns:
            Dictionary with cache statistics
        """
        with self.lock:
            stats = self.stats.copy()
            stats["hit_ratio"] = (
                stats["hits"] / max(1, stats["hits"] + stats["misses"])
            )
            return stats
    
    def close(self) -> None:
        """Close the cache and release resources."""
        self._stop_event.set()
        self._cleanup_thread.join(timeout=1.0)
    
    def _evict_one(self) -> bool:
        """
        Evict one item from the cache using LRU policy.
        
        Returns:
            True if an item was evicted
        """
        if not self.cache:
            return False
        
        # Find the least recently used item
        lru_key = min(
            self.cache,
            key=lambda k: self.cache[k].last_accessed
        )
        
        # Remove the item
        entry = self.cache.pop(lru_key)
        self.size_bytes -= entry.size_bytes
        self.stats["evictions"] += 1
        return True
    
    def _cleanup_loop(self) -> None:
        """Background thread for cleaning up expired entries."""
        while not self._stop_event.is_set():
            try:
                # Sleep first to give the cache time to be used
                time.sleep(self.cleanup_interval)
                
                # Run cleanup
                self._cleanup_expired()
            except Exception as e:
                self.logger.error(f"Error in cache cleanup: {str(e)}")
    
    def _cleanup_expired(self) -> None:
        """Remove expired entries from the cache."""
        now = time.time()
        removed_count = 0
        removed_bytes = 0
        
        with self.lock:
            # Identify expired entries
            expired_keys = [
                key for key, entry in self.cache.items()
                if entry.is_expired()
            ]
            
            # Remove expired entries
            for key in expired_keys:
                entry = self.cache.pop(key)
                self.size_bytes -= entry.size_bytes
                removed_count += 1
                removed_bytes += entry.size_bytes
            
            # Update stats
            self.stats["expirations"] += removed_count
            self.stats["size_bytes"] = self.size_bytes
            self.stats["items"] = len(self.cache)
            self.stats["last_cleanup"] = now
            
            if removed_count > 0:
                self.logger.debug(
                    f"Removed {removed_count} expired entries "
                    f"({removed_bytes / 1024:.1f} KB)"
                )


class DiskCache(CacheBackend[str, Any]):
    """
    Disk-based cache backend.
    
    This class provides a persistent cache that stores data on disk,
    allowing cache persistence between application runs.
    """
    
    def __init__(
        self,
        cache_dir: str = ".cache",
        max_size_bytes: int = 1024 * 1024 * 1024,  # 1GB
        max_items: int = 100000,
        default_ttl: Optional[float] = 86400.0,  # 1 day
        cleanup_interval: float = 300.0,  # 5 minutes
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the disk cache.
        
        Args:
            cache_dir: Directory for cache files
            max_size_bytes: Maximum cache size in bytes
            max_items: Maximum number of items in the cache
            default_ttl: Default time-to-live in seconds
            cleanup_interval: Interval for cleanup in seconds
            logger: Logger instance
        """
        super().__init__(logger)
        self.cache_dir = os.path.abspath(cache_dir)
        self.max_size_bytes = max_size_bytes
        self.max_items = max_items
        self.default_ttl = default_ttl
        self.cleanup_interval = cleanup_interval
        
        # Create cache directory if it doesn't exist
        os.makedirs(self.cache_dir, exist_ok=True)
        
        # Metadata database
        self.db_path = os.path.join(self.cache_dir, "metadata.db")
        self.db = None
        self._init_database()
        
        # Cache statistics
        self.stats = {
            "hits": 0,
            "misses": 0,
            "insertions": 0,
            "evictions": 0,
            "expirations": 0,
            "size_bytes": 0,
            "items": 0,
            "last_cleanup": 0.0
        }
        
        # Load metadata and update stats
        self._load_metadata()
        
        # Start cleanup thread
        self._stop_event = threading.Event()
        self._cleanup_thread = threading.Thread(
            target=self._cleanup_loop,
            daemon=True,
            name="disk-cache-cleanup"
        )
        self._cleanup_thread.start()
    
    def _init_database(self) -> None:
        """Initialize the metadata database."""
        self.db = sqlite3.connect(self.db_path, check_same_thread=False)
        self.db.row_factory = sqlite3.Row
        cursor = self.db.cursor()
        
        # Create tables if they don't exist
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS cache_entries (
                key TEXT PRIMARY KEY,
                created_at REAL,
                expires_at REAL,
                size_bytes INTEGER,
                source TEXT,
                hit_count INTEGER,
                last_accessed REAL
            )
        """)
        
        # Create indexes
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_expires_at ON cache_entries(expires_at)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_last_accessed ON cache_entries(last_accessed)")
        
        self.db.commit()
    
    def _load_metadata(self) -> None:
        """Load metadata from the database and update stats."""
        cursor = self.db.cursor()
        
        # Get total size and count
        cursor.execute("SELECT COUNT(*) as count, SUM(size_bytes) as total_size FROM cache_entries")
        row = cursor.fetchone()
        
        self.stats["items"] = row["count"] or 0
        self.stats["size_bytes"] = row["total_size"] or 0
        
        # Check for orphaned files
        self._cleanup_orphaned_files()
    
    def _cleanup_orphaned_files(self) -> None:
        """Remove cache files that don't have corresponding metadata."""
        cursor = self.db.cursor()
        
        # Get all keys from the database
        cursor.execute("SELECT key FROM cache_entries")
        db_keys = {row["key"] for row in cursor.fetchall()}
        
        # Get all files in the cache directory
        cache_files = [
            f for f in os.listdir(self.cache_dir)
            if f.endswith(".cache") and os.path.isfile(os.path.join(self.cache_dir, f))
        ]
        
        # Find orphaned files
        for file_name in cache_files:
            key = file_name[:-6]  # Remove .cache suffix
            if key not in db_keys:
                file_path = os.path.join(self.cache_dir, file_name)
                try:
                    os.remove(file_path)
                    self.logger.debug(f"Removed orphaned cache file: {file_path}")
                except Exception as e:
                    self.logger.error(f"Error removing orphaned cache file: {str(e)}")
    
    def _get_cache_path(self, key: str) -> str:
        """
        Get the file path for a cache key.
        
        Args:
            key: Cache key
            
        Returns:
            Absolute path to the cache file
        """
        # Hash the key to create a safe filename
        hashed_key = hashlib.md5(key.encode("utf-8")).hexdigest()
        return os.path.join(self.cache_dir, f"{hashed_key}.cache")
    
    def get(self, key: str) -> Optional[CacheEntry[str, Any]]:
        """
        Get a value from the cache.
        
        Args:
            key: Cache key
            
        Returns:
            Cache entry or None if not found
        """
        # Get metadata
        cursor = self.db.cursor()
        cursor.execute(
            "SELECT * FROM cache_entries WHERE key = ?",
            (key,)
        )
        row = cursor.fetchone()
        
        if not row:
            self.stats["misses"] += 1
            return None
        
        # Check if expired
        if row["expires_at"] and time.time() > row["expires_at"]:
            # Remove expired entry
            self.delete(key)
            self.stats["expirations"] += 1
            self.stats["misses"] += 1
            return None
        
        # Get value from file
        cache_path = self._get_cache_path(key)
        try:
            with open(cache_path, "rb") as f:
                value = pickle.load(f)
        except Exception as e:
            self.logger.error(f"Error loading cache file: {str(e)}")
            self.delete(key)
            self.stats["misses"] += 1
            return None
        
        # Update hit count and last accessed
        cursor.execute(
            "UPDATE cache_entries SET hit_count = hit_count + 1, last_accessed = ? WHERE key = ?",
            (time.time(), key)
        )
        self.db.commit()
        
        # Create cache entry
        entry = CacheEntry(
            key=key,
            value=value,
            created_at=row["created_at"],
            expires_at=row["expires_at"],
            size_bytes=row["size_bytes"],
            source=row["source"],
            hit_count=row["hit_count"] + 1,
            last_accessed=time.time()
        )
        
        self.stats["hits"] += 1
        return entry
    
    def put(self, entry: CacheEntry[str, Any]) -> bool:
        """
        Store a value in the cache.
        
        Args:
            entry: Cache entry to store
            
        Returns:
            True if the value was stored successfully
        """
        # Set expiry if not provided
        if entry.expires_at is None and self.default_ttl is not None:
            entry.expires_at = time.time() + self.default_ttl
        
        # Check if we need to evict items
        cursor = self.db.cursor()
        cursor.execute("SELECT COUNT(*) as count, SUM(size_bytes) as total_size FROM cache_entries")
        row = cursor.fetchone()
        
        current_count = row["count"] or 0
        current_size = row["total_size"] or 0
        
        while (current_count >= self.max_items or 
               current_size + entry.size_bytes > self.max_size_bytes):
            if not self._evict_one():
                # Can't evict any more items
                self.logger.warning(
                    f"Cache full, can't store item of size {entry.size_bytes} bytes"
                )
                return False
            
            # Refresh stats
            cursor.execute("SELECT COUNT(*) as count, SUM(size_bytes) as total_size FROM cache_entries")
            row = cursor.fetchone()
            current_count = row["count"] or 0
            current_size = row["total_size"] or 0
        
        # Save value to file
        cache_path = self._get_cache_path(entry.key)
        try:
            with open(cache_path, "wb") as f:
                pickle.dump(entry.value, f)
        except Exception as e:
            self.logger.error(f"Error saving cache file: {str(e)}")
            return False
        
        # Update or insert metadata
        cursor.execute(
            """
            INSERT OR REPLACE INTO cache_entries
            (key, created_at, expires_at, size_bytes, source, hit_count, last_accessed)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                entry.key,
                entry.created_at,
                entry.expires_at,
                entry.size_bytes,
                entry.source,
                entry.hit_count,
                entry.last_accessed
            )
        )
        self.db.commit()
        
        # Update stats
        self.stats["insertions"] += 1
        self.stats["items"] = current_count + 1
        self.stats["size_bytes"] = current_size + entry.size_bytes
        
        return True
    
    def delete(self, key: str) -> bool:
        """
        Delete a value from the cache.
        
        Args:
            key: Cache key
            
        Returns:
            True if the value was deleted successfully
        """
        # Get metadata
        cursor = self.db.cursor()
        cursor.execute(
            "SELECT size_bytes FROM cache_entries WHERE key = ?",
            (key,)
        )
        row = cursor.fetchone()
        
        if not row:
            return False
        
        # Delete file
        cache_path = self._get_cache_path(key)
        try:
            if os.path.exists(cache_path):
                os.remove(cache_path)
        except Exception as e:
            self.logger.error(f"Error removing cache file: {str(e)}")
        
        # Delete metadata
        cursor.execute(
            "DELETE FROM cache_entries WHERE key = ?",
            (key,)
        )
        self.db.commit()
        
        # Update stats
        self.stats["items"] -= 1
        self.stats["size_bytes"] -= row["size_bytes"]
        
        return True
    
    def clear(self) -> bool:
        """
        Clear all values from the cache.
        
        Returns:
            True if the cache was cleared successfully
        """
        # Delete all cache files
        cache_files = [
            f for f in os.listdir(self.cache_dir)
            if f.endswith(".cache") and os.path.isfile(os.path.join(self.cache_dir, f))
        ]
        
        for file_name in cache_files:
            try:
                os.remove(os.path.join(self.cache_dir, file_name))
            except Exception as e:
                self.logger.error(f"Error removing cache file: {str(e)}")
        
        # Clear metadata
        cursor = self.db.cursor()
        cursor.execute("DELETE FROM cache_entries")
        self.db.commit()
        
        # Update stats
        self.stats["items"] = 0
        self.stats["size_bytes"] = 0
        
        return True
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about the cache.
        
        Returns:
            Dictionary with cache statistics
        """
        cursor = self.db.cursor()
        
        # Get total size and count
        cursor.execute("SELECT COUNT(*) as count, SUM(size_bytes) as total_size FROM cache_entries")
        row = cursor.fetchone()
        
        self.stats["items"] = row["count"] or 0
        self.stats["size_bytes"] = row["total_size"] or 0
        
        stats = self.stats.copy()
        stats["hit_ratio"] = (
            stats["hits"] / max(1, stats["hits"] + stats["misses"])
        )
        
        # Get additional stats
        cursor.execute("SELECT COUNT(*) as count FROM cache_entries WHERE expires_at < ?", (time.time(),))
        stats["expired_items"] = cursor.fetchone()["count"] or 0
        
        cursor.execute("SELECT AVG(hit_count) as avg_hits FROM cache_entries")
        stats["avg_hits"] = cursor.fetchone()["avg_hits"] or 0
        
        return stats
    
    def close(self) -> None:
        """Close the cache and release resources."""
        self._stop_event.set()
        self._cleanup_thread.join(timeout=1.0)
        
        if self.db:
            self.db.close()
            self.db = None
    
    def _evict_one(self) -> bool:
        """
        Evict one item from the cache using LRU policy.
        
        Returns:
            True if an item was evicted
        """
        cursor = self.db.cursor()
        
        # Find the least recently used item
        cursor.execute(
            "SELECT key FROM cache_entries ORDER BY last_accessed ASC LIMIT 1"
        )
        row = cursor.fetchone()
        
        if not row:
            return False
        
        # Delete the item
        return self.delete(row["key"])
    
    def _cleanup_loop(self) -> None:
        """Background thread for cleaning up expired entries."""
        while not self._stop_event.is_set():
            try:
                # Sleep first to give the cache time to be used
                time.sleep(self.cleanup_interval)
                
                # Run cleanup
                self._cleanup_expired()
            except Exception as e:
                self.logger.error(f"Error in cache cleanup: {str(e)}")
    
    def _cleanup_expired(self) -> None:
        """Remove expired entries from the cache."""
        now = time.time()
        cursor = self.db.cursor()
        
        # Get expired entries
        cursor.execute(
            "SELECT key FROM cache_entries WHERE expires_at < ? AND expires_at IS NOT NULL",
            (now,)
        )
        expired_keys = [row["key"] for row in cursor.fetchall()]
        
        # Remove expired entries
        removed_count = 0
        for key in expired_keys:
            if self.delete(key):
                removed_count += 1
        
        # Update stats
        self.stats["expirations"] += removed_count
        self.stats["last_cleanup"] = now
        
        if removed_count > 0:
            self.logger.debug(f"Removed {removed_count} expired cache entries")


class CacheManager:
    """
    Manages multiple cache backends with different strategies.
    
    This class provides a unified interface for working with
    multiple cache backends (e.g., memory, disk) with different
    policies based on data types and sources.
    """
    
    def __init__(
        self,
        memory_cache: Optional[MemoryCache] = None,
        disk_cache: Optional[DiskCache] = None,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the cache manager.
        
        Args:
            memory_cache: Memory cache backend
            disk_cache: Disk cache backend
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
        # Create default caches if not provided
        self.memory_cache = memory_cache or MemoryCache(logger=self.logger)
        self.disk_cache = disk_cache or DiskCache(logger=self.logger)
        
        # Cache policy configuration
        self.policies = {
            "default": {
                "backend": "memory",
                "ttl": 3600.0  # 1 hour
            },
            "pubchem": {
                "backend": "disk",
                "ttl": 86400.0 * 7  # 7 days
            },
            "chembl": {
                "backend": "disk",
                "ttl": 86400.0 * 30  # 30 days
            },
            "transient": {
                "backend": "memory",
                "ttl": 300.0  # 5 minutes
            }
        }
    
    def get_backend(self, policy: str) -> CacheBackend:
        """
        Get the appropriate cache backend for a policy.
        
        Args:
            policy: Cache policy name
            
        Returns:
            Cache backend
        """
        if policy not in self.policies:
            policy = "default"
        
        backend_name = self.policies[policy]["backend"]
        
        if backend_name == "memory":
            return self.memory_cache
        elif backend_name == "disk":
            return self.disk_cache
        else:
            return self.memory_cache
    
    def get_ttl(self, policy: str) -> Optional[float]:
        """
        Get the TTL for a policy.
        
        Args:
            policy: Cache policy name
            
        Returns:
            TTL in seconds or None for no expiration
        """
        if policy not in self.policies:
            policy = "default"
        
        return self.policies[policy].get("ttl")
    
    def get(
        self,
        key: str,
        policy: str = "default"
    ) -> Optional[Any]:
        """
        Get a value from the cache.
        
        Args:
            key: Cache key
            policy: Cache policy name
            
        Returns:
            Cached value or None if not found
        """
        backend = self.get_backend(policy)
        entry = backend.get(key)
        
        if entry is None:
            return None
        
        return entry.value
    
    def put(
        self,
        key: str,
        value: Any,
        policy: str = "default",
        ttl: Optional[float] = None,
        source: str = ""
    ) -> bool:
        """
        Store a value in the cache.
        
        Args:
            key: Cache key
            value: Value to cache
            policy: Cache policy name
            ttl: Time-to-live in seconds (overrides policy TTL)
            source: Source of the data
            
        Returns:
            True if the value was stored successfully
        """
        backend = self.get_backend(policy)
        
        # Use provided TTL or get from policy
        if ttl is None:
            ttl = self.get_ttl(policy)
        
        # Calculate expiry time
        expires_at = None
        if ttl is not None:
            expires_at = time.time() + ttl
        
        # Create cache entry
        entry = CacheEntry(
            key=key,
            value=value,
            expires_at=expires_at,
            source=source
        )
        
        return backend.put(entry)
    
    def delete(
        self,
        key: str,
        policy: str = "default"
    ) -> bool:
        """
        Delete a value from the cache.
        
        Args:
            key: Cache key
            policy: Cache policy name
            
        Returns:
            True if the value was deleted successfully
        """
        backend = self.get_backend(policy)
        return backend.delete(key)
    
    def clear(
        self,
        policy: Optional[str] = None
    ) -> bool:
        """
        Clear values from the cache.
        
        Args:
            policy: Cache policy name (None to clear all caches)
            
        Returns:
            True if the cache was cleared successfully
        """
        if policy is None:
            # Clear all caches
            self.memory_cache.clear()
            self.disk_cache.clear()
            return True
        
        backend = self.get_backend(policy)
        return backend.clear()
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about all caches.
        
        Returns:
            Dictionary with cache statistics
        """
        return {
            "memory": self.memory_cache.get_stats(),
            "disk": self.disk_cache.get_stats(),
            "policies": self.policies
        }
    
    def close(self) -> None:
        """Close all caches and release resources."""
        self.memory_cache.close()
        self.disk_cache.close()


class AsyncCacheManager:
    """
    Asynchronous wrapper for the cache manager.
    
    This class provides async/await support for cache operations
    by running them in a thread pool.
    """
    
    def __init__(
        self,
        cache_manager: Optional[CacheManager] = None,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the async cache manager.
        
        Args:
            cache_manager: Cache manager
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        self.cache_manager = cache_manager or CacheManager(logger=self.logger)
    
    async def get(
        self,
        key: str,
        policy: str = "default"
    ) -> Optional[Any]:
        """
        Get a value from the cache asynchronously.
        
        Args:
            key: Cache key
            policy: Cache policy name
            
        Returns:
            Cached value or None if not found
        """
        return await asyncio.to_thread(
            self.cache_manager.get,
            key,
            policy
        )
    
    async def put(
        self,
        key: str,
        value: Any,
        policy: str = "default",
        ttl: Optional[float] = None,
        source: str = ""
    ) -> bool:
        """
        Store a value in the cache asynchronously.
        
        Args:
            key: Cache key
            value: Value to cache
            policy: Cache policy name
            ttl: Time-to-live in seconds (overrides policy TTL)
            source: Source of the data
            
        Returns:
            True if the value was stored successfully
        """
        return await asyncio.to_thread(
            self.cache_manager.put,
            key,
            value,
            policy,
            ttl,
            source
        )
    
    async def delete(
        self,
        key: str,
        policy: str = "default"
    ) -> bool:
        """
        Delete a value from the cache asynchronously.
        
        Args:
            key: Cache key
            policy: Cache policy name
            
        Returns:
            True if the value was deleted successfully
        """
        return await asyncio.to_thread(
            self.cache_manager.delete,
            key,
            policy
        )
    
    async def clear(
        self,
        policy: Optional[str] = None
    ) -> bool:
        """
        Clear values from the cache asynchronously.
        
        Args:
            policy: Cache policy name (None to clear all caches)
            
        Returns:
            True if the cache was cleared successfully
        """
        return await asyncio.to_thread(
            self.cache_manager.clear,
            policy
        )
    
    async def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about all caches asynchronously.
        
        Returns:
            Dictionary with cache statistics
        """
        return await asyncio.to_thread(
            self.cache_manager.get_stats
        )
    
    def close(self) -> None:
        """Close all caches and release resources."""
        self.cache_manager.close()


# Higher-level API functions

def create_cache_key(
    prefix: str,
    parts: List[Any],
    max_length: int = 200
) -> str:
    """
    Create a deterministic cache key from multiple parts.
    
    Args:
        prefix: Key prefix for namespace
        parts: List of key parts (will be converted to strings)
        max_length: Maximum key length
        
    Returns:
        Cache key string
    """
    # Convert all parts to strings
    str_parts = [str(part) for part in parts]
    
    # Join parts
    key = f"{prefix}:{':'.join(str_parts)}"
    
    # If the key is too long, hash the excess
    if len(key) > max_length:
        excess = key[max_length - 40:]
        hashed = hashlib.md5(excess.encode("utf-8")).hexdigest()
        key = f"{key[:max_length - 40]}{hashed}"
    
    return key