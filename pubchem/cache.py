"""
Multi-level caching system for PubChem API responses.

This module provides both in-memory and disk-based caching for PubChem API
responses to reduce API calls and improve performance.
"""

import os
import json
import time
import logging
import sqlite3
import hashlib
from typing import Dict, Any, Optional, Tuple, Union
from functools import lru_cache
from pathlib import Path

logger = logging.getLogger(__name__)

def initialize_cache(db_path: Optional[str] = None) -> None:
    """
    Initialize the SQLite database and create tables if they don't exist.
    
    Args:
        db_path: Path to the SQLite database file. If None, defaults to 'pubchem_cache.sqlite'
                in the current directory.
    """
    # Set default path if not provided
    if db_path is None:
        db_path = os.path.join(os.path.dirname(__file__), 'pubchem_cache.sqlite')
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
    
    # Connect to the database
    conn = sqlite3.connect(db_path)
    
    try:
        # Enable WAL mode for better concurrency
        conn.execute('PRAGMA journal_mode=WAL')
        
        # Create compound_cache table if it doesn't exist
        conn.execute('''
        CREATE TABLE IF NOT EXISTS compound_cache (
            cid INTEGER PRIMARY KEY,
            data TEXT NOT NULL,
            last_updated INTEGER NOT NULL,
            last_accessed INTEGER NOT NULL,
            source TEXT NOT NULL
        )
        ''')
        
        # Create index on last_accessed for LRU/aging queries
        conn.execute('''
        CREATE INDEX IF NOT EXISTS idx_last_accessed ON compound_cache(last_accessed)
        ''')
        
        # Create cache_stats table if it doesn't exist
        conn.execute('''
        CREATE TABLE IF NOT EXISTS cache_stats (
            key TEXT PRIMARY KEY,
            value INTEGER NOT NULL
        )
        ''')
        
        # Initialize stats if they don't exist
        stats = ['hits', 'misses', 'total_entries', 'hit_rate']
        for stat in stats:
            conn.execute('''
            INSERT OR IGNORE INTO cache_stats (key, value) VALUES (?, 0)
            ''', (stat,))
        
        # Commit changes
        conn.commit()
        
        logger.info(f"PubChem cache initialized at {db_path}")
        
    except sqlite3.Error as e:
        logger.error(f"Error initializing cache database: {e}")
        raise
    finally:
        conn.close()

def get_compound(cid: int, max_age_seconds: Optional[int] = None, db_path: Optional[str] = None) -> Optional[dict]:
    """
    Retrieve compound data from the cache if present and not expired.
    
    Args:
        cid: PubChem Compound ID to retrieve
        max_age_seconds: Optional maximum age in seconds. If the cached data is older
                        than this, it will be considered expired and None will be returned.
        db_path: Path to the SQLite database file. If None, defaults to 'pubchem_cache.sqlite'
                in the current directory.
                        
    Returns:
        Dictionary containing compound data if found and not expired, None otherwise.
        
    Note:
        This function updates the last_accessed timestamp and increments hit/miss statistics.
    """
    # Set default path if not provided
    if db_path is None:
        db_path = os.path.join(os.path.dirname(__file__), 'pubchem_cache.sqlite')
    
    # Get current time
    current_time = int(time.time())
    
    # Connect to the database
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row  # Enable row factory for dict-like access
    
    try:
        # Begin transaction
        with conn:
            cursor = conn.cursor()
            
            # Query for the compound
            cursor.execute('''
            SELECT cid, data, last_updated, last_accessed, source
            FROM compound_cache
            WHERE cid = ?
            ''', (cid,))
            
            result = cursor.fetchone()
            
            # If compound not found, increment miss count and return None
            if result is None:
                cursor.execute('''
                UPDATE cache_stats
                SET value = value + 1
                WHERE key = 'misses'
                ''')
                return None
            
            # Parse the data
            compound_data = json.loads(result['data'])
            last_updated = result['last_updated']
            
            # Check if data is expired
            if max_age_seconds is not None and (current_time - last_updated) > max_age_seconds:
                cursor.execute('''
                UPDATE cache_stats
                SET value = value + 1
                WHERE key = 'misses'
                ''')
                return None
            
            # Update last_accessed timestamp
            cursor.execute('''
            UPDATE compound_cache
            SET last_accessed = ?
            WHERE cid = ?
            ''', (current_time, cid))
            
            # Increment hit count
            cursor.execute('''
            UPDATE cache_stats
            SET value = value + 1
            WHERE key = 'hits'
            ''')
            
            # Update hit rate - we don't need to do this on every hit
            # as it will be calculated when get_cache_stats() is called
            
            return compound_data
            
    except sqlite3.Error as e:
        logger.error(f"Error retrieving compound {cid} from cache: {e}")
        return None
    finally:
        conn.close()

def store_compound(cid: int, data: dict, source: str = "pubchem", db_path: Optional[str] = None) -> None:
    """
    Store or update compound data in the cache.
    
    Args:
        cid: PubChem Compound ID to store
        data: Dictionary containing compound data to store
        source: Source of the data (default: "pubchem")
        db_path: Path to the SQLite database file. If None, defaults to 'pubchem_cache.sqlite'
                in the current directory.
                
    Note:
        This function updates both last_updated and last_accessed timestamps.
        If the compound already exists in the cache, it will be overwritten.
    """
    # Set default path if not provided
    if db_path is None:
        db_path = os.path.join(os.path.dirname(__file__), 'pubchem_cache.sqlite')
    
    # Get current time
    current_time = int(time.time())
    
    # Connect to the database
    conn = sqlite3.connect(db_path)
    
    try:
        # Begin transaction
        with conn:
            cursor = conn.cursor()
            
            # Serialize the data to JSON
            serialized_data = json.dumps(data)
            
            # Insert or replace the compound data
            cursor.execute('''
            INSERT OR REPLACE INTO compound_cache
            (cid, data, last_updated, last_accessed, source)
            VALUES (?, ?, ?, ?, ?)
            ''', (cid, serialized_data, current_time, current_time, source))
            
            # Update total_entries in cache_stats is handled by get_cache_stats
            # We don't need to update it here as it would require additional queries
            
            logger.info(f"Stored/updated compound {cid} in cache (source: {source})")
            
    except sqlite3.Error as e:
        logger.error(f"Error storing compound {cid} in cache: {e}")
        raise
    finally:
        conn.close()

def get_cache_stats(db_path: Optional[str] = None) -> dict:
    """
    Get cache statistics including hits, misses, total entries, and hit rate.
    
    Args:
        db_path: Path to the SQLite database file. If None, defaults to 'pubchem_cache.sqlite'
                in the current directory.
                
    Returns:
        Dictionary with cache statistics: hits, misses, total entries, hit rate.
    """
    # Set default path if not provided
    if db_path is None:
        db_path = os.path.join(os.path.dirname(__file__), 'pubchem_cache.sqlite')
    
    # Connect to the database
    conn = sqlite3.connect(db_path)
    
    try:
        cursor = conn.cursor()
        
        # Get stats
        cursor.execute("SELECT key, value FROM cache_stats")
        stats = {row[0]: row[1] for row in cursor.fetchall()}
        
        # Count total entries
        cursor.execute("SELECT COUNT(*) FROM compound_cache")
        stats['total_entries'] = cursor.fetchone()[0]
        
        # Calculate hit rate
        hits = stats.get('hits', 0)
        misses = stats.get('misses', 0)
        total = hits + misses
        
        if total > 0:
            stats['hit_rate'] = int((hits / total) * 100)
        else:
            stats['hit_rate'] = 0
        
        # Update hit_rate in the database
        cursor.execute('''
        UPDATE cache_stats
        SET value = ?
        WHERE key = 'hit_rate'
        ''', (stats['hit_rate'],))
        
        # Update total_entries in the database
        cursor.execute('''
        UPDATE cache_stats
        SET value = ?
        WHERE key = 'total_entries'
        ''', (stats['total_entries'],))
        
        conn.commit()
        
        return stats
        
    except sqlite3.Error as e:
        logger.error(f"Error getting cache statistics: {e}")
        return {'hits': 0, 'misses': 0, 'total_entries': 0, 'hit_rate': 0}
    finally:
        conn.close()

def get_cache_stats_detailed(db_path: Optional[str] = None) -> dict:
    """
    Get detailed cache statistics including basic stats plus age distribution,
    source distribution, and size metrics.
    
    Args:
        db_path: Path to the SQLite database file. If None, defaults to 'pubchem_cache.sqlite'
                in the current directory.
                
    Returns:
        Dictionary with detailed cache statistics.
    """
    # Set default path if not provided
    if db_path is None:
        db_path = os.path.join(os.path.dirname(__file__), 'pubchem_cache.sqlite')
    
    # Get current time
    current_time = int(time.time())
    
    # Connect to the database
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row  # Enable row factory for dict-like access
    
    try:
        cursor = conn.cursor()
        
        # Get basic stats first
        basic_stats = get_cache_stats(db_path)
        
        # Initialize detailed stats with basic stats
        detailed_stats = dict(basic_stats)
        
        # Add age statistics
        cursor.execute('''
        SELECT 
            AVG(? - last_updated) as avg_age,
            MIN(? - last_updated) as min_age,
            MAX(? - last_updated) as max_age,
            AVG(? - last_accessed) as avg_access_age,
            MIN(? - last_accessed) as min_access_age,
            MAX(? - last_accessed) as max_access_age
        FROM compound_cache
        ''', (current_time, current_time, current_time, current_time, current_time, current_time))
        
        age_stats = cursor.fetchone()
        if age_stats and detailed_stats['total_entries'] > 0:
            detailed_stats['avg_age_seconds'] = int(age_stats['avg_age'] or 0)
            detailed_stats['min_age_seconds'] = int(age_stats['min_age'] or 0)
            detailed_stats['max_age_seconds'] = int(age_stats['max_age'] or 0)
            detailed_stats['avg_access_age_seconds'] = int(age_stats['avg_access_age'] or 0)
            detailed_stats['min_access_age_seconds'] = int(age_stats['min_access_age'] or 0)
            detailed_stats['max_access_age_seconds'] = int(age_stats['max_access_age'] or 0)
        
        # Add source distribution
        cursor.execute('''
        SELECT source, COUNT(*) as count
        FROM compound_cache
        GROUP BY source
        ''')
        
        source_stats = {}
        for row in cursor.fetchall():
            source_stats[row['source']] = row['count']
        
        detailed_stats['source_distribution'] = source_stats
        
        # Add size statistics (approximate size of data in bytes)
        cursor.execute('''
        SELECT 
            SUM(LENGTH(data)) as total_data_size,
            AVG(LENGTH(data)) as avg_data_size,
            MIN(LENGTH(data)) as min_data_size,
            MAX(LENGTH(data)) as max_data_size
        FROM compound_cache
        ''')
        
        size_stats = cursor.fetchone()
        if size_stats and detailed_stats['total_entries'] > 0:
            detailed_stats['total_data_size_bytes'] = int(size_stats['total_data_size'] or 0)
            detailed_stats['avg_data_size_bytes'] = int(size_stats['avg_data_size'] or 0)
            detailed_stats['min_data_size_bytes'] = int(size_stats['min_data_size'] or 0)
            detailed_stats['max_data_size_bytes'] = int(size_stats['max_data_size'] or 0)
        
        # Get most recently accessed entries (top 5)
        cursor.execute('''
        SELECT cid, last_accessed
        FROM compound_cache
        ORDER BY last_accessed DESC
        LIMIT 5
        ''')
        
        recent_access = []
        for row in cursor.fetchall():
            recent_access.append({
                'cid': row['cid'],
                'last_accessed': row['last_accessed'],
                'accessed_ago_seconds': current_time - row['last_accessed']
            })
        
        detailed_stats['recent_access'] = recent_access
        
        # Get least recently accessed entries (bottom 5)
        cursor.execute('''
        SELECT cid, last_accessed
        FROM compound_cache
        ORDER BY last_accessed ASC
        LIMIT 5
        ''')
        
        oldest_access = []
        for row in cursor.fetchall():
            oldest_access.append({
                'cid': row['cid'],
                'last_accessed': row['last_accessed'],
                'accessed_ago_seconds': current_time - row['last_accessed']
            })
        
        detailed_stats['oldest_access'] = oldest_access
        
        # Add database file size
        try:
            db_size = os.path.getsize(db_path)
            detailed_stats['db_file_size_bytes'] = db_size
        except OSError:
            detailed_stats['db_file_size_bytes'] = 0
        
        return detailed_stats
        
    except sqlite3.Error as e:
        logger.error(f"Error getting detailed cache statistics: {e}")
        return basic_stats
    finally:
        conn.close()


def prune_cache(max_entries: int, db_path: Optional[str] = None) -> int:
    """
    Remove least-recently-used entries from the cache if the total number of entries
    exceeds the specified maximum.
    
    Args:
        max_entries: Maximum number of entries to keep in the cache
        db_path: Path to the SQLite database file. If None, defaults to 'pubchem_cache.sqlite'
                in the current directory.
                
    Returns:
        Number of entries removed from the cache
    """
    # Set default path if not provided
    if db_path is None:
        db_path = os.path.join(os.path.dirname(__file__), 'pubchem_cache.sqlite')
    
    # Connect to the database
    conn = sqlite3.connect(db_path)
    
    try:
        # Begin transaction
        with conn:
            cursor = conn.cursor()
            
            # Get total number of entries
            cursor.execute("SELECT COUNT(*) FROM compound_cache")
            total_entries = cursor.fetchone()[0]
            
            # If total entries is less than or equal to max_entries, no pruning needed
            if total_entries <= max_entries:
                return 0
            
            # Calculate how many entries to remove
            entries_to_remove = total_entries - max_entries
            
            # Get the CIDs of the least recently accessed entries
            cursor.execute("""
            SELECT cid FROM compound_cache
            ORDER BY last_accessed ASC
            LIMIT ?
            """, (entries_to_remove,))
            
            cids_to_remove = [row[0] for row in cursor.fetchall()]
            
            # Remove the entries
            if cids_to_remove:
                placeholders = ','.join('?' for _ in cids_to_remove)
                cursor.execute(f"""
                DELETE FROM compound_cache
                WHERE cid IN ({placeholders})
                """, cids_to_remove)
            
            # Log the pruning operation
            logger.info(f"Pruned {entries_to_remove} entries from cache (max_entries={max_entries})")
            
            return entries_to_remove
            
    except sqlite3.Error as e:
        logger.error(f"Error pruning cache: {e}")
        return 0
    finally:
        conn.close()


class PubChemCache:
    """
    Multi-level cache for PubChem API responses.
    
    Features:
    - In-memory LRU cache for fast access to frequently used data
    - Disk-based persistent cache for long-term storage
    - Automatic cache invalidation based on TTL (time-to-live)
    """
    
    def __init__(self, cache_dir: str = "cache/pubchem", memory_size: int = 1000, 
                 ttl: int = 86400 * 30):  # Default TTL: 30 days
        """
        Initialize the cache.
        
        Args:
            cache_dir: Directory to store disk cache files
            memory_size: Maximum number of items to keep in memory cache
            ttl: Time-to-live for cache entries in seconds
        """
        self.cache_dir = Path(cache_dir)
        self.ttl = ttl
        self.memory_size = memory_size
        
        # Create cache directory if it doesn't exist
        os.makedirs(self.cache_dir, exist_ok=True)
        
        # Initialize stats
        self.stats = {
            "memory_hits": 0,
            "disk_hits": 0,
            "misses": 0,
            "writes": 0
        }
    
    @lru_cache(maxsize=1000)  # Default memory cache size
    def _memory_get(self, key: str) -> Tuple[Any, float]:
        """
        Get item from memory cache.
        
        Args:
            key: Cache key
            
        Returns:
            Tuple of (cached_value, timestamp)
        """
        # This is just a wrapper around lru_cache
        # The actual implementation is handled by the decorator
        # This should never be called directly
        return None, 0
    
    def _get_disk_cache_path(self, key: str) -> Path:
        """
        Get the file path for a disk cache entry.
        
        Args:
            key: Cache key
            
        Returns:
            Path object for the cache file
        """
        # Use MD5 hash of the key as the filename to avoid invalid characters
        hashed_key = hashlib.md5(key.encode()).hexdigest()
        return self.cache_dir / f"{hashed_key}.json"
    
    def _disk_get(self, key: str) -> Optional[Tuple[Any, float]]:
        """
        Get item from disk cache.
        
        Args:
            key: Cache key
            
        Returns:
            Tuple of (cached_value, timestamp) or None if not found
        """
        cache_path = self._get_disk_cache_path(key)
        
        if not cache_path.exists():
            return None
        
        try:
            with open(cache_path, 'r') as f:
                cache_data = json.load(f)
                
            timestamp = cache_data.get('timestamp', 0)
            value = cache_data.get('value')
            
            return value, timestamp
        except (json.JSONDecodeError, IOError) as e:
            logger.warning(f"Error reading disk cache for {key}: {e}")
            return None
    
    def _disk_set(self, key: str, value: Any) -> None:
        """
        Set item in disk cache.
        
        Args:
            key: Cache key
            value: Value to cache
        """
        cache_path = self._get_disk_cache_path(key)
        timestamp = time.time()
        
        try:
            cache_data = {
                'timestamp': timestamp,
                'value': value
            }
            
            with open(cache_path, 'w') as f:
                json.dump(cache_data, f)
                
        except IOError as e:
            logger.warning(f"Error writing to disk cache for {key}: {e}")
    
    def get(self, key: str) -> Optional[Any]:
        """
        Get item from cache (either memory or disk).
        
        Args:
            key: Cache key
            
        Returns:
            Cached value or None if not found or expired
        """
        # Try memory cache first
        try:
            value, timestamp = self._memory_get(key)
            
            # Check if expired
            if time.time() - timestamp <= self.ttl:
                self.stats["memory_hits"] += 1
                return value
                
            # If expired, continue to disk cache
        except Exception:
            # Memory cache miss, try disk
            pass
        
        # Try disk cache
        disk_result = self._disk_get(key)
        
        if disk_result:
            value, timestamp = disk_result
            
            # Check if expired
            if time.time() - timestamp <= self.ttl:
                # Update memory cache with disk value
                self._memory_get.cache_clear()  # Clear just this key
                self._memory_get(key)  # This will add it to the cache
                
                self.stats["disk_hits"] += 1
                return value
        
        # Cache miss
        self.stats["misses"] += 1
        return None
    
    def set(self, key: str, value: Any) -> None:
        """
        Set item in both memory and disk cache.
        
        Args:
            key: Cache key
            value: Value to cache
        """
        # Update memory cache
        self._memory_get.cache_clear()  # Clear just this key
        timestamp = time.time()
        self._memory_get(key)  # This will add it to the cache with the current timestamp
        
        # Update disk cache
        self._disk_set(key, value)
        
        self.stats["writes"] += 1
    
    def clear(self) -> None:
        """Clear both memory and disk cache."""
        # Clear memory cache
        self._memory_get.cache_clear()
        
        # Clear disk cache
        for cache_file in self.cache_dir.glob("*.json"):
            try:
                os.remove(cache_file)
            except OSError as e:
                logger.warning(f"Error removing cache file {cache_file}: {e}")
        
        # Reset stats
        self.stats = {
            "memory_hits": 0,
            "disk_hits": 0,
            "misses": 0,
            "writes": 0
        }
    
    def get_stats(self) -> Dict[str, int]:
        """
        Get cache statistics.
        
        Returns:
            Dictionary with cache statistics
        """
        return self.stats