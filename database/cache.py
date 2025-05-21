#!/usr/bin/env python3
"""
Redis-based caching system for database queries.

This module provides a caching layer for database queries to improve performance
for frequently accessed data. It uses Redis as the caching backend and supports
various caching strategies and invalidation mechanisms.
"""

import os
import json
import time
import hashlib
import logging
import functools
from typing import Any, Dict, List, Optional, Union, Callable, Tuple
from datetime import datetime, timedelta

import redis
from redis.exceptions import RedisError

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Default cache configuration
DEFAULT_CACHE_TTL = 3600  # 1 hour in seconds
DEFAULT_CACHE_PREFIX = 'cryoprotect:'
DEFAULT_CACHE_ENABLED = True

# Initialize Redis connection
_redis_client = None

def get_redis_client() -> redis.Redis:
    """Get Redis client singleton."""
    global _redis_client
    
    if _redis_client is None:
        # Get Redis connection details from environment variables
        redis_url = os.environ.get('REDIS_URL')
        redis_host = os.environ.get('REDIS_HOST', 'localhost')
        redis_port = int(os.environ.get('REDIS_PORT', 6379))
        redis_db = int(os.environ.get('REDIS_DB', 0))
        redis_password = os.environ.get('REDIS_PASSWORD')
        
        try:
            if redis_url:
                _redis_client = redis.from_url(redis_url)
                logger.info(f"Connected to Redis using URL: {redis_url.split('@')[-1]}")
            else:
                _redis_client = redis.Redis(
                    host=redis_host,
                    port=redis_port,
                    db=redis_db,
                    password=redis_password,
                    decode_responses=False,  # We'll handle decoding ourselves
                    socket_timeout=5,        # 5 second timeout for operations
                    socket_connect_timeout=5 # 5 second timeout for connections
                )
                logger.info(f"Connected to Redis at {redis_host}:{redis_port} (DB: {redis_db})")
            
            # Test connection
            _redis_client.ping()
        except RedisError as e:
            logger.warning(f"Failed to connect to Redis: {e}")
            
            # Create a dummy implementation that doesn't cache anything
            class DummyRedis:
                def get(self, *args, **kwargs):
                    return None
                def set(self, *args, **kwargs):
                    pass
                def delete(self, *args, **kwargs):
                    pass
                def exists(self, *args, **kwargs):
                    return False
                def expire(self, *args, **kwargs):
                    pass
                def ping(self, *args, **kwargs):
                    return True
                def flushdb(self, *args, **kwargs):
                    pass
                def keys(self, *args, **kwargs):
                    return []
                def scan_iter(self, *args, **kwargs):
                    return []
                def pipeline(self, *args, **kwargs):
                    return self
                def execute(self, *args, **kwargs):
                    return []
            
            _redis_client = DummyRedis()
            logger.warning("Using dummy Redis client (no caching)")
    
    return _redis_client

def generate_cache_key(namespace: str, *args, **kwargs) -> str:
    """
    Generate a cache key from the namespace and arguments.
    
    Args:
        namespace: The namespace for the cache key
        *args: Positional arguments to include in the key
        **kwargs: Keyword arguments to include in the key
        
    Returns:
        A cache key string
    """
    # Create a string representation of the arguments
    key_parts = [str(arg) for arg in args]
    key_parts.extend(f"{k}={v}" for k, v in sorted(kwargs.items()))
    
    # Join the parts and create a hash
    key_str = '|'.join(key_parts)
    key_hash = hashlib.md5(key_str.encode('utf-8')).hexdigest()
    
    # Create the final key with namespace and hash
    cache_key = f"{DEFAULT_CACHE_PREFIX}{namespace}:{key_hash}"
    
    return cache_key

def cache_get(key: str) -> Optional[Any]:
    """
    Get a value from the cache.
    
    Args:
        key: The cache key
        
    Returns:
        The cached value or None if not found
    """
    if not DEFAULT_CACHE_ENABLED:
        return None
    
    client = get_redis_client()
    try:
        # Get the value from Redis
        cached_data = client.get(key)
        
        if cached_data:
            # Deserialize the JSON data
            try:
                return json.loads(cached_data)
            except json.JSONDecodeError:
                logger.warning(f"Failed to deserialize cached data for key: {key}")
                return None
        
        return None
    except RedisError as e:
        logger.warning(f"Redis error getting key {key}: {e}")
        return None

def cache_set(key: str, value: Any, ttl: int = DEFAULT_CACHE_TTL) -> bool:
    """
    Set a value in the cache.
    
    Args:
        key: The cache key
        value: The value to cache (must be JSON serializable)
        ttl: Time to live in seconds
        
    Returns:
        True if successful, False otherwise
    """
    if not DEFAULT_CACHE_ENABLED:
        return False
    
    client = get_redis_client()
    try:
        # Serialize the value to JSON
        serialized_value = json.dumps(value)
        
        # Set the value in Redis with TTL
        client.set(key, serialized_value, ex=ttl)
        return True
    except (RedisError, TypeError, json.JSONDecodeError) as e:
        logger.warning(f"Redis error setting key {key}: {e}")
        return False

def cache_delete(key: str) -> bool:
    """
    Delete a value from the cache.
    
    Args:
        key: The cache key
        
    Returns:
        True if successful, False otherwise
    """
    if not DEFAULT_CACHE_ENABLED:
        return False
    
    client = get_redis_client()
    try:
        client.delete(key)
        return True
    except RedisError as e:
        logger.warning(f"Redis error deleting key {key}: {e}")
        return False

def cache_invalidate_pattern(pattern: str) -> int:
    """
    Invalidate all cache keys matching a pattern.
    
    Args:
        pattern: The pattern to match (e.g., 'molecules:*')
        
    Returns:
        Number of keys invalidated
    """
    if not DEFAULT_CACHE_ENABLED:
        return 0
    
    client = get_redis_client()
    try:
        # Find all keys matching the pattern
        key_pattern = f"{DEFAULT_CACHE_PREFIX}{pattern}"
        matching_keys = list(client.scan_iter(match=key_pattern))
        
        if matching_keys:
            # Delete all matching keys
            client.delete(*matching_keys)
            
        return len(matching_keys)
    except RedisError as e:
        logger.warning(f"Redis error invalidating pattern {pattern}: {e}")
        return 0

def cached(namespace: str, ttl: int = DEFAULT_CACHE_TTL):
    """
    Decorator for caching function results.
    
    Args:
        namespace: The namespace for the cache key
        ttl: Time to live in seconds
        
    Returns:
        Decorated function
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if not DEFAULT_CACHE_ENABLED:
                return func(*args, **kwargs)
            
            # Generate a cache key
            cache_key = generate_cache_key(namespace, func.__name__, *args, **kwargs)
            
            # Check if the result is already cached
            cached_result = cache_get(cache_key)
            if cached_result is not None:
                logger.debug(f"Cache hit for {func.__name__} with key {cache_key}")
                return cached_result
            
            # If not cached, call the function
            logger.debug(f"Cache miss for {func.__name__} with key {cache_key}")
            result = func(*args, **kwargs)
            
            # Cache the result
            cache_set(cache_key, result, ttl=ttl)
            
            return result
        return wrapper
    return decorator

def clear_cache() -> bool:
    """
    Clear the entire cache.
    
    Returns:
        True if successful, False otherwise
    """
    if not DEFAULT_CACHE_ENABLED:
        return False
    
    client = get_redis_client()
    try:
        # Delete all keys with our prefix
        pattern = f"{DEFAULT_CACHE_PREFIX}*"
        keys = list(client.scan_iter(match=pattern))
        
        if keys:
            client.delete(*keys)
            
        logger.info(f"Cleared {len(keys)} keys from cache")
        return True
    except RedisError as e:
        logger.warning(f"Redis error clearing cache: {e}")
        return False

def get_cache_stats() -> Dict[str, Any]:
    """
    Get cache statistics.
    
    Returns:
        Dictionary with cache statistics
    """
    if not DEFAULT_CACHE_ENABLED:
        return {
            'enabled': False,
            'keys': 0,
            'memory_used': 0,
            'hit_rate': 0,
            'miss_rate': 0
        }
    
    client = get_redis_client()
    try:
        # Get all keys with our prefix
        pattern = f"{DEFAULT_CACHE_PREFIX}*"
        keys = list(client.scan_iter(match=pattern))
        
        # Get info from Redis
        info = client.info()
        
        # Return statistics
        return {
            'enabled': True,
            'keys': len(keys),
            'memory_used': info.get('used_memory_human', 'unknown'),
            'hit_rate': info.get('keyspace_hits', 0),
            'miss_rate': info.get('keyspace_misses', 0)
        }
    except RedisError as e:
        logger.warning(f"Redis error getting cache stats: {e}")
        return {
            'enabled': False,
            'error': str(e)
        }

# Cache decorator for specific database operations
def cache_query(ttl: int = DEFAULT_CACHE_TTL):
    """
    Decorator for caching database query results.
    
    Args:
        ttl: Time to live in seconds
        
    Returns:
        Decorated function
    """
    return cached('query', ttl=ttl)

def cache_molecule(ttl: int = DEFAULT_CACHE_TTL):
    """
    Decorator for caching molecule data.
    
    Args:
        ttl: Time to live in seconds
        
    Returns:
        Decorated function
    """
    return cached('molecule', ttl=ttl)

def cache_property(ttl: int = DEFAULT_CACHE_TTL):
    """
    Decorator for caching property data.
    
    Args:
        ttl: Time to live in seconds
        
    Returns:
        Decorated function
    """
    return cached('property', ttl=ttl)

# Invalidation functions for specific types
def invalidate_molecule(molecule_id: str) -> int:
    """
    Invalidate all cache keys related to a molecule.
    
    Args:
        molecule_id: The molecule ID
        
    Returns:
        Number of keys invalidated
    """
    pattern = f"molecule:*{molecule_id}*"
    count = cache_invalidate_pattern(pattern)
    
    # Also invalidate query cache that might include this molecule
    query_pattern = f"query:*"
    count += cache_invalidate_pattern(query_pattern)
    
    return count

def invalidate_property(property_type_id: str = None) -> int:
    """
    Invalidate all cache keys related to properties.
    
    Args:
        property_type_id: The property type ID (optional)
        
    Returns:
        Number of keys invalidated
    """
    if property_type_id:
        pattern = f"property:*{property_type_id}*"
    else:
        pattern = f"property:*"
        
    count = cache_invalidate_pattern(pattern)
    
    # Also invalidate query cache that might include properties
    query_pattern = f"query:*"
    count += cache_invalidate_pattern(query_pattern)
    
    return count