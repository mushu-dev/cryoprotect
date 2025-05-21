# Caching Mechanism

The caching mechanism provides a way to store and retrieve data efficiently to improve the performance of the unified molecular importer. It reduces redundant API calls, database operations, and computational tasks by storing the results of expensive operations for later reuse.

## Features

- **Multiple Cache Backends**: In-memory and disk-based persistence
- **Smart Cache Policies**: Different caching policies based on data source and type
- **Automatic Expiration**: Time-based expiration with configurable TTLs (Time-To-Live)
- **Cache Statistics**: Comprehensive statistics for monitoring and optimization
- **Resource Management**: Size limits, item count limits, and LRU (Least Recently Used) eviction
- **Asynchronous Support**: Full async/await support for non-blocking operations
- **Thread Safety**: All operations are thread-safe for concurrent access

## Architecture

The caching system consists of several key components:

1. **CacheEntry**: Data structure for storing cache entries with metadata
2. **CacheBackend**: Abstract interface for cache storage backends
3. **MemoryCache**: In-memory implementation of the cache backend
4. **DiskCache**: Persistent disk-based implementation of the cache backend
5. **CacheManager**: Unified interface for working with multiple cache backends
6. **AsyncCacheManager**: Asynchronous wrapper for the cache manager

## Usage

### Basic Usage

```python
from unified_importer.core.caching import CacheManager

# Create a cache manager
cache_manager = CacheManager()

# Store a value with default policy (memory cache, 1 hour TTL)
cache_manager.put("example:key", {"data": "value"})

# Store a value with PubChem policy (disk cache, 7 days TTL)
cache_manager.put("pubchem:cid:12345", compound_data, policy="pubchem", source="pubchem")

# Retrieve a value
value = cache_manager.get("example:key")
compound = cache_manager.get("pubchem:cid:12345", policy="pubchem")

# Delete a value
cache_manager.delete("example:key")

# Clear the cache
cache_manager.clear()

# Get cache statistics
stats = cache_manager.get_stats()
```

### Asynchronous Usage

```python
from unified_importer.core.caching import AsyncCacheManager

# Create an async cache manager
cache_manager = AsyncCacheManager()

async def example():
    # Store a value
    await cache_manager.put("example:key", {"data": "value"})
    
    # Retrieve a value
    value = await cache_manager.get("example:key")
    
    # Delete a value
    await cache_manager.delete("example:key")
    
    # Clear the cache
    await cache_manager.clear()
    
    # Get cache statistics
    stats = await cache_manager.get_stats()
```

### Creating Cache Keys

```python
from unified_importer.core.caching import create_cache_key

# Create a cache key from multiple parts
key = create_cache_key("pubchem", ["compound", 12345, "properties"])
# Result: "pubchem:compound:12345:properties"

# Long keys will be automatically hashed
key = create_cache_key("query", ["very_long_search_term_with_many_parameters"], max_length=50)
```

## Cache Policies

The caching system supports different policies for different types of data:

| Policy | Backend | TTL | Use Case |
|--------|---------|-----|----------|
| default | memory | 1 hour | General-purpose, temporary data |
| pubchem | disk | 7 days | PubChem API results, which change infrequently |
| chembl | disk | 30 days | ChEMBL API results, which change very infrequently |
| transient | memory | 5 minutes | Short-lived, frequently changing data |

## Configuration Options

### MemoryCache

| Parameter | Default | Description |
|-----------|---------|-------------|
| max_size_bytes | 100MB | Maximum size of the cache in bytes |
| max_items | 10000 | Maximum number of items in the cache |
| default_ttl | 3600.0 (1 hour) | Default time-to-live in seconds |
| cleanup_interval | 60.0 (1 minute) | Interval for cleaning up expired entries |

### DiskCache

| Parameter | Default | Description |
|-----------|---------|-------------|
| cache_dir | `.cache` | Directory for cache files |
| max_size_bytes | 1GB | Maximum size of the cache in bytes |
| max_items | 100000 | Maximum number of items in the cache |
| default_ttl | 86400.0 (1 day) | Default time-to-live in seconds |
| cleanup_interval | 300.0 (5 minutes) | Interval for cleaning up expired entries |

## Integration with Molecular Importer

The caching mechanism is integrated with the unified molecular importer in several key areas:

### 1. API Request Caching

When fetching data from external APIs like PubChem and ChEMBL, the importer caches the responses to reduce the number of API calls:

```python
# Example of caching API requests
async def fetch_compound(self, identifier: str) -> Optional[Dict[str, Any]]:
    # Create cache key
    cache_key = create_cache_key("pubchem", ["compound", identifier])
    
    # Try to get from cache
    cached_data = await self.cache_manager.get(cache_key, policy="pubchem")
    if cached_data is not None:
        return cached_data
    
    # Fetch from API
    response = await self.api_client.get_compound(identifier)
    
    # Cache the result
    await self.cache_manager.put(
        cache_key,
        response,
        policy="pubchem",
        source="pubchem"
    )
    
    return response
```

### 2. Transformation Caching

Expensive molecule transformations and calculations are cached to avoid redundant computation:

```python
# Example of caching transformations
async def transform_molecule(self, molecule_data: Dict[str, Any]) -> Dict[str, Any]:
    # Create cache key from molecule identifier
    molecule_id = molecule_data.get("id") or molecule_data.get("cid")
    cache_key = create_cache_key("transform", ["molecule", molecule_id])
    
    # Try to get from cache
    cached_data = await self.cache_manager.get(cache_key)
    if cached_data is not None:
        return cached_data
    
    # Perform transformation
    transformed = await self._perform_transformation(molecule_data)
    
    # Cache the result
    await self.cache_manager.put(cache_key, transformed)
    
    return transformed
```

### 3. Search Results Caching

Search results from APIs are cached to improve performance of repeated searches:

```python
# Example of caching search results
async def search_compounds(self, query: str, max_results: Optional[int] = None) -> List[str]:
    # Create cache key
    cache_key = create_cache_key("search", ["pubchem", query, str(max_results)])
    
    # Try to get from cache
    cached_results = await self.cache_manager.get(cache_key, policy="transient")
    if cached_results is not None:
        return cached_results
    
    # Perform search
    results = await self._perform_search(query, max_results)
    
    # Cache the results
    await self.cache_manager.put(
        cache_key,
        results,
        policy="transient"
    )
    
    return results
```

## Monitoring and Statistics

The caching system provides detailed statistics for monitoring and optimization:

```python
# Get cache statistics
stats = cache_manager.get_stats()

# Example output
{
    "memory": {
        "hits": 1250,
        "misses": 350,
        "hit_ratio": 0.78,
        "insertions": 1600,
        "evictions": 200,
        "expirations": 150,
        "size_bytes": 2500000,
        "items": 1250,
        "last_cleanup": 1621234567.89
    },
    "disk": {
        "hits": 9500,
        "misses": 500,
        "hit_ratio": 0.95,
        "insertions": 10000,
        "evictions": 0,
        "expirations": 500,
        "size_bytes": 50000000,
        "items": 9500,
        "last_cleanup": 1621234567.89,
        "expired_items": 0,
        "avg_hits": 2.5
    },
    "policies": {
        "default": { "backend": "memory", "ttl": 3600.0 },
        "pubchem": { "backend": "disk", "ttl": 604800.0 },
        "chembl": { "backend": "disk", "ttl": 2592000.0 },
        "transient": { "backend": "memory", "ttl": 300.0 }
    }
}
```

## Best Practices

1. **Use Appropriate Policies**: Choose the right cache policy based on the data source and expected change frequency.
2. **Create Specific Cache Keys**: Use the `create_cache_key` function with a prefix and specific parts to ensure unique keys.
3. **Set Appropriate TTLs**: Choose TTLs based on how frequently the data changes.
4. **Monitor Statistics**: Regularly check cache statistics to identify performance issues and optimization opportunities.
5. **Handle Cache Misses Gracefully**: Always be prepared to fetch the data from the original source when it's not in the cache.
6. **Use Transactions with Database Caching**: When caching database operations, ensure consistency by using transactions.
7. **Consider Memory Constraints**: Set appropriate size limits to prevent excessive memory usage.
8. **Use Asynchronous API**: When working with async code, use the `AsyncCacheManager` to avoid blocking the event loop.

## Limitations and Considerations

1. **Disk Cache Overhead**: Disk caching involves file I/O operations, which can be slower than memory caching.
2. **Memory Usage**: Memory caching can consume significant RAM for large datasets.
3. **Cache Invalidation**: Determining when to invalidate cache entries can be challenging.
4. **Stale Data**: Cached data can become stale if the source data changes before the cache entry expires.
5. **Thread Safety**: While the cache is thread-safe, concurrent operations can still lead to race conditions in the application code.

## Implementation Details

### Cache Entry Serialization

Cache entries are serialized using Python's `pickle` module, which allows for efficient storage and retrieval of complex data structures. However, this comes with some limitations:

- Pickled objects may not be compatible across different Python versions
- Custom classes may not be unpicklable unless they're available in the import path
- Some objects (like file handles or database connections) cannot be pickled

### LRU Eviction

Both the memory and disk caches use an LRU (Least Recently Used) eviction policy. When the cache reaches its size or item limit, the least recently accessed items are removed first. This is implemented by tracking the last access time of each cache entry.

### Thread Safety

All cache operations are protected by locks to ensure thread safety. The memory cache uses a `threading.RLock` for synchronization, while the disk cache uses SQLite's built-in concurrency controls.

### Automatic Cleanup

Both cache backends run background threads that periodically clean up expired entries. The cleanup interval is configurable and defaults to 1 minute for the memory cache and 5 minutes for the disk cache.