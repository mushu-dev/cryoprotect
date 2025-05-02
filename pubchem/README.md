# ResilientPubChemClient

A production-ready client for the PubChem API with advanced resilience features for CryoProtect v2.

## Features

- **Adaptive Rate Limiting**: Automatically adjusts request rates based on weekdays vs weekends to respect PubChem API limits
- **Multi-level Caching**: Combines in-memory and disk-based caching for optimal performance
- **Exponential Backoff Retry Logic**: Intelligently retries failed requests with increasing delays
- **Circuit Breaker Pattern**: Prevents cascading failures by temporarily disabling API calls after repeated failures
- **Fallback Mechanism**: Uses cached data when the API is unavailable
- **Weekend Job Scheduler**: Runs bulk operations during weekends when API limits are less restrictive

## Installation

The package is part of the CryoProtect v2 project and doesn't require separate installation.

## Usage

### Basic Usage

```python
from pubchem import ResilientPubChemClient

# Create a client with default settings
client = ResilientPubChemClient()

# Get molecule properties
properties = client.get_molecule_properties(2244)
print(properties)

# Get compound synonyms
synonyms = client.get_compound_synonyms(2244)
print(synonyms["Synonyms"])

# Search for compounds
results = client.search_compounds("glucose")
print(results["CIDs"])
```

### Advanced Usage

```python
# Create a client with custom settings
client = ResilientPubChemClient(
    cache_dir="custom_cache",
    weekday_requests_per_second=1.0,  # Slower on weekdays
    weekend_requests_per_second=5.0,  # Faster on weekends
    max_retries=3,
    failure_threshold=5,
    recovery_timeout=60,
    cache_ttl=86400 * 7,  # 7 days
    memory_cache_size=500,
    enable_scheduler=True
)

# Schedule a job to prefetch molecule properties
cids = [2244, 5793, 5288, 702]
job_id = client.prefetch_molecule_properties(cids)

# Check job status
status = client.get_job_status(job_id)
print(status)

# Get statistics
cache_stats = client.get_cache_stats()
rate_limiter_stats = client.get_rate_limiter_stats()
circuit_breaker_stats = client.get_circuit_breaker_stats()
```

## Components

### ResilientPubChemClient

The main client class that integrates all resilience features.

### PubChemCache

Multi-level caching system with in-memory LRU cache and disk-based persistent storage.

### AdaptiveRateLimiter

Rate limiter that adjusts request rates based on the day of the week.

### WeekendJobScheduler

Scheduler for running bulk operations during weekends when API limits are less restrictive.

### CircuitBreaker

Circuit breaker implementation to prevent repeated failures.

## Testing

The package includes comprehensive unit tests for all components:

```bash
# Run all tests
python -m unittest discover pubchem

# Run specific test file
python -m unittest pubchem.test_client
```

## Error Handling

The client includes robust error handling:

- API errors are caught and logged
- Circuit breaker prevents repeated failures
- Fallback to cached data when the API is unavailable
- Exponential backoff for transient errors

## Performance Considerations

- Use the cache to minimize API calls
- Schedule bulk operations for weekends
- Monitor rate limiter and circuit breaker statistics
- Adjust cache TTL based on data freshness requirements