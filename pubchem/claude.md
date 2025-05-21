# PubChem Integration Module Reference

## Overview

The PubChem module provides a comprehensive client for integrating with the PubChem database of chemical substances, supporting large-scale data retrieval operations with advanced resilience features.

## Key Components

### Client Architecture
- **Core Client**: `ResilientPubChemClient` in `client.py`
- **Simplified Client**: `PubChemClient` in `simple_client.py` for basic operations
- **Rate Limiting**: Advanced adaptive rate limiting in `rate_limiter.py`
- **Caching**: SQLite-backed persistent cache in `cache.py`
- **Fallback Handling**: RDKit-based property calculation fallback in `rdkit_fallback.py`

### Resilience Features
- **Adaptive Rate Limiting**: Dynamically adjusts request frequency
- **Exponential Backoff**: Increases wait times after failures
- **Request Batching**: Optimizes multiple requests
- **Cache Management**: Time-based cache invalidation
- **Error Recovery**: Multiple retry strategies

## Worker Pool System

The `worker_pool.py` module implements a robust parallel processing system:

- **Concurrent Workers**: Process multiple compounds simultaneously
- **Work Queue**: Priority-based task scheduling
- **Result Collection**: Aggregated result handling
- **Error Management**: Graceful handling of worker failures
- **Progress Tracking**: Real-time operation monitoring

## Data Integration Process

1. **CID Resolution**: Convert identifiers to PubChem CIDs
2. **Basic Properties**: Retrieve core compound properties
3. **Detailed Information**: Fetch extended compound data
4. **Synonym Collection**: Gather alternative names
5. **Structure Data**: Retrieve structural representations
6. **Property Calculation**: Compute or retrieve chemical properties

## Scheduler System

The `scheduler.py` provides intelligent task scheduling:

- **Rate-Aware Scheduling**: Respects PubChem's rate limits
- **Priority Queuing**: Critical tasks processed first
- **Batched Execution**: Optimizes similar requests
- **Throttling**: Prevents overwhelming the API
- **Fairness**: Ensures all operations make progress

## Configuration Options

- **Retry Settings**: Customizable retry attempts and strategies
- **Cache Configuration**: Cache lifetime and size parameters
- **Rate Limit Parameters**: Fine-tuning of request rates
- **Timeout Values**: Connection and read timeout adjustments
- **Batch Sizes**: Control of request batching

## Best Practices

1. **Use CIDs**: Always work with PubChem CIDs when possible
2. **Batch Similar Requests**: Group similar operations
3. **Enable Caching**: Always use the cache for repeated operations
4. **Implement Checkpointing**: For large operations
5. **Use Worker Pools**: For parallel processing needs
6. **Handle Partial Data**: Not all compounds have complete information

## Common Pitfalls

1. **Rate Limit Blocks**: PubChem can block aggressive requests
2. **Incomplete Data**: Some compounds lack certain properties
3. **Large Result Sets**: Some operations return massive datasets
4. **Timeout Issues**: Complex queries may exceed default timeouts
5. **Identity Resolution**: Name/SMILES to CID mapping can be ambiguous
6. **API Changes**: PubChem occasionally modifies response formats

## Usage Examples

```python
# Basic client usage
from pubchem import ResilientPubChemClient

client = ResilientPubChemClient()
compound = client.get_compound_by_cid(2244)

# Batch operation
compounds = client.get_compounds_by_cids([2244, 3672])

# With caching enabled
client.enable_caching()
# First call hits API
aspirin = client.get_compound_by_name("aspirin")
# Second call uses cache
aspirin_again = client.get_compound_by_name("aspirin")

# Using the worker pool for parallel operations
from pubchem.worker_pool import WorkerPool

pool = WorkerPool(workers=5)
for cid in range(1, 100):
    pool.submit_task(client.get_compound_by_cid, cid)
results = pool.get_all_results()

# Simplified client for basic needs
from pubchem.simple_client import PubChemClient

simple_client = PubChemClient()
aspirin_cid = simple_client.get_cid_by_name("aspirin")
```