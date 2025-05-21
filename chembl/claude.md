# ChEMBL Integration Module Reference

## Overview

The ChEMBL module provides a robust client for integrating with the ChEMBL public database of bioactive molecules, featuring advanced functionality for reliable data retrieval in a production environment.

## Key Components

### Client Architecture
- **Core Client**: `ResilientChEMBLClient` in `client.py`
- **Rate Limiting**: Adaptive rate limiting with backoff in `rate_limiter.py`
- **Caching**: Multi-level caching system in `cache.py`
- **Error Handling**: Comprehensive error management in `error_handler.py`

### Resilience Features
- **Circuit Breaking**: Prevents cascading failures during API disruptions
- **Adaptive Backoff**: Dynamically adjusts request rate based on API response
- **Request Retry**: Intelligent retry logic for transient failures
- **Cache Fallback**: Uses cached data when API is unavailable
- **Batch Processing**: Optimized batch requests to minimize API calls

## Data Integration Workflow

1. **Search Terms Generation**: `search_terms.py` for optimized search strategies
2. **Initial Query**: Basic information retrieval with pagination
3. **Detail Enrichment**: Fetching detailed compound information
4. **Property Calculation**: Chemical property extraction and standardization
5. **Reference Registration**: Association with CryoProtect reference system
6. **Checkpoint Management**: Progress tracking for resumable operations

## Worker System

The module includes a parallel worker system (`worker.py`) that enables:

- **Concurrent Processing**: Parallel API requests within rate limits
- **Work Queue Management**: Prioritized task scheduling
- **Progress Tracking**: Real-time progress monitoring
- **Fault Isolation**: Worker failures don't affect entire system

## Configuration Options

- **Rate Limits**: Configurable requests per second/minute
- **Cache Size**: Adjustable memory and disk cache sizes
- **Retry Settings**: Customizable retry attempts and backoff
- **Timeout Values**: Connection and read timeout configuration
- **Circuit Breaker**: Failure threshold and recovery settings

## Best Practices

1. **Prefer Batch Operations**: Use batch queries over individual requests
2. **Leverage Caching**: Enable caching for repeated queries
3. **Monitor Rate Limits**: Stay below ChEMBL's published rate limits
4. **Use Checkpoints**: Always implement checkpoint-based resumption
5. **Handle Partial Results**: Some compounds may not have complete data

## Common Pitfalls

1. **API Changes**: ChEMBL API occasionally changes response formats
2. **Rate Limiting**: Aggressive querying can trigger IP blocks
3. **Incomplete Data**: Not all compounds have all properties
4. **Large Result Sets**: Some queries return extremely large datasets
5. **Timeout Issues**: Some complex queries take longer than default timeouts
6. **Reference Resolution**: ChEMBL to CID mapping can be ambiguous

## Usage Examples

```python
# Basic client usage
from chembl import ResilientChEMBLClient

client = ResilientChEMBLClient()
results = client.search_compounds("aspirin")

# Batch operation
compounds = client.get_compounds_batch(["CHEMBL25", "CHEMBL2"])

# With caching
client.enable_caching()
# First call hits API
compounds = client.get_compound_by_chembl_id("CHEMBL25")
# Second call uses cache
compounds = client.get_compound_by_chembl_id("CHEMBL25")

# Worker pool for parallel processing
from chembl.worker import ChEMBLWorkerPool

pool = ChEMBLWorkerPool(workers=5)
pool.submit_task(client.get_compound_by_chembl_id, "CHEMBL25")
pool.submit_task(client.get_compound_by_chembl_id, "CHEMBL2")
results = pool.get_results()
```