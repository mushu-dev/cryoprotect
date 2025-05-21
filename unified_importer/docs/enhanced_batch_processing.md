# Enhanced Batch Processing

The enhanced batch processing system optimizes the processing of large datasets in the unified molecular importer. It provides adaptive batch sizing, parallel processing capabilities, and detailed metrics to improve performance and reliability.

## Features

- **Adaptive Batch Sizing**: Automatically adjusts batch sizes based on processing times, memory usage, and error rates
- **Parallel Processing**: Efficiently processes batches in parallel using multiple worker threads
- **Resource Monitoring**: Tracks memory usage and processing times to optimize resource utilization
- **Error Handling**: Gracefully handles errors in batch processing without failing the entire operation
- **Chemical Data Optimizations**: Specialized strategies for chemical data sources like ChEMBL and PubChem
- **Asynchronous Support**: Full async/await support for efficient I/O handling

## Architecture

The enhanced batch processing system consists of several key components:

1. **AdaptiveBatchProcessor**: Core class that manages batch processing with adaptive sizing
2. **BatchProcessingStrategy**: Strategy pattern implementation for batch size determination
3. **ChemicalDataBatchStrategy**: Specialized strategy for chemical data processing
4. **BatchProcessingManager**: Integration layer that connects batch processing with the importer
5. **Batch Transformation Functions**: Specialized functions for molecule and property transformations

## Usage

### Basic Usage

```python
from unified_importer.core.enhanced_batch_processing import AdaptiveBatchProcessor

# Define a function to process batches
def process_batch(items):
    results = []
    for item in items:
        # Process item
        results.append(processed_item)
    return results

# Create a processor
processor = AdaptiveBatchProcessor(
    process_func=process_batch,
    initial_batch_size=100,
    min_batch_size=10,
    max_batch_size=500,
    max_workers=4
)

# Process items
results = processor.process_items_parallel(items)
```

### Using the BatchProcessingManager

```python
from unified_importer.core.batch_integration import BatchProcessingManager

# Create a manager with database integration
manager = BatchProcessingManager(
    db=database_operations,
    max_workers=4,
    adaptive_batch_sizing=True,
    parallel_processing=True
)

# Process items
results = manager.process_items(
    items=molecules,
    operation="transform_molecules",
    data_source="chembl",
    process_func=molecule_transform_function,
    initial_batch_size=50
)

# Get processing statistics
stats = manager.get_stats()
print(f"Processed {stats['total_items_processed']} items in {stats['total_duration']:.2f} seconds")
```

### Asynchronous Processing

```python
async def process_data():
    # Create async processor
    processor = AdaptiveBatchProcessor(
        process_func=process_batch,
        initial_batch_size=100
    )
    
    # Process items asynchronously
    results = await processor.process_items_async(items)
    
    # Or with parallel processing
    results = await processor.process_items_parallel_async(items)
```

## Configuration Options

### AdaptiveBatchProcessor

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| process_func | Callable | - | Function to process a batch of items |
| initial_batch_size | int | 100 | Initial size of each batch |
| min_batch_size | int | 10 | Minimum batch size |
| max_batch_size | int | 1000 | Maximum batch size |
| max_workers | int | CPU count | Maximum worker threads/processes |
| strategy | BatchProcessingStrategy | None | Batch processing strategy |
| logger | Logger | None | Logger instance |

### BatchProcessingStrategy

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| initial_batch_size | int | 100 | Initial size of each batch |
| min_batch_size | int | 10 | Minimum batch size |
| max_batch_size | int | 1000 | Maximum batch size |
| target_duration | float | 5.0 | Target duration for batch processing in seconds |
| memory_limit | int | 1GB | Memory usage limit in bytes |
| logger | Logger | None | Logger instance |

### ChemicalDataBatchStrategy

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| data_source | str | "generic" | Source of the chemical data ("chembl", "pubchem", "generic") |
| initial_batch_size | int | 100 | Initial size of each batch |
| min_batch_size | int | 10 | Minimum batch size |
| max_batch_size | int | 1000 | Maximum batch size |
| target_duration | float | 5.0 | Target duration for batch processing in seconds |
| memory_limit | int | 1GB | Memory usage limit in bytes |
| logger | Logger | None | Logger instance |

### BatchProcessingManager

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| db | EnhancedDatabaseOperations | - | Database operations instance |
| max_workers | int | CPU count | Maximum worker threads/processes |
| adaptive_batch_sizing | bool | True | Whether to use adaptive batch sizing |
| parallel_processing | bool | True | Whether to use parallel processing |
| logger | Logger | None | Logger instance |

## Performance Considerations

### Optimal Batch Size

The optimal batch size depends on several factors:

- **Memory Usage**: Larger batch sizes consume more memory
- **Processing Time**: Larger batches may improve throughput but increase latency
- **Error Rate**: Smaller batches limit the impact of errors
- **Data Complexity**: Complex molecules may require smaller batches

The adaptive batch sizing system automatically adjusts based on these factors, but you can provide initial guidance through configuration.

### Parallel Processing

Parallel processing can significantly improve throughput, but consider these factors:

- **CPU Bound vs. I/O Bound**: Use more workers for I/O bound tasks, fewer for CPU bound
- **Resource Contention**: Too many workers can lead to diminishing returns
- **Memory Usage**: Each worker requires additional memory
- **Database Connection Limits**: Consider database connection limits when using parallel processing

## Monitoring and Statistics

The batch processing system provides detailed statistics to help monitor performance:

```python
stats = processor.get_stats()
```

Key metrics include:

- **Total Items Processed**: Number of items processed
- **Batch Count**: Number of batches processed
- **Processing Rate**: Items processed per second
- **Error Rate**: Percentage of items that failed processing
- **Memory Usage**: Memory consumed during processing
- **Batch Size Evolution**: How batch sizes adapted over time

## Integration with Database Operations

The batch processing system integrates with the enhanced database operations:

```python
async def batch_db_insert(db, table, records):
    async with db.transaction_async() as tx_id:
        inserted_ids = await db.batch_insert(
            table=table,
            records=records,
            transaction_id=tx_id
        )
        return inserted_ids
```

This integration provides:

- **Transactional Safety**: Ensures all records in a batch are inserted atomically
- **Connection Pooling**: Efficiently uses database connections
- **Error Handling**: Gracefully handles database errors
- **Performance Metrics**: Tracks database operation performance

## Best Practices

1. **Start with Conservative Batch Sizes**: Begin with smaller batch sizes (50-100) and let the system adapt
2. **Monitor Memory Usage**: Watch for excessive memory consumption, especially with complex molecules
3. **Balance Workers**: Set `max_workers` based on your system capabilities and workload characteristics
4. **Use Appropriate Strategies**: Use `ChemicalDataBatchStrategy` for chemical data sources
5. **Handle Errors Gracefully**: Ensure your processing functions handle errors without crashing
6. **Transaction Management**: Use transactions for database operations to ensure data integrity
7. **Regular Monitoring**: Periodically check processing statistics to identify performance issues