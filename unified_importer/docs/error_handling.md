# Error Handling System Documentation

## Overview

The CryoProtect Unified Importer Error Handling System provides a comprehensive framework for managing errors in a consistent, recoverable, and informative way. The system classifies errors, implements intelligent recovery strategies, and provides detailed context for debugging and analysis.

## Key Components

### Error Classification

Errors are classified along multiple dimensions:

#### Error Categories (`ErrorCategory` enum)

- **NETWORK**: Connection issues, DNS failures, etc.
- **DATABASE**: Database connectivity or query problems
- **VALIDATION**: Data validation failures
- **RESOURCE**: Resource unavailability (memory, disk space, etc.)
- **AUTHORIZATION**: Permission and authentication errors
- **CONFIGURATION**: Configuration errors, missing settings
- **EXTERNAL_SERVICE**: Failures in external services or APIs
- **TIMEOUT**: Operations that exceed time limits
- **DATA_FORMAT**: Issues with data structure or format
- **UNKNOWN**: Unclassified errors

#### Error Severity (`ErrorSeverity` enum)

- **CRITICAL**: System-wide failures requiring immediate attention
- **HIGH**: Severe issues affecting major functionality
- **MEDIUM**: Significant issues affecting some functionality
- **LOW**: Minor issues with limited impact
- **INFO**: Informational issues with negligible impact

#### Recovery Strategies (`RecoveryStrategy` enum)

- **RETRY**: Attempt the operation again immediately
- **FALLBACK**: Use an alternative implementation or data source
- **SKIP**: Skip the problematic item and continue
- **ABORT**: Abort the current operation entirely
- **LOG_ONLY**: Record the error but continue normally
- **DELAYED_RETRY**: Retry after a delay with exponential backoff
- **CIRCUIT_BREAKER**: Prevent cascading failures by stopping attempts temporarily

### Core Classes

#### `ErrorContext`

Captures detailed information about the context in which an error occurred:

```python
context = ErrorContext(
    component="ChEMBL Importer",
    operation="fetch_molecule",
    data={"molecule_id": "CHEMBL123"},
    metadata={"attempt": 2, "timestamp": time.time()}
)
```

#### `ClassifiedError`

Combines an exception with its classification and context:

```python
classified_error = ClassifiedError(
    error=error,
    category=ErrorCategory.VALIDATION,
    severity=ErrorSeverity.MEDIUM,
    recovery_strategy=RecoveryStrategy.SKIP,
    context=context
)
```

#### `ErrorClassifier`

Analyzes exceptions and assigns categories, severity levels, and recovery strategies:

```python
classifier = ErrorClassifier(logger=logger)
classified = classifier.classify(error, context)
```

The classifier includes built-in rules for common error types and can be extended with custom rules:

```python
classifier.register_rule(
    error_class=CustomDatabaseError,
    category=ErrorCategory.DATABASE,
    severity=ErrorSeverity.HIGH,
    recovery_strategy=RecoveryStrategy.RETRY
)
```

#### `RetryManager`

Implements retry logic with exponential backoff and jitter:

```python
retry_manager = RetryManager(
    max_attempts=3,
    initial_delay=0.1,
    max_delay=10.0,
    backoff_factor=2,
    logger=logger
)

result = retry_manager.execute(
    function_to_retry,
    component="Component Name",
    operation="Operation Name"
)
```

The RetryManager supports:
- Configurable retry attempts
- Exponential backoff with jitter for distributed systems
- Non-retryable error types
- Detailed logging of retry attempts
- Circuit breaker pattern for persistent failures

#### `ValidationErrorHandler`

Manages validation errors with customizable strategies:

```python
handler = ValidationErrorHandler(logger=logger)

result = handler.handle_validation(
    validation_function, input_data,
    error_strategy=RecoveryStrategy.FALLBACK,
    component="Data Validator",
    operation="validate_structure",
    fallback_func=lambda x: default_conversion(x),
    default_value=None
)
```

The ValidationErrorHandler supports multiple error strategies:
- SKIP: Return a default value and continue
- FALLBACK: Use a fallback function to transform invalid data
- ABORT: Raise the validation error
- LOG_ONLY: Log the error but return the original value

#### `ErrorManager`

Provides a unified interface for all error handling functionality:

```python
error_manager = ErrorManager(logger=logger)

# Classify an error
classified = error_manager.classify_error(error, context)

# Retry a function
result = error_manager.retry(
    function, component="Component", operation="Operation",
    max_attempts=5, initial_delay=0.1
)

# Handle validation
validated = error_manager.validate(
    validation_function, data,
    component="Validator", operation="validate_data",
    error_strategy=RecoveryStrategy.SKIP, default_value=None
)

# Process items in bulk with error handling
results = error_manager.process_bulk(
    process_function, items,
    component="Bulk Processor", operation="process_items",
    error_strategy=RecoveryStrategy.SKIP,
    error_callback=lambda err, item: record_error(err, item)
)
```

## Circuit Breaker Pattern

The error handling system implements the Circuit Breaker pattern to prevent cascading failures. When a component fails repeatedly:

1. The circuit "opens" after a threshold of failures
2. Subsequent requests fail fast without attempting the operation
3. After a timeout, the circuit enters a "half-open" state
4. A single success restores the circuit to "closed"
5. Further failures immediately open the circuit again

```python
error_manager.configure_circuit_breaker(
    failure_threshold=5,     # Number of failures before opening
    reset_timeout=60.0,      # Seconds before attempting half-open
    half_open_max_calls=1    # Max calls in half-open state
)
```

## Bulk Operation Error Handling

The system provides specialized support for handling errors in bulk operations:

```python
results = error_manager.process_bulk(
    process_function,   # Function to process each item
    items,              # List of items to process
    component="Bulk Processor",
    operation="process_batch",
    error_strategy=RecoveryStrategy.SKIP,
    parallel=True,      # Process in parallel
    max_workers=4,      # Number of parallel workers
    error_callback=lambda err, item: record_error(err, item)
)
```

The bulk processor supports:
- Parallel processing of items
- Consistent error handling across all items
- Error callbacks for custom error recording
- Different error strategies per operation

## Best Practices

### 1. Use Context Consistently

Always provide detailed context information for errors:

```python
context = ErrorContext(
    component="ChEMBL Importer",
    operation="fetch_molecule",
    data={"molecule_id": id, "attempt": attempt},
    metadata={"source": source_name, "timestamp": time.time()}
)
```

### 2. Classify Errors Appropriately

Register custom rules for application-specific exceptions:

```python
classifier.register_rule(
    error_class=MoleculeFormatError,
    category=ErrorCategory.DATA_FORMAT,
    severity=ErrorSeverity.MEDIUM,
    recovery_strategy=RecoveryStrategy.FALLBACK
)
```

### 3. Use Appropriate Recovery Strategies

Choose recovery strategies based on error type and context:

- **RETRY**: For transient network or database errors
- **FALLBACK**: When alternative data sources are available
- **SKIP**: For individual data item validation errors
- **ABORT**: For critical system configuration errors
- **CIRCUIT_BREAKER**: For dependent service failures

### 4. Configure Retries Appropriately

Adjust retry parameters based on the operation:

```python
# For quick operations with potentially transient failures
retry_manager.execute(quick_operation, max_attempts=5, initial_delay=0.1)

# For slow external API calls
retry_manager.execute(api_operation, max_attempts=3, initial_delay=1.0)
```

### 5. Implement Validation Error Handling

Use appropriate strategies for validation errors:

```python
# For critical data integrity
handler.handle_validation(validate_id, id, 
                          error_strategy=RecoveryStrategy.ABORT)

# For non-critical fields
handler.handle_validation(validate_description, description,
                          error_strategy=RecoveryStrategy.SKIP,
                          default_value="")
```

## Integration with Logging System

The error handling system integrates with the application logging system:

```python
logger = logging.getLogger("unified_importer")
error_manager = ErrorManager(logger=logger)
```

Errors are logged with detailed context information:

```
ERROR [ChEMBL Importer] Error during fetch_molecule: ConnectionError: 
Failed to connect to ChEMBL API (attempt 2/5)
Context: {'molecule_id': 'CHEMBL123', 'source': 'ChEMBL', 'timestamp': 1620000000.0}
Recovery strategy: RETRY with backoff delay of 0.4s
```

## Integration with Monitoring System

The error handling system collects error metrics for monitoring:

```python
error_metrics = error_manager.get_error_metrics()

# Example metrics:
# {
#   'error_counts_by_category': {'NETWORK': 12, 'DATABASE': 3, ...},
#   'error_counts_by_component': {'ChEMBL Importer': 8, 'PubChem Importer': 7, ...},
#   'recovery_success_rate': 0.85,
#   'retry_success_counts': {'attempt_1': 15, 'attempt_2': 8, 'attempt_3': 3, ...},
#   'circuit_breaker_open_count': 2,
#   'avg_recovery_time': 0.385
# }
```

These metrics can be exported to monitoring systems for alerting and dashboards.

## Example Usage

### Basic Error Handling

```python
from unified_importer.core.error_handling import ErrorManager, RecoveryStrategy, ErrorContext

error_manager = ErrorManager(logger=logger)

try:
    result = process_data(data)
except Exception as e:
    context = ErrorContext(
        component="Data Processor",
        operation="process_data",
        data={"data_id": data.id}
    )
    classified_error = error_manager.classify_error(e, context)
    
    if classified_error.recovery_strategy == RecoveryStrategy.RETRY:
        result = error_manager.retry(
            lambda: process_data(data),
            component="Data Processor",
            operation="process_data"
        )
    elif classified_error.recovery_strategy == RecoveryStrategy.SKIP:
        result = default_result
        logger.warning(f"Skipped processing data {data.id}: {e}")
    else:
        raise classified_error.error
```

### Simplified Error Handling with ErrorManager

```python
from unified_importer.core.error_handling import ErrorManager, RecoveryStrategy

error_manager = ErrorManager(logger=logger)

# With automatic retry for appropriate errors
result = error_manager.execute_with_handling(
    lambda: process_data(data),
    component="Data Processor",
    operation="process_data",
    default_value=default_result
)

# For validation operations
validated_data = error_manager.validate(
    validate_molecule_data, molecule_data,
    component="Molecule Validator",
    operation="validate_molecule",
    error_strategy=RecoveryStrategy.SKIP,
    default_value=None
)

# For bulk operations
results = error_manager.process_bulk(
    process_molecule, molecules,
    component="Molecule Processor",
    operation="process_molecules",
    error_strategy=RecoveryStrategy.SKIP,
    parallel=True,
    max_workers=4
)
```

## Conclusion

The CryoProtect Unified Importer Error Handling System provides a comprehensive approach to managing errors, implementing recovery strategies, and collecting detailed diagnostics. By using this system consistently throughout the application, we achieve more robust error handling, better user experiences, and more maintainable code.

For implementation details, see the `unified_importer/core/error_handling.py` module and the corresponding unit tests in `unified_importer/tests/test_error_handling.py`.