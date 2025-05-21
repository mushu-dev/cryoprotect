# Enhanced Retry Mechanism Documentation

## Overview

The Enhanced Retry Mechanism is a sophisticated error recovery system that provides robust handling of transient failures in the CryoProtect Unified Importer. It implements industry-standard reliability patterns like exponential backoff, jitter, and circuit breakers to enhance the system's resilience against common failures when interacting with external services and databases.

## Key Features

- **Exponential Backoff**: Gradually increases delay between retries to reduce load during service degradation
- **Jitter**: Adds randomness to retry delays to prevent retry storms in distributed systems
- **Circuit Breaker Pattern**: Prevents cascading failures by temporarily disabling operations after repeated failures
- **Custom Retry Policies**: Configurable policies for different error types and contexts
- **Persistent Circuit State**: Maintains circuit breaker state across application restarts
- **Comprehensive Metrics**: Collects detailed metrics on retry attempts and success rates
- **Decorators**: Easy-to-use decorators for adding retry capability to functions
- **Async Support**: Full support for asynchronous operations with `async/await`

## Core Components

### RetryConfig

Configuration for retry behavior with the following parameters:

```python
class RetryConfig:
    max_attempts: int = 3                  # Maximum number of retry attempts
    initial_delay: float = 0.1             # Initial delay before first retry (seconds)
    max_delay: float = 60.0                # Maximum delay between retries (seconds)
    backoff_factor: float = 2.0            # Exponential backoff multiplier
    jitter_factor: float = 0.1             # Random jitter factor (0.0-1.0) to add to delay
    non_retryable_errors: Set[type]        # Exception types that shouldn't be retried
    retryable_error_categories: Set[ErrorCategory]  # Error categories eligible for retry
```

### CircuitBreakerConfig

Configuration for the circuit breaker pattern:

```python
class CircuitBreakerConfig:
    failure_threshold: int = 5             # Number of failures before opening the circuit
    reset_timeout: float = 60.0            # Seconds before attempting half-open state
    half_open_max_calls: int = 1           # Maximum number of calls in half-open state
    success_threshold: int = 1             # Successes needed to close circuit from half-open
    open_timeout_increment_factor: float = 1.5  # Factor to increase timeout on repeated failures
    max_open_timeout: float = 300.0        # Maximum timeout in seconds (5 minutes)
```

### CircuitBreakerState

Tracks the state of a circuit breaker:

```python
class CircuitBreakerState:
    state: CircuitState                # Current circuit state (CLOSED, OPEN, HALF_OPEN)
    failure_count: int                 # Count of consecutive failures
    success_count: int                 # Count of consecutive successes (for half-open)
    last_failure_time: float           # Time of last failure
    open_until: float                  # Time when circuit will try half-open
    current_timeout: float             # Current timeout duration
    half_open_calls: int               # Number of calls in half-open state
    service_identifier: str            # Identifier for the service
    operation_identifier: str          # Identifier for the operation
```

### CircuitBreakerRegistry

Manages circuit breakers across the application:

```python
class CircuitBreakerRegistry:
    def get_breaker(self, service: str, operation: str) -> CircuitBreakerState
    def report_success(self, service: str, operation: str) -> None
    def report_failure(self, service: str, operation: str, config: CircuitBreakerConfig) -> None
```

### EnhancedRetryManager

Main interface for the retry mechanism:

```python
class EnhancedRetryManager:
    def execute(self, func, component, operation, args=None, kwargs=None, 
                config=None, use_circuit_breaker=True, circuit_breaker_config=None)
                
    def execute_async(self, func, component, operation, args=None, kwargs=None,
                      config=None, use_circuit_breaker=True, circuit_breaker_config=None)
                      
    def retry(self, func, component, operation, *args, max_attempts=None, 
              initial_delay=None, max_delay=None, backoff_factor=None, 
              jitter_factor=None, use_circuit_breaker=True, **kwargs)
              
    def retry_decorator(self, component, operation, max_attempts=None, 
                        initial_delay=None, max_delay=None, backoff_factor=None,
                        jitter_factor=None, use_circuit_breaker=True)
                        
    def configure_from_dict(self, config_dict)
    
    def get_circuit_breaker_status() -> Dict[str, Dict[str, Any]]
    
    def reset_circuit_breaker(self, service: str, operation: str) -> None
```

## Circuit Breaker Pattern

The circuit breaker pattern prevents cascading failures by temporarily disabling operations after repeated failures. It has three states:

### Closed State (Normal Operation)

- All requests are allowed to flow through
- Failures are counted
- When failure count exceeds the failure threshold, transition to Open state

### Open State (Failure Mode)

- All requests are blocked without attempting execution
- After a timeout period, transition to Half-Open state
- Each consecutive failure increases the timeout exponentially 

### Half-Open State (Recovery Mode)

- Limited requests are allowed through to test if the system has recovered
- If requests succeed, transition back to Closed state
- If requests fail, transition back to Open state with increased timeout

## Usage Examples

### Basic Usage

```python
from unified_importer.core.error_handling.retry_enhancement import EnhancedRetryManager

# Create a retry manager
retry_manager = EnhancedRetryManager()

# Execute a function with retry
try:
    result = retry_manager.execute(
        fetch_data, "ExternalAPI", "fetch_molecules",
        args=(molecule_ids,),
        kwargs={"include_properties": True}
    )
    process_data(result)
except Exception as e:
    handle_error(e)
```

### Async Usage

```python
import asyncio
from unified_importer.core.error_handling.retry_enhancement import EnhancedRetryManager

# Create a retry manager
retry_manager = EnhancedRetryManager()

# Execute an async function with retry
async def fetch_data_with_retry():
    try:
        result = await retry_manager.execute_async(
            fetch_data_async, "ExternalAPI", "fetch_molecules_async",
            args=(molecule_ids,),
            kwargs={"include_properties": True}
        )
        return result
    except Exception as e:
        handle_error(e)
        return None

# Run the async function
asyncio.run(fetch_data_with_retry())
```

### Using the Decorator

```python
from unified_importer.core.error_handling.retry_enhancement import EnhancedRetryManager

# Create a retry manager
retry_manager = EnhancedRetryManager()

# Define a function with retry capability
@retry_manager.retry_decorator(
    component="DataSource",
    operation="fetch_molecule_data",
    max_attempts=5,
    initial_delay=1.0,
    backoff_factor=2.0
)
def fetch_molecule_data(molecule_id):
    # This function will automatically retry on failure
    response = requests.get(f"https://api.example.com/molecules/{molecule_id}")
    response.raise_for_status()
    return response.json()

# Call the function - retries are handled automatically
try:
    data = fetch_molecule_data("CHEMBL123")
    process_data(data)
except Exception as e:
    # This exception means all retries failed
    handle_ultimate_failure(e)
```

### Using Simplified Retry

```python
from unified_importer.core.error_handling.retry_enhancement import EnhancedRetryManager

# Create a retry manager
retry_manager = EnhancedRetryManager()

# Use the simplified retry interface
try:
    # Will retry automatically with default settings
    result = retry_manager.retry(
        fetch_data, "ExternalAPI", "fetch_molecules", 
        molecule_ids, include_properties=True
    )
    process_data(result)
except Exception as e:
    handle_ultimate_failure(e)
```

## Custom Configuration

```python
from unified_importer.core.error_handling.retry_enhancement import (
    EnhancedRetryManager, RetryConfig, CircuitBreakerConfig
)
from unified_importer.core.error_handling.error_classification import ErrorCategory

# Create custom retry config
retry_config = RetryConfig(
    max_attempts=5,
    initial_delay=0.5,
    max_delay=30.0,
    backoff_factor=1.5,
    jitter_factor=0.2,
    retryable_error_categories={
        ErrorCategory.NETWORK,
        ErrorCategory.TIMEOUT,
        ErrorCategory.DATABASE
    },
    non_retryable_errors={
        ValueError,
        TypeError,
        KeyError
    }
)

# Create custom circuit breaker config
circuit_breaker_config = CircuitBreakerConfig(
    failure_threshold=3,
    reset_timeout=60.0,
    half_open_max_calls=2,
    success_threshold=2,
    open_timeout_increment_factor=2.0,
    max_open_timeout=300.0
)

# Create retry manager with custom configs
retry_manager = EnhancedRetryManager(
    default_config=retry_config,
    circuit_breaker_config=circuit_breaker_config
)

# Use in execution
result = retry_manager.execute(
    fetch_data, "ExternalAPI", "fetch_molecules",
    args=(molecule_ids,),
    kwargs={"include_properties": True},
    config=retry_config,  # Can also override per call
    use_circuit_breaker=True
)
```

## Configuration from Dictionary

```python
retry_manager = EnhancedRetryManager()

# Configure from dictionary (e.g., from config file)
retry_manager.configure_from_dict({
    'retry': {
        'max_attempts': 5,
        'initial_delay': 0.5,
        'max_delay': 30.0,
        'backoff_factor': 1.5,
        'jitter_factor': 0.2
    },
    'circuit_breaker': {
        'failure_threshold': 3,
        'reset_timeout': 60.0,
        'half_open_max_calls': 2,
        'success_threshold': 2,
        'open_timeout_increment_factor': 2.0,
        'max_open_timeout': 300.0
    }
})
```

## Circuit Breaker Monitoring and Reset

```python
# Get status of all circuit breakers
status = retry_manager.get_circuit_breaker_status()

# Print circuit breaker status
for key, state in status.items():
    service, operation = key.split(':')
    print(f"Circuit {service}:{operation} - State: {state['state']}")
    if state['state'] == 'OPEN':
        open_until = time.strftime('%H:%M:%S', time.localtime(state['open_until']))
        print(f"  Open until: {open_until}")
    print(f"  Failure count: {state['failure_count']}")

# Reset a circuit breaker manually
retry_manager.reset_circuit_breaker("ExternalAPI", "fetch_molecules")
```

## Integration with Error Classification

The Enhanced Retry Mechanism integrates seamlessly with the Error Classification System:

```python
from unified_importer.core.error_handling.error_classification import ErrorClassifier
from unified_importer.core.error_handling.retry_enhancement import EnhancedRetryManager

# Create error classifier
error_classifier = ErrorClassifier()

# Create retry manager with the classifier
retry_manager = EnhancedRetryManager(error_classifier=error_classifier)

# Execution will use the classifier to determine if errors are retryable
result = retry_manager.execute(
    fetch_data, "ExternalAPI", "fetch_molecules",
    args=(molecule_ids,)
)
```

## Best Practices

### 1. Configure Appropriate Retry Policies

Adjust retry parameters based on the operation:
- For quick operations with potentially transient failures, use short initial delays and more attempts
- For slow external API calls, use longer initial delays and fewer attempts

```python
# Quick database operations
retry_manager.retry(
    db_operation, "Database", "query_molecules",
    max_attempts=5, initial_delay=0.1, backoff_factor=1.5
)

# Slow external API calls
retry_manager.retry(
    api_operation, "ExternalAPI", "fetch_molecules",
    max_attempts=3, initial_delay=1.0, backoff_factor=2.0
)
```

### 2. Use Circuit Breakers for External Dependencies

Circuit breakers are essential for preventing cascading failures when external services degrade:

```python
retry_manager.execute(
    external_service_call, "ExternalAPI", "fetch_data",
    use_circuit_breaker=True,
    circuit_breaker_config=CircuitBreakerConfig(
        failure_threshold=3,
        reset_timeout=60.0
    )
)
```

### 3. Group Related Operations Under the Same Circuit Breaker

Use consistent component and operation identifiers to group related operations:

```python
# These share the same circuit breaker
retry_manager.retry(fetch_molecules, "ChEMBL", "fetch_data")
retry_manager.retry(search_molecules, "ChEMBL", "fetch_data")

# This has a separate circuit breaker
retry_manager.retry(fetch_properties, "ChEMBL", "fetch_properties")
```

### 4. Monitor Circuit Breaker Status

Regularly check the status of circuit breakers to detect service degradation:

```python
# Get all circuit breaker statuses
status = retry_manager.get_circuit_breaker_status()

# Log any open circuits
for key, state in status.items():
    if state['state'] == 'OPEN':
        logger.warning(f"Circuit {key} is OPEN until {time.ctime(state['open_until'])}")
```

### 5. Use Persistent Circuit Breaker State

Configure the registry to persist circuit breaker state across application restarts:

```python
# Create a registry with persistence
registry = CircuitBreakerRegistry(state_file_path="/path/to/circuit_state.json")

# Create retry manager with the registry  
retry_manager = EnhancedRetryManager(circuit_breaker_registry=registry)
```

## Tuning Advice

### Exponential Backoff

Adjust backoff parameters based on the expected recovery time of the service:
- For services with quick recovery (e.g., database connections), use smaller backoff factors (1.2-1.5)
- For services with slower recovery (e.g., external APIs with rate limiting), use larger backoff factors (2.0-3.0)

### Jitter

Jitter prevents "thundering herd" problems in distributed systems:
- Use jitter factors of 0.1-0.2 for most services
- For heavily loaded services, increase jitter factor to 0.2-0.3

### Circuit Breakers

Tune circuit breaker parameters based on service characteristics:
- For critical services, use higher failure thresholds (5-10)
- For non-critical services, use lower failure thresholds (2-3)
- For services with slow recovery, use longer reset timeouts (60-300 seconds)
- For services with quick recovery, use shorter reset timeouts (10-30 seconds)

## Conclusion

The Enhanced Retry Mechanism provides a comprehensive solution for handling transient failures in distributed systems. By configuring appropriate retry policies and leveraging circuit breakers, you can build resilient applications that gracefully handle service degradation and recovery.

For more information, see:
- [Error Handling System](error_handling.md) - The complete error handling framework
- [Circuit Breaker Pattern](https://martinfowler.com/bliki/CircuitBreaker.html) - Martin Fowler's explanation of the pattern
- [Retries with Exponential Backoff](https://aws.amazon.com/blogs/architecture/exponential-backoff-and-jitter/) - AWS Architecture Blog