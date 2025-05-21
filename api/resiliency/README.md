# CryoProtect Resiliency Patterns

This directory contains the resiliency patterns implemented for the CryoProtect API v2. These patterns help ensure robust operation in production environments by handling transient failures, preventing cascading failures, and providing better user experience during service degradation.

## Core Components

### 1. Retry Mechanism (`retry.py`)

The retry mechanism automatically retries failed operations with exponential backoff to handle transient failures.

```python
from api.resiliency.retry import retry_with_backoff

@retry_with_backoff(max_retries=3, base_delay=0.5)
def external_api_call():
    # This function will retry up to 3 times with exponential backoff
    response = requests.get("https://api.example.com/data")
    response.raise_for_status()
    return response.json()
```

Key features:
- Configurable retry attempts, delay, backoff factor, and jitter
- Support for specific exception types
- Built-in logging and observability
- Specialized decorators for common operations

### 2. Circuit Breaker (`circuit_breaker.py`)

The circuit breaker pattern prevents cascading failures by temporarily stopping calls to failing services.

```python
from api.resiliency.circuit_breaker import circuit_breaker

@circuit_breaker(name="database", failure_threshold=3, recovery_timeout=30)
def database_operation():
    # Circuit opens after 3 consecutive failures
    # Recovery is attempted after 30 seconds
    return db.query("SELECT * FROM data")
```

Key features:
- Three states: CLOSED (normal), OPEN (failing), HALF-OPEN (testing recovery)
- Configurable failure thresholds and recovery timeouts
- Thread-safe implementation
- Optional fallback function for when the circuit is open

### 3. Timeout Management (`timeout.py`)

The timeout mechanism enforces time limits on operations to prevent resource blockage.

```python
from api.resiliency.timeout import with_timeout

@with_timeout(seconds=5)
def time_sensitive_operation():
    # This function will raise TimeoutError if it doesn't complete within 5 seconds
    process_data()
```

Key features:
- Thread-safe implementation for multi-threaded applications
- Configurable timeout thresholds
- Integration with logging and observability

### 4. Service Health Tracking (`service_health.py`)

The service health tracker monitors the health of dependent services based on success/failure rates and response times.

```python
from api.resiliency.service_health import track_service_health, get_service_health

@track_service_health("external_api")
def call_external_api():
    # This function's success/failure will be tracked
    return requests.get("https://api.example.com/data").json()

# Check service health
health = get_service_health("external_api")
if health.is_healthy():
    # Service is healthy
    pass
elif health.is_degraded():
    # Service is degraded
    pass
else:
    # Service is unhealthy
    pass
```

Key features:
- Real-time health status based on success rates and response times
- Automatic detection of degraded or failing services
- Configurable thresholds for degradation and failure

## Combined Usage

These patterns can be combined for maximum resiliency:

```python
from api.resiliency import retry_with_backoff, circuit_breaker, with_timeout, track_service_health

@retry_with_backoff(max_retries=3)
@circuit_breaker(name="database")
@with_timeout(seconds=5)
@track_service_health("database")
def robust_database_operation():
    # This function has comprehensive resiliency
    return db.query("SELECT * FROM data")
```

The execution order (from innermost to outermost) is important:
1. **Timeout** (innermost) - Prevents operations from hanging
2. **Service Health** - Tracks metrics for the operation
3. **Circuit Breaker** - Prevents calls when the service is failing
4. **Retry** (outermost) - Retries on transient failures

## Demo

For a complete demonstration of these patterns, run:

```bash
python examples/resiliency_demo.py
```

This script demonstrates each pattern individually and in combination, with simulated failures to show how the patterns respond.

## Unit Tests

Unit tests for these patterns are in the `tests/unit/` directory:

- `test_retry_mechanism.py`: Tests for the retry mechanism
- `test_circuit_breaker.py`: Tests for the circuit breaker pattern
- `test_timeout.py`: Tests for the timeout mechanism
- `test_service_health.py`: Tests for the service health tracker

## Integration with Observability

All resiliency patterns integrate with the observability system to provide metrics and logs for monitoring and debugging. If the observability system is not available, they fall back to basic logging.

### Metrics

- **Retry Attempts**: Count of retry attempts by function and exception type
- **Circuit Breaker State**: Current state of each circuit breaker
- **Circuit Breaker Failures/Successes**: Count of failures and successes by circuit
- **Timeout Errors**: Count of timeout errors by function
- **Service Health Status**: Current health status of each service

## Documentation

For detailed information about each pattern, refer to the docstrings in the source files. For a high-level overview of the production readiness plan, see the [Production Readiness Plan](../../PRODUCTION_READINESS_PLAN.md) document.