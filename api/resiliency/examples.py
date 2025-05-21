"""
CryoProtect v2 API - Resiliency Patterns Examples

This module provides examples of how to use the resiliency patterns implemented
in the CryoProtect API.

Examples include:
- Retry mechanism with exponential backoff
- Circuit breaker pattern
- Timeout handling
- Service health tracking
- Combined patterns for maximum resiliency

These examples can be used as templates for implementing resiliency in your own
code, whether it's for API calls, database operations, or any other operations
that may fail transiently.
"""

import time
import random
import logging
import requests
from typing import Dict, Any, List, Optional

# Import resiliency patterns
from .retry import retry_with_backoff, retry_network_operations, retry_database_operations
from .circuit_breaker import circuit_breaker, CircuitBreakerState
from .timeout import with_timeout, TimeoutError
from .service_health import track_service_health, get_service_health

# Get logger
logger = logging.getLogger(__name__)

# ===== Retry Pattern Examples =====

@retry_with_backoff(max_retries=3, base_delay=0.1)
def example_basic_retry():
    """
    Basic retry example.
    
    This function will retry up to 3 times with exponential backoff
    if it raises any exception.
    """
    # Simulate a transient error with 80% probability
    if random.random() < 0.8:
        raise ConnectionError("Simulated connection error")
    
    return "Operation completed successfully"

@retry_with_backoff(
    max_retries=5,
    exceptions=(ConnectionError, TimeoutError, requests.RequestException),
    base_delay=0.5,
    max_delay=30,
    backoff_factor=2,
    jitter=0.1
)
def example_advanced_retry(url: str) -> Dict[str, Any]:
    """
    Advanced retry example with HTTP requests.
    
    This function will retry up to 5 times with exponential backoff
    if it raises ConnectionError, TimeoutError, or RequestException.
    
    Args:
        url: URL to request
        
    Returns:
        dict: JSON response
    """
    # Make an HTTP request with a timeout
    response = requests.get(url, timeout=5)
    
    # Raise an exception for any 4xx/5xx errors
    response.raise_for_status()
    
    # Return the JSON response
    return response.json()

@retry_network_operations()
def example_network_retry():
    """
    Network operation retry example.
    
    This function uses a specialized retry decorator for network operations,
    which already includes common network-related exceptions and appropriate
    backoff settings.
    """
    # Make an HTTP request
    response = requests.get("https://api.example.com/data")
    response.raise_for_status()
    return response.json()

@retry_database_operations()
def example_database_retry():
    """
    Database operation retry example.
    
    This function uses a specialized retry decorator for database operations,
    which already includes common database-related exceptions and appropriate
    backoff settings.
    """
    # Simulate a database query
    # In a real application, this would be a database query
    if random.random() < 0.3:
        # Simulate a transient database error
        raise ConnectionError("Database connection error")
    
    return [{"id": 1, "name": "Example"}]

# ===== Circuit Breaker Pattern Examples =====

@circuit_breaker(name="example_basic", failure_threshold=3)
def example_basic_circuit_breaker():
    """
    Basic circuit breaker example.
    
    This function will trip the circuit breaker after 3 consecutive failures,
    and will attempt recovery after 60 seconds (default).
    """
    # Simulate a failure with 80% probability
    if random.random() < 0.8:
        raise ConnectionError("Simulated connection error")
    
    return "Operation completed successfully"

@circuit_breaker(
    name="example_advanced",
    failure_threshold=5,
    recovery_timeout=30,
    exceptions=(ConnectionError, TimeoutError),
    fallback=lambda: {"error": "Service unavailable, please try again later"}
)
def example_advanced_circuit_breaker():
    """
    Advanced circuit breaker example with fallback.
    
    This function will trip the circuit breaker after 5 consecutive failures,
    and will attempt recovery after 30 seconds. It uses a fallback function
    to provide a default response when the circuit is open.
    """
    # Simulate a failure with 80% probability
    if random.random() < 0.8:
        raise ConnectionError("Simulated connection error")
    
    return {"data": "Operation completed successfully"}

def check_circuit_state():
    """
    Check the state of all circuit breakers.
    
    Returns:
        dict: Circuit breaker states
    """
    from .circuit_breaker import CircuitBreaker
    
    # Get all circuit breakers
    circuit_breakers = CircuitBreaker.get_all_circuit_breakers()
    
    # Build a dictionary of circuit breaker states
    states = {}
    for name, circuit in circuit_breakers.items():
        states[name] = {
            "state": circuit.get_state().name,
            "metrics": circuit.get_metrics()
        }
    
    return states

# ===== Timeout Pattern Examples =====

@with_timeout(seconds=5)
def example_basic_timeout():
    """
    Basic timeout example.
    
    This function will raise TimeoutError if it doesn't complete within 5 seconds.
    """
    # Simulate a long-running operation
    time.sleep(random.uniform(1, 10))
    
    return "Operation completed successfully"

@with_timeout(seconds=3, exception_class=requests.Timeout)
def example_custom_timeout_exception():
    """
    Timeout example with custom exception.
    
    This function will raise requests.Timeout if it doesn't complete within 3 seconds.
    """
    # Simulate a long-running operation
    time.sleep(random.uniform(1, 5))
    
    return "Operation completed successfully"

# ===== Service Health Tracking Examples =====

@track_service_health("example_service")
def example_service_health():
    """
    Service health tracking example.
    
    This function's health metrics will be tracked under the name "example_service".
    """
    # Simulate variable response times
    response_time = random.uniform(0.1, 2.0)
    time.sleep(response_time)
    
    # Simulate occasional failures
    if random.random() < 0.2:
        raise ConnectionError("Simulated connection error")
    
    return "Operation completed successfully"

def check_service_health():
    """
    Check the health of all services.
    
    Returns:
        dict: Service health metrics
    """
    from .service_health import ServiceHealth
    
    # Get all service health trackers
    service_health_trackers = ServiceHealth.get_all_service_health()
    
    # Build a dictionary of service health metrics
    health_metrics = {}
    for name, tracker in service_health_trackers.items():
        health_metrics[name] = tracker.get_metrics()
    
    return health_metrics

# ===== Combined Patterns Examples =====

@retry_with_backoff(max_retries=3)
@circuit_breaker(name="combined")
@with_timeout(seconds=5)
@track_service_health("combined_service")
def example_combined_patterns():
    """
    Combined resiliency patterns example.
    
    This function combines all resiliency patterns:
    - Retry with exponential backoff (3 attempts)
    - Circuit breaker (open after 5 consecutive failures, default)
    - Timeout (5 seconds)
    - Service health tracking (under "combined_service")
    
    The execution order is:
    1. Timeout (innermost) - Prevents operations from hanging
    2. Service health - Tracks metrics for the operation
    3. Circuit breaker - Prevents calls when the service is failing
    4. Retry (outermost) - Retries on transient failures
    
    This order ensures that the timeout is always enforced first, then service
    health is tracked, then the circuit breaker checks if the operation should
    be attempted, and finally retry handles any transient failures.
    """
    # Simulate a complex operation with various failure modes
    r = random.random()
    
    if r < 0.4:
        # Simulate a timeout
        time.sleep(10)
        return "This should never be reached due to timeout"
    
    elif r < 0.7:
        # Simulate a transient error
        raise ConnectionError("Simulated connection error")
    
    # Simulate successful operation with variable response time
    time.sleep(random.uniform(0.1, 1.0))
    
    return "Operation completed successfully"

# ===== Real-world Examples =====

@retry_network_operations(max_retries=3)
@circuit_breaker(name="external_api")
@with_timeout(seconds=10)
@track_service_health("external_api")
def fetch_from_external_api(url: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Fetch data from an external API with full resiliency.
    
    Args:
        url: URL to request
        params: Optional query parameters
        
    Returns:
        dict: JSON response
    """
    # Make a request to the external API
    response = requests.get(url, params=params, timeout=9)
    
    # Raise an exception for any 4xx/5xx errors
    response.raise_for_status()
    
    # Return the JSON response
    return response.json()

@retry_database_operations(max_retries=3)
@circuit_breaker(name="database", failure_threshold=3, recovery_timeout=60)
@with_timeout(seconds=5)
@track_service_health("database")
def get_database_records(query: str, params: Optional[List[Any]] = None) -> List[Dict[str, Any]]:
    """
    Get records from the database with full resiliency.
    
    Args:
        query: SQL query
        params: Optional query parameters
        
    Returns:
        list: Query results
    """
    # In a real application, this would be a database query
    # using a database client like psycopg2, sqlite3, etc.
    
    # Simulate a database query with various failure modes
    r = random.random()
    
    if r < 0.1:
        # Simulate a timeout
        time.sleep(10)
        return "This should never be reached due to timeout"
    
    elif r < 0.2:
        # Simulate a transient error
        raise ConnectionError("Database connection error")
    
    # Simulate successful operation with variable response time
    time.sleep(random.uniform(0.05, 0.5))
    
    # Return simulated database records
    return [{"id": i, "name": f"Record {i}"} for i in range(10)]