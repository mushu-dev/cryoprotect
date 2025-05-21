"""
CryoProtect v2 API - Resiliency Module

This module provides resiliency patterns for making the API more robust against failures:

1. Retry Mechanism:
   - Decorator-based approach with configurable retry attempts
   - Exponential backoff with jitter for optimal retry timing
   - Customizable retry conditions based on exception types
   - Automatic logging and observability integration

2. Circuit Breaker:
   - Prevents cascading failures by stopping repeated failing calls
   - Automatic recovery with half-open state
   - Configurable thresholds and cooldown periods
   - Integration with monitoring systems

3. Timeout Management:
   - Enforce time limits on operations to prevent blocked resources
   - Configurable timeout thresholds
   - Integration with the retry mechanism

4. Service Health Tracking:
   - Monitor dependent service health
   - Automatic service degradation detection
   - Circuit breaker integration

Usage:
    from api.resiliency import retry_with_backoff, circuit_breaker, with_timeout
    
    @retry_with_backoff(max_retries=3, exceptions=(ConnectionError, TimeoutError))
    def external_api_call():
        # This function will retry up to 3 times if it raises ConnectionError or TimeoutError
        pass
        
    @circuit_breaker(failure_threshold=5, recovery_timeout=30)
    def database_operation():
        # This function will stop being called after 5 consecutive failures
        # and will recover after 30 seconds
        pass
        
    @with_timeout(seconds=5)
    def time_sensitive_operation():
        # This function will raise TimeoutError if it doesn't complete within 5 seconds
        pass
"""

# Import exposed components
from .retry import retry_with_backoff
from .circuit_breaker import circuit_breaker, CircuitBreakerState
from .timeout import with_timeout