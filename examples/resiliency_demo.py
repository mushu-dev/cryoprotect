#!/usr/bin/env python3
"""
CryoProtect API - Resiliency Patterns Demo

This script demonstrates how to use the resiliency patterns implemented
in the CryoProtect API, including:

1. Retry Mechanism with Exponential Backoff
2. Circuit Breaker Pattern
3. Timeout Management
4. Service Health Tracking
5. Combined Patterns for Maximum Resiliency

Run this script to see the resiliency patterns in action.
"""

import time
import random
import logging
import sys
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[logging.StreamHandler()]
)

# Import resiliency patterns
sys.path.append('/home/mushu/Projects/cryoprotect')
from api.resiliency.retry import retry_with_backoff
from api.resiliency.circuit_breaker import circuit_breaker, CircuitBreakerState
from api.resiliency.timeout import with_timeout, TimeoutError
from api.resiliency.service_health import track_service_health, get_service_health

# Get logger
logger = logging.getLogger(__name__)

# ===== Simulated Services =====

class UnreliableService:
    """
    Simulates an unreliable service with various failure modes.
    """
    
    def __init__(self, name, success_rate=0.7, avg_response_time=0.5, timeout_rate=0.1):
        """
        Initialize the unreliable service.
        
        Args:
            name: Service name
            success_rate: Probability of success (0.0-1.0)
            avg_response_time: Average response time in seconds
            timeout_rate: Probability of a timeout (0.0-1.0)
        """
        self.name = name
        self.success_rate = success_rate
        self.avg_response_time = avg_response_time
        self.timeout_rate = timeout_rate
        self.call_count = 0
    
    def call(self):
        """
        Call the unreliable service.
        
        Returns:
            dict: Response data
            
        Raises:
            ConnectionError: If the service fails
            TimeoutError: If the service times out
        """
        # Increment call count
        self.call_count += 1
        
        # Simulate variable response time
        response_time = random.expovariate(1.0 / self.avg_response_time)
        
        # Check for timeout
        if random.random() < self.timeout_rate:
            # Sleep for a long time to simulate a timeout
            logger.debug(f"Service {self.name} is timing out")
            time.sleep(10)
            return {"error": "Timeout"}
        
        # Sleep to simulate response time
        time.sleep(response_time)
        
        # Check for failure
        if random.random() > self.success_rate:
            logger.debug(f"Service {self.name} is failing")
            raise ConnectionError(f"Service {self.name} failed")
        
        # Return success
        logger.debug(f"Service {self.name} succeeded in {response_time:.2f}s")
        return {
            "data": f"Response from {self.name}",
            "timestamp": datetime.now().isoformat(),
            "call_count": self.call_count
        }

# ===== Resiliency Patterns =====

# Create unreliable services
database_service = UnreliableService("database", success_rate=0.6, avg_response_time=0.3)
api_service = UnreliableService("api", success_rate=0.5, avg_response_time=0.7, timeout_rate=0.2)

# Retry only
@retry_with_backoff(max_retries=3, base_delay=0.5, backoff_factor=2, jitter=0.1)
def call_with_retry(service):
    """Call a service with retry."""
    return service.call()

# Circuit breaker only
@circuit_breaker(name="service_circuit", failure_threshold=3, recovery_timeout=5)
def call_with_circuit_breaker(service):
    """Call a service with circuit breaker."""
    return service.call()

# Timeout only
@with_timeout(seconds=2)
def call_with_timeout(service):
    """Call a service with timeout."""
    return service.call()

# Service health tracking only
@track_service_health("service_health")
def call_with_health_tracking(service):
    """Call a service with health tracking."""
    return service.call()

# Combined patterns
@retry_with_backoff(max_retries=3, base_delay=0.5, backoff_factor=2, jitter=0.1)
@circuit_breaker(name="combined_circuit", failure_threshold=3, recovery_timeout=5)
@with_timeout(seconds=2)
@track_service_health("combined_service")
def call_with_all_patterns(service):
    """Call a service with all resiliency patterns."""
    return service.call()

def demo_retry():
    """Demonstrate the retry pattern."""
    logger.info("===== RETRY PATTERN DEMO =====")
    
    success_count = 0
    failure_count = 0
    
    for i in range(10):
        try:
            logger.info(f"Call {i+1}: Calling service with retry...")
            result = call_with_retry(api_service)
            logger.info(f"Call {i+1}: Success! Result: {result}")
            success_count += 1
        except Exception as e:
            logger.error(f"Call {i+1}: Failed after retries: {str(e)}")
            failure_count += 1
    
    logger.info(f"Retry demo complete: {success_count} successes, {failure_count} failures")

def demo_circuit_breaker():
    """Demonstrate the circuit breaker pattern."""
    logger.info("===== CIRCUIT BREAKER PATTERN DEMO =====")
    
    # Reset the service for demonstration
    api_service.success_rate = 0.3
    
    success_count = 0
    failure_count = 0
    circuit_open_count = 0
    
    for i in range(15):
        try:
            logger.info(f"Call {i+1}: Calling service with circuit breaker...")
            result = call_with_circuit_breaker(api_service)
            logger.info(f"Call {i+1}: Success! Result: {result}")
            success_count += 1
        except Exception as e:
            if "circuit breaker" in str(e).lower():
                logger.warning(f"Call {i+1}: Circuit open: {str(e)}")
                circuit_open_count += 1
            else:
                logger.error(f"Call {i+1}: Service error: {str(e)}")
                failure_count += 1
    
    # Reset service to normal
    api_service.success_rate = 0.7
    
    logger.info("Waiting for circuit recovery...")
    time.sleep(6)
    
    # Try again after recovery
    for i in range(5):
        try:
            logger.info(f"Recovery call {i+1}: Calling service after circuit recovery...")
            result = call_with_circuit_breaker(api_service)
            logger.info(f"Recovery call {i+1}: Success! Result: {result}")
            success_count += 1
        except Exception as e:
            if "circuit breaker" in str(e).lower():
                logger.warning(f"Recovery call {i+1}: Circuit still open: {str(e)}")
                circuit_open_count += 1
            else:
                logger.error(f"Recovery call {i+1}: Service error: {str(e)}")
                failure_count += 1
    
    logger.info(
        f"Circuit breaker demo complete: {success_count} successes, "
        f"{failure_count} failures, {circuit_open_count} circuit open"
    )

def demo_timeout():
    """Demonstrate the timeout pattern."""
    logger.info("===== TIMEOUT PATTERN DEMO =====")
    
    # Increase timeout rate for demonstration
    api_service.timeout_rate = 0.5
    
    success_count = 0
    timeout_count = 0
    
    for i in range(10):
        try:
            logger.info(f"Call {i+1}: Calling service with timeout...")
            result = call_with_timeout(api_service)
            logger.info(f"Call {i+1}: Success! Result: {result}")
            success_count += 1
        except TimeoutError as e:
            logger.warning(f"Call {i+1}: Timeout: {str(e)}")
            timeout_count += 1
        except Exception as e:
            logger.error(f"Call {i+1}: Other error: {str(e)}")
    
    # Reset service to normal
    api_service.timeout_rate = 0.2
    
    logger.info(f"Timeout demo complete: {success_count} successes, {timeout_count} timeouts")

def demo_service_health():
    """Demonstrate the service health tracking pattern."""
    logger.info("===== SERVICE HEALTH TRACKING DEMO =====")
    
    # Call the service multiple times to build up health metrics
    for i in range(20):
        try:
            logger.info(f"Call {i+1}: Calling service with health tracking...")
            result = call_with_health_tracking(database_service)
            logger.info(f"Call {i+1}: Success!")
        except Exception as e:
            logger.error(f"Call {i+1}: Error: {str(e)}")
    
    # Get service health metrics
    health = get_service_health("service_health")
    metrics = health.get_metrics()
    
    logger.info(f"Service health metrics: {metrics}")
    logger.info(f"Service is healthy: {health.is_healthy()}")
    logger.info(f"Service is degraded: {health.is_degraded()}")
    logger.info(f"Service is unhealthy: {health.is_unhealthy()}")

def demo_combined_patterns():
    """Demonstrate combined resiliency patterns."""
    logger.info("===== COMBINED PATTERNS DEMO =====")
    
    success_count = 0
    failure_count = 0
    timeout_count = 0
    circuit_open_count = 0
    
    for i in range(20):
        try:
            logger.info(f"Call {i+1}: Calling service with all patterns...")
            result = call_with_all_patterns(api_service)
            logger.info(f"Call {i+1}: Success! Result: {result}")
            success_count += 1
        except TimeoutError as e:
            logger.warning(f"Call {i+1}: Timeout: {str(e)}")
            timeout_count += 1
        except Exception as e:
            if "circuit breaker" in str(e).lower():
                logger.warning(f"Call {i+1}: Circuit open: {str(e)}")
                circuit_open_count += 1
            else:
                logger.error(f"Call {i+1}: Failed: {str(e)}")
                failure_count += 1
    
    # Wait for circuit recovery
    logger.info("Waiting for circuit recovery...")
    time.sleep(6)
    
    # Try again after recovery
    for i in range(5):
        try:
            logger.info(f"Recovery call {i+1}: Calling service after circuit recovery...")
            result = call_with_all_patterns(api_service)
            logger.info(f"Recovery call {i+1}: Success! Result: {result}")
            success_count += 1
        except Exception as e:
            logger.error(f"Recovery call {i+1}: Failed: {str(e)}")
            failure_count += 1
    
    # Get service health metrics
    health = get_service_health("combined_service")
    metrics = health.get_metrics()
    
    logger.info(f"Service health metrics: {metrics}")
    logger.info(
        f"Combined patterns demo complete: {success_count} successes, "
        f"{failure_count} failures, {timeout_count} timeouts, "
        f"{circuit_open_count} circuit open"
    )

def main():
    """Run the resiliency patterns demo."""
    logger.info("Starting resiliency patterns demo...")
    
    # Demonstrate each pattern
    demo_retry()
    demo_circuit_breaker()
    demo_timeout()
    demo_service_health()
    demo_combined_patterns()
    
    logger.info("Resiliency patterns demo complete!")

if __name__ == "__main__":
    main()