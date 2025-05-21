"""
CryoProtect v2 API - Service Health Tracking

This module provides a mechanism to track the health of dependent services and 
automatically adapt to degraded or failing services.

Key features:
- Tracks service health based on success and failure rates
- Automatically detects degraded performance
- Integrates with circuit breaker and retry mechanisms
- Provides health status for monitoring and alerting

Usage:
    from api.resiliency.service_health import (
        track_service_health,
        get_service_health,
        ServiceHealth
    )
    
    # Track service health
    @track_service_health("database")
    def database_operation():
        # This function's success/failure will be tracked
        pass
    
    # Check service health status
    health_status = get_service_health("database")
    if health_status.is_healthy():
        # Service is healthy
        pass
    elif health_status.is_degraded():
        # Service is degraded
        pass
    else:
        # Service is unhealthy
        pass

Benefits:
- Provides real-time health status of dependent services
- Enables adaptive behavior based on service health
- Supports automatic degradation detection
- Integrates with monitoring and alerting systems
"""

import time
import enum
import functools
import threading
import logging
from typing import Dict, Any, Optional, Callable, List
import statistics
from collections import deque

# Get logger
logger = logging.getLogger(__name__)

# Import observability if available
try:
    from ..observability import (
        log_with_context, 
        report_error,
        get_correlation_id,
        get_request_context
    )
    from monitoring.prometheus_metrics import (
        SERVICE_HEALTH_STATUS,
        SERVICE_RESPONSE_TIME
    )
    HAS_OBSERVABILITY = True
except ImportError:
    # Fallback to basic logging if observability is not available
    def log_with_context(logger, level, message, context=None, exc_info=False):
        getattr(logger, level)(message, exc_info=exc_info)
    
    def report_error(error, context=None):
        logger.error(f"Error: {str(error)}", exc_info=error)
    
    def get_correlation_id():
        return None
    
    def get_request_context():
        return {}
    
    # Dummy metrics
    class DummyGauge:
        def __init__(self):
            pass
        
        def labels(self, **kwargs):
            return self
        
        def set(self, value):
            pass
    
    SERVICE_HEALTH_STATUS = DummyGauge()
    SERVICE_RESPONSE_TIME = DummyGauge()
    HAS_OBSERVABILITY = False

class HealthStatus(enum.Enum):
    """Service health status."""
    HEALTHY = 0      # Service is operating normally
    DEGRADED = 1     # Service is experiencing some issues but still operational
    UNHEALTHY = 2    # Service is experiencing severe issues or is unavailable

class PerformanceStats:
    """
    Performance statistics for a service.
    
    This class tracks response times, error rates, and success rates
    for a service over a sliding window of time.
    """
    
    def __init__(self, window_size=100, window_time=300):
        """
        Initialize performance statistics.
        
        Args:
            window_size: Maximum number of data points to keep (default: 100)
            window_time: Maximum age of data points in seconds (default: 300)
        """
        self.window_size = window_size
        self.window_time = window_time
        self.response_times = deque(maxlen=window_size)
        self.response_times_timestamps = deque(maxlen=window_size)
        self.success_count = 0
        self.failure_count = 0
        self.total_count = 0
        self.lock = threading.RLock()
    
    def record_success(self, response_time):
        """
        Record a successful operation.
        
        Args:
            response_time: Response time in seconds
        """
        with self.lock:
            # Clean up old data
            self._clean_old_data()
            
            # Record response time
            self.response_times.append(response_time)
            self.response_times_timestamps.append(time.time())
            
            # Update counters
            self.success_count += 1
            self.total_count += 1
    
    def record_failure(self):
        """Record a failed operation."""
        with self.lock:
            # Update counters
            self.failure_count += 1
            self.total_count += 1
    
    def _clean_old_data(self):
        """Remove data points older than window_time."""
        if not self.response_times_timestamps:
            return
        
        current_time = time.time()
        while (self.response_times_timestamps and 
               current_time - self.response_times_timestamps[0] > self.window_time):
            self.response_times.popleft()
            self.response_times_timestamps.popleft()
    
    def get_average_response_time(self):
        """
        Get the average response time.
        
        Returns:
            float: Average response time in seconds, or None if no data
        """
        with self.lock:
            # Clean up old data
            self._clean_old_data()
            
            if not self.response_times:
                return None
            
            return statistics.mean(self.response_times)
    
    def get_p95_response_time(self):
        """
        Get the 95th percentile response time.
        
        Returns:
            float: 95th percentile response time in seconds, or None if no data
        """
        with self.lock:
            # Clean up old data
            self._clean_old_data()
            
            if not self.response_times or len(self.response_times) < 10:
                return None
            
            sorted_times = sorted(self.response_times)
            index = int(len(sorted_times) * 0.95)
            return sorted_times[index]
    
    def get_error_rate(self):
        """
        Get the error rate.
        
        Returns:
            float: Error rate (0.0-1.0), or None if no data
        """
        with self.lock:
            if self.total_count == 0:
                return None
            
            return self.failure_count / self.total_count
    
    def get_success_rate(self):
        """
        Get the success rate.
        
        Returns:
            float: Success rate (0.0-1.0), or None if no data
        """
        with self.lock:
            if self.total_count == 0:
                return None
            
            return self.success_count / self.total_count
    
    def reset(self):
        """Reset all statistics."""
        with self.lock:
            self.response_times.clear()
            self.response_times_timestamps.clear()
            self.success_count = 0
            self.failure_count = 0
            self.total_count = 0

class ServiceHealth:
    """
    Service health status and metrics.
    
    This class tracks the health of a service based on performance statistics
    and thresholds for degradation and failure.
    """
    
    # Class-level registry of service health trackers
    _registry = {}
    _registry_lock = threading.RLock()
    
    @classmethod
    def get_service_health(cls, name):
        """
        Get or create a service health tracker by name.
        
        Args:
            name: Name of the service
            
        Returns:
            ServiceHealth instance
        """
        with cls._registry_lock:
            if name not in cls._registry:
                cls._registry[name] = ServiceHealth(name)
            return cls._registry[name]
    
    @classmethod
    def get_all_service_health(cls):
        """
        Get all registered service health trackers.
        
        Returns:
            Dict of service health instances by name
        """
        with cls._registry_lock:
            return cls._registry.copy()
    
    def __init__(self, name):
        """
        Initialize a new service health tracker.
        
        Args:
            name: Name of the service for monitoring
        """
        self.name = name
        self.status = HealthStatus.HEALTHY
        self.last_status_change = time.time()
        self.performance_stats = PerformanceStats()
        self.lock = threading.RLock()
        
        # Default thresholds (can be overridden)
        self.degraded_response_time_threshold = 1.0  # seconds
        self.unhealthy_response_time_threshold = 3.0  # seconds
        self.degraded_error_rate_threshold = 0.05  # 5%
        self.unhealthy_error_rate_threshold = 0.2   # 20%
        
        # Update metrics
        if HAS_OBSERVABILITY:
            SERVICE_HEALTH_STATUS.labels(service=self.name).set(self.status.value)
    
    def on_success(self, response_time):
        """
        Record a successful operation.
        
        Args:
            response_time: Response time in seconds
        """
        with self.lock:
            # Record success in performance stats
            self.performance_stats.record_success(response_time)
            
            # Update status based on new stats
            self._update_status()
            
            # Update metrics
            if HAS_OBSERVABILITY:
                SERVICE_HEALTH_STATUS.labels(service=self.name).set(self.status.value)
                SERVICE_RESPONSE_TIME.labels(service=self.name).set(response_time)
    
    def on_failure(self):
        """Record a failed operation."""
        with self.lock:
            # Record failure in performance stats
            self.performance_stats.record_failure()
            
            # Update status based on new stats
            self._update_status()
            
            # Update metrics
            if HAS_OBSERVABILITY:
                SERVICE_HEALTH_STATUS.labels(service=self.name).set(self.status.value)
    
    def _update_status(self):
        """Update service health status based on performance statistics."""
        prev_status = self.status
        
        # Get current statistics
        avg_response_time = self.performance_stats.get_average_response_time()
        p95_response_time = self.performance_stats.get_p95_response_time()
        error_rate = self.performance_stats.get_error_rate()
        
        # Check if service is unhealthy
        if (error_rate is not None and error_rate >= self.unhealthy_error_rate_threshold or
                p95_response_time is not None and p95_response_time >= self.unhealthy_response_time_threshold):
            self.status = HealthStatus.UNHEALTHY
        
        # Check if service is degraded
        elif (error_rate is not None and error_rate >= self.degraded_error_rate_threshold or
                avg_response_time is not None and avg_response_time >= self.degraded_response_time_threshold):
            self.status = HealthStatus.DEGRADED
        
        # Otherwise, service is healthy
        else:
            self.status = HealthStatus.HEALTHY
        
        # Log status change
        if self.status != prev_status:
            self.last_status_change = time.time()
            
            log_with_context(
                logger, 'warning' if self.status != HealthStatus.HEALTHY else 'info',
                f"Service '{self.name}' health status changed from {prev_status.name} to {self.status.name}",
                context={
                    'event_type': 'service_health_change',
                    'service': {
                        'name': self.name,
                        'status': self.status.name,
                        'previous_status': prev_status.name,
                        'avg_response_time': avg_response_time,
                        'p95_response_time': p95_response_time,
                        'error_rate': error_rate,
                        'success_rate': self.performance_stats.get_success_rate()
                    },
                    'thresholds': {
                        'degraded_response_time': self.degraded_response_time_threshold,
                        'unhealthy_response_time': self.unhealthy_response_time_threshold,
                        'degraded_error_rate': self.degraded_error_rate_threshold,
                        'unhealthy_error_rate': self.unhealthy_error_rate_threshold
                    }
                }
            )
    
    def is_healthy(self):
        """
        Check if the service is healthy.
        
        Returns:
            bool: True if the service is healthy, False otherwise
        """
        with self.lock:
            return self.status == HealthStatus.HEALTHY
    
    def is_degraded(self):
        """
        Check if the service is degraded.
        
        Returns:
            bool: True if the service is degraded, False otherwise
        """
        with self.lock:
            return self.status == HealthStatus.DEGRADED
    
    def is_unhealthy(self):
        """
        Check if the service is unhealthy.
        
        Returns:
            bool: True if the service is unhealthy, False otherwise
        """
        with self.lock:
            return self.status == HealthStatus.UNHEALTHY
    
    def get_status(self):
        """
        Get the current health status.
        
        Returns:
            HealthStatus: Current health status
        """
        with self.lock:
            return self.status
    
    def get_metrics(self):
        """
        Get service health metrics.
        
        Returns:
            dict: Service health metrics
        """
        with self.lock:
            # Get current statistics
            avg_response_time = self.performance_stats.get_average_response_time()
            p95_response_time = self.performance_stats.get_p95_response_time()
            error_rate = self.performance_stats.get_error_rate()
            success_rate = self.performance_stats.get_success_rate()
            
            return {
                'name': self.name,
                'status': self.status.name,
                'last_status_change': self.last_status_change,
                'statistics': {
                    'average_response_time': avg_response_time,
                    'p95_response_time': p95_response_time,
                    'error_rate': error_rate,
                    'success_rate': success_rate,
                    'success_count': self.performance_stats.success_count,
                    'failure_count': self.performance_stats.failure_count,
                    'total_count': self.performance_stats.total_count
                },
                'thresholds': {
                    'degraded_response_time': self.degraded_response_time_threshold,
                    'unhealthy_response_time': self.unhealthy_response_time_threshold,
                    'degraded_error_rate': self.degraded_error_rate_threshold,
                    'unhealthy_error_rate': self.unhealthy_error_rate_threshold
                }
            }
    
    def reset(self):
        """Reset service health to healthy status."""
        with self.lock:
            self.status = HealthStatus.HEALTHY
            self.last_status_change = time.time()
            self.performance_stats.reset()
            
            # Log reset
            log_with_context(
                logger, 'info',
                f"Service '{self.name}' health status reset to HEALTHY",
                context={
                    'event_type': 'service_health_reset',
                    'service': {
                        'name': self.name,
                        'status': self.status.name
                    }
                }
            )
            
            # Update metrics
            if HAS_OBSERVABILITY:
                SERVICE_HEALTH_STATUS.labels(service=self.name).set(self.status.value)

def track_service_health(service_name):
    """
    Decorator that tracks service health for a function.
    
    Args:
        service_name: Name of the service to track
        
    Returns:
        Decorated function
        
    The decorator will track the success, failure, and response time of each function call
    and update the service health status accordingly.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Get service health tracker
            service_health = ServiceHealth.get_service_health(service_name)
            
            # Record start time
            start_time = time.time()
            
            try:
                # Call the function
                result = func(*args, **kwargs)
                
                # Record success and response time
                response_time = time.time() - start_time
                service_health.on_success(response_time)
                
                return result
            
            except Exception as e:
                # Record failure
                service_health.on_failure()
                
                # Re-raise the exception
                raise
        
        return wrapper
    
    return decorator

def get_service_health(service_name):
    """
    Get the health status of a service.
    
    Args:
        service_name: Name of the service
        
    Returns:
        ServiceHealth: Service health tracker
    """
    return ServiceHealth.get_service_health(service_name)

def get_all_service_health():
    """
    Get the health status of all registered services.
    
    Returns:
        Dict of service health trackers by name
    """
    return ServiceHealth.get_all_service_health()