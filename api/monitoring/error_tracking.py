#!/usr/bin/env python3
"""
CryoProtect v2 - Centralized Error Tracking System

This module provides a comprehensive error tracking system with:
- Error aggregation and deduplication
- Context-aware error collection
- Error categorization, prioritization, and trending
- Alert thresholds for different error types
- Integration with structured logging system
- Real-time error metrics and dashboards
- Historical error analysis and reporting

Usage:
    from api.monitoring.error_tracking import (
        ErrorTracker,
        track_error,
        get_error_stats,
        ErrorCategory,
        ErrorSeverity
    )
"""

import os
import sys
import time
import json
import uuid
import enum
import copy
import threading
import logging
import traceback
import datetime
import hashlib
from functools import wraps
from typing import Dict, List, Any, Set, Optional, Union, Callable, Type, TypeVar, Tuple, NamedTuple

# Import structured logging
from logging_enhanced import log_with_context, get_logger

# Import Prometheus metrics if available
try:
    from monitoring.prometheus_metrics import ERROR_COUNT
    PROMETHEUS_AVAILABLE = True
except ImportError:
    PROMETHEUS_AVAILABLE = False

# Set up logger
logger = get_logger(__name__)

# Type variables
T = TypeVar('T')  # Return type for decorated functions


class ErrorSeverity(enum.Enum):
    """Error severity levels."""
    DEBUG = 0
    INFO = 1
    WARNING = 2
    ERROR = 3
    CRITICAL = 4


class ErrorCategory(enum.Enum):
    """Categories for error classification."""
    # Infrastructure errors
    NETWORK = "network"                  # Network connectivity issues
    DATABASE = "database"                # Database access issues
    SYSTEM = "system"                    # System-level issues (OS, hardware)
    RESOURCE = "resource"                # Resource exhaustion (memory, disk, CPU)
    
    # API errors
    API_CLIENT = "api_client"            # Client-side API issues
    API_SERVER = "api_server"            # Server-side API issues
    API_GATEWAY = "api_gateway"          # API gateway issues
    
    # External service errors
    SERVICE_DEPENDENCY = "service_dependency"  # External service dependency issues
    SERVICE_TIMEOUT = "service_timeout"        # External service timeout
    SERVICE_UNAVAILABLE = "service_unavailable"  # External service unavailable
    
    # Data errors
    DATA_VALIDATION = "data_validation"  # Data validation failures
    DATA_INTEGRITY = "data_integrity"    # Data integrity issues
    DATA_CORRUPTION = "data_corruption"  # Data corruption
    
    # Security errors
    AUTHENTICATION = "authentication"    # Authentication failures
    AUTHORIZATION = "authorization"      # Authorization failures
    SECURITY_POLICY = "security_policy"  # Security policy violations
    
    # Application errors
    BUSINESS_LOGIC = "business_logic"    # Business logic errors
    CONFIGURATION = "configuration"      # Configuration errors
    CODE = "code"                        # Code errors (bugs)
    
    # User-facing errors
    USER_INPUT = "user_input"            # User input errors
    USER_EXPERIENCE = "user_experience"  # User experience issues
    
    # Observability errors
    MONITORING = "monitoring"            # Monitoring system issues
    LOGGING = "logging"                  # Logging system issues
    METRICS = "metrics"                  # Metrics collection issues
    
    # Unknown/uncategorized
    UNKNOWN = "unknown"                  # Uncategorized errors


class ErrorFingerprint:
    """
    Creates a unique fingerprint for an error to enable deduplication.
    
    The fingerprint is based on the error type, message, and stacktrace,
    but ignores variable data like timestamps, request IDs, etc.
    """
    
    @staticmethod
    def generate(
        error_type: str,
        message: str,
        traceback_str: str,
        context: Optional[Dict[str, Any]] = None
    ) -> str:
        """
        Generate a fingerprint for an error.
        
        Args:
            error_type: Type of the error
            message: Error message
            traceback_str: Error traceback as string
            context: Error context
            
        Returns:
            Fingerprint string (SHA-256 hash)
        """
        # Normalize error message (remove variable parts)
        normalized_message = ErrorFingerprint._normalize_message(message)
        
        # Extract relevant parts of traceback (filenames, line numbers, function names)
        normalized_traceback = ErrorFingerprint._normalize_traceback(traceback_str)
        
        # Create fingerprint input
        fingerprint_input = f"{error_type}:{normalized_message}:{normalized_traceback}"
        
        # Generate SHA-256 hash
        return hashlib.sha256(fingerprint_input.encode('utf-8')).hexdigest()
    
    @staticmethod
    def _normalize_message(message: str) -> str:
        """
        Normalize error message by removing variable parts.
        
        Args:
            message: Error message
            
        Returns:
            Normalized message
        """
        # Remove UUIDs
        import re
        uuid_pattern = re.compile(
            r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}',
            re.IGNORECASE
        )
        message = uuid_pattern.sub('<uuid>', message)
        
        # Remove timestamps
        timestamp_pattern = re.compile(
            r'\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:Z|[+-]\d{2}:\d{2})?'
        )
        message = timestamp_pattern.sub('<timestamp>', message)
        
        # Remove IP addresses
        ip_pattern = re.compile(r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3}')
        message = ip_pattern.sub('<ip>', message)
        
        # Remove numbers that are likely to be variable
        # (but keep small numbers that might be part of error codes)
        number_pattern = re.compile(r'\b\d{3,}\b')
        message = number_pattern.sub('<number>', message)
        
        return message
    
    @staticmethod
    def _normalize_traceback(traceback_str: str) -> str:
        """
        Extract relevant parts from traceback.
        
        Args:
            traceback_str: Error traceback as string
            
        Returns:
            Normalized traceback
        """
        # Extract filenames, line numbers, and function names
        import re
        
        frames = []
        
        # Match "File "path/to/file.py", line X, in function_name"
        frame_pattern = re.compile(r'File "([^"]+)", line (\d+), in (\w+)')
        
        for match in frame_pattern.finditer(traceback_str):
            filename = match.group(1)
            # Get only the filename, not the full path
            filename = os.path.basename(filename)
            line = match.group(2)
            function = match.group(3)
            frames.append(f"{filename}:{line}:{function}")
        
        # Join frames to create normalized traceback
        return "|".join(frames)


class TrackedError:
    """
    Represents a tracked error instance with all relevant information.
    """
    
    def __init__(
        self,
        exception: Optional[Exception] = None,
        message: Optional[str] = None,
        category: ErrorCategory = ErrorCategory.UNKNOWN,
        severity: ErrorSeverity = ErrorSeverity.ERROR,
        fingerprint: Optional[str] = None,
        context: Optional[Dict[str, Any]] = None,
        correlation_id: Optional[str] = None,
        traceback_str: Optional[str] = None,
        timestamp: Optional[float] = None
    ):
        """
        Initialize a tracked error.
        
        Args:
            exception: The exception object
            message: Error message (required if exception is None)
            category: Error category
            severity: Error severity
            fingerprint: Error fingerprint (generated if None)
            context: Error context
            correlation_id: Correlation ID for distributed tracing
            traceback_str: Traceback as string
            timestamp: Error timestamp (defaults to current time)
        """
        # Set error information
        self.timestamp = timestamp or time.time()
        
        # Extract information from exception
        if exception is not None:
            self.exception_type = type(exception).__name__
            self.message = str(exception)
            if traceback_str is None:
                self.traceback = ''.join(traceback.format_exception(
                    type(exception),
                    exception,
                    exception.__traceback__
                ))
            else:
                self.traceback = traceback_str
        else:
            if message is None:
                raise ValueError("Either exception or message must be provided")
            self.exception_type = "Manual"
            self.message = message
            self.traceback = traceback_str or ""
        
        # Set categorization
        self.category = category
        self.severity = severity
        
        # Set context
        self.context = context or {}
        
        # Set correlation ID
        self.correlation_id = correlation_id or self.context.get('correlation_id', str(uuid.uuid4()))
        
        # Set error ID and fingerprint
        self.error_id = str(uuid.uuid4())
        
        # Generate fingerprint if not provided
        if fingerprint is None:
            self.fingerprint = ErrorFingerprint.generate(
                self.exception_type,
                self.message,
                self.traceback,
                self.context
            )
        else:
            self.fingerprint = fingerprint
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the tracked error to a dictionary.
        
        Returns:
            Dictionary representation of the error
        """
        return {
            'error_id': self.error_id,
            'fingerprint': self.fingerprint,
            'timestamp': self.timestamp,
            'exception_type': self.exception_type,
            'message': self.message,
            'category': self.category.value,
            'severity': self.severity.name,
            'correlation_id': self.correlation_id,
            'context': self.context,
            'traceback': self.traceback
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'TrackedError':
        """
        Create a TrackedError from a dictionary.
        
        Args:
            data: Dictionary representation of the error
            
        Returns:
            TrackedError instance
        """
        # Convert category and severity from string to enum
        category = ErrorCategory(data['category'])
        severity = ErrorSeverity[data['severity']]
        
        return cls(
            message=data['message'],
            category=category,
            severity=severity,
            fingerprint=data['fingerprint'],
            context=data['context'],
            correlation_id=data['correlation_id'],
            traceback_str=data['traceback'],
            timestamp=data['timestamp']
        )
    
    def to_json(self) -> str:
        """
        Convert the tracked error to a JSON string.
        
        Returns:
            JSON string representation of the error
        """
        return json.dumps(self.to_dict(), default=str)
    
    @classmethod
    def from_json(cls, json_str: str) -> 'TrackedError':
        """
        Create a TrackedError from a JSON string.
        
        Args:
            json_str: JSON string representation of the error
            
        Returns:
            TrackedError instance
        """
        data = json.loads(json_str)
        return cls.from_dict(data)


class ErrorGroup:
    """
    Represents a group of similar errors (based on fingerprint).
    """
    
    def __init__(self, first_error: TrackedError):
        """
        Initialize an error group with the first error.
        
        Args:
            first_error: First error in the group
        """
        self.fingerprint = first_error.fingerprint
        self.exception_type = first_error.exception_type
        self.category = first_error.category
        self.first_seen = first_error.timestamp
        self.last_seen = first_error.timestamp
        self.count = 1
        self.errors: List[TrackedError] = [first_error]
        self.sample_error_ids: List[str] = [first_error.error_id]
        self.max_samples = 10  # Maximum number of sample errors to keep
        self.severities: Dict[str, int] = {first_error.severity.name: 1}
        self.status = "active"  # active, resolved, muted
        self.resolution_data: Optional[Dict[str, Any]] = None
    
    def add_error(self, error: TrackedError) -> None:
        """
        Add an error to the group.
        
        Args:
            error: Error to add
        """
        # Update group statistics
        self.last_seen = error.timestamp
        self.count += 1
        
        # Update severities
        severity_name = error.severity.name
        if severity_name in self.severities:
            self.severities[severity_name] += 1
        else:
            self.severities[severity_name] = 1
        
        # Keep the error in the group (up to max_samples)
        self.errors.append(error)
        if len(self.errors) > self.max_samples:
            self.errors.pop(0)  # Remove oldest error
        
        # Keep track of sample error IDs
        self.sample_error_ids.append(error.error_id)
        if len(self.sample_error_ids) > self.max_samples:
            self.sample_error_ids.pop(0)  # Remove oldest error ID
    
    def get_sample_errors(self) -> List[TrackedError]:
        """
        Get sample errors from the group.
        
        Returns:
            List of sample errors
        """
        return list(self.errors)
    
    def get_dominant_severity(self) -> ErrorSeverity:
        """
        Get the dominant severity for the group.
        
        Returns:
            Dominant error severity
        """
        # Find the severity with the highest count
        dominant_severity_name = max(self.severities.items(), key=lambda x: x[1])[0]
        return ErrorSeverity[dominant_severity_name]
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the error group to a dictionary.
        
        Returns:
            Dictionary representation of the error group
        """
        return {
            'fingerprint': self.fingerprint,
            'exception_type': self.exception_type,
            'category': self.category.value,
            'first_seen': self.first_seen,
            'last_seen': self.last_seen,
            'count': self.count,
            'sample_error_ids': self.sample_error_ids,
            'severities': self.severities,
            'status': self.status,
            'resolution_data': self.resolution_data
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any], errors: Optional[List[TrackedError]] = None) -> 'ErrorGroup':
        """
        Create an ErrorGroup from a dictionary.
        
        Args:
            data: Dictionary representation of the error group
            errors: List of errors to include in the group
            
        Returns:
            ErrorGroup instance
        """
        if errors is None or not errors:
            # Create a dummy first error if none provided
            dummy_error = TrackedError(
                message=f"Placeholder for {data['exception_type']}",
                category=ErrorCategory(data['category']),
                fingerprint=data['fingerprint'],
                timestamp=data['first_seen']
            )
            group = cls(dummy_error)
        else:
            group = cls(errors[0])
            for error in errors[1:]:
                group.add_error(error)
        
        # Override properties from data
        group.fingerprint = data['fingerprint']
        group.exception_type = data['exception_type']
        group.category = ErrorCategory(data['category'])
        group.first_seen = data['first_seen']
        group.last_seen = data['last_seen']
        group.count = data['count']
        group.sample_error_ids = data['sample_error_ids']
        group.severities = data['severities']
        group.status = data['status']
        group.resolution_data = data['resolution_data']
        
        return group


class ErrorStats:
    """
    Statistics about tracked errors.
    """
    
    def __init__(self, time_window: float = 3600):
        """
        Initialize error statistics.
        
        Args:
            time_window: Time window for statistics in seconds (default: 1 hour)
        """
        self.time_window = time_window
        self.reset()
    
    def reset(self) -> None:
        """Reset statistics."""
        self.total_errors = 0
        self.unique_fingerprints = 0
        self.total_errors_by_category: Dict[str, int] = {}
        self.total_errors_by_severity: Dict[str, int] = {}
        self.errors_over_time: List[Tuple[float, int]] = []
        self.top_error_groups: List[Dict[str, Any]] = []
        self.new_error_groups: List[Dict[str, Any]] = []
        self.resolved_error_groups: List[Dict[str, Any]] = []
        self.recent_errors_per_min: int = 0
        self.error_rate_change: float = 0
    
    def update(self, tracker: 'ErrorTracker') -> None:
        """
        Update statistics from a tracker.
        
        Args:
            tracker: Error tracker to get statistics from
        """
        with tracker.lock:
            # Reset statistics
            self.reset()
            
            # Calculate time threshold for recent errors
            now = time.time()
            cutoff_time = now - self.time_window
            
            # Count total errors and unique fingerprints
            self.total_errors = sum(group.count for group in tracker.error_groups.values())
            self.unique_fingerprints = len(tracker.error_groups)
            
            # Count errors by category and severity
            for group in tracker.error_groups.values():
                category = group.category.value
                if category not in self.total_errors_by_category:
                    self.total_errors_by_category[category] = 0
                self.total_errors_by_category[category] += group.count
                
                for severity, count in group.severities.items():
                    if severity not in self.total_errors_by_severity:
                        self.total_errors_by_severity[severity] = 0
                    self.total_errors_by_severity[severity] += count
            
            # Get errors over time (1-minute intervals)
            error_times = []
            for group in tracker.error_groups.values():
                for error in group.get_sample_errors():
                    if error.timestamp >= cutoff_time:
                        error_times.append(error.timestamp)
            
            if error_times:
                # Group errors by minute
                minute_bins: Dict[int, int] = {}
                for timestamp in error_times:
                    minute = int(timestamp / 60)
                    if minute not in minute_bins:
                        minute_bins[minute] = 0
                    minute_bins[minute] += 1
                
                # Convert to list of (timestamp, count) tuples
                self.errors_over_time = [(min_time * 60, count) for min_time, count in minute_bins.items()]
                self.errors_over_time.sort()
                
                # Calculate recent errors per minute
                recent_window_start = now - 300  # Last 5 minutes
                recent_errors = sum(count for timestamp, count in self.errors_over_time 
                                   if timestamp >= recent_window_start)
                self.recent_errors_per_min = recent_errors / 5
                
                # Calculate error rate change (comparing last 5 min to previous 5 min)
                prev_window_start = recent_window_start - 300
                prev_errors = sum(count for timestamp, count in self.errors_over_time 
                                 if timestamp >= prev_window_start and timestamp < recent_window_start)
                
                if prev_errors > 0:
                    self.error_rate_change = (recent_errors - prev_errors) / prev_errors
                else:
                    self.error_rate_change = float('inf') if recent_errors > 0 else 0
            
            # Get top error groups
            sorted_groups = sorted(
                tracker.error_groups.values(),
                key=lambda g: g.count,
                reverse=True
            )
            
            self.top_error_groups = [group.to_dict() for group in sorted_groups[:10]]
            
            # Get new error groups (first seen in the time window)
            new_groups = [
                group for group in tracker.error_groups.values()
                if group.first_seen >= cutoff_time
            ]
            
            self.new_error_groups = [group.to_dict() for group in new_groups]
            
            # Get recently resolved error groups
            resolved_groups = [
                group for group in tracker.error_groups.values()
                if group.status == 'resolved' and 
                (group.resolution_data and group.resolution_data.get('resolved_at', 0) >= cutoff_time)
            ]
            
            self.resolved_error_groups = [group.to_dict() for group in resolved_groups]
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the statistics to a dictionary.
        
        Returns:
            Dictionary representation of the statistics
        """
        return {
            'time_window': self.time_window,
            'total_errors': self.total_errors,
            'unique_fingerprints': self.unique_fingerprints,
            'total_errors_by_category': self.total_errors_by_category,
            'total_errors_by_severity': self.total_errors_by_severity,
            'errors_over_time': self.errors_over_time,
            'top_error_groups': self.top_error_groups,
            'new_error_groups': self.new_error_groups,
            'resolved_error_groups': self.resolved_error_groups,
            'recent_errors_per_min': self.recent_errors_per_min,
            'error_rate_change': self.error_rate_change
        }


class AlertConfig:
    """
    Configuration for error alerts.
    """
    
    def __init__(
        self,
        error_count_threshold: int = 10,
        error_rate_threshold: float = 0.5,  # 50% increase
        time_window: float = 300,  # 5 minutes
        cooldown_period: float = 3600,  # 1 hour
        severity_thresholds: Optional[Dict[str, int]] = None,
        category_thresholds: Optional[Dict[str, int]] = None
    ):
        """
        Initialize alert configuration.
        
        Args:
            error_count_threshold: Minimum error count to trigger an alert
            error_rate_threshold: Minimum error rate change to trigger an alert
            time_window: Time window for error rate calculation in seconds
            cooldown_period: Cooldown period between alerts for the same error group
            severity_thresholds: Custom thresholds by severity
            category_thresholds: Custom thresholds by category
        """
        self.error_count_threshold = error_count_threshold
        self.error_rate_threshold = error_rate_threshold
        self.time_window = time_window
        self.cooldown_period = cooldown_period
        
        # Default severity thresholds
        self.severity_thresholds = {
            'DEBUG': 100,
            'INFO': 50,
            'WARNING': 20,
            'ERROR': 10,
            'CRITICAL': 1
        }
        
        # Apply custom severity thresholds
        if severity_thresholds:
            self.severity_thresholds.update(severity_thresholds)
        
        # Default category thresholds
        self.category_thresholds = {}
        for category in ErrorCategory:
            self.category_thresholds[category.value] = error_count_threshold
        
        # Apply custom category thresholds
        if category_thresholds:
            self.category_thresholds.update(category_thresholds)
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the alert configuration to a dictionary.
        
        Returns:
            Dictionary representation of the alert configuration
        """
        return {
            'error_count_threshold': self.error_count_threshold,
            'error_rate_threshold': self.error_rate_threshold,
            'time_window': self.time_window,
            'cooldown_period': self.cooldown_period,
            'severity_thresholds': self.severity_thresholds,
            'category_thresholds': self.category_thresholds
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'AlertConfig':
        """
        Create an AlertConfig from a dictionary.
        
        Args:
            data: Dictionary representation of the alert configuration
            
        Returns:
            AlertConfig instance
        """
        return cls(
            error_count_threshold=data['error_count_threshold'],
            error_rate_threshold=data['error_rate_threshold'],
            time_window=data['time_window'],
            cooldown_period=data['cooldown_period'],
            severity_thresholds=data['severity_thresholds'],
            category_thresholds=data['category_thresholds']
        )


class Alert:
    """
    Represents an error alert.
    """
    
    def __init__(
        self,
        error_group: ErrorGroup,
        alert_reason: str,
        timestamp: Optional[float] = None
    ):
        """
        Initialize an alert.
        
        Args:
            error_group: Error group that triggered the alert
            alert_reason: Reason for the alert
            timestamp: Alert timestamp (defaults to current time)
        """
        self.alert_id = str(uuid.uuid4())
        self.fingerprint = error_group.fingerprint
        self.error_count = error_group.count
        self.exception_type = error_group.exception_type
        self.category = error_group.category
        self.dominant_severity = error_group.get_dominant_severity()
        self.alert_reason = alert_reason
        self.timestamp = timestamp or time.time()
        self.sample_error_ids = error_group.sample_error_ids.copy()
        self.status = "open"  # open, acknowledged, resolved
        self.acknowledged_at: Optional[float] = None
        self.resolved_at: Optional[float] = None
    
    def acknowledge(self) -> None:
        """Acknowledge the alert."""
        self.status = "acknowledged"
        self.acknowledged_at = time.time()
    
    def resolve(self) -> None:
        """Resolve the alert."""
        self.status = "resolved"
        self.resolved_at = time.time()
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the alert to a dictionary.
        
        Returns:
            Dictionary representation of the alert
        """
        return {
            'alert_id': self.alert_id,
            'fingerprint': self.fingerprint,
            'error_count': self.error_count,
            'exception_type': self.exception_type,
            'category': self.category.value,
            'dominant_severity': self.dominant_severity.name,
            'alert_reason': self.alert_reason,
            'timestamp': self.timestamp,
            'sample_error_ids': self.sample_error_ids,
            'status': self.status,
            'acknowledged_at': self.acknowledged_at,
            'resolved_at': self.resolved_at
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Alert':
        """
        Create an Alert from a dictionary.
        
        Args:
            data: Dictionary representation of the alert
            
        Returns:
            Alert instance
        """
        # Create a dummy error group for the alert
        dummy_group = ErrorGroup(
            TrackedError(
                message=f"Placeholder for {data['exception_type']}",
                category=ErrorCategory(data['category']),
                fingerprint=data['fingerprint']
            )
        )
        dummy_group.count = data['error_count']
        dummy_group.sample_error_ids = data['sample_error_ids']
        
        # Create alert
        alert = cls(
            error_group=dummy_group,
            alert_reason=data['alert_reason'],
            timestamp=data['timestamp']
        )
        
        # Set alert properties
        alert.alert_id = data['alert_id']
        alert.status = data['status']
        alert.acknowledged_at = data['acknowledged_at']
        alert.resolved_at = data['resolved_at']
        
        return alert


class ErrorTracker:
    """
    Central error tracking system that collects, aggregates, and analyzes errors.
    """
    
    _instance = None
    
    @classmethod
    def get_instance(cls) -> 'ErrorTracker':
        """
        Get or create the singleton instance of ErrorTracker.
        
        Returns:
            ErrorTracker instance
        """
        if cls._instance is None:
            cls._instance = ErrorTracker()
        return cls._instance
    
    def __init__(self):
        """Initialize error tracker."""
        # Ensure singleton pattern
        if ErrorTracker._instance is not None:
            raise RuntimeError("ErrorTracker is a singleton. Use get_instance() instead.")
        
        ErrorTracker._instance = self
        
        # Error storage
        self.errors: Dict[str, TrackedError] = {}  # error_id -> TrackedError
        self.error_groups: Dict[str, ErrorGroup] = {}  # fingerprint -> ErrorGroup
        
        # Alert system
        self.alerts: Dict[str, Alert] = {}  # alert_id -> Alert
        self.last_alert_times: Dict[str, float] = {}  # fingerprint -> timestamp
        self.alert_config = AlertConfig()
        
        # Alert callbacks
        self.alert_callbacks: List[Callable[[Alert], None]] = []
        
        # Statistics
        self.stats = ErrorStats()
        
        # Storage limits
        self.max_stored_errors = 10000
        self.max_stored_groups = 1000
        self.max_stored_alerts = 1000
        
        # Persistence
        self.persistence_enabled = False
        self.persistence_path = "monitoring/error_tracking"
        
        # Thread safety
        self.lock = threading.RLock()
        
        # Create logger
        self.logger = get_logger(__name__)
    
    def track_error(
        self,
        error: Union[Exception, str],
        category: Optional[ErrorCategory] = None,
        severity: Optional[ErrorSeverity] = None,
        context: Optional[Dict[str, Any]] = None,
        correlation_id: Optional[str] = None
    ) -> TrackedError:
        """
        Track an error.
        
        Args:
            error: Exception or error message
            category: Error category (auto-detected if None)
            severity: Error severity (auto-detected if None)
            context: Error context
            correlation_id: Correlation ID for distributed tracing
            
        Returns:
            Tracked error
        """
        with self.lock:
            # Normalize error
            if isinstance(error, str):
                exception = None
                message = error
                if category is None:
                    category = ErrorCategory.UNKNOWN
                if severity is None:
                    severity = ErrorSeverity.ERROR
            else:
                exception = error
                message = None
                if category is None:
                    category = self._detect_category(exception)
                if severity is None:
                    severity = self._detect_severity(exception)
            
            # Create tracked error
            tracked_error = TrackedError(
                exception=exception,
                message=message,
                category=category,
                severity=severity,
                context=context,
                correlation_id=correlation_id
            )
            
            # Store the error
            self.errors[tracked_error.error_id] = tracked_error
            
            # Add to or create error group
            if tracked_error.fingerprint in self.error_groups:
                # Add to existing group
                group = self.error_groups[tracked_error.fingerprint]
                group.add_error(tracked_error)
            else:
                # Create new group
                group = ErrorGroup(tracked_error)
                self.error_groups[tracked_error.fingerprint] = group
            
            # Check for alerts
            self._check_for_alerts(group)
            
            # Update Prometheus metrics if available
            if PROMETHEUS_AVAILABLE:
                ERROR_COUNT.labels(
                    type=tracked_error.exception_type,
                    endpoint=tracked_error.context.get('endpoint', 'unknown')
                ).inc()
            
            # Log the error
            self._log_error(tracked_error)
            
            # Enforce storage limits
            self._enforce_storage_limits()
            
            # Save to persistent storage if enabled
            if self.persistence_enabled:
                self._save_error(tracked_error)
            
            return tracked_error
    
    def _detect_category(self, exception: Exception) -> ErrorCategory:
        """
        Detect the error category based on the exception type.
        
        Args:
            exception: Exception to categorize
            
        Returns:
            Error category
        """
        import socket
        import requests
        import sqlite3
        import json
        import ssl
        
        exception_type = type(exception)
        exception_name = exception_type.__name__
        message = str(exception).lower()
        
        # Network errors
        if (
            isinstance(exception, (socket.error, socket.gaierror, socket.timeout)) or
            exception_name in ('ConnectionError', 'ConnectionRefusedError', 'ConnectionResetError')
        ):
            return ErrorCategory.NETWORK
        
        # API errors
        if (
            exception_name in ('HTTPError', 'RequestException', 'Timeout', 'ReadTimeout', 'ConnectTimeout') or
            ('api' in message and ('failed' in message or 'error' in message))
        ):
            return ErrorCategory.API_CLIENT
        
        # Database errors
        if (
            isinstance(exception, sqlite3.Error) or
            exception_name in ('DatabaseError', 'OperationalError', 'IntegrityError')
        ):
            return ErrorCategory.DATABASE
        
        # Service dependency errors
        if (
            'service' in message and ('unavailable' in message or 'failed' in message) or
            exception_name in ('ServiceUnavailableError', 'ServiceError')
        ):
            return ErrorCategory.SERVICE_DEPENDENCY
        
        # Timeout errors
        if (
            exception_name in ('TimeoutError', 'DeadlineExceededError') or
            'timeout' in message or 'timed out' in message
        ):
            return ErrorCategory.SERVICE_TIMEOUT
        
        # Data validation errors
        if (
            exception_name in ('ValidationError', 'ValueError', 'TypeError') or
            ('validation' in message or 'invalid' in message)
        ):
            return ErrorCategory.DATA_VALIDATION
        
        # Data integrity errors
        if (
            exception_name in ('IntegrityError', 'DataError') or
            ('integrity' in message or 'constraint' in message)
        ):
            return ErrorCategory.DATA_INTEGRITY
        
        # Authentication/Authorization errors
        if (
            exception_name in ('AuthenticationError', 'NotAuthenticated') or
            'authentication' in message or 'not authenticated' in message
        ):
            return ErrorCategory.AUTHENTICATION
        
        if (
            exception_name in ('PermissionError', 'NotAuthorized', 'Forbidden') or
            'permission' in message or 'not authorized' in message or 'forbidden' in message
        ):
            return ErrorCategory.AUTHORIZATION
        
        # Configuration errors
        if (
            exception_name in ('ConfigurationError', 'ImproperlyConfigured') or
            'configuration' in message or 'config' in message
        ):
            return ErrorCategory.CONFIGURATION
        
        # Default to unknown category
        return ErrorCategory.UNKNOWN
    
    def _detect_severity(self, exception: Exception) -> ErrorSeverity:
        """
        Detect the error severity based on the exception type.
        
        Args:
            exception: Exception to assess
            
        Returns:
            Error severity
        """
        exception_type = type(exception)
        exception_name = exception_type.__name__
        message = str(exception).lower()
        
        # Critical errors
        if (
            exception_name in ('SystemExit', 'KeyboardInterrupt', 'MemoryError', 'SystemError') or
            'fatal' in message or 'critical' in message
        ):
            return ErrorSeverity.CRITICAL
        
        # Error level exceptions (most common)
        if (
            exception_name in ('Exception', 'RuntimeError', 'ValueError', 'IOError') or
            'error' in message or 'exception' in message or 'failed' in message
        ):
            return ErrorSeverity.ERROR
        
        # Warning level exceptions
        if (
            exception_name in ('Warning', 'DeprecationWarning', 'PendingDeprecationWarning') or
            'warning' in message or 'deprecated' in message
        ):
            return ErrorSeverity.WARNING
        
        # Default to ERROR severity
        return ErrorSeverity.ERROR
    
    def _log_error(self, error: TrackedError) -> None:
        """
        Log a tracked error.
        
        Args:
            error: Error to log
        """
        # Map severity to log level
        severity_to_level = {
            ErrorSeverity.DEBUG: logging.DEBUG,
            ErrorSeverity.INFO: logging.INFO,
            ErrorSeverity.WARNING: logging.WARNING,
            ErrorSeverity.ERROR: logging.ERROR,
            ErrorSeverity.CRITICAL: logging.CRITICAL
        }
        
        log_level = severity_to_level.get(error.severity, logging.ERROR)
        
        # Create log message
        log_message = f"{error.category.value.upper()}: {error.exception_type}: {error.message}"
        
        # Log using structured logging
        log_with_context(
            self.logger,
            logging.getLevelName(log_level).lower(),
            log_message,
            context={
                'error_tracking': {
                    'error_id': error.error_id,
                    'fingerprint': error.fingerprint,
                    'category': error.category.value,
                    'severity': error.severity.name
                },
                'correlation_id': error.correlation_id,
                **error.context
            }
        )
    
    def _check_for_alerts(self, group: ErrorGroup) -> None:
        """
        Check if an error group should trigger an alert.
        
        Args:
            group: Error group to check
        """
        now = time.time()
        fingerprint = group.fingerprint
        
        # Check cooldown period
        if (
            fingerprint in self.last_alert_times and
            now - self.last_alert_times[fingerprint] < self.alert_config.cooldown_period
        ):
            return
        
        # Check count threshold by category
        category_threshold = self.alert_config.category_thresholds.get(
            group.category.value,
            self.alert_config.error_count_threshold
        )
        
        # Check count threshold by severity
        dominant_severity = group.get_dominant_severity()
        severity_threshold = self.alert_config.severity_thresholds.get(
            dominant_severity.name,
            self.alert_config.error_count_threshold
        )
        
        # Use the lower of the two thresholds
        threshold = min(category_threshold, severity_threshold)
        
        # Check if the error count exceeds the threshold
        if group.count >= threshold:
            # Calculate the error rate (errors per minute)
            time_window = self.alert_config.time_window
            cutoff_time = now - time_window
            
            errors_in_window = sum(1 for error in group.get_sample_errors() if error.timestamp >= cutoff_time)
            error_rate = errors_in_window / (time_window / 60)  # errors per minute
            
            # Create alert
            alert_reason = f"Error count {group.count} exceeds threshold {threshold}"
            alert = Alert(group, alert_reason)
            
            # Add the alert
            self.alerts[alert.alert_id] = alert
            
            # Update last alert time
            self.last_alert_times[fingerprint] = now
            
            # Notify alert callbacks
            for callback in self.alert_callbacks:
                try:
                    callback(alert)
                except Exception as e:
                    self.logger.error(f"Error in alert callback: {str(e)}")
    
    def _enforce_storage_limits(self) -> None:
        """Enforce storage limits for errors, groups, and alerts."""
        # Limit errors
        if len(self.errors) > self.max_stored_errors:
            # Sort errors by timestamp (oldest first)
            sorted_errors = sorted(
                self.errors.items(),
                key=lambda x: x[1].timestamp
            )
            
            # Remove oldest errors
            to_remove = sorted_errors[:len(sorted_errors) - self.max_stored_errors]
            for error_id, _ in to_remove:
                del self.errors[error_id]
        
        # Limit error groups
        if len(self.error_groups) > self.max_stored_groups:
            # Sort groups by last_seen (oldest first)
            sorted_groups = sorted(
                self.error_groups.items(),
                key=lambda x: x[1].last_seen
            )
            
            # Remove oldest groups
            to_remove = sorted_groups[:len(sorted_groups) - self.max_stored_groups]
            for fingerprint, _ in to_remove:
                del self.error_groups[fingerprint]
        
        # Limit alerts
        if len(self.alerts) > self.max_stored_alerts:
            # Sort alerts by timestamp (oldest first)
            sorted_alerts = sorted(
                self.alerts.items(),
                key=lambda x: x[1].timestamp
            )
            
            # Remove oldest alerts
            to_remove = sorted_alerts[:len(sorted_alerts) - self.max_stored_alerts]
            for alert_id, _ in to_remove:
                del self.alerts[alert_id]
    
    def get_error(self, error_id: str) -> Optional[TrackedError]:
        """
        Get a tracked error by ID.
        
        Args:
            error_id: Error ID
            
        Returns:
            TrackedError instance or None if not found
        """
        with self.lock:
            return self.errors.get(error_id)
    
    def get_error_group(self, fingerprint: str) -> Optional[ErrorGroup]:
        """
        Get an error group by fingerprint.
        
        Args:
            fingerprint: Error fingerprint
            
        Returns:
            ErrorGroup instance or None if not found
        """
        with self.lock:
            return self.error_groups.get(fingerprint)
    
    def get_alert(self, alert_id: str) -> Optional[Alert]:
        """
        Get an alert by ID.
        
        Args:
            alert_id: Alert ID
            
        Returns:
            Alert instance or None if not found
        """
        with self.lock:
            return self.alerts.get(alert_id)
    
    def get_stats(self, time_window: Optional[float] = None) -> Dict[str, Any]:
        """
        Get error statistics.
        
        Args:
            time_window: Time window for statistics in seconds
            
        Returns:
            Dictionary with error statistics
        """
        with self.lock:
            if time_window is not None:
                stats = ErrorStats(time_window)
            else:
                stats = self.stats
            
            stats.update(self)
            return stats.to_dict()
    
    def add_alert_callback(self, callback: Callable[[Alert], None]) -> None:
        """
        Add a callback function for alerts.
        
        Args:
            callback: Function to call when an alert is triggered
        """
        with self.lock:
            self.alert_callbacks.append(callback)
    
    def resolve_error_group(
        self,
        fingerprint: str,
        resolution_data: Optional[Dict[str, Any]] = None
    ) -> bool:
        """
        Mark an error group as resolved.
        
        Args:
            fingerprint: Error fingerprint
            resolution_data: Optional resolution data
            
        Returns:
            True if the group was found and resolved, False otherwise
        """
        with self.lock:
            if fingerprint in self.error_groups:
                group = self.error_groups[fingerprint]
                group.status = "resolved"
                group.resolution_data = resolution_data or {
                    "resolved_at": time.time(),
                    "resolution_type": "manual"
                }
                
                # Also resolve any open alerts for this group
                for alert in self.alerts.values():
                    if alert.fingerprint == fingerprint and alert.status != "resolved":
                        alert.resolve()
                
                return True
            
            return False
    
    def mute_error_group(
        self,
        fingerprint: str,
        mute_data: Optional[Dict[str, Any]] = None
    ) -> bool:
        """
        Mute an error group (stop alerting but keep tracking).
        
        Args:
            fingerprint: Error fingerprint
            mute_data: Optional muting data
            
        Returns:
            True if the group was found and muted, False otherwise
        """
        with self.lock:
            if fingerprint in self.error_groups:
                group = self.error_groups[fingerprint]
                group.status = "muted"
                
                if mute_data is None:
                    mute_data = {
                        "muted_at": time.time(),
                        "mute_type": "manual"
                    }
                
                if group.resolution_data is None:
                    group.resolution_data = {}
                
                group.resolution_data.update(mute_data)
                
                return True
            
            return False
    
    def acknowledge_alert(self, alert_id: str) -> bool:
        """
        Acknowledge an alert.
        
        Args:
            alert_id: Alert ID
            
        Returns:
            True if the alert was found and acknowledged, False otherwise
        """
        with self.lock:
            if alert_id in self.alerts:
                alert = self.alerts[alert_id]
                if alert.status == "open":
                    alert.acknowledge()
                    return True
            
            return False
    
    def resolve_alert(self, alert_id: str) -> bool:
        """
        Resolve an alert.
        
        Args:
            alert_id: Alert ID
            
        Returns:
            True if the alert was found and resolved, False otherwise
        """
        with self.lock:
            if alert_id in self.alerts:
                alert = self.alerts[alert_id]
                if alert.status != "resolved":
                    alert.resolve()
                    return True
            
            return False
    
    def configure_alerts(self, config: AlertConfig) -> None:
        """
        Configure the alert system.
        
        Args:
            config: Alert configuration
        """
        with self.lock:
            self.alert_config = config
    
    def _save_error(self, error: TrackedError) -> None:
        """
        Save an error to persistent storage.
        
        Args:
            error: Error to save
        """
        try:
            # Ensure directory exists
            os.makedirs(self.persistence_path, exist_ok=True)
            
            # Save error
            error_path = os.path.join(self.persistence_path, f"error_{error.error_id}.json")
            with open(error_path, 'w') as f:
                f.write(error.to_json())
        except Exception as e:
            self.logger.error(f"Error saving tracked error: {str(e)}")
    
    def enable_persistence(self, path: Optional[str] = None) -> None:
        """
        Enable persistence of error data.
        
        Args:
            path: Directory path for persistence
        """
        with self.lock:
            self.persistence_enabled = True
            if path:
                self.persistence_path = path
            
            # Ensure directory exists
            os.makedirs(self.persistence_path, exist_ok=True)
    
    def disable_persistence(self) -> None:
        """Disable persistence of error data."""
        with self.lock:
            self.persistence_enabled = False
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the tracker to a dictionary.
        
        Returns:
            Dictionary representation of the tracker
        """
        with self.lock:
            return {
                'errors': {error_id: error.to_dict() for error_id, error in self.errors.items()},
                'error_groups': {fingerprint: group.to_dict() for fingerprint, group in self.error_groups.items()},
                'alerts': {alert_id: alert.to_dict() for alert_id, alert in self.alerts.items()},
                'last_alert_times': self.last_alert_times,
                'alert_config': self.alert_config.to_dict(),
                'max_stored_errors': self.max_stored_errors,
                'max_stored_groups': self.max_stored_groups,
                'max_stored_alerts': self.max_stored_alerts,
                'persistence_enabled': self.persistence_enabled,
                'persistence_path': self.persistence_path
            }


# Global error tracker instance
tracker = ErrorTracker.get_instance()


def track_error(
    error: Union[Exception, str],
    category: Optional[ErrorCategory] = None,
    severity: Optional[ErrorSeverity] = None,
    context: Optional[Dict[str, Any]] = None,
    correlation_id: Optional[str] = None
) -> TrackedError:
    """
    Track an error using the global error tracker.
    
    Args:
        error: Exception or error message
        category: Error category (auto-detected if None)
        severity: Error severity (auto-detected if None)
        context: Error context
        correlation_id: Correlation ID for distributed tracing
        
    Returns:
        Tracked error
    """
    return tracker.track_error(
        error=error,
        category=category,
        severity=severity,
        context=context,
        correlation_id=correlation_id
    )


def get_error_stats(time_window: Optional[float] = None) -> Dict[str, Any]:
    """
    Get error statistics from the global error tracker.
    
    Args:
        time_window: Time window for statistics in seconds
        
    Returns:
        Dictionary with error statistics
    """
    return tracker.get_stats(time_window)


def configure_alerts(config: AlertConfig) -> None:
    """
    Configure the alert system of the global error tracker.
    
    Args:
        config: Alert configuration
    """
    tracker.configure_alerts(config)


def add_alert_callback(callback: Callable[[Alert], None]) -> None:
    """
    Add a callback function for alerts to the global error tracker.
    
    Args:
        callback: Function to call when an alert is triggered
    """
    tracker.add_alert_callback(callback)


def resolve_error_group(
    fingerprint: str,
    resolution_data: Optional[Dict[str, Any]] = None
) -> bool:
    """
    Mark an error group as resolved in the global error tracker.
    
    Args:
        fingerprint: Error fingerprint
        resolution_data: Optional resolution data
        
    Returns:
        True if the group was found and resolved, False otherwise
    """
    return tracker.resolve_error_group(fingerprint, resolution_data)


def enable_persistence(path: Optional[str] = None) -> None:
    """
    Enable persistence of error data in the global error tracker.
    
    Args:
        path: Directory path for persistence
    """
    tracker.enable_persistence(path)


def with_error_tracking(
    category: Optional[ErrorCategory] = None,
    severity: Optional[ErrorSeverity] = None
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator to track errors in a function.
    
    Args:
        category: Error category (auto-detected if None)
        severity: Error severity (auto-detected if None)
        
    Returns:
        Decorator function
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def wrapper(*args, **kwargs) -> T:
            try:
                return func(*args, **kwargs)
            except Exception as e:
                # Create context with function information
                context = {
                    'function': func.__name__,
                    'module': func.__module__,
                    'args': str(args),
                    'kwargs': str(kwargs)
                }
                
                # Track the error
                track_error(
                    error=e,
                    category=category,
                    severity=severity,
                    context=context
                )
                
                # Re-raise the exception
                raise
        
        return wrapper
    
    return decorator