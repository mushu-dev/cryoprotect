#!/usr/bin/env python3
"""
CryoProtect v2 - Performance Monitoring System

This module provides comprehensive performance monitoring capabilities:
- Application-level performance tracking
- Database query performance monitoring
- Resource usage tracking (CPU, memory, disk)
- Response time analysis and profiling
- Performance bottleneck detection
- Long-running operation tracking
- Integration with Prometheus metrics
- Real-time performance dashboards

Usage:
    from api.monitoring.performance import (
        track_execution_time,
        track_resource_usage,
        monitor_query,
        start_performance_tracking,
        get_performance_stats,
        PerformanceTracker
    )
"""

import os
import gc
import time
import json
import uuid
import enum
import cProfile
import pstats
import io
import functools
import threading
import logging
import psutil
import traceback
from typing import Dict, List, Any, Set, Optional, Union, Callable, Type, TypeVar, Tuple

# Import structured logging
from logging_enhanced import log_with_context, get_logger

# Import Prometheus metrics if available
try:
    from monitoring.prometheus_metrics import (
        REQUEST_LATENCY,
        DB_QUERY_LATENCY,
        CPU_USAGE,
        MEMORY_USAGE,
        DISK_USAGE
    )
    PROMETHEUS_AVAILABLE = True
except ImportError:
    PROMETHEUS_AVAILABLE = False

# Set up logger
logger = get_logger(__name__)

# Type variables
T = TypeVar('T')  # Return type for decorated functions


class PerformanceCategory(enum.Enum):
    """Performance monitoring categories."""
    API = "api"                  # API endpoint performance
    DATABASE = "database"        # Database operation performance
    BACKGROUND = "background"    # Background task performance
    RESOURCE = "resource"        # Resource usage
    COMPUTATION = "computation"  # Computation/algorithm performance
    EXTERNAL = "external"        # External service performance
    SYSTEM = "system"            # System-level performance
    APPLICATION = "application"  # Application-level performance


class ResourceType(enum.Enum):
    """Resource types for monitoring."""
    CPU = "cpu"          # CPU usage
    MEMORY = "memory"    # Memory usage
    DISK = "disk"        # Disk usage
    NETWORK = "network"  # Network usage
    THREAD = "thread"    # Thread usage
    CONNECTION = "connection"  # Connection pool usage


class PerformanceThreshold:
    """
    Performance threshold configuration for generating alerts.
    """
    
    def __init__(
        self,
        category: PerformanceCategory,
        warning_threshold: float,
        critical_threshold: float,
        name: Optional[str] = None
    ):
        """
        Initialize performance threshold.
        
        Args:
            category: Performance category
            warning_threshold: Warning threshold in seconds
            critical_threshold: Critical threshold in seconds
            name: Optional name for the threshold (e.g., specific operation)
        """
        self.category = category
        self.warning_threshold = warning_threshold
        self.critical_threshold = critical_threshold
        self.name = name
    
    def check(self, duration: float) -> Optional[str]:
        """
        Check if a duration exceeds thresholds.
        
        Args:
            duration: Duration to check in seconds
            
        Returns:
            Level of threshold exceeded ('critical', 'warning', or None)
        """
        if duration >= self.critical_threshold:
            return 'critical'
        elif duration >= self.warning_threshold:
            return 'warning'
        return None


class ExecutionContext:
    """
    Context information for performance tracking.
    """
    
    def __init__(
        self,
        category: PerformanceCategory,
        name: str,
        context: Optional[Dict[str, Any]] = None,
        correlation_id: Optional[str] = None
    ):
        """
        Initialize execution context.
        
        Args:
            category: Performance category
            name: Operation name
            context: Additional context information
            correlation_id: Correlation ID for distributed tracing
        """
        self.category = category
        self.name = name
        self.context = context or {}
        self.correlation_id = correlation_id or str(uuid.uuid4())
        self.start_time = time.time()
        self.end_time: Optional[float] = None
        self.duration: Optional[float] = None
        self.cpu_usage_start = psutil.Process(os.getpid()).cpu_percent(interval=None)
        self.memory_usage_start = psutil.Process(os.getpid()).memory_info().rss
        self.thread_id = threading.get_ident()
        self.process_id = os.getpid()
        
        # Store current stack trace
        self.stack_trace = ''.join(traceback.format_stack(limit=15))
        
        # For database operations
        if category == PerformanceCategory.DATABASE:
            self.query = context.get('query', '') if context else ''
            self.params = context.get('params', None) if context else None
            self.db_name = context.get('db_name', '') if context else ''
        
        # Profile data (if profiling is enabled)
        self.profile_data = None
    
    def complete(self) -> None:
        """Complete the execution context and calculate the duration."""
        self.end_time = time.time()
        self.duration = self.end_time - self.start_time
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the execution context to a dictionary.
        
        Returns:
            Dictionary representation of the execution context
        """
        result = {
            'category': self.category.value,
            'name': self.name,
            'correlation_id': self.correlation_id,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'duration': self.duration,
            'thread_id': self.thread_id,
            'process_id': self.process_id,
            'context': self.context
        }
        
        if self.category == PerformanceCategory.DATABASE:
            result.update({
                'query': self.query,
                'db_name': self.db_name
            })
            
            # Only include params if they can be serialized
            if self.params:
                try:
                    json.dumps({'params': self.params})
                    result['params'] = self.params
                except (TypeError, OverflowError):
                    result['params'] = str(self.params)
        
        if self.profile_data:
            result['profile_data'] = self.profile_data
        
        return result


class PerformanceTracker:
    """
    Performance tracking system that monitors and analyzes execution times.
    """
    
    _instance = None
    
    @classmethod
    def get_instance(cls) -> 'PerformanceTracker':
        """
        Get or create the singleton instance of PerformanceTracker.
        
        Returns:
            PerformanceTracker instance
        """
        if cls._instance is None:
            cls._instance = PerformanceTracker()
        return cls._instance
    
    def __init__(self):
        """Initialize performance tracker."""
        # Ensure singleton pattern
        if PerformanceTracker._instance is not None:
            raise RuntimeError("PerformanceTracker is a singleton. Use get_instance() instead.")
        
        PerformanceTracker._instance = self
        
        # Performance data storage
        self.executions: Dict[str, ExecutionContext] = {}  # correlation_id -> ExecutionContext
        self.active_executions: Dict[str, ExecutionContext] = {}  # correlation_id -> ExecutionContext
        self.recent_executions: List[ExecutionContext] = []
        
        # Thresholds
        self.thresholds: Dict[Tuple[PerformanceCategory, Optional[str]], PerformanceThreshold] = {}
        
        # Statistics
        self.stats_by_category: Dict[str, Dict[str, Any]] = {}
        self.stats_by_name: Dict[str, Dict[str, Any]] = {}
        self.bottlenecks: List[Dict[str, Any]] = []
        
        # Resource usage tracking
        self.resource_snapshots: List[Dict[str, Any]] = []
        self.resource_tracking_interval = 60  # seconds
        self.last_resource_snapshot = 0
        
        # Performance threshold callbacks
        self.threshold_callbacks: List[Callable[[ExecutionContext, str], None]] = []
        
        # Profiling
        self.profiling_enabled = False
        self.profile_threshold = 1.0  # seconds
        
        # Storage limits
        self.max_recent_executions = 1000
        self.max_resource_snapshots = 100
        
        # Thread safety
        self.lock = threading.RLock()
        
        # Create logger
        self.logger = get_logger(__name__)
        
        # Default thresholds
        self._set_default_thresholds()
    
    def _set_default_thresholds(self) -> None:
        """Set default performance thresholds."""
        # API thresholds
        self.add_threshold(PerformanceCategory.API, 1.0, 3.0)
        
        # Database thresholds
        self.add_threshold(PerformanceCategory.DATABASE, 0.5, 2.0)
        
        # Background task thresholds
        self.add_threshold(PerformanceCategory.BACKGROUND, 10.0, 60.0)
        
        # Computation thresholds
        self.add_threshold(PerformanceCategory.COMPUTATION, 2.0, 10.0)
        
        # External service thresholds
        self.add_threshold(PerformanceCategory.EXTERNAL, 2.0, 10.0)
    
    def start_tracking(
        self,
        category: PerformanceCategory,
        name: str,
        context: Optional[Dict[str, Any]] = None,
        correlation_id: Optional[str] = None,
        enable_profiling: bool = False
    ) -> ExecutionContext:
        """
        Start tracking an operation's performance.
        
        Args:
            category: Performance category
            name: Operation name
            context: Additional context information
            correlation_id: Correlation ID for distributed tracing
            enable_profiling: Whether to enable profiling for this execution
            
        Returns:
            ExecutionContext for the operation
        """
        with self.lock:
            # Create execution context
            exec_context = ExecutionContext(
                category=category,
                name=name,
                context=context,
                correlation_id=correlation_id
            )
            
            # Store in active executions
            self.active_executions[exec_context.correlation_id] = exec_context
            
            # Start profiling if enabled
            if self.profiling_enabled and enable_profiling:
                self._start_profiling(exec_context.correlation_id)
            
            return exec_context
    
    def stop_tracking(
        self,
        correlation_id: str,
        additional_context: Optional[Dict[str, Any]] = None
    ) -> Optional[ExecutionContext]:
        """
        Stop tracking an operation's performance.
        
        Args:
            correlation_id: Correlation ID of the execution
            additional_context: Additional context to add to the execution
            
        Returns:
            Completed ExecutionContext or None if not found
        """
        with self.lock:
            if correlation_id not in self.active_executions:
                return None
            
            # Get execution context
            exec_context = self.active_executions[correlation_id]
            
            # Complete the execution
            exec_context.complete()
            
            # Add additional context
            if additional_context:
                exec_context.context.update(additional_context)
            
            # Stop profiling if enabled
            if self.profiling_enabled:
                self._stop_profiling(correlation_id, exec_context)
            
            # Check thresholds
            self._check_thresholds(exec_context)
            
            # Update Prometheus metrics if available
            if PROMETHEUS_AVAILABLE:
                # Update metrics based on category
                if exec_context.category == PerformanceCategory.API:
                    REQUEST_LATENCY.labels(
                        method=exec_context.context.get('method', 'unknown'),
                        endpoint=exec_context.name
                    ).observe(exec_context.duration)
                elif exec_context.category == PerformanceCategory.DATABASE:
                    DB_QUERY_LATENCY.labels(
                        operation=exec_context.name,
                        table=exec_context.context.get('table', 'unknown')
                    ).observe(exec_context.duration)
            
            # Remove from active executions
            del self.active_executions[correlation_id]
            
            # Add to executions and recent executions
            self.executions[correlation_id] = exec_context
            self.recent_executions.append(exec_context)
            
            # Enforce storage limits
            while len(self.recent_executions) > self.max_recent_executions:
                self.recent_executions.pop(0)
            
            # Update statistics
            self._update_statistics(exec_context)
            
            # Log slow operations
            self._log_slow_operation(exec_context)
            
            return exec_context
    
    def _start_profiling(self, correlation_id: str) -> None:
        """
        Start profiling an operation.
        
        Args:
            correlation_id: Correlation ID of the execution
        """
        # Not implemented - would use cProfile
        pass
    
    def _stop_profiling(self, correlation_id: str, exec_context: ExecutionContext) -> None:
        """
        Stop profiling an operation and capture the results.
        
        Args:
            correlation_id: Correlation ID of the execution
            exec_context: Execution context to update with profile data
        """
        # Not implemented - would use cProfile
        pass
    
    def _check_thresholds(self, exec_context: ExecutionContext) -> None:
        """
        Check if the execution time exceeds any thresholds.
        
        Args:
            exec_context: Completed execution context
        """
        category = exec_context.category
        name = exec_context.name
        duration = exec_context.duration
        
        # Check category + name specific threshold
        threshold = self.thresholds.get((category, name))
        
        # Fall back to category threshold if no specific threshold
        if threshold is None:
            threshold = self.thresholds.get((category, None))
        
        # Check threshold
        if threshold:
            level = threshold.check(duration)
            if level:
                # Log threshold violation
                self._log_threshold_violation(exec_context, threshold, level)
                
                # Notify callbacks
                for callback in self.threshold_callbacks:
                    try:
                        callback(exec_context, level)
                    except Exception as e:
                        self.logger.error(f"Error in threshold callback: {str(e)}")
    
    def _log_threshold_violation(
        self,
        exec_context: ExecutionContext,
        threshold: PerformanceThreshold,
        level: str
    ) -> None:
        """
        Log a performance threshold violation.
        
        Args:
            exec_context: Execution context that violated the threshold
            threshold: Threshold that was violated
            level: Violation level ('critical' or 'warning')
        """
        log_level = logging.CRITICAL if level == 'critical' else logging.WARNING
        
        message = (
            f"Performance {level}: {exec_context.category.value} operation '{exec_context.name}' "
            f"took {exec_context.duration:.3f}s (threshold: {getattr(threshold, f'{level}_threshold'):.3f}s)"
        )
        
        log_with_context(
            self.logger,
            logging.getLevelName(log_level).lower(),
            message,
            context={
                'performance_tracking': {
                    'category': exec_context.category.value,
                    'name': exec_context.name,
                    'duration': exec_context.duration,
                    'threshold_level': level,
                    'threshold_value': getattr(threshold, f'{level}_threshold'),
                    'correlation_id': exec_context.correlation_id
                },
                **exec_context.context
            }
        )
    
    def _log_slow_operation(self, exec_context: ExecutionContext) -> None:
        """
        Log details about slow operations (regardless of threshold).
        
        Args:
            exec_context: Completed execution context
        """
        # Define slow operation thresholds by category
        slow_thresholds = {
            PerformanceCategory.API: 0.5,
            PerformanceCategory.DATABASE: 0.1,
            PerformanceCategory.BACKGROUND: 5.0,
            PerformanceCategory.COMPUTATION: 1.0,
            PerformanceCategory.EXTERNAL: 1.0,
            PerformanceCategory.SYSTEM: 0.5,
            PerformanceCategory.APPLICATION: 0.5,
        }
        
        threshold = slow_thresholds.get(exec_context.category, 1.0)
        
        if exec_context.duration >= threshold:
            # Log at INFO level (not a threshold violation, just information)
            message = (
                f"Slow operation: {exec_context.category.value} '{exec_context.name}' "
                f"took {exec_context.duration:.3f}s"
            )
            
            context = {
                'performance_tracking': {
                    'category': exec_context.category.value,
                    'name': exec_context.name,
                    'duration': exec_context.duration,
                    'correlation_id': exec_context.correlation_id
                }
            }
            
            # Add database-specific context
            if exec_context.category == PerformanceCategory.DATABASE:
                context['performance_tracking']['query'] = exec_context.query
                context['performance_tracking']['db_name'] = exec_context.db_name
            
            log_with_context(
                self.logger,
                'info',
                message,
                context={**context, **exec_context.context}
            )
    
    def _update_statistics(self, exec_context: ExecutionContext) -> None:
        """
        Update performance statistics with a completed execution.
        
        Args:
            exec_context: Completed execution context
        """
        category = exec_context.category.value
        name = exec_context.name
        duration = exec_context.duration
        
        # Update category statistics
        if category not in self.stats_by_category:
            self.stats_by_category[category] = {
                'count': 0,
                'total_time': 0,
                'min_time': float('inf'),
                'max_time': 0,
                'avg_time': 0,
                'recent_times': []
            }
        
        cat_stats = self.stats_by_category[category]
        cat_stats['count'] += 1
        cat_stats['total_time'] += duration
        cat_stats['min_time'] = min(cat_stats['min_time'], duration)
        cat_stats['max_time'] = max(cat_stats['max_time'], duration)
        cat_stats['avg_time'] = cat_stats['total_time'] / cat_stats['count']
        cat_stats['recent_times'].append(duration)
        
        # Keep only the last 100 executions for moving average
        if len(cat_stats['recent_times']) > 100:
            cat_stats['recent_times'].pop(0)
        
        # Update name statistics
        key = f"{category}.{name}"
        if key not in self.stats_by_name:
            self.stats_by_name[key] = {
                'category': category,
                'name': name,
                'count': 0,
                'total_time': 0,
                'min_time': float('inf'),
                'max_time': 0,
                'avg_time': 0,
                'recent_times': []
            }
        
        name_stats = self.stats_by_name[key]
        name_stats['count'] += 1
        name_stats['total_time'] += duration
        name_stats['min_time'] = min(name_stats['min_time'], duration)
        name_stats['max_time'] = max(name_stats['max_time'], duration)
        name_stats['avg_time'] = name_stats['total_time'] / name_stats['count']
        name_stats['recent_times'].append(duration)
        
        # Keep only the last 100 executions for moving average
        if len(name_stats['recent_times']) > 100:
            name_stats['recent_times'].pop(0)
        
        # Update bottlenecks (top 10 slowest operations by average time)
        # Only consider operations with at least 5 executions
        if name_stats['count'] >= 5:
            # Update the bottlenecks list
            self.bottlenecks = sorted(
                [stats for stats in self.stats_by_name.values() if stats['count'] >= 5],
                key=lambda x: x['avg_time'],
                reverse=True
            )[:10]
    
    def add_threshold(
        self,
        category: PerformanceCategory,
        warning_threshold: float,
        critical_threshold: float,
        name: Optional[str] = None
    ) -> None:
        """
        Add a performance threshold.
        
        Args:
            category: Performance category
            warning_threshold: Warning threshold in seconds
            critical_threshold: Critical threshold in seconds
            name: Optional name for specific operation thresholds
        """
        with self.lock:
            threshold = PerformanceThreshold(
                category=category,
                warning_threshold=warning_threshold,
                critical_threshold=critical_threshold,
                name=name
            )
            
            self.thresholds[(category, name)] = threshold
    
    def add_threshold_callback(self, callback: Callable[[ExecutionContext, str], None]) -> None:
        """
        Add a callback function for threshold violations.
        
        Args:
            callback: Function to call when a threshold is violated
        """
        with self.lock:
            self.threshold_callbacks.append(callback)
    
    def track_resource_usage(self, force: bool = False) -> None:
        """
        Track system resource usage.
        
        Args:
            force: Force tracking even if interval hasn't elapsed
        """
        with self.lock:
            now = time.time()
            
            # Check if it's time to take a snapshot
            if not force and now - self.last_resource_snapshot < self.resource_tracking_interval:
                return
            
            self.last_resource_snapshot = now
            
            # Get resource usage
            process = psutil.Process(os.getpid())
            
            # CPU usage
            cpu_percent = process.cpu_percent(interval=0.1)
            
            # Memory usage
            memory_info = process.memory_info()
            memory_percent = process.memory_percent()
            
            # Disk usage
            disk_io = process.io_counters() if hasattr(process, 'io_counters') else None
            
            # Thread count
            thread_count = threading.active_count()
            
            # System-wide resource usage
            system_cpu = psutil.cpu_percent(interval=0.1)
            system_memory = psutil.virtual_memory()
            system_disk = psutil.disk_usage('/')
            
            # Create snapshot
            snapshot = {
                'timestamp': now,
                'process': {
                    'cpu_percent': cpu_percent,
                    'memory_rss': memory_info.rss,
                    'memory_vms': memory_info.vms,
                    'memory_percent': memory_percent,
                    'thread_count': thread_count
                },
                'system': {
                    'cpu_percent': system_cpu,
                    'memory_total': system_memory.total,
                    'memory_available': system_memory.available,
                    'memory_percent': system_memory.percent,
                    'disk_total': system_disk.total,
                    'disk_free': system_disk.free,
                    'disk_percent': system_disk.percent
                }
            }
            
            if disk_io:
                snapshot['process']['disk_read_bytes'] = disk_io.read_bytes
                snapshot['process']['disk_write_bytes'] = disk_io.write_bytes
            
            # Add snapshot to list
            self.resource_snapshots.append(snapshot)
            
            # Enforce storage limits
            if len(self.resource_snapshots) > self.max_resource_snapshots:
                self.resource_snapshots.pop(0)
            
            # Update Prometheus metrics if available
            if PROMETHEUS_AVAILABLE:
                CPU_USAGE.set(cpu_percent)
                MEMORY_USAGE.labels(type='rss').set(memory_info.rss)
                MEMORY_USAGE.labels(type='vms').set(memory_info.vms)
                DISK_USAGE.labels(type='total').set(system_disk.total)
                DISK_USAGE.labels(type='free').set(system_disk.free)
                DISK_USAGE.labels(type='used').set(system_disk.used)
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get performance statistics.
        
        Returns:
            Dictionary with performance statistics
        """
        with self.lock:
            # Update resource usage
            self.track_resource_usage(force=True)
            
            # Calculate current resource usage
            if self.resource_snapshots:
                latest_snapshot = self.resource_snapshots[-1]
                current_resources = {
                    'cpu_percent': latest_snapshot['process']['cpu_percent'],
                    'memory_rss': latest_snapshot['process']['memory_rss'],
                    'memory_percent': latest_snapshot['process']['memory_percent'],
                    'thread_count': latest_snapshot['process']['thread_count']
                }
            else:
                current_resources = {}
            
            # Calculate resource usage trends
            resource_trends = self._calculate_resource_trends()
            
            # Get active executions
            active = [
                {
                    'category': exec_context.category.value,
                    'name': exec_context.name,
                    'correlation_id': exec_context.correlation_id,
                    'start_time': exec_context.start_time,
                    'elapsed': time.time() - exec_context.start_time
                }
                for exec_context in self.active_executions.values()
            ]
            
            # Sort active executions by elapsed time (descending)
            active.sort(key=lambda x: x['elapsed'], reverse=True)
            
            # Prepare recent executions
            recent = [
                {
                    'category': exec_context.category.value,
                    'name': exec_context.name,
                    'correlation_id': exec_context.correlation_id,
                    'start_time': exec_context.start_time,
                    'duration': exec_context.duration,
                    'end_time': exec_context.end_time
                }
                for exec_context in reversed(self.recent_executions)
            ][:100]  # Only return the 100 most recent
            
            # Prepare bottlenecks
            bottlenecks = [
                {
                    'category': item['category'],
                    'name': item['name'],
                    'avg_time': item['avg_time'],
                    'max_time': item['max_time'],
                    'count': item['count']
                }
                for item in self.bottlenecks
            ]
            
            # Prepare statistics by category
            by_category = {}
            for category, stats in self.stats_by_category.items():
                by_category[category] = {
                    'count': stats['count'],
                    'total_time': stats['total_time'],
                    'min_time': stats['min_time'],
                    'max_time': stats['max_time'],
                    'avg_time': stats['avg_time'],
                    'p95_time': self._calculate_percentile(stats['recent_times'], 95)
                }
            
            # Prepare overall statistics
            overall = {
                'total_executions': sum(stats['count'] for stats in self.stats_by_category.values()),
                'active_executions': len(self.active_executions),
                'total_time': sum(stats['total_time'] for stats in self.stats_by_category.values())
            }
            
            return {
                'overall': overall,
                'current_resources': current_resources,
                'resource_trends': resource_trends,
                'by_category': by_category,
                'bottlenecks': bottlenecks,
                'active_executions': active,
                'recent_executions': recent
            }
    
    def _calculate_percentile(self, values: List[float], percentile: float) -> float:
        """
        Calculate the percentile of a list of values.
        
        Args:
            values: List of values
            percentile: Percentile to calculate (0-100)
            
        Returns:
            Percentile value
        """
        if not values:
            return 0
        
        # Sort values
        sorted_values = sorted(values)
        
        # Calculate index
        k = (len(sorted_values) - 1) * (percentile / 100.0)
        f = int(k)
        c = math.ceil(k)
        
        if f == c:
            return sorted_values[int(k)]
        else:
            # Interpolate
            d0 = sorted_values[f] * (c - k)
            d1 = sorted_values[c] * (k - f)
            return d0 + d1
    
    def _calculate_resource_trends(self) -> Dict[str, Any]:
        """
        Calculate resource usage trends.
        
        Returns:
            Dictionary with resource usage trends
        """
        if len(self.resource_snapshots) < 2:
            return {}
        
        # Get the oldest and newest snapshots
        oldest = self.resource_snapshots[0]
        newest = self.resource_snapshots[-1]
        
        # Calculate time difference
        time_diff = newest['timestamp'] - oldest['timestamp']
        if time_diff <= 0:
            return {}
        
        # Calculate trends
        cpu_trend = (newest['process']['cpu_percent'] - oldest['process']['cpu_percent']) / time_diff
        memory_trend = (newest['process']['memory_rss'] - oldest['process']['memory_rss']) / time_diff
        
        return {
            'cpu_percent_per_minute': cpu_trend * 60,
            'memory_bytes_per_minute': memory_trend * 60
        }
    
    def enable_profiling(self, threshold: float = 1.0) -> None:
        """
        Enable profiling for slow operations.
        
        Args:
            threshold: Minimum duration in seconds to trigger profiling
        """
        with self.lock:
            self.profiling_enabled = True
            self.profile_threshold = threshold
    
    def disable_profiling(self) -> None:
        """Disable profiling."""
        with self.lock:
            self.profiling_enabled = False
    
    def get_execution(self, correlation_id: str) -> Optional[ExecutionContext]:
        """
        Get an execution context by correlation ID.
        
        Args:
            correlation_id: Correlation ID
            
        Returns:
            ExecutionContext or None if not found
        """
        with self.lock:
            return self.executions.get(correlation_id)
    
    def clear_stats(self) -> None:
        """Clear all performance statistics."""
        with self.lock:
            # Clear statistics
            self.stats_by_category.clear()
            self.stats_by_name.clear()
            self.bottlenecks.clear()
            
            # Keep active executions and resource snapshots
            # but clear completed executions
            self.executions.clear()
            self.recent_executions.clear()


# Global performance tracker instance
tracker = PerformanceTracker.get_instance()


def track_execution_time(
    category: PerformanceCategory,
    name: Optional[str] = None,
    context: Optional[Dict[str, Any]] = None,
    correlation_id: Optional[str] = None,
    enable_profiling: bool = False
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator to track execution time of a function.
    
    Args:
        category: Performance category
        name: Operation name (defaults to function name)
        context: Additional context
        correlation_id: Correlation ID for distributed tracing
        enable_profiling: Whether to enable profiling for this execution
        
    Returns:
        Decorator function
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        nonlocal name
        
        # Use function name if name not provided
        if name is None:
            name = func.__name__
        
        @functools.wraps(func)
        def wrapper(*args, **kwargs) -> T:
            # Create context with function information
            func_context = {
                'function': func.__name__,
                'module': func.__module__
            }
            
            if context:
                func_context.update(context)
            
            # Start tracking
            exec_context = tracker.start_tracking(
                category=category,
                name=name,
                context=func_context,
                correlation_id=correlation_id,
                enable_profiling=enable_profiling
            )
            
            try:
                # Execute function
                result = func(*args, **kwargs)
                
                # Additional context for successful execution
                additional_context = {'success': True}
                
                return result
            except Exception as e:
                # Additional context for failed execution
                additional_context = {
                    'success': False,
                    'error': str(e),
                    'error_type': type(e).__name__
                }
                
                # Re-raise the exception
                raise
            finally:
                # Stop tracking
                tracker.stop_tracking(
                    correlation_id=exec_context.correlation_id,
                    additional_context=additional_context
                )
        
        return wrapper
    
    return decorator


def track_database_query(
    query: str,
    db_name: str = 'default',
    params: Optional[Any] = None,
    context: Optional[Dict[str, Any]] = None,
    correlation_id: Optional[str] = None
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator to track database query execution time.
    
    Args:
        query: SQL query or operation description
        db_name: Database name
        params: Query parameters
        context: Additional context
        correlation_id: Correlation ID for distributed tracing
        
    Returns:
        Decorator function
    """
    # Create query context
    query_context = {
        'query': query,
        'db_name': db_name,
        'params': params
    }
    
    if context:
        query_context.update(context)
    
    # Use track_execution_time with DATABASE category
    return track_execution_time(
        category=PerformanceCategory.DATABASE,
        name=query[:50] + ('...' if len(query) > 50 else ''),  # Use truncated query as name
        context=query_context,
        correlation_id=correlation_id
    )


class monitor_query:
    """
    Context manager to track database query execution time.
    
    Usage:
        with monitor_query("SELECT * FROM users", db_name="users_db"):
            # Execute query
            cursor.execute("SELECT * FROM users")
            results = cursor.fetchall()
    """
    
    def __init__(
        self,
        query: str,
        db_name: str = 'default',
        params: Optional[Any] = None,
        context: Optional[Dict[str, Any]] = None,
        correlation_id: Optional[str] = None
    ):
        """
        Initialize the context manager.
        
        Args:
            query: SQL query or operation description
            db_name: Database name
            params: Query parameters
            context: Additional context
            correlation_id: Correlation ID for distributed tracing
        """
        self.query = query
        self.db_name = db_name
        self.params = params
        self.context = context or {}
        self.correlation_id = correlation_id
        
        # Create query context
        self.query_context = {
            'query': query,
            'db_name': db_name,
            'params': params
        }
        
        if context:
            self.query_context.update(context)
    
    def __enter__(self) -> str:
        """
        Start tracking the query.
        
        Returns:
            Correlation ID for the tracking
        """
        # Start tracking
        self.exec_context = tracker.start_tracking(
            category=PerformanceCategory.DATABASE,
            name=self.query[:50] + ('...' if len(self.query) > 50 else ''),
            context=self.query_context,
            correlation_id=self.correlation_id
        )
        
        return self.exec_context.correlation_id
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """
        Stop tracking the query.
        
        Args:
            exc_type: Exception type
            exc_val: Exception value
            exc_tb: Exception traceback
        """
        # Additional context for execution result
        additional_context = {'success': exc_type is None}
        
        if exc_type is not None:
            additional_context.update({
                'error': str(exc_val),
                'error_type': exc_type.__name__
            })
        
        # Stop tracking
        tracker.stop_tracking(
            correlation_id=self.exec_context.correlation_id,
            additional_context=additional_context
        )


def track_resource_usage() -> None:
    """Track system resource usage using the global performance tracker."""
    tracker.track_resource_usage()


def enable_profiling(threshold: float = 1.0) -> None:
    """
    Enable profiling for slow operations.
    
    Args:
        threshold: Minimum duration in seconds to trigger profiling
    """
    tracker.enable_profiling(threshold)


def add_performance_threshold(
    category: PerformanceCategory,
    warning_threshold: float,
    critical_threshold: float,
    name: Optional[str] = None
) -> None:
    """
    Add a performance threshold to the global tracker.
    
    Args:
        category: Performance category
        warning_threshold: Warning threshold in seconds
        critical_threshold: Critical threshold in seconds
        name: Optional name for specific operation thresholds
    """
    tracker.add_threshold(
        category=category,
        warning_threshold=warning_threshold,
        critical_threshold=critical_threshold,
        name=name
    )


def get_performance_stats() -> Dict[str, Any]:
    """
    Get performance statistics from the global tracker.
    
    Returns:
        Dictionary with performance statistics
    """
    return tracker.get_stats()


def start_performance_tracking(
    category: PerformanceCategory,
    name: str,
    context: Optional[Dict[str, Any]] = None,
    correlation_id: Optional[str] = None
) -> str:
    """
    Start tracking an operation's performance.
    
    Args:
        category: Performance category
        name: Operation name
        context: Additional context
        correlation_id: Correlation ID for distributed tracing
        
    Returns:
        Correlation ID for the tracking
    """
    exec_context = tracker.start_tracking(
        category=category,
        name=name,
        context=context,
        correlation_id=correlation_id
    )
    
    return exec_context.correlation_id


def stop_performance_tracking(
    correlation_id: str,
    additional_context: Optional[Dict[str, Any]] = None
) -> None:
    """
    Stop tracking an operation's performance.
    
    Args:
        correlation_id: Correlation ID returned from start_performance_tracking
        additional_context: Additional context to add to the execution
    """
    tracker.stop_tracking(
        correlation_id=correlation_id,
        additional_context=additional_context
    )


import math  # Required for _calculate_percentile method