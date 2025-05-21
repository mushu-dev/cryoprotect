#!/usr/bin/env python3
"""
Unified Monitoring Module for CryoProtect v2

This module provides a comprehensive monitoring solution that combines:
1. Database connection health monitoring
2. API endpoint health checking
3. Performance metrics collection
4. Resource usage tracking
5. Observability and distributed tracing
6. Error aggregation and reporting
7. Web dashboard for visualization

Usage:
    from unified_monitoring import MonitoringService
    
    # Initialize with application
    monitoring = MonitoringService(app)
    
    # Or initialize manually
    monitoring = MonitoringService()
    monitoring.init_app(app)
    
    # Start all monitoring services
    monitoring.start()
    
    # Use specific monitoring features
    with monitoring.track_operation("database_query"):
        # Your code here
        pass
"""

import os
import time
import json
import uuid
import logging
import threading
import traceback
import socket
import platform
import psutil
from typing import Dict, Any, List, Optional, Union, Callable, Tuple, Set
from datetime import datetime, timedelta
from contextlib import contextmanager
from functools import wraps
import inspect

# For web dashboard
try:
    from flask import Flask, request, g, jsonify, render_template_string
except ImportError:
    Flask = None
    request = None
    g = None
    jsonify = None
    render_template_string = None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('unified_monitoring')

# Constants for header names
CORRELATION_ID_HEADER = 'X-Correlation-ID'
REQUEST_ID_HEADER = 'X-Request-ID'
TRACE_ID_HEADER = 'X-Trace-ID'
SPAN_ID_HEADER = 'X-Span-ID'
TIMING_HEADER = 'X-Response-Time-Ms'

class MonitoringService:
    """
    Unified monitoring service that integrates multiple monitoring capabilities.
    
    This class provides:
    - Database connection health monitoring
    - API endpoint health checking
    - Performance metrics collection
    - Resource usage tracking
    - Observability and distributed tracing
    - Error aggregation and reporting
    - Web dashboard for visualization
    """
    
    _instance = None
    _lock = threading.Lock()
    
    @classmethod
    def get_instance(cls) -> 'MonitoringService':
        """
        Get singleton instance of MonitoringService.
        
        Returns:
            MonitoringService: Singleton instance
        """
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = MonitoringService()
        return cls._instance
    
    def __init__(self, app=None):
        """
        Initialize monitoring service.
        
        Args:
            app: Optional Flask application to initialize with
        """
        if MonitoringService._instance is not None:
            raise RuntimeError("MonitoringService is a singleton. Use get_instance() instead.")
            
        MonitoringService._instance = self
        
        # Application reference
        self.app = None
        
        # Monitoring state
        self.monitoring_active = False
        self.monitoring_threads = {}
        
        # Initialize component states
        self._init_component_states()
        
        # Initialize metrics collectors
        self._init_metrics_collectors()
        
        # Initialize with app if provided
        if app is not None:
            self.init_app(app)
        
        logger.info("MonitoringService initialized")
    
    def _init_component_states(self):
        """Initialize monitoring component states."""
        # Health monitoring state
        self.health_status = {
            "database": {
                "status": "Unknown",
                "last_check_time": 0,
                "is_healthy": False,
                "active_connection_type": None,
                "connection_types": {},
                "recent_failures": [],
                "recent_successes": [],
                "check_frequency_seconds": 60
            },
            "api": {
                "status": "Unknown",
                "endpoints": {},
                "last_check_time": 0,
                "overall_health_percentage": 0,
                "check_frequency_seconds": 60
            },
            "system": {
                "status": "Unknown",
                "cpu_usage": 0,
                "memory_usage": 0,
                "disk_usage": 0,
                "last_check_time": 0,
                "check_frequency_seconds": 30
            }
        }
        
        # Observability state
        self.tracing_enabled = False
        self.current_traces = {}
        self.active_operations = {}
        
        # Error aggregation state
        self.errors = []
        self.error_counts = {}
        self.max_error_history = 1000
        
        # Alert state
        self.alerts = []
        self.alert_handlers = []
        self.alert_silence_periods = {}
        
        # Dashboard state
        self.dashboard_enabled = False
        self.dashboard_port = 5001
        
        # Monitoring configuration
        self.config = {
            "monitoring_dir": "monitoring",
            "metrics_retention_days": 7,
            "health_check_intervals": {
                "database": 60,
                "api": 60,
                "system": 30
            },
            "alert_thresholds": {
                "database_failures": 3,
                "api_errors": 3,
                "memory_usage": 90,
                "cpu_usage": 90,
                "disk_usage": 90
            }
        }
    
    def _init_metrics_collectors(self):
        """Initialize metrics collectors."""
        # Performance metrics collectors
        self.performance_metrics = {
            "database": PerformanceMetrics("database_operations"),
            "api": PerformanceMetrics("api_operations"),
            "system": PerformanceMetrics("system_resources")
        }
        
        # Progress trackers by operation
        self.progress_trackers = {}
        
        # Resource usage metrics
        self.resource_metrics = {
            "cpu_history": [],
            "memory_history": [],
            "disk_history": [],
            "network_history": [],
            "timestamps": []
        }
        
        # Operation timing data
        self.operation_timings = {}
    
    def init_app(self, app):
        """
        Initialize monitoring with a Flask application.
        
        Args:
            app: Flask application
        """
        if not app or not hasattr(app, 'before_request'):
            logger.warning("Not a valid Flask application")
            return
        
        self.app = app
        
        # Register middleware for request tracking
        self._register_middleware(app)
        
        # Add dashboard routes if Flask is available
        if Flask is not None:
            self._register_dashboard_routes(app)
        
        logger.info("Initialized monitoring with Flask application")
    
    def _register_middleware(self, app):
        """
        Register middleware with Flask application.
        
        Args:
            app: Flask application
        """
        @app.before_request
        def before_request():
            """Middleware to run before each request."""
            # Start request timing
            g.start_time = time.time()
            
            # Generate request ID if not already present
            g.request_id = request.headers.get(REQUEST_ID_HEADER) or str(uuid.uuid4())
            
            # Use existing correlation ID from header or generate a new one
            g.correlation_id = request.headers.get(CORRELATION_ID_HEADER) or str(uuid.uuid4())
            
            # Generate trace ID if not already present
            g.trace_id = request.headers.get(TRACE_ID_HEADER) or str(uuid.uuid4())
            
            # Generate span ID for this request
            g.span_id = str(uuid.uuid4())
            
            # Record trace start in monitoring service
            if self.tracing_enabled:
                self.current_traces[g.trace_id] = {
                    "start_time": time.time(),
                    "correlation_id": g.correlation_id,
                    "request_id": g.request_id,
                    "spans": {
                        g.span_id: {
                            "name": request.path,
                            "start_time": time.time(),
                            "end_time": None,
                            "parent_span_id": None
                        }
                    }
                }
        
        @app.after_request
        def after_request(response):
            """Middleware to run after each request."""
            # Calculate request duration
            if hasattr(g, 'start_time'):
                duration = time.time() - g.start_time
                
                # Add timing header to response
                if duration is not None:
                    response.headers[TIMING_HEADER] = str(int(duration * 1000))
                    
                    # Record API operation performance
                    self.performance_metrics["api"].record_operation(
                        operation_type="http_request",
                        success=response.status_code < 400,
                        execution_time=duration
                    )
            
            # Add tracing headers to response
            if hasattr(g, 'correlation_id'):
                response.headers[CORRELATION_ID_HEADER] = g.correlation_id
            if hasattr(g, 'request_id'):
                response.headers[REQUEST_ID_HEADER] = g.request_id
            if hasattr(g, 'trace_id'):
                response.headers[TRACE_ID_HEADER] = g.trace_id
            if hasattr(g, 'span_id'):
                response.headers[SPAN_ID_HEADER] = g.span_id
            
            # Record trace completion
            if self.tracing_enabled and hasattr(g, 'trace_id') and g.trace_id in self.current_traces:
                trace = self.current_traces[g.trace_id]
                if hasattr(g, 'span_id') and g.span_id in trace["spans"]:
                    trace["spans"][g.span_id]["end_time"] = time.time()
                    trace["spans"][g.span_id]["status_code"] = response.status_code
            
            return response
        
        @app.errorhandler(Exception)
        def handle_exception(e):
            """Catch and log all exceptions."""
            # Log the error
            self.record_error("api_exception", str(e), traceback.format_exc())
            
            # Re-raise the exception (Flask will handle it)
            raise e
    
    def _register_dashboard_routes(self, app):
        """
        Register dashboard routes with Flask application.
        
        Args:
            app: Flask application
        """
        # Import dashboard template
        from api_monitoring_dashboard import DASHBOARD_TEMPLATE
        
        @app.route('/monitoring')
        def monitoring_dashboard():
            """Render the monitoring dashboard."""
            data = self.get_monitoring_data()
            return render_template_string(DASHBOARD_TEMPLATE, data=data)
        
        @app.route('/monitoring/api/status')
        def monitoring_api_status():
            """Return the current monitoring data as JSON."""
            return jsonify(self.get_monitoring_data())
        
        @app.route('/monitoring/api/health')
        def monitoring_api_health():
            """Return a simplified health status."""
            health_data = {
                "status": "healthy",
                "database": self.health_status["database"]["status"],
                "api": self.health_status["api"]["status"],
                "system": self.health_status["system"]["status"],
                "timestamp": datetime.now().isoformat()
            }
            
            # Set overall status based on component statuses
            if any(status in ["Error", "Unhealthy"] for status in [
                health_data["database"], health_data["api"], health_data["system"]
            ]):
                health_data["status"] = "unhealthy"
            elif any(status in ["Warning"] for status in [
                health_data["database"], health_data["api"], health_data["system"]
            ]):
                health_data["status"] = "warning"
            
            return jsonify(health_data)
    
    def start(self):
        """Start all monitoring services."""
        self._ensure_monitoring_directory()
        
        # Start monitoring threads if not already running
        if not self.monitoring_active:
            self.monitoring_active = True
            
            # Start database health monitoring
            self._start_monitoring_thread(
                "database_health",
                self._database_health_check_thread,
                self.config["health_check_intervals"]["database"]
            )
            
            # Start API health monitoring
            self._start_monitoring_thread(
                "api_health",
                self._api_health_check_thread,
                self.config["health_check_intervals"]["api"]
            )
            
            # Start system resource monitoring
            self._start_monitoring_thread(
                "system_health",
                self._system_health_check_thread,
                self.config["health_check_intervals"]["system"]
            )
            
            # Start metrics persistence thread
            self._start_monitoring_thread(
                "metrics_persistence",
                self._metrics_persistence_thread,
                300  # 5 minutes
            )
            
            # Enable tracing
            self.tracing_enabled = True
            
            logger.info("Started all monitoring services")
    
    def stop(self):
        """Stop all monitoring services."""
        if self.monitoring_active:
            self.monitoring_active = False
            
            # Stop all monitoring threads
            for thread_name, thread_info in list(self.monitoring_threads.items()):
                thread_info["stop_event"].set()
                thread_info["thread"].join(timeout=5)
                if thread_info["thread"].is_alive():
                    logger.warning(f"Thread {thread_name} did not terminate gracefully")
                else:
                    logger.info(f"Thread {thread_name} stopped")
                    del self.monitoring_threads[thread_name]
            
            # Disable tracing
            self.tracing_enabled = False
            
            logger.info("Stopped all monitoring services")
    
    def _start_monitoring_thread(self, thread_name, thread_func, interval_seconds):
        """
        Start a monitoring thread.
        
        Args:
            thread_name: Name of the thread
            thread_func: Thread function
            interval_seconds: Interval between checks in seconds
        """
        if thread_name in self.monitoring_threads:
            logger.warning(f"Thread {thread_name} is already running")
            return
        
        stop_event = threading.Event()
        
        thread = threading.Thread(
            target=self._thread_wrapper,
            args=(thread_func, stop_event, interval_seconds),
            name=f"monitoring_{thread_name}",
            daemon=True
        )
        
        thread.start()
        
        self.monitoring_threads[thread_name] = {
            "thread": thread,
            "stop_event": stop_event,
            "interval_seconds": interval_seconds,
            "last_executed": time.time()
        }
        
        logger.info(f"Started {thread_name} thread with interval {interval_seconds}s")
    
    def _thread_wrapper(self, thread_func, stop_event, interval_seconds):
        """
        Wrapper for monitoring threads.
        
        Args:
            thread_func: Thread function
            stop_event: Stop event
            interval_seconds: Interval between checks in seconds
        """
        while not stop_event.is_set():
            try:
                thread_func()
            except Exception as e:
                logger.error(f"Error in monitoring thread: {str(e)}")
                traceback.print_exc()
            
            # Sleep until next interval or stop event is set
            stop_event.wait(interval_seconds)
    
    def _database_health_check_thread(self):
        """Thread function for database health checking."""
        # Try to import connection utilities
        try:
            from db_connection_utils import ConnectionManager
            connection_manager = ConnectionManager.get_instance()
        except (ImportError, AttributeError):
            connection_manager = None
            logger.warning("ConnectionManager not available, cannot check database health")
        
        # Skip if no connection manager
        if not connection_manager:
            self.health_status["database"]["status"] = "Unknown"
            self.health_status["database"]["last_check_time"] = time.time()
            return
        
        start_time = time.time()
        
        # Get connection manager stats
        conn_stats = connection_manager.get_stats() if hasattr(connection_manager, 'get_stats') else {}
        
        # Check circuit breaker states
        circuit_breakers = {}
        if hasattr(connection_manager, 'circuit_breakers'):
            for name, breaker in connection_manager.circuit_breakers.items():
                circuit_breakers[name] = breaker.get_state()
        
        # Check active connection pool
        active_pool = connection_manager.active_pool if hasattr(connection_manager, 'active_pool') else None
        
        # Try to get a test connection
        connection = None
        is_healthy = False
        error_message = None
        
        try:
            connection = connection_manager.get_connection()
            if connection:
                # Test the connection with a simple query
                if hasattr(connection, 'execute_query'):
                    result = connection.execute_query("SELECT 1 as test")
                    is_healthy = result and len(result) > 0 and result[0].get('test') == 1
                else:
                    # Assume connection is healthy if we got one
                    is_healthy = True
                
                self.performance_metrics["database"].record_operation(
                    operation_type="connection",
                    success=is_healthy,
                    execution_time=time.time() - start_time
                )
                
                if is_healthy:
                    self.health_status["database"]["recent_successes"].append({
                        "timestamp": time.time(),
                        "connection_type": active_pool,
                        "latency": time.time() - start_time
                    })
                    # Keep only last 10 successes
                    if len(self.health_status["database"]["recent_successes"]) > 10:
                        self.health_status["database"]["recent_successes"] = self.health_status["database"]["recent_successes"][-10:]
                else:
                    error_message = "Connection test query failed"
                    self.health_status["database"]["recent_failures"].append({
                        "timestamp": time.time(),
                        "connection_type": active_pool,
                        "error": error_message
                    })
                    # Keep only last 10 failures
                    if len(self.health_status["database"]["recent_failures"]) > 10:
                        self.health_status["database"]["recent_failures"] = self.health_status["database"]["recent_failures"][-10:]
        except Exception as e:
            is_healthy = False
            error_message = str(e)
            
            self.performance_metrics["database"].record_operation(
                operation_type="connection",
                success=False,
                execution_time=time.time() - start_time
            )
            
            self.health_status["database"]["recent_failures"].append({
                "timestamp": time.time(),
                "connection_type": active_pool,
                "error": error_message
            })
            # Keep only last 10 failures
            if len(self.health_status["database"]["recent_failures"]) > 10:
                self.health_status["database"]["recent_failures"] = self.health_status["database"]["recent_failures"][-10:]
        finally:
            if connection and hasattr(connection, 'close'):
                connection.close()
        
        # Update health status
        self.health_status["database"].update({
            "last_check_time": time.time(),
            "is_healthy": is_healthy,
            "status": "OK" if is_healthy else "Error",
            "active_connection_type": active_pool,
            "connection_types": circuit_breakers,
            "check_duration": time.time() - start_time,
            "error_message": error_message
        })
        
        # Check for persistent failures and trigger alerts
        if not is_healthy:
            recent_failures = self.health_status["database"]["recent_failures"]
            if len(recent_failures) >= self.config["alert_thresholds"]["database_failures"]:
                # Check if all failures are within the last 5 minutes
                recent_failure_timestamps = [f["timestamp"] for f in recent_failures[-3:]]
                if all(time.time() - ts < 300 for ts in recent_failure_timestamps):
                    self._trigger_alert(
                        "database_connection_failure",
                        f"Database connection failures: {error_message}",
                        "database"
                    )
        
        # Log health status
        if is_healthy:
            logger.info(f"Database health check passed (type: {active_pool}, "
                       f"latency: {self.health_status['database']['check_duration']:.3f}s)")
        else:
            logger.warning(f"Database health check failed (type: {active_pool}, "
                          f"error: {error_message}")
    
    def _api_health_check_thread(self):
        """Thread function for API health checking."""
        # Skip if no app
        if not self.app:
            self.health_status["api"]["status"] = "Unknown"
            self.health_status["api"]["last_check_time"] = time.time()
            return
        
        # Define endpoints to monitor (can be configured)
        endpoints = self._get_api_endpoints_to_monitor()
        
        total_response_time = 0
        error_count = 0
        endpoint_count = len(endpoints)
        
        for endpoint_info in endpoints:
            endpoint = endpoint_info["endpoint"]
            method = endpoint_info["method"]
            data = endpoint_info.get("data")
            
            status, message, response_time = self._check_api_endpoint(endpoint, method, data)
            
            if status == "Error":
                error_count += 1
                self.record_error(
                    "api_endpoint_error",
                    f"Error on {method} {endpoint}: {message}",
                    None,
                    {"endpoint": endpoint, "method": method}
                )
            
            total_response_time += response_time
            
            self.health_status["api"]["endpoints"][endpoint] = {
                "method": method,
                "status": status,
                "message": message,
                "response_time": response_time,
                "last_checked": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        
        # Update performance metrics
        avg_response_time = total_response_time / endpoint_count if endpoint_count > 0 else 0
        error_rate = (error_count / endpoint_count) * 100 if endpoint_count > 0 else 0
        
        # Update API health status
        if error_count == 0:
            api_status = "OK"
        elif error_count < endpoint_count / 2:
            api_status = "Warning"
        else:
            api_status = "Error"
        
        health_percentage = 100 - error_rate
        
        self.health_status["api"].update({
            "status": api_status,
            "last_check_time": time.time(),
            "overall_health_percentage": health_percentage,
            "avg_response_time": avg_response_time,
            "error_rate": error_rate
        })
        
        # Check for persistent failures and trigger alerts
        if error_count > 0 and error_count >= self.config["alert_thresholds"]["api_errors"]:
            self._trigger_alert(
                "api_health_degraded",
                f"API health check failed for {error_count}/{endpoint_count} endpoints",
                "api"
            )
        
        logger.info(f"API health check completed: {health_percentage:.1f}% healthy, "
                   f"avg response time: {avg_response_time:.3f}s")
    
    def _get_api_endpoints_to_monitor(self):
        """
        Get list of API endpoints to monitor.
        
        Returns:
            List of endpoint definitions
        """
        # Default endpoints to monitor
        endpoints = [
            {"endpoint": "/health", "method": "GET"},
            {"endpoint": "/api/v1/molecules", "method": "GET"}
        ]
        
        # Try to get endpoints from app rules
        if self.app and hasattr(self.app, 'url_map'):
            for rule in self.app.url_map.iter_rules():
                # Skip static files and non-GET endpoints
                if rule.endpoint == 'static' or 'GET' not in rule.methods:
                    continue
                
                # Add endpoint if it's not already in the list
                if not any(e["endpoint"] == rule.rule for e in endpoints):
                    endpoints.append({
                        "endpoint": rule.rule,
                        "method": "GET"
                    })
        
        return endpoints
    
    def _check_api_endpoint(self, endpoint, method="GET", data=None):
        """
        Check if an API endpoint is working.
        
        Args:
            endpoint: Endpoint to check
            method: HTTP method
            data: Request data
            
        Returns:
            Tuple of (status, message, response_time)
        """
        if not self.app:
            return "Unknown", "No Flask app available", 0
        
        try:
            # Use Flask test client
            client = self.app.test_client()
            
            start_time = time.time()
            
            with self.app.app_context():
                if method == "GET":
                    response = client.get(endpoint)
                elif method == "POST":
                    response = client.post(endpoint, json=data)
                elif method == "PUT":
                    response = client.put(endpoint, json=data)
                elif method == "DELETE":
                    response = client.delete(endpoint)
                else:
                    return "Error", f"Unsupported method: {method}", 0
                
                response_time = int((time.time() - start_time) * 1000)  # in milliseconds
                
                if response.status_code < 300:
                    status = "OK"
                elif response.status_code < 500:
                    status = "Warning"
                else:
                    status = "Error"
                
                return status, f"Status code: {response.status_code}", response_time
        except Exception as e:
            return "Error", str(e), 0
    
    def _system_health_check_thread(self):
        """Thread function for system resource monitoring."""
        try:
            # Get system resource usage
            cpu_percent = psutil.cpu_percent(interval=1)
            memory = psutil.virtual_memory()
            disk = psutil.disk_usage('/')
            
            memory_percent = memory.percent
            disk_percent = disk.percent
            
            # Determine system status based on resource usage
            if cpu_percent > 90 or memory_percent > 90 or disk_percent > 90:
                system_status = "Error"
            elif cpu_percent > 75 or memory_percent > 75 or disk_percent > 75:
                system_status = "Warning"
            else:
                system_status = "OK"
            
            # Update system health status
            self.health_status["system"].update({
                "status": system_status,
                "last_check_time": time.time(),
                "cpu_usage": cpu_percent,
                "memory_usage": memory_percent,
                "disk_usage": disk_percent,
                "cpu_info": {
                    "cores": psutil.cpu_count(logical=False),
                    "threads": psutil.cpu_count(logical=True)
                },
                "memory_info": {
                    "total": memory.total,
                    "available": memory.available,
                    "used": memory.used
                },
                "disk_info": {
                    "total": disk.total,
                    "free": disk.free,
                    "used": disk.used
                },
                "process_info": {
                    "pid": os.getpid(),
                    "cpu_percent": psutil.Process(os.getpid()).cpu_percent(),
                    "memory_percent": psutil.Process(os.getpid()).memory_percent(),
                    "threads": len(psutil.Process(os.getpid()).threads())
                }
            })
            
            # Update resource metrics history
            timestamp = time.time()
            self.resource_metrics["cpu_history"].append(cpu_percent)
            self.resource_metrics["memory_history"].append(memory_percent)
            self.resource_metrics["disk_history"].append(disk_percent)
            self.resource_metrics["timestamps"].append(timestamp)
            
            # Try to get network stats
            try:
                net_io = psutil.net_io_counters()
                self.resource_metrics["network_history"].append({
                    "bytes_sent": net_io.bytes_sent,
                    "bytes_recv": net_io.bytes_recv
                })
            except Exception:
                # Network stats may not be available
                self.resource_metrics["network_history"].append(None)
            
            # Keep only the last 100 data points
            for key in ["cpu_history", "memory_history", "disk_history", "network_history", "timestamps"]:
                if len(self.resource_metrics[key]) > 100:
                    self.resource_metrics[key] = self.resource_metrics[key][-100:]
            
            # Trigger alerts for high resource usage
            if cpu_percent > self.config["alert_thresholds"]["cpu_usage"]:
                self._trigger_alert(
                    "high_cpu_usage",
                    f"High CPU usage: {cpu_percent:.1f}%",
                    "system"
                )
            
            if memory_percent > self.config["alert_thresholds"]["memory_usage"]:
                self._trigger_alert(
                    "high_memory_usage",
                    f"High memory usage: {memory_percent:.1f}%",
                    "system"
                )
            
            if disk_percent > self.config["alert_thresholds"]["disk_usage"]:
                self._trigger_alert(
                    "high_disk_usage",
                    f"High disk usage: {disk_percent:.1f}%",
                    "system"
                )
            
            logger.info(f"System health check: CPU: {cpu_percent:.1f}%, "
                       f"Memory: {memory_percent:.1f}%, Disk: {disk_percent:.1f}%")
            
        except Exception as e:
            logger.error(f"Error in system health check: {str(e)}")
            self.health_status["system"].update({
                "status": "Unknown",
                "last_check_time": time.time(),
                "error": str(e)
            })
    
    def _metrics_persistence_thread(self):
        """Thread function for metrics persistence."""
        self._ensure_monitoring_directory()
        
        # Save performance metrics
        for name, collector in self.performance_metrics.items():
            file_path = os.path.join(self.config["monitoring_dir"], f"{name}_metrics.json")
            try:
                with open(file_path, 'w') as f:
                    json.dump(collector.get_metrics(), f, indent=2)
            except Exception as e:
                logger.warning(f"Could not save metrics to {file_path}: {str(e)}")
        
        # Save health status
        health_file_path = os.path.join(self.config["monitoring_dir"], "health_status.json")
        try:
            with open(health_file_path, 'w') as f:
                json.dump(self.health_status, f, indent=2)
        except Exception as e:
            logger.warning(f"Could not save health status to {health_file_path}: {str(e)}")
        
        # Save resource metrics
        resource_file_path = os.path.join(self.config["monitoring_dir"], "resource_metrics.json")
        try:
            with open(resource_file_path, 'w') as f:
                json.dump(self.resource_metrics, f, indent=2)
        except Exception as e:
            logger.warning(f"Could not save resource metrics to {resource_file_path}: {str(e)}")
        
        # Save errors
        errors_file_path = os.path.join(self.config["monitoring_dir"], "errors.json")
        try:
            with open(errors_file_path, 'w') as f:
                json.dump({
                    "errors": self.errors[-100:],  # Keep only the last 100 errors
                    "error_counts": self.error_counts
                }, f, indent=2)
        except Exception as e:
            logger.warning(f"Could not save errors to {errors_file_path}: {str(e)}")
        
        # Save alerts
        alerts_file_path = os.path.join(self.config["monitoring_dir"], "alerts.json")
        try:
            with open(alerts_file_path, 'w') as f:
                json.dump({
                    "alerts": self.alerts[-100:],  # Keep only the last 100 alerts
                }, f, indent=2)
        except Exception as e:
            logger.warning(f"Could not save alerts to {alerts_file_path}: {str(e)}")
        
        logger.info("Persisted monitoring metrics and status to disk")
    
    def _ensure_monitoring_directory(self):
        """Create monitoring directory if it doesn't exist."""
        monitoring_dir = self.config["monitoring_dir"]
        if not os.path.exists(monitoring_dir):
            try:
                os.makedirs(monitoring_dir)
                logger.info(f"Created monitoring directory: {monitoring_dir}")
            except Exception as e:
                logger.warning(f"Could not create monitoring directory: {str(e)}")
    
    def record_error(self, error_type, message, stack_trace=None, context=None):
        """
        Record an error.
        
        Args:
            error_type: Type of error
            message: Error message
            stack_trace: Optional stack trace
            context: Optional error context
        """
        timestamp = time.time()
        
        # Create error record
        error = {
            "error_type": error_type,
            "message": message,
            "timestamp": timestamp,
            "formatted_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "stack_trace": stack_trace
        }
        
        # Add correlation ID if available
        if request and g and hasattr(g, 'correlation_id'):
            error["correlation_id"] = g.correlation_id
        
        # Add trace ID if available
        if request and g and hasattr(g, 'trace_id'):
            error["trace_id"] = g.trace_id
        
        # Add additional context if provided
        if context:
            error["context"] = context
        
        # Add to errors list
        self.errors.append(error)
        
        # Keep errors list within limit
        if len(self.errors) > self.max_error_history:
            self.errors = self.errors[-self.max_error_history:]
        
        # Update error counts
        if error_type not in self.error_counts:
            self.error_counts[error_type] = 0
        self.error_counts[error_type] += 1
        
        # Log error
        logger.error(f"{error_type}: {message}")
        
        return error
    
    def _trigger_alert(self, alert_type, message, source, severity="error", context=None):
        """
        Trigger an alert.
        
        Args:
            alert_type: Type of alert
            message: Alert message
            source: Source of the alert (e.g., "database", "api")
            severity: Alert severity (e.g., "info", "warning", "error")
            context: Optional alert context
        """
        # Check if alert is silenced
        if alert_type in self.alert_silence_periods:
            silence_until = self.alert_silence_periods[alert_type]
            if time.time() < silence_until:
                logger.info(f"Alert {alert_type} is silenced until {datetime.fromtimestamp(silence_until)}")
                return
        
        timestamp = time.time()
        
        # Create alert record
        alert = {
            "alert_type": alert_type,
            "message": message,
            "source": source,
            "severity": severity,
            "timestamp": timestamp,
            "formatted_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        # Add additional context if provided
        if context:
            alert["context"] = context
        
        # Add to alerts list
        self.alerts.append(alert)
        
        # Keep alerts list within limit
        if len(self.alerts) > 100:
            self.alerts = self.alerts[-100:]
        
        # Log alert
        logger.warning(f"ALERT - {source} - {severity}: {message}")
        
        # Call alert handlers
        for handler in self.alert_handlers:
            try:
                handler(alert)
            except Exception as e:
                logger.error(f"Error in alert handler: {str(e)}")
        
        return alert
    
    def add_alert_handler(self, handler):
        """
        Add an alert handler.
        
        Args:
            handler: Handler function that takes an alert object
        """
        if callable(handler):
            self.alert_handlers.append(handler)
            return True
        return False
    
    def silence_alert(self, alert_type, duration_seconds=3600):
        """
        Silence an alert for a specified duration.
        
        Args:
            alert_type: Type of alert to silence
            duration_seconds: Duration to silence the alert in seconds
        """
        self.alert_silence_periods[alert_type] = time.time() + duration_seconds
        logger.info(f"Alert {alert_type} silenced for {duration_seconds} seconds")
    
    def unsilence_alert(self, alert_type):
        """
        Unsilence an alert.
        
        Args:
            alert_type: Type of alert to unsilence
        """
        if alert_type in self.alert_silence_periods:
            del self.alert_silence_periods[alert_type]
            logger.info(f"Alert {alert_type} unsilenced")
    
    def get_monitoring_data(self):
        """
        Get current monitoring data.
        
        Returns:
            Dictionary with monitoring data
        """
        return {
            "health": self.health_status,
            "performance": {
                "database": self.performance_metrics["database"].get_metrics(),
                "api": self.performance_metrics["api"].get_metrics(),
                "system": self.performance_metrics["system"].get_metrics()
            },
            "resources": self.resource_metrics,
            "errors": self.errors[-10:],  # Last 10 errors
            "alerts": self.alerts[-10:],  # Last 10 alerts
            "system_info": {
                "hostname": socket.gethostname(),
                "platform": platform.platform(),
                "python_version": platform.python_version(),
                "start_time": psutil.boot_time()
            },
            "last_updated": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
    
    @contextmanager
    def track_operation(self, operation_type, items_count=0, context=None):
        """
        Context manager for tracking operation performance.
        
        Args:
            operation_type: Type of operation
            items_count: Number of items processed
            context: Optional operation context
            
        Yields:
            None
        """
        start_time = time.time()
        operation_id = str(uuid.uuid4())
        success = True
        
        # Generate span ID if tracing is enabled
        span_id = str(uuid.uuid4()) if self.tracing_enabled else None
        
        # Get trace ID from request context or generate a new one
        trace_id = None
        if request and g and hasattr(g, 'trace_id'):
            trace_id = g.trace_id
        elif self.tracing_enabled:
            trace_id = str(uuid.uuid4())
        
        # Get correlation ID from request context or generate a new one
        correlation_id = None
        if request and g and hasattr(g, 'correlation_id'):
            correlation_id = g.correlation_id
        
        # Track active operation
        self.active_operations[operation_id] = {
            "operation_type": operation_type,
            "start_time": start_time,
            "trace_id": trace_id,
            "span_id": span_id,
            "correlation_id": correlation_id,
            "context": context
        }
        
        # If tracing is enabled and we have a trace ID, record the span
        if self.tracing_enabled and trace_id:
            # Create trace if it doesn't exist
            if trace_id not in self.current_traces:
                self.current_traces[trace_id] = {
                    "start_time": start_time,
                    "correlation_id": correlation_id,
                    "spans": {}
                }
            
            # Add span to trace
            self.current_traces[trace_id]["spans"][span_id] = {
                "name": operation_type,
                "start_time": start_time,
                "end_time": None,
                "parent_span_id": getattr(g, 'span_id', None) if request and g else None
            }
        
        try:
            yield
        except Exception as e:
            success = False
            
            # Record error
            self.record_error(
                f"{operation_type}_error",
                str(e),
                traceback.format_exc(),
                context
            )
            
            raise
        finally:
            end_time = time.time()
            execution_time = end_time - start_time
            
            # Record operation performance
            category = "database" if operation_type.startswith(("db_", "database_", "query_")) else "api"
            self.performance_metrics[category].record_operation(
                operation_type=operation_type,
                success=success,
                execution_time=execution_time,
                items_processed=items_count
            )
            
            # Update operation timings for this type
            if operation_type not in self.operation_timings:
                self.operation_timings[operation_type] = {
                    "count": 0,
                    "success_count": 0,
                    "total_time": 0,
                    "min_time": float('inf'),
                    "max_time": 0,
                    "avg_time": 0,
                    "last_execution_time": None
                }
            
            timings = self.operation_timings[operation_type]
            timings["count"] += 1
            if success:
                timings["success_count"] += 1
            timings["total_time"] += execution_time
            timings["min_time"] = min(timings["min_time"], execution_time)
            timings["max_time"] = max(timings["max_time"], execution_time)
            timings["avg_time"] = timings["total_time"] / timings["count"]
            timings["last_execution_time"] = end_time
            
            # If tracing is enabled and we have a trace ID, update the span
            if self.tracing_enabled and trace_id and trace_id in self.current_traces:
                trace = self.current_traces[trace_id]
                if span_id and span_id in trace["spans"]:
                    trace["spans"][span_id]["end_time"] = end_time
                    trace["spans"][span_id]["duration"] = execution_time
                    trace["spans"][span_id]["success"] = success
            
            # Clean up active operation
            if operation_id in self.active_operations:
                del self.active_operations[operation_id]
    
    def track_progress(self, name, total_items, checkpoint_file=None):
        """
        Create a progress tracker for batch operations.
        
        Args:
            name: Name of the operation being tracked
            total_items: Total number of items to process
            checkpoint_file: Optional path to checkpoint file for resumable operations
            
        Returns:
            ProgressTracker instance
        """
        self._ensure_monitoring_directory()
        
        # Create checkpoint directory if needed
        if checkpoint_file:
            checkpoint_dir = os.path.dirname(checkpoint_file)
            if checkpoint_dir and not os.path.exists(checkpoint_dir):
                try:
                    os.makedirs(checkpoint_dir)
                    logger.info(f"Created checkpoint directory: {checkpoint_dir}")
                except Exception as e:
                    logger.warning(f"Could not create checkpoint directory: {str(e)}")
        
        # Create progress tracker
        tracker = ProgressTracker(name, total_items, checkpoint_file)
        
        # Store in progress trackers
        self.progress_trackers[name] = tracker
        
        return tracker
    
    def get_progress_trackers(self):
        """
        Get all progress trackers.
        
        Returns:
            Dictionary of progress trackers
        """
        return {name: tracker.get_status() for name, tracker in self.progress_trackers.items()}


class PerformanceMetrics:
    """
    Collects and stores performance metrics for operations.
    
    This class provides methods to track various performance metrics such as
    execution time, memory usage, and throughput. It also provides
    methods to generate reports and alerts based on these metrics.
    """
    
    def __init__(self, name, metrics_file=None):
        """
        Initialize performance metrics collector.
        
        Args:
            name: Name of the metrics collector (for identification)
            metrics_file: Optional path to file for persisting metrics
        """
        self.name = name
        self.metrics_file = metrics_file
        self.start_time = time.time()
        self.lock = threading.RLock()
        self.metrics = {
            "name": name,
            "start_time": self.start_time,
            "operations": {
                "total": 0,
                "successful": 0,
                "failed": 0
            },
            "timing": {
                "total_execution_time": 0,
                "min_execution_time": float('inf'),
                "max_execution_time": 0,
                "avg_execution_time": 0
            },
            "throughput": {
                "items_processed": 0,
                "items_per_second": 0,
                "bytes_processed": 0,
                "bytes_per_second": 0
            },
            "memory": {
                "peak_memory_usage": 0,
                "current_memory_usage": 0
            },
            "database": {
                "connection_attempts": 0,
                "connection_failures": 0,
                "query_count": 0,
                "transaction_count": 0,
                "rollback_count": 0
            },
            "custom_metrics": {},
            "history": []
        }
        
        # Load existing metrics if file exists
        if metrics_file and os.path.exists(metrics_file):
            try:
                with open(metrics_file, 'r') as f:
                    saved_metrics = json.load(f)
                    # Merge saved metrics with current metrics
                    self._merge_metrics(saved_metrics)
                    logger.info(f"Loaded metrics from {metrics_file}")
            except Exception as e:
                logger.warning(f"Could not load metrics from {metrics_file}: {str(e)}")
    
    def _merge_metrics(self, saved_metrics):
        """
        Merge saved metrics with current metrics.
        
        Args:
            saved_metrics: Dictionary of saved metrics
        """
        with self.lock:
            # Preserve history
            if "history" in saved_metrics:
                self.metrics["history"] = saved_metrics["history"]
            
            # Preserve operation counts
            if "operations" in saved_metrics:
                for key in self.metrics["operations"]:
                    if key in saved_metrics["operations"]:
                        self.metrics["operations"][key] = saved_metrics["operations"][key]
            
            # Preserve throughput
            if "throughput" in saved_metrics:
                for key in self.metrics["throughput"]:
                    if key in saved_metrics["throughput"]:
                        self.metrics["throughput"][key] = saved_metrics["throughput"][key]
            
            # Preserve database metrics
            if "database" in saved_metrics:
                for key in self.metrics["database"]:
                    if key in saved_metrics["database"]:
                        self.metrics["database"][key] = saved_metrics["database"][key]
            
            # Preserve custom metrics
            if "custom_metrics" in saved_metrics:
                self.metrics["custom_metrics"] = saved_metrics["custom_metrics"]
    
    def record_operation(self, operation_type, success, execution_time,
                        items_processed=0, bytes_processed=0):
        """
        Record an operation.
        
        Args:
            operation_type: Type of operation (e.g., "query", "transaction")
            success: Whether the operation was successful
            execution_time: Time taken to execute the operation (in seconds)
            items_processed: Number of items processed in the operation
            bytes_processed: Number of bytes processed in the operation
        """
        with self.lock:
            # Update operation counts
            self.metrics["operations"]["total"] += 1
            if success:
                self.metrics["operations"]["successful"] += 1
            else:
                self.metrics["operations"]["failed"] += 1
            
            # Update timing metrics
            self.metrics["timing"]["total_execution_time"] += execution_time
            self.metrics["timing"]["min_execution_time"] = min(
                self.metrics["timing"]["min_execution_time"], execution_time
            ) if self.metrics["timing"]["min_execution_time"] != float('inf') else execution_time
            self.metrics["timing"]["max_execution_time"] = max(
                self.metrics["timing"]["max_execution_time"], execution_time
            )
            self.metrics["timing"]["avg_execution_time"] = (
                self.metrics["timing"]["total_execution_time"] / 
                self.metrics["operations"]["total"]
            )
            
            # Update throughput metrics
            self.metrics["throughput"]["items_processed"] += items_processed
            self.metrics["throughput"]["bytes_processed"] += bytes_processed
            
            elapsed = time.time() - self.start_time
            if elapsed > 0:
                self.metrics["throughput"]["items_per_second"] = (
                    self.metrics["throughput"]["items_processed"] / elapsed
                )
                self.metrics["throughput"]["bytes_per_second"] = (
                    self.metrics["throughput"]["bytes_processed"] / elapsed
                )
            
            # Update memory metrics
            process = psutil.Process(os.getpid())
            current_memory = process.memory_info().rss
            self.metrics["memory"]["current_memory_usage"] = current_memory
            self.metrics["memory"]["peak_memory_usage"] = max(
                self.metrics["memory"]["peak_memory_usage"], current_memory
            )
            
            # Update database metrics based on operation type
            if operation_type == "connection":
                self.metrics["database"]["connection_attempts"] += 1
                if not success:
                    self.metrics["database"]["connection_failures"] += 1
            elif operation_type == "query":
                self.metrics["database"]["query_count"] += 1
            elif operation_type == "transaction":
                self.metrics["database"]["transaction_count"] += 1
            elif operation_type == "rollback":
                self.metrics["database"]["rollback_count"] += 1
            
            # Add to history (with limited size to prevent memory issues)
            history_entry = {
                "timestamp": time.time(),
                "operation_type": operation_type,
                "success": success,
                "execution_time": execution_time,
                "items_processed": items_processed
            }
            
            self.metrics["history"].append(history_entry)
            if len(self.metrics["history"]) > 1000:  # Keep last 1000 operations
                self.metrics["history"] = self.metrics["history"][-1000:]
            
            # Save metrics to file if specified
            if self.metrics_file:
                try:
                    with open(self.metrics_file, 'w') as f:
                        json.dump(self.metrics, f, indent=2)
                except Exception as e:
                    logger.warning(f"Could not save metrics to {self.metrics_file}: {str(e)}")
    
    def record_custom_metric(self, metric_name, value):
        """
        Record a custom metric.
        
        Args:
            metric_name: Name of the metric
            value: Value of the metric
        """
        with self.lock:
            self.metrics["custom_metrics"][metric_name] = value
            
            # Save metrics to file if specified
            if self.metrics_file:
                try:
                    with open(self.metrics_file, 'w') as f:
                        json.dump(self.metrics, f, indent=2)
                except Exception as e:
                    logger.warning(f"Could not save metrics to {self.metrics_file}: {str(e)}")
    
    def get_metrics(self):
        """
        Get current metrics.
        
        Returns:
            Dictionary of current metrics
        """
        with self.lock:
            # Update elapsed time and throughput calculations
            elapsed = time.time() - self.start_time
            self.metrics["elapsed_time"] = elapsed
            
            if elapsed > 0:
                self.metrics["throughput"]["items_per_second"] = (
                    self.metrics["throughput"]["items_processed"] / elapsed
                )
                self.metrics["throughput"]["bytes_per_second"] = (
                    self.metrics["throughput"]["bytes_processed"] / elapsed
                )
            
            # Update memory metrics
            process = psutil.Process(os.getpid())
            self.metrics["memory"]["current_memory_usage"] = process.memory_info().rss
            
            return self.metrics.copy()
    
    def generate_report(self):
        """
        Generate a human-readable report of the metrics.
        
        Returns:
            String containing the report
        """
        metrics = self.get_metrics()
        elapsed = metrics["elapsed_time"]
        
        report = [
            f"Performance Report for {self.name}",
            f"=======================================",
            f"Duration: {timedelta(seconds=int(elapsed))}",
            f"",
            f"Operations:",
            f"  Total: {metrics['operations']['total']}",
            f"  Successful: {metrics['operations']['successful']}",
            f"  Failed: {metrics['operations']['failed']}",
            f"  Success Rate: {metrics['operations']['successful'] / metrics['operations']['total'] * 100:.2f}% (if total > 0)",
            f"",
            f"Timing:",
            f"  Average Execution Time: {metrics['timing']['avg_execution_time'] * 1000:.2f} ms",
            f"  Min Execution Time: {metrics['timing']['min_execution_time'] * 1000:.2f} ms (if not inf)",
            f"  Max Execution Time: {metrics['timing']['max_execution_time'] * 1000:.2f} ms",
            f"",
            f"Throughput:",
            f"  Items Processed: {metrics['throughput']['items_processed']}",
            f"  Items Per Second: {metrics['throughput']['items_per_second']:.2f}",
            f"  Bytes Processed: {metrics['throughput']['bytes_processed']}",
            f"  Bytes Per Second: {metrics['throughput']['bytes_per_second']:.2f}",
            f"",
            f"Memory:",
            f"  Current Memory Usage: {metrics['memory']['current_memory_usage'] / (1024 * 1024):.2f} MB",
            f"  Peak Memory Usage: {metrics['memory']['peak_memory_usage'] / (1024 * 1024):.2f} MB",
            f"",
            f"Database:",
            f"  Connection Attempts: {metrics['database']['connection_attempts']}",
            f"  Connection Failures: {metrics['database']['connection_failures']}",
            f"  Query Count: {metrics['database']['query_count']}",
            f"  Transaction Count: {metrics['database']['transaction_count']}",
            f"  Rollback Count: {metrics['database']['rollback_count']}",
            f"",
            f"Custom Metrics:"
        ]
        
        for name, value in metrics["custom_metrics"].items():
            report.append(f"  {name}: {value}")
        
        return "\n".join(report)
    
    def reset(self):
        """Reset metrics to initial state."""
        with self.lock:
            self.start_time = time.time()
            self.metrics = {
                "name": self.name,
                "start_time": self.start_time,
                "operations": {
                    "total": 0,
                    "successful": 0,
                    "failed": 0
                },
                "timing": {
                    "total_execution_time": 0,
                    "min_execution_time": float('inf'),
                    "max_execution_time": 0,
                    "avg_execution_time": 0
                },
                "throughput": {
                    "items_processed": 0,
                    "items_per_second": 0,
                    "bytes_processed": 0,
                    "bytes_per_second": 0
                },
                "memory": {
                    "peak_memory_usage": 0,
                    "current_memory_usage": 0
                },
                "database": {
                    "connection_attempts": 0,
                    "connection_failures": 0,
                    "query_count": 0,
                    "transaction_count": 0,
                    "rollback_count": 0
                },
                "custom_metrics": {},
                "history": []
            }
            
            # Save reset metrics to file if specified
            if self.metrics_file:
                try:
                    with open(self.metrics_file, 'w') as f:
                        json.dump(self.metrics, f, indent=2)
                except Exception as e:
                    logger.warning(f"Could not save reset metrics to {self.metrics_file}: {str(e)}")


class ProgressTracker:
    """
    Tracks progress of batch operations with time estimation.
    
    This class provides methods to track progress of batch operations,
    estimate completion time, and generate progress reports.
    """
    
    def __init__(self, name, total_items, checkpoint_file=None):
        """
        Initialize progress tracker.
        
        Args:
            name: Name of the operation being tracked
            total_items: Total number of items to process
            checkpoint_file: Optional path to checkpoint file for resumable operations
        """
        self.name = name
        self.total_items = total_items
        self.checkpoint_file = checkpoint_file
        self.start_time = time.time()
        self.last_update_time = self.start_time
        self.processed_items = 0
        self.successful_items = 0
        self.failed_items = 0
        self.current_batch = 0
        self.total_batches = 0
        self.batch_sizes = []
        self.batch_times = []
        self.lock = threading.RLock()
        
        # Generate a unique ID for this tracking session
        self.tracking_id = str(uuid.uuid4())
        
        # Load checkpoint if exists
        if checkpoint_file and os.path.exists(checkpoint_file):
            try:
                with open(checkpoint_file, 'r') as f:
                    checkpoint = json.load(f)
                    if "processed" in checkpoint:
                        self.processed_items = checkpoint["processed"]
                        if "successful" in checkpoint:
                            self.successful_items = checkpoint["successful"]
                        if "failed" in checkpoint:
                            self.failed_items = checkpoint["failed"]
                        if "last_batch" in checkpoint:
                            self.current_batch = checkpoint["last_batch"]
                        logger.info(f"Loaded progress from checkpoint: {self.processed_items}/{self.total_items} items processed")
            except Exception as e:
                logger.warning(f"Could not load checkpoint from {checkpoint_file}: {str(e)}")
        
        # Initialize metrics
        self.metrics = PerformanceMetrics(f"progress_{name}")
        
        logger.info(f"Progress tracker initialized for {name}: {self.processed_items}/{total_items} items")
    
    def update(self, items_processed, successful=None, failed=None,
              batch_size=None, save_checkpoint=True):
        """
        Update progress.
        
        Args:
            items_processed: Number of items processed in this update
            successful: Number of items successfully processed (default: all)
            failed: Number of items that failed processing (default: 0)
            batch_size: Size of the current batch (for batch operations)
            save_checkpoint: Whether to save checkpoint (if checkpoint_file is set)
            
        Returns:
            Dictionary with current progress status
        """
        with self.lock:
            current_time = time.time()
            
            # Set defaults
            if successful is None:
                successful = items_processed
            if failed is None:
                failed = 0
            
            # Update counters
            self.processed_items += items_processed
            self.successful_items += successful
            self.failed_items += failed
            
            # Ensure we don't exceed total items
            self.processed_items = min(self.processed_items, self.total_items)
            
            # Update batch information
            if batch_size is not None:
                self.current_batch += 1
                self.batch_sizes.append(batch_size)
                self.batch_times.append(current_time - self.last_update_time)
            
            # Calculate progress percentage
            progress_pct = (self.processed_items / self.total_items) * 100 if self.total_items > 0 else 0
            
            # Calculate elapsed time
            elapsed = current_time - self.start_time
            
            # Calculate items per second
            items_per_sec = self.processed_items / elapsed if elapsed > 0 else 0
            
            # Calculate estimated time remaining
            remaining_items = self.total_items - self.processed_items
            est_remaining_time = remaining_items / items_per_sec if items_per_sec > 0 else 0
            
            # Calculate estimated completion time
            est_completion = datetime.now() + timedelta(seconds=est_remaining_time)
            
            # Calculate average batch processing time (from last 10 batches)
            recent_batch_times = self.batch_times[-10:] if self.batch_times else [0]
            avg_batch_time = sum(recent_batch_times) / len(recent_batch_times)
            
            # Calculate average batch size (from last 10 batches)
            recent_batch_sizes = self.batch_sizes[-10:] if self.batch_sizes else [0]
            avg_batch_size = sum(recent_batch_sizes) / len(recent_batch_sizes) if recent_batch_sizes else 0
            
            # Calculate items per batch
            items_per_batch = avg_batch_size
            
            # Calculate batches per second
            batches_per_sec = 1 / avg_batch_time if avg_batch_time > 0 else 0
            
            # Update metrics
            self.metrics.record_operation(
                operation_type="batch",
                success=(failed == 0),
                execution_time=current_time - self.last_update_time,
                items_processed=items_processed
            )
            
            # Save checkpoint if requested
            if save_checkpoint and self.checkpoint_file:
                checkpoint = {
                    "position": self.processed_items,
                    "processed": self.processed_items,
                    "successful": self.successful_items,
                    "failed": self.failed_items,
                    "last_updated": datetime.now().isoformat(),
                    "last_batch": self.current_batch
                }
                
                try:
                    with open(self.checkpoint_file, 'w') as f:
                        json.dump(checkpoint, f)
                except Exception as e:
                    logger.warning(f"Could not save checkpoint to {self.checkpoint_file}: {str(e)}")
            
            # Update last update time
            self.last_update_time = current_time
            
            # Create progress status
            status = {
                "name": self.name,
                "tracking_id": self.tracking_id,
                "total_items": self.total_items,
                "processed_items": self.processed_items,
                "successful_items": self.successful_items,
                "failed_items": self.failed_items,
                "progress_percentage": progress_pct,
                "elapsed_time": elapsed,
                "items_per_second": items_per_sec,
                "estimated_remaining_time": est_remaining_time,
                "estimated_completion": est_completion.isoformat(),
                "current_batch": self.current_batch,
                "avg_batch_time": avg_batch_time,
                "avg_batch_size": avg_batch_size,
                "items_per_batch": items_per_batch,
                "batches_per_second": batches_per_sec,
                "last_update_time": current_time
            }
            
            # Log progress
            logger.info(f"Progress [{self.name}]: {self.processed_items}/{self.total_items} items "
                       f"({progress_pct:.1f}%) at {items_per_sec:.1f} items/sec, "
                       f"ETA: {est_completion.strftime('%Y-%m-%d %H:%M:%S')}")
            
            return status
    
    def get_status(self):
        """
        Get current progress status.
        
        Returns:
            Dictionary with current progress status
        """
        with self.lock:
            current_time = time.time()
            
            # Calculate progress percentage
            progress_pct = (self.processed_items / self.total_items) * 100 if self.total_items > 0 else 0
            
            # Calculate elapsed time
            elapsed = current_time - self.start_time
            
            # Calculate items per second
            items_per_sec = self.processed_items / elapsed if elapsed > 0 else 0
            
            # Calculate estimated time remaining
            remaining_items = self.total_items - self.processed_items
            est_remaining_time = remaining_items / items_per_sec if items_per_sec > 0 else 0
            
            # Calculate estimated completion time
            est_completion = datetime.now() + timedelta(seconds=est_remaining_time)
            
            # Calculate average batch processing time (from last 10 batches)
            recent_batch_times = self.batch_times[-10:] if self.batch_times else [0]
            avg_batch_time = sum(recent_batch_times) / len(recent_batch_times)
            
            # Calculate average batch size (from last 10 batches)
            recent_batch_sizes = self.batch_sizes[-10:] if self.batch_sizes else [0]
            avg_batch_size = sum(recent_batch_sizes) / len(recent_batch_sizes) if recent_batch_sizes else 0
            
            # Calculate items per batch
            items_per_batch = avg_batch_size
            
            # Calculate batches per second
            batches_per_sec = 1 / avg_batch_time if avg_batch_time > 0 else 0
            
            # Create progress status
            return {
                "name": self.name,
                "tracking_id": self.tracking_id,
                "total_items": self.total_items,
                "processed_items": self.processed_items,
                "successful_items": self.successful_items,
                "failed_items": self.failed_items,
                "progress_percentage": progress_pct,
                "elapsed_time": elapsed,
                "items_per_second": items_per_sec,
                "estimated_remaining_time": est_remaining_time,
                "estimated_completion": est_completion.isoformat(),
                "current_batch": self.current_batch,
                "avg_batch_time": avg_batch_time,
                "avg_batch_size": avg_batch_size,
                "items_per_batch": items_per_batch,
                "batches_per_second": batches_per_sec
            }
    
    def generate_report(self):
        """
        Generate a human-readable progress report.
        
        Returns:
            String containing the progress report
        """
        status = self.get_status()
        
        report = [
            f"Progress Report for {self.name}",
            f"=======================================",
            f"Status: {status['processed_items']}/{status['total_items']} items processed "
            f"({status['progress_percentage']:.1f}%)",
            f"",
            f"Time:",
            f"  Elapsed: {timedelta(seconds=int(status['elapsed_time']))}",
            f"  Estimated Remaining: {timedelta(seconds=int(status['estimated_remaining_time']))}",
            f"  Estimated Completion: {datetime.fromisoformat(status['estimated_completion']).strftime('%Y-%m-%d %H:%M:%S')}",
            f"",
            f"Performance:",
            f"  Items Per Second: {status['items_per_second']:.1f}",
            f"  Batches Per Second: {status['batches_per_second']:.2f}",
            f"  Average Batch Size: {status['avg_batch_size']:.1f} items",
            f"  Average Batch Time: {status['avg_batch_time'] * 1000:.1f} ms",
            f"",
            f"Results:",
            f"  Successful: {status['successful_items']} items",
            f"  Failed: {status['failed_items']} items"
        ]
        
        if status['processed_items'] > 0:
            report.append(f"  Success Rate: {status['successful_items'] / status['processed_items'] * 100:.1f}%")
        
        return "\n".join(report)
    
    def reset(self):
        """Reset progress tracker to initial state."""
        with self.lock:
            self.start_time = time.time()
            self.last_update_time = self.start_time
            self.processed_items = 0
            self.successful_items = 0
            self.failed_items = 0
            self.current_batch = 0
            self.batch_sizes = []
            self.batch_times = []
            
            # Generate a new tracking ID
            self.tracking_id = str(uuid.uuid4())
            
            # Reset metrics
            self.metrics.reset()
            
            # Save reset checkpoint if specified
            if self.checkpoint_file:
                checkpoint = {
                    "position": 0,
                    "processed": 0,
                    "successful": 0,
                    "failed": 0,
                    "last_updated": datetime.now().isoformat(),
                    "last_batch": 0
                }
                
                try:
                    with open(self.checkpoint_file, 'w') as f:
                        json.dump(checkpoint, f)
                except Exception as e:
                    logger.warning(f"Could not save reset checkpoint to {self.checkpoint_file}: {str(e)}")
            
            logger.info(f"Progress tracker reset for {self.name}")


# Utility functions

def create_monitoring_service(app=None):
    """
    Create or get singleton monitoring service.
    
    Args:
        app: Optional Flask application to initialize with
        
    Returns:
        MonitoringService instance
    """
    if app:
        # Create a new monitoring service with the app
        try:
            return MonitoringService(app)
        except RuntimeError:
            # Instance already exists, get it and init with app
            monitor = MonitoringService.get_instance()
            monitor.init_app(app)
            return monitor
    else:
        # Just get the instance
        return MonitoringService.get_instance()


def start_monitoring(app=None, dashboard_port=5001):
    """
    Start all monitoring services.
    
    Args:
        app: Optional Flask application to initialize with
        dashboard_port: Port for the monitoring dashboard
    """
    monitor = create_monitoring_service(app)
    monitor.config["dashboard_port"] = dashboard_port
    monitor.start()
    
    logger.info(f"Started monitoring services with dashboard on port {dashboard_port}")
    
    return monitor


def track_operation(operation_type, items_count=0, context=None):
    """
    Decorator for tracking operation performance.
    
    Args:
        operation_type: Type of operation
        items_count: Number of items processed
        context: Optional operation context
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            monitor = MonitoringService.get_instance()
            
            with monitor.track_operation(operation_type, items_count, context):
                return func(*args, **kwargs)
                
        return wrapper
    
    return decorator


def track_database_operation(items_count=0):
    """
    Decorator for tracking database operation performance.
    
    Args:
        items_count: Number of items processed
    """
    def decorator(func):
        func_name = func.__name__
        
        @wraps(func)
        def wrapper(*args, **kwargs):
            module_name = func.__module__
            operation_type = f"db_{module_name}.{func_name}"
            
            monitor = MonitoringService.get_instance()
            
            with monitor.track_operation(operation_type, items_count):
                return func(*args, **kwargs)
                
        return wrapper
    
    return decorator


def track_api_operation(name=None):
    """
    Decorator for tracking API operation performance.
    
    Args:
        name: Optional name for the operation
    """
    def decorator(func):
        func_name = func.__name__
        
        @wraps(func)
        def wrapper(*args, **kwargs):
            operation_type = name or f"api_{func_name}"
            
            monitor = MonitoringService.get_instance()
            
            with monitor.track_operation(operation_type):
                return func(*args, **kwargs)
                
        return wrapper
    
    return decorator


# Examples of usage

def example_flask_app():
    """Example of using MonitoringService with Flask."""
    if Flask is None:
        logger.warning("Flask is not available, skipping example")
        return
    
    app = Flask(__name__)
    
    # Initialize monitoring
    monitor = start_monitoring(app)
    
    # Example route with tracking
    @app.route('/example')
    @track_api_operation()
    def example_route():
        return jsonify({"status": "ok"})
    
    # Example database function with tracking
    @track_database_operation(items_count=10)
    def example_db_function():
        # Simulate database operation
        time.sleep(0.1)
        return {"result": "success"}
    
    # Example of manual tracking
    with monitor.track_operation("custom_operation", items_count=5):
        # Your code here
        time.sleep(0.1)
    
    return app


# If this module is executed directly, run the example
if __name__ == "__main__":
    # Ensure this script is not accidentally run in production
    if os.environ.get("FLASK_ENV") == "production":
        print("This script should not be executed directly in production environment")
        sys.exit(1)
    
    # Run example Flask app if available
    if Flask is not None:
        app = example_flask_app()
        app.run(debug=True, port=5000)
    else:
        # Run non-Flask example
        monitor = start_monitoring()
        
        # Example of tracking progress
        progress = monitor.track_progress("example_operation", 100)
        
        for i in range(0, 100, 10):
            # Simulate processing
            time.sleep(0.5)
            
            # Update progress
            progress.update(10)
        
        # Print report
        print(monitor.get_monitoring_data())
        print(progress.generate_report())
        
        # Clean up
        monitor.stop()