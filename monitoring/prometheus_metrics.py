#!/usr/bin/env python3
"""
CryoProtect v2 - Prometheus Metrics Collection

This module provides Prometheus metrics collection for the CryoProtect API.
It defines custom metrics and exposes them via a /metrics endpoint for Prometheus scraping.

Metrics collected:
- API request count (by endpoint, method, status code)
- API request latency (by endpoint, method)
- Database query count and latency
- Error count (by type)
- System resource usage (CPU, memory, disk)
"""

import time
import logging
import psutil
from flask import request, Response, Blueprint
from prometheus_client import Counter, Histogram, Gauge, Summary, generate_latest, REGISTRY, CONTENT_TYPE_LATEST

# Set up logging
logger = logging.getLogger(__name__)

# Create a Blueprint for the metrics endpoint
metrics_bp = Blueprint('metrics', __name__)

# Define Prometheus metrics
# API metrics
REQUEST_COUNT = Counter(
    'cryoprotect_http_requests_total',
    'Total number of HTTP requests',
    ['method', 'endpoint', 'status']
)

REQUEST_LATENCY = Histogram(
    'cryoprotect_http_request_duration_seconds',
    'HTTP request latency in seconds',
    ['method', 'endpoint'],
    buckets=(0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0, float('inf'))
)

# Database metrics
DB_QUERY_COUNT = Counter(
    'cryoprotect_db_queries_total',
    'Total number of database queries',
    ['operation', 'table']
)

DB_QUERY_LATENCY = Histogram(
    'cryoprotect_db_query_duration_seconds',
    'Database query latency in seconds',
    ['operation', 'table'],
    buckets=(0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, float('inf'))
)

# Error metrics
ERROR_COUNT = Counter(
    'cryoprotect_errors_total',
    'Total number of errors',
    ['type', 'endpoint']
)

# System resource metrics
CPU_USAGE = Gauge(
    'cryoprotect_cpu_usage_percent',
    'CPU usage in percent'
)

MEMORY_USAGE = Gauge(
    'cryoprotect_memory_usage_bytes',
    'Memory usage in bytes',
    ['type']
)

DISK_USAGE = Gauge(
    'cryoprotect_disk_usage_bytes',
    'Disk usage in bytes',
    ['type']
)

# Connection pool metrics
CONNECTION_POOL_SIZE = Gauge(
    'cryoprotect_db_connection_pool_size',
    'Database connection pool size',
    ['state']
)

# Custom business metrics
MOLECULE_COUNT = Gauge(
    'cryoprotect_molecules_total',
    'Total number of molecules in the system'
)

MIXTURE_COUNT = Gauge(
    'cryoprotect_mixtures_total',
    'Total number of mixtures in the system'
)

PREDICTION_COUNT = Gauge(
    'cryoprotect_predictions_total',
    'Total number of predictions in the system'
)

EXPERIMENT_COUNT = Gauge(
    'cryoprotect_experiments_total',
    'Total number of experiments in the system'
)

# Request context metrics
REQUEST_IN_PROGRESS = Gauge(
    'cryoprotect_http_requests_in_progress',
    'Number of HTTP requests currently in progress',
    ['method', 'endpoint']
)

# Define a decorator to track request metrics
def track_request_metrics():
    """
    Decorator to track request metrics.
    
    This decorator should be applied to Flask request handlers to track
    request count, latency, and other metrics.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            # Extract endpoint and method
            endpoint = request.endpoint
            method = request.method
            
            # Track request in progress
            REQUEST_IN_PROGRESS.labels(method=method, endpoint=endpoint).inc()
            
            # Track request latency
            start_time = time.time()
            
            try:
                # Execute the request handler
                response = func(*args, **kwargs)
                
                # Get status code
                status_code = response.status_code
                
                # Track request count
                REQUEST_COUNT.labels(method=method, endpoint=endpoint, status=status_code).inc()
                
                return response
            except Exception as e:
                # Track error
                ERROR_COUNT.labels(type=type(e).__name__, endpoint=endpoint).inc()
                
                # Re-raise the exception
                raise
            finally:
                # Track request latency
                latency = time.time() - start_time
                REQUEST_LATENCY.labels(method=method, endpoint=endpoint).observe(latency)
                
                # Decrement request in progress
                REQUEST_IN_PROGRESS.labels(method=method, endpoint=endpoint).dec()
        
        # Preserve function metadata
        wrapper.__name__ = func.__name__
        wrapper.__doc__ = func.__doc__
        
        return wrapper
    
    return decorator

# Database query tracking context manager
class DatabaseQueryMetrics:
    """
    Context manager to track database query metrics.
    
    Usage:
        with DatabaseQueryMetrics('select', 'molecules'):
            # Execute database query
    """
    def __init__(self, operation, table):
        self.operation = operation
        self.table = table
        self.start_time = None
    
    def __enter__(self):
        self.start_time = time.time()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        latency = time.time() - self.start_time
        
        # Track query count
        DB_QUERY_COUNT.labels(operation=self.operation, table=self.table).inc()
        
        # Track query latency
        DB_QUERY_LATENCY.labels(operation=self.operation, table=self.table).observe(latency)
        
        # Track error if any
        if exc_type is not None:
            ERROR_COUNT.labels(type=exc_type.__name__, endpoint='database').inc()

# Function to update system resource metrics
def update_system_metrics():
    """Update system resource metrics."""
    try:
        # Update CPU usage
        cpu_percent = psutil.cpu_percent(interval=0.1)
        CPU_USAGE.set(cpu_percent)
        
        # Update memory usage
        memory = psutil.virtual_memory()
        MEMORY_USAGE.labels(type='total').set(memory.total)
        MEMORY_USAGE.labels(type='available').set(memory.available)
        MEMORY_USAGE.labels(type='used').set(memory.used)
        
        # Update disk usage
        disk = psutil.disk_usage('/')
        DISK_USAGE.labels(type='total').set(disk.total)
        DISK_USAGE.labels(type='free').set(disk.free)
        DISK_USAGE.labels(type='used').set(disk.used)
    except Exception as e:
        logger.error(f"Error updating system metrics: {str(e)}")

# Function to update connection pool metrics
def update_connection_pool_metrics(active, idle, max_connections):
    """
    Update connection pool metrics.
    
    Args:
        active: Number of active connections
        idle: Number of idle connections
        max_connections: Maximum number of connections
    """
    CONNECTION_POOL_SIZE.labels(state='active').set(active)
    CONNECTION_POOL_SIZE.labels(state='idle').set(idle)
    CONNECTION_POOL_SIZE.labels(state='max').set(max_connections)

# Function to update business metrics
def update_business_metrics(molecule_count, mixture_count, prediction_count, experiment_count):
    """
    Update business metrics.
    
    Args:
        molecule_count: Number of molecules
        mixture_count: Number of mixtures
        prediction_count: Number of predictions
        experiment_count: Number of experiments
    """
    MOLECULE_COUNT.set(molecule_count)
    MIXTURE_COUNT.set(mixture_count)
    PREDICTION_COUNT.set(prediction_count)
    EXPERIMENT_COUNT.set(experiment_count)

# Define the metrics endpoint
@metrics_bp.route('/metrics')
def metrics():
    """
    Expose Prometheus metrics.
    
    This endpoint is scraped by Prometheus to collect metrics.
    """
    # Update system metrics before generating the response
    update_system_metrics()
    
    # Generate and return metrics
    return Response(generate_latest(REGISTRY), mimetype=CONTENT_TYPE_LATEST)

def init_app(app):
    """
    Initialize the metrics collection for the Flask app.
    
    Args:
        app: Flask application instance
    """
    # Register the metrics blueprint
    app.register_blueprint(metrics_bp)
    
    # Log initialization
    logger.info("Prometheus metrics collection initialized")