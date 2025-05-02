#!/usr/bin/env python3
"""
CryoProtect v2 - Monitoring Middleware

This module provides middleware for integrating monitoring with the Flask application.
It includes middleware for tracking request metrics, database queries, and other metrics.
"""

import time
import logging
from flask import request, g
from functools import wraps

from monitoring.prometheus_metrics import (
    REQUEST_COUNT, REQUEST_LATENCY, ERROR_COUNT, REQUEST_IN_PROGRESS,
    update_system_metrics, update_connection_pool_metrics, update_business_metrics
)

# Set up logging
logger = logging.getLogger(__name__)

class PrometheusMiddleware:
    """
    Middleware for tracking Prometheus metrics in Flask applications.
    
    This middleware tracks request count, latency, and other metrics for all endpoints.
    """
    
    def __init__(self, app=None):
        self.app = app
        if app is not None:
            self.init_app(app)
    
    def init_app(self, app):
        """
        Initialize the middleware with a Flask application.
        
        Args:
            app: Flask application instance
        """
        # Register before_request handler
        @app.before_request
        def before_request():
            # Store request start time (using datetime for consistency with logging_enhanced.py)
            from datetime import datetime
            g.start_time = datetime.utcnow()
            
            # Track request in progress
            endpoint = request.endpoint or 'unknown'
            method = request.method
            REQUEST_IN_PROGRESS.labels(method=method, endpoint=endpoint).inc()
        
        # Register after_request handler
        @app.after_request
        def after_request(response):
            # Skip metrics endpoint to avoid circular tracking
            if request.path == '/metrics':
                return response
            
            # Calculate request duration
            if hasattr(g, 'start_time'):
                # Use datetime for consistency with logging_enhanced.py
                from datetime import datetime
                duration = (datetime.utcnow() - g.start_time).total_seconds()
                
                # Get endpoint and method
                endpoint = request.endpoint or 'unknown'
                method = request.method
                
                # Track request count
                REQUEST_COUNT.labels(
                    method=method,
                    endpoint=endpoint,
                    status=response.status_code
                ).inc()
                
                # Track request latency
                REQUEST_LATENCY.labels(
                    method=method,
                    endpoint=endpoint
                ).observe(duration)
                
                # Decrement request in progress
                REQUEST_IN_PROGRESS.labels(method=method, endpoint=endpoint).dec()
            
            return response
        
        # Register error handler
        @app.errorhandler(Exception)
        def handle_exception(error):
            # Track error
            endpoint = request.endpoint or 'unknown'
            ERROR_COUNT.labels(
                type=type(error).__name__,
                endpoint=endpoint
            ).inc()
            
            # Let the default error handler handle the exception
            raise error
        
        # Log initialization
        logger.info("Prometheus middleware initialized")

def track_database_metrics(func):
    """
    Decorator to track database metrics.
    
    This decorator should be applied to database query functions to track
    query count, latency, and other metrics.
    
    Args:
        func: Function to decorate
        
    Returns:
        Decorated function
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        from monitoring.prometheus_metrics import DB_QUERY_COUNT, DB_QUERY_LATENCY, ERROR_COUNT
        
        # Extract operation and table from function name or arguments
        operation = kwargs.get('operation', func.__name__.split('_')[0])
        table = kwargs.get('table', 'unknown')
        
        # Track query latency
        start_time = time.time()
        
        try:
            # Execute the query
            result = func(*args, **kwargs)
            
            return result
        except Exception as e:
            # Track error
            ERROR_COUNT.labels(type=type(e).__name__, endpoint='database').inc()
            
            # Re-raise the exception
            raise
        finally:
            # Track query count and latency
            latency = time.time() - start_time
            DB_QUERY_COUNT.labels(operation=operation, table=table).inc()
            DB_QUERY_LATENCY.labels(operation=operation, table=table).observe(latency)
    
    return wrapper

def update_connection_pool_status(active, idle, max_connections):
    """
    Update connection pool status metrics.
    
    Args:
        active: Number of active connections
        idle: Number of idle connections
        max_connections: Maximum number of connections
    """
    update_connection_pool_metrics(active, idle, max_connections)

def update_entity_counts(supabase):
    """
    Update entity count metrics.
    
    Args:
        supabase: Supabase client instance
    """
    try:
        # Get molecule count
        molecule_response = supabase.from_("molecules").select("id", count="exact").execute()
        molecule_count = molecule_response.count if hasattr(molecule_response, 'count') else 0
        
        # Get mixture count
        mixture_response = supabase.from_("mixtures").select("id", count="exact").execute()
        mixture_count = mixture_response.count if hasattr(mixture_response, 'count') else 0
        
        # Get prediction count
        prediction_response = supabase.from_("predictions").select("id", count="exact").execute()
        prediction_count = prediction_response.count if hasattr(prediction_response, 'count') else 0
        
        # Get experiment count
        experiment_response = supabase.from_("experiments").select("id", count="exact").execute()
        experiment_count = experiment_response.count if hasattr(experiment_response, 'count') else 0
        
        # Update metrics
        update_business_metrics(
            molecule_count=molecule_count,
            mixture_count=mixture_count,
            prediction_count=prediction_count,
            experiment_count=experiment_count
        )
    except Exception as e:
        logger.error(f"Error updating entity counts: {str(e)}")