"""
CryoProtect v2 - Monitoring Package

This package provides monitoring and metrics collection for the CryoProtect API.
It includes Prometheus metrics collection, middleware for tracking requests,
and utilities for monitoring system resources.
"""

from monitoring.prometheus_metrics import init_app as init_metrics
from monitoring.middleware import PrometheusMiddleware

__all__ = ['init_metrics', 'PrometheusMiddleware']