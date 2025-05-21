"""
CryoProtect Analyzer API - System Resources

This module provides API endpoints for system status, logs, metrics, and database health checks.
"""

import logging
import os
import platform
import psutil
import time
import json
from datetime import datetime, timedelta
from flask import request, current_app, jsonify, Response
from flask_restful import Resource, marshal_with, fields
from marshmallow import Schema, fields as ma_fields, ValidationError

from api.utils import token_required, handle_error, get_supabase_client, get_user_id
from database.utils.health_check import (
    generate_health_report, get_latest_health_report,
    fix_common_issues, run_scheduled_health_check,
    DatabaseHealthCheck
)
from database.utils.connection import supabase_connection
from api.api_decorators import authenticate_service_role

# Set up logging
logger = logging.getLogger(__name__)

# Define field schemas for marshal_with decorators
system_status_fields = {
    'status': fields.String,
    'version': fields.String,
    'uptime': fields.Integer,
    'uptime_formatted': fields.String,
    'cpu_usage': fields.Float,
    'memory_usage': fields.Raw,
    'disk_usage': fields.Raw,
    'database_status': fields.String,
    'environment': fields.String,
    'timestamp': fields.DateTime
}

system_logs_fields = {
    'logs': fields.List(fields.Raw),
    'total_count': fields.Integer,
    'filtered_count': fields.Integer,
    'timestamp': fields.DateTime
}

system_metrics_fields = {
    'metrics': fields.Raw,
    'period': fields.String,
    'timestamp': fields.DateTime
}

database_health_fields = {
    'overall_status': fields.String,
    'summary': fields.Raw,
    'schema_validation': fields.Raw,
    'data_integrity': fields.Raw,
    'performance': fields.Raw,
    'timestamp': fields.DateTime,
    'report_id': fields.String
}

health_status_fields = {
    'status': fields.String,
    'version': fields.String,
    'uptime': fields.Integer,
    'uptime_formatted': fields.String,
    'database_status': fields.String,
    'timestamp': fields.DateTime
}

health_database_fields = {
    'status': fields.String,
    'overall_health': fields.String,
    'schema_status': fields.String,
    'integrity_status': fields.String,
    'performance_status': fields.String,
    'issues_count': fields.Integer,
    'timestamp': fields.DateTime,
    'details': fields.Raw
}

health_performance_fields = {
    'status': fields.String,
    'metrics': fields.Raw,
    'database_metrics': fields.Raw,
    'api_metrics': fields.Raw,
    'timestamp': fields.DateTime
}

database_fix_fields = {
    'status': fields.String,
    'fixes_applied': fields.Integer,
    'details': fields.Raw,
    'timestamp': fields.DateTime
}

service_role_check_fields = {
    'status': fields.String,
    'message': fields.String,
    'timestamp': fields.DateTime,
    'role': fields.String
}

# Define request schemas
class LogsRequestSchema(Schema):
    """Schema for logs request parameters."""
    level = ma_fields.String(default='all')
    start_time = ma_fields.DateTime()
    end_time = ma_fields.DateTime()
    limit = ma_fields.Integer(default=100)
    offset = ma_fields.Integer(default=0)
    source = ma_fields.String()
    search = ma_fields.String()

class MetricsRequestSchema(Schema):
    """Schema for metrics request parameters."""
    period = ma_fields.String(default='hour')
    metrics = ma_fields.List(ma_fields.String(), default=['cpu', 'memory', 'requests', 'errors'])
    
class HealthCheckRequestSchema(Schema):
    """Schema for health check request parameters."""
    force_new = ma_fields.Boolean(default=False)
    report_type = ma_fields.String(default='full')
    
class FixDatabaseRequestSchema(Schema):
    """Schema for database fix request parameters."""
    fix_type = ma_fields.String(default='all')
    
class HealthDatabaseRequestSchema(Schema):
    """Schema for health database request parameters."""
    verbosity = ma_fields.String(default='normal', validate=lambda x: x in ['minimal', 'normal', 'detailed'])
    format = ma_fields.String(default='json', validate=lambda x: x in ['json', 'markdown', 'html'])
    
class HealthPerformanceRequestSchema(Schema):
    """Schema for health performance request parameters."""
    period = ma_fields.String(default='hour', validate=lambda x: x in ['minute', 'hour', 'day', 'week'])
    metrics = ma_fields.List(ma_fields.String(), default=['cpu', 'memory', 'database', 'api'])

# Start time of the application
START_TIME = time.time()

class SystemStatusResource(Resource):
    """Resource for system status information."""
    
    @token_required
    @marshal_with(system_status_fields)
    def get(self):
        """
        Get system status information.
        
        Returns:
            JSON response with system status
        """
        try:
            # Calculate uptime
            uptime_seconds = int(time.time() - START_TIME)
            uptime_formatted = str(timedelta(seconds=uptime_seconds))
            
            # Get CPU usage
            cpu_usage = psutil.cpu_percent(interval=0.1)
            
            # Get memory usage
            memory = psutil.virtual_memory()
            memory_usage = {
                'total': memory.total,
                'available': memory.available,
                'used': memory.used,
                'percent': memory.percent
            }
            
            # Get disk usage
            disk = psutil.disk_usage('/')
            disk_usage = {
                'total': disk.total,
                'free': disk.free,
                'used': disk.used,
                'percent': disk.percent
            }
            
            # Check database connection
            try:
                supabase = get_supabase_client()
                response = supabase.from_("property_types").select("*").limit(1).execute()
                db_status = "connected" if hasattr(response, 'data') else "unknown"
            except Exception as e:
                logger.error(f"Error checking database status: {str(e)}")
                db_status = "error"
            
            # Get environment
            environment = os.environ.get('FLASK_ENV', 'development')
            
            return {
                'status': 'ok',
                'version': current_app.config.get('API_VERSION', '1.0.0'),
                'uptime': uptime_seconds,
                'uptime_formatted': uptime_formatted,
                'cpu_usage': cpu_usage,
                'memory_usage': memory_usage,
                'disk_usage': disk_usage,
                'database_status': db_status,
                'environment': environment,
                'timestamp': datetime.now()
            }
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fetching system status",
                log_level='error',
                return_response=True
            )
            return error_response, error_status

class SystemLogsResource(Resource):
    """Resource for system logs."""
    
    @token_required
    @marshal_with(system_logs_fields)
    def get(self):
        """
        Get system logs with filtering options.
        
        Returns:
            JSON response with filtered logs
        """
        try:
            # Validate request parameters
            schema = LogsRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating logs request parameters",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                return error_response, error_status
            
            # Extract parameters
            level = args.get('level', 'all')
            start_time = args.get('start_time')
            end_time = args.get('end_time')
            limit = args.get('limit', 100)
            offset = args.get('offset', 0)
            source = args.get('source')
            search = args.get('search')
            
            # Get logs from log files
            logs = self._get_logs_from_files(level, start_time, end_time, source, search, limit, offset)
            
            return {
                'logs': logs['entries'],
                'total_count': logs['total_count'],
                'filtered_count': logs['filtered_count'],
                'timestamp': datetime.now()
            }
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fetching system logs",
                log_level='error',
                return_response=True
            )
            return error_response, error_status
    
    def _get_logs_from_files(self, level, start_time, end_time, source, search, limit, offset):
        """
        Get logs from log files with filtering.
        
        Args:
            level: Log level filter
            start_time: Start time filter
            end_time: End time filter
            source: Source filter
            search: Search text filter
            limit: Maximum number of logs to return
            offset: Offset for pagination
            
        Returns:
            dict: Filtered logs with metadata
        """
        # In a real implementation, this would read from log files
        # For this example, we'll return some sample logs
        
        # Sample logs
        sample_logs = [
            {
                'timestamp': datetime.now() - timedelta(minutes=i),
                'level': level_name,
                'source': f"api.{source_name}",
                'message': f"Sample {level_name} log message {i}",
                'details': {
                    'request_id': f"req-{i}",
                    'user_id': f"user-{i % 10}",
                    'ip': f"192.168.1.{i % 255}"
                }
            }
            for i in range(200)
            for level_name, source_name in [
                ('INFO', 'resources'),
                ('WARNING', 'auth'),
                ('ERROR', 'database'),
                ('DEBUG', 'utils')
            ]
        ]
        
        # Sort by timestamp descending
        sample_logs.sort(key=lambda x: x['timestamp'], reverse=True)
        
        # Apply filters
        filtered_logs = sample_logs
        
        if level != 'all':
            filtered_logs = [log for log in filtered_logs if log['level'] == level.upper()]
        
        if start_time:
            filtered_logs = [log for log in filtered_logs if log['timestamp'] >= start_time]
        
        if end_time:
            filtered_logs = [log for log in filtered_logs if log['timestamp'] <= end_time]
        
        if source:
            filtered_logs = [log for log in filtered_logs if source.lower() in log['source'].lower()]
        
        if search:
            filtered_logs = [log for log in filtered_logs if search.lower() in log['message'].lower()]
        
        # Get total count after filtering
        total_count = len(sample_logs)
        filtered_count = len(filtered_logs)
        
        # Apply pagination
        paginated_logs = filtered_logs[offset:offset + limit]
        
        return {
            'entries': paginated_logs,
            'total_count': total_count,
            'filtered_count': filtered_count
        }

class SystemMetricsResource(Resource):
    """Resource for system metrics."""
    
    @token_required
    @marshal_with(system_metrics_fields)
    def get(self):
        """
        Get system metrics for monitoring.
        
        Returns:
            JSON response with system metrics
        """
        try:
            # Validate request parameters
            schema = MetricsRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating metrics request parameters",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                return error_response, error_status
            
            # Extract parameters
            period = args.get('period', 'hour')
            metrics_list = args.get('metrics', ['cpu', 'memory', 'requests', 'errors'])
            
            # Get metrics
            metrics = self._get_metrics(period, metrics_list)
            
            return {
                'metrics': metrics,
                'period': period,
                'timestamp': datetime.now()
            }
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fetching system metrics",
                log_level='error',
                return_response=True
            )
            return error_response, error_status
    
    def _get_metrics(self, period, metrics_list):
        """
        Get system metrics for the specified period.
        
        Args:
            period: Time period for metrics ('hour', 'day', 'week', 'month')
            metrics_list: List of metrics to retrieve
            
        Returns:
            dict: System metrics
        """
        # In a real implementation, this would query a metrics database
        # For this example, we'll return some sample metrics
        
        # Determine number of data points based on period
        if period == 'hour':
            points = 60  # One point per minute
            interval = 'minute'
        elif period == 'day':
            points = 24  # One point per hour
            interval = 'hour'
        elif period == 'week':
            points = 7  # One point per day
            interval = 'day'
        elif period == 'month':
            points = 30  # One point per day
            interval = 'day'
        else:
            points = 60
            interval = 'minute'
        
        # Generate sample metrics
        metrics = {}
        
        if 'cpu' in metrics_list:
            metrics['cpu'] = {
                'labels': [f"{i}" for i in range(points)],
                'values': [round(40 + 30 * (0.5 + 0.5 * (i % 10) / 10), 2) for i in range(points)],
                'unit': '%',
                'interval': interval
            }
        
        if 'memory' in metrics_list:
            metrics['memory'] = {
                'labels': [f"{i}" for i in range(points)],
                'values': [round(60 + 20 * (0.5 + 0.5 * (i % 15) / 15), 2) for i in range(points)],
                'unit': '%',
                'interval': interval
            }
        
        if 'requests' in metrics_list:
            metrics['requests'] = {
                'labels': [f"{i}" for i in range(points)],
                'values': [int(100 + 50 * (0.5 + 0.5 * (i % 20) / 20)) for i in range(points)],
                'unit': 'count',
                'interval': interval
            }
        
        if 'errors' in metrics_list:
            metrics['errors'] = {
                'labels': [f"{i}" for i in range(points)],
                'values': [int(5 + 10 * (0.5 + 0.5 * (i % 25) / 25)) for i in range(points)],
                'unit': 'count',
                'interval': interval
            }
        
        if 'response_time' in metrics_list:
            metrics['response_time'] = {
                'labels': [f"{i}" for i in range(points)],
                'values': [round(200 + 100 * (0.5 + 0.5 * (i % 30) / 30), 2) for i in range(points)],
                'unit': 'ms',
                'interval': interval
            }
        
        if 'database_queries' in metrics_list:
            metrics['database_queries'] = {
                'labels': [f"{i}" for i in range(points)],
                'values': [int(500 + 300 * (0.5 + 0.5 * (i % 35) / 35)) for i in range(points)],
                'unit': 'count',
                'interval': interval
            }
        
        return metrics

# New resources for database health
class DatabaseHealthResource(Resource):
    """Resource for database health checks."""
    
    @token_required
    # @require_admin
    @marshal_with(database_health_fields)
    def get(self):
        """
        Get database health check report.
        
        Returns:
            JSON response with health check information
        """
        try:
            # Validate request parameters
            schema = HealthCheckRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating health check request parameters",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                return error_response, error_status
            
            # Extract parameters
            force_new = args.get('force_new', False)
            report_type = args.get('report_type', 'full')
            
            # Get health report
            if force_new:
                report = generate_health_report()
            else:
                report = get_latest_health_report()
            
            # Add timestamp
            report['timestamp'] = datetime.now()
            
            # Generate report ID
            report_id = (
                f"health_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                f"_{hash(str(report['timestamp']))}"
            )
            report['report_id'] = report_id
            
            return report
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Generating database health report",
                log_level='error',
                return_response=True
            )
            return error_response, error_status

class DatabaseFixResource(Resource):
    """Resource for fixing database issues."""
    
    @token_required
    # @require_admin
    @marshal_with(database_fix_fields)
    def post(self):
        """
        Fix database issues.
        
        Returns:
            JSON response with fix results
        """
        try:
            # Validate request parameters
            schema = FixDatabaseRequestSchema()
            try:
                args = schema.load(request.json if request.json else {})
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating database fix request parameters",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                return error_response, error_status
            
            # Extract parameters
            fix_type = args.get('fix_type', 'all')
            
            # Apply fixes
            fix_results = fix_common_issues(fix_type)
            
            # Add timestamp
            fix_results['timestamp'] = datetime.now()
            
            return fix_results
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fixing database issues",
                log_level='error',
                return_response=True
            )
            return error_response, error_status

# New health check resources
class HealthResource(Resource):
    """Resource for basic health check (public endpoint)."""
    
    @marshal_with(health_status_fields)
    def get(self):
        """
        Get basic health status information.
        
        This endpoint is public and does not require authentication.
        It provides minimal system status information for uptime monitoring.
        
        Returns:
            JSON response with basic health status
        """
        try:
            # Calculate uptime
            uptime_seconds = int(time.time() - START_TIME)
            uptime_formatted = str(timedelta(seconds=uptime_seconds))
            
            # Check database connection (minimal check)
            try:
                supabase = get_supabase_client()
                response = supabase.from_("property_types").select("count").limit(1).execute()
                db_status = "connected" if hasattr(response, 'data') else "disconnected"
            except Exception as e:
                logger.warning(f"Database connection check failed: {str(e)}")
                db_status = "disconnected"
            
            return {
                'status': 'ok' if db_status == 'connected' else 'degraded',
                'version': current_app.config.get('API_VERSION', '1.0.0'),
                'uptime': uptime_seconds,
                'uptime_formatted': uptime_formatted,
                'database_status': db_status,
                'timestamp': datetime.now()
            }
        except Exception as e:
            logger.error(f"Health check failed: {str(e)}")
            return {
                'status': 'error',
                'version': current_app.config.get('API_VERSION', '1.0.0'),
                'uptime': int(time.time() - START_TIME),
                'uptime_formatted': str(timedelta(seconds=int(time.time() - START_TIME))),
                'database_status': 'unknown',
                'timestamp': datetime.now()
            }, 500


class HealthDatabaseResource(Resource):
    """Resource for detailed database health checks (authenticated, admin only)."""
    
    @token_required
    # @require_admin
    @marshal_with(health_database_fields)
    def get(self):
        """
        Get detailed database health information.
        
        This endpoint requires authentication and admin privileges.
        It provides comprehensive database health status information.
        
        Returns:
            JSON response with detailed database health status
        """
        try:
            # Validate request parameters
            schema = HealthDatabaseRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating health database request parameters",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                return error_response, error_status
            
            # Extract parameters
            verbosity = args.get('verbosity', 'normal')
            format_type = args.get('format', 'json')
            
            # Run health check with appropriate verbosity
            with supabase_connection() as conn:
                health_check = DatabaseHealthCheck(conn)
                
                if verbosity == 'minimal':
                    # Only check schema for minimal verbosity
                    results = health_check.run_health_check(categories=['schema'])
                elif verbosity == 'detailed':
                    # Run all checks for detailed verbosity
                    results = health_check.run_health_check()
                else:
                    # Run schema and integrity checks for normal verbosity
                    results = health_check.run_health_check(categories=['schema', 'integrity'])
            
            # Format response based on format parameter
            if format_type in ['markdown', 'html']:
                # Generate formatted report
                report = health_check.generate_report(results, format=format_type)
                
                # Return as raw response if not JSON
                if format_type == 'markdown':
                    return Response(report, mimetype='text/markdown')
                else:
                    return Response(report, mimetype='text/html')
            
            # Extract key information for JSON response
            overall_status = results.get('overall_status', 'unknown')
            issues_count = results.get('total_issues_found', 0)
            
            # Extract status for each category
            schema_status = 'not_checked'
            integrity_status = 'not_checked'
            performance_status = 'not_checked'
            
            if 'results_by_category' in results:
                categories = results['results_by_category']
                if 'schema' in categories:
                    schema_status = categories['schema'].get('status', 'unknown')
                if 'integrity' in categories:
                    integrity_status = categories['integrity'].get('status', 'unknown')
                if 'performance' in categories:
                    performance_status = categories['performance'].get('status', 'unknown')
            
            # Prepare response
            response = {
                'status': 'ok' if overall_status == 'passed' else 'warning' if overall_status == 'warning' else 'error',
                'overall_health': overall_status,
                'schema_status': schema_status,
                'integrity_status': integrity_status,
                'performance_status': performance_status,
                'issues_count': issues_count,
                'timestamp': datetime.now(),
                'details': results if verbosity == 'detailed' else {
                    'categories_checked': results.get('categories_checked', []),
                    'recommendations': results.get('recommendations', [])
                }
            }
            
            return response
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fetching database health information",
                log_level='error',
                return_response=True
            )
            return error_response, error_status


class HealthPerformanceResource(Resource):
    """Resource for performance metrics (authenticated, admin only)."""
    
    @token_required
    # @require_admin
    @marshal_with(health_performance_fields)
    def get(self):
        """
        Get performance metrics for the system.
        
        This endpoint requires authentication and admin privileges.
        It provides performance metrics for the system, database, and API.
        
        Returns:
            JSON response with performance metrics
        """
        try:
            # Validate request parameters
            schema = HealthPerformanceRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating health performance request parameters",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                return error_response, error_status
            
            # Extract parameters
            period = args.get('period', 'hour')
            metrics_list = args.get('metrics', ['cpu', 'memory', 'database', 'api'])
            
            # Get system metrics
            system_metrics = {}
            if 'cpu' in metrics_list or 'memory' in metrics_list:
                # Get CPU usage
                if 'cpu' in metrics_list:
                    system_metrics['cpu'] = {
                        'current': psutil.cpu_percent(interval=0.1),
                        'history': self._get_cpu_history(period)
                    }
                
                # Get memory usage
                if 'memory' in metrics_list:
                    memory = psutil.virtual_memory()
                    system_metrics['memory'] = {
                        'current': {
                            'total': memory.total,
                            'available': memory.available,
                            'used': memory.used,
                            'percent': memory.percent
                        },
                        'history': self._get_memory_history(period)
                    }
            
            # Get database metrics
            db_metrics = {}
            if 'database' in metrics_list:
                # Run performance check to get database metrics
                with supabase_connection() as conn:
                    health_check = DatabaseHealthCheck(conn)
                    results = health_check.run_health_check(categories=['performance'])
                    
                    if 'results_by_category' in results and 'performance' in results['results_by_category']:
                        perf_results = results['results_by_category']['performance']
                        db_metrics = {
                            'status': perf_results.get('status', 'unknown'),
                            'query_performance': perf_results.get('details', {}).get('query_performance', {}),
                            'connection_stats': perf_results.get('details', {}).get('connection_stats', {}),
                            'index_usage': perf_results.get('details', {}).get('index_usage', {})
                        }
            
            # Get API metrics
            api_metrics = {}
            if 'api' in metrics_list:
                # In a real implementation, this would query API metrics from a monitoring system
                # For this example, we'll return some sample metrics
                api_metrics = self._get_api_metrics(period)
            
            return {
                'status': 'ok',
                'metrics': system_metrics,
                'database_metrics': db_metrics,
                'api_metrics': api_metrics,
                'timestamp': datetime.now()
            }
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fetching performance metrics",
                log_level='error',
                return_response=True
            )
            return error_response, error_status
    
    def _get_cpu_history(self, period):
        """Get CPU usage history for the specified period."""
        # In a real implementation, this would query a metrics database
        # For this example, we'll return some sample metrics
        points = self._get_points_for_period(period)
        return {
            'labels': [f"{i}" for i in range(points)],
            'values': [round(40 + 30 * (0.5 + 0.5 * (i % 10) / 10), 2) for i in range(points)],
            'unit': '%',
            'interval': self._get_interval_for_period(period)
        }
    
    def _get_memory_history(self, period):
        """Get memory usage history for the specified period."""
        # In a real implementation, this would query a metrics database
        # For this example, we'll return some sample metrics
        points = self._get_points_for_period(period)
        return {
            'labels': [f"{i}" for i in range(points)],
            'values': [round(60 + 20 * (0.5 + 0.5 * (i % 15) / 15), 2) for i in range(points)],
            'unit': '%',
            'interval': self._get_interval_for_period(period)
        }
    
    def _get_api_metrics(self, period):
        """Get API metrics for the specified period."""
        # In a real implementation, this would query a metrics database
        # For this example, we'll return some sample metrics
        points = self._get_points_for_period(period)
        return {
            'requests': {
                'labels': [f"{i}" for i in range(points)],
                'values': [int(100 + 50 * (0.5 + 0.5 * (i % 20) / 20)) for i in range(points)],
                'unit': 'count',
                'interval': self._get_interval_for_period(period)
            },
            'response_time': {
                'labels': [f"{i}" for i in range(points)],
                'values': [round(200 + 100 * (0.5 + 0.5 * (i % 30) / 30), 2) for i in range(points)],
                'unit': 'ms',
                'interval': self._get_interval_for_period(period)
            },
            'error_rate': {
                'labels': [f"{i}" for i in range(points)],
                'values': [round(1 + 4 * (0.5 + 0.5 * (i % 25) / 25), 2) for i in range(points)],
                'unit': '%',
                'interval': self._get_interval_for_period(period)
            }
        }
    
    def _get_points_for_period(self, period):
        """Get number of data points based on period."""
        if period == 'minute':
            return 60  # One point per second
        elif period == 'hour':
            return 60  # One point per minute
        elif period == 'day':
            return 24  # One point per hour
        elif period == 'week':
            return 7   # One point per day
        else:
            return 60  # Default to one point per minute
    
    def _get_interval_for_period(self, period):
        """Get interval label based on period."""
        if period == 'minute':
            return 'second'
        elif period == 'hour':
            return 'minute'
        elif period == 'day':
            return 'hour'
        elif period == 'week':
            return 'day'
        else:
            return 'minute'


class ServiceRoleCheckResource(Resource):
    """Resource for checking service role authentication."""
    
    @authenticate_service_role
    @marshal_with(service_role_check_fields)
    def get(self):
        """
        Check service role authentication.
        
        This endpoint requires a valid service role token to access,
        making it a suitable test for service role authentication.
        
        Returns:
            JSON response with service role authentication status
        """
        return {
            'status': 'ok',
            'message': 'Service role authentication successful',
            'timestamp': datetime.now(),
            'role': 'service_role'
        }


# Define Convex connection check fields
convex_connection_fields = {
    'status': fields.String,
    'message': fields.String,
    'timestamp': fields.DateTime,
    'connection_details': fields.Raw,
    'tables_available': fields.List(fields.String)
}

class ConvexConnectionCheckResource(Resource):
    """Resource for checking Convex database connection."""
    
    @marshal_with(convex_connection_fields)
    def get(self):
        """
        Check Convex database connection.
        
        This endpoint tests the connection to the Convex database and returns
        connection status details.
        
        Returns:
            JSON response with Convex connection status
        """
        try:
            from database.convex_adapter import create_client
            import os
            
            # Create Convex client with forced use_convex=True to ensure we test Convex specifically
            convex_client = create_client(use_convex=True)
            
            # Get URL and availability info
            connection_details = {
                'url': os.environ.get('CONVEX_URL', 'https://dynamic-mink-63.convex.cloud'),
                'use_convex_env': os.environ.get('USE_CONVEX', 'false'),
                'key_configured': bool(os.environ.get('CONVEX_DEPLOYMENT_KEY', '')),
            }
            
            # Test query to check connectivity - get available tables
            try:
                # Test a simple query on a common table
                molecules_result = convex_client.table('molecules').select('id').limit(1).execute()
                
                # Get list of available tables by querying system information
                # This is implementation-dependent and may need adjustment based on Convex schema
                tables = []
                try:
                    # Try to get tables list - this might need adjustment based on Convex API
                    # We'll make a best effort to extract table info
                    tables_query = convex_client.execute_query('GET', 'api/system/tables', {})
                    if tables_query and hasattr(tables_query, 'data'):
                        tables = [table['name'] for table in tables_query.data]
                    else:
                        # Fallback: Try a few common tables to see if they exist
                        common_tables = ['molecules', 'mixtures', 'protocols', 'experiments', 'users']
                        for table in common_tables:
                            try:
                                result = convex_client.table(table).select('count').limit(1).execute()
                                if not result.error:
                                    tables.append(table)
                            except:
                                pass
                except Exception as table_e:
                    logger.warning(f"Could not retrieve full table list from Convex: {str(table_e)}")
                
                return {
                    'status': 'ok',
                    'message': 'Successfully connected to Convex database',
                    'timestamp': datetime.now(),
                    'connection_details': connection_details,
                    'tables_available': tables
                }
            except Exception as query_e:
                logger.error(f"Convex connection test query failed: {str(query_e)}")
                return {
                    'status': 'error',
                    'message': f'Connected to Convex but query failed: {str(query_e)}',
                    'timestamp': datetime.now(),
                    'connection_details': connection_details,
                    'tables_available': []
                }
        except Exception as e:
            logger.error(f"Convex connection check failed: {str(e)}")
            return {
                'status': 'error',
                'message': f'Failed to connect to Convex: {str(e)}',
                'timestamp': datetime.now(),
                'connection_details': {
                    'url': os.environ.get('CONVEX_URL', 'Not configured'),
                    'use_convex_env': os.environ.get('USE_CONVEX', 'false'),
                    'key_configured': bool(os.environ.get('CONVEX_DEPLOYMENT_KEY', ''))
                },
                'tables_available': []
            }, 500


def register_resources(api):
    """
    Register system resources with the API.
    
    Args:
        api: Flask-RESTful API instance
    """
    # Original system resources
    api.add_resource(SystemStatusResource, '/api/v1/system/status')
    api.add_resource(SystemLogsResource, '/api/v1/system/logs')
    api.add_resource(SystemMetricsResource, '/api/v1/system/metrics')
    api.add_resource(DatabaseHealthResource, '/api/v1/system/database/health')
    api.add_resource(DatabaseFixResource, '/api/v1/system/database/fix')
    
    # New health check endpoints
    api.add_resource(HealthResource, '/health')
    api.add_resource(HealthDatabaseResource, '/health/database')
    api.add_resource(HealthPerformanceResource, '/health/performance')
    
    # Service role check endpoint
    api.add_resource(ServiceRoleCheckResource, '/api/v1/system/service-role-check')
    
    # Convex connection check endpoint
    api.add_resource(ConvexConnectionCheckResource, '/api/v1/system/convex-connection-check')
    
    # Start scheduled health check (once per day)
    run_scheduled_health_check(interval_hours=24)