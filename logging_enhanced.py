#!/usr/bin/env python3
"""
Enhanced Logging System for CryoProtect v2

This module provides a production-ready logging system with:
- JSON-structured logging format
- Correlation IDs for request tracking
- Contextual information
- Log rotation and retention policies
- ELK stack integration (Elasticsearch, Logstash, Kibana)
"""

import os
import sys
import json
import uuid
import logging
import logging.handlers
import datetime
import socket
import traceback
from typing import Dict, Any, Optional, Union
from flask import Flask
from pythonjsonlogger import jsonlogger
import ecs_logging


# Default configuration values
DEFAULT_CONFIG = {
    'LOG_LEVEL': 'INFO',
    'LOG_FORMAT': '%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    'LOG_FILE': 'logs/cryoprotect.log',
    'LOG_TO_FILE': '1',
    'LOG_TO_CONSOLE': '1',
    'LOG_TO_ELK': '0',
    'LOG_ROTATION_SIZE': '10485760',  # 10MB
    'LOG_BACKUP_COUNT': '10',
    'LOG_JSON_FORMAT': '1',
    'ELASTICSEARCH_HOST': 'elasticsearch:9200',
    'ELASTICSEARCH_INDEX': 'cryoprotect-logs',
    'LOG_CORRELATION_ID_HEADER': 'X-Correlation-ID',
    'LOG_RETENTION_DAYS': '30',
    'LOG_INCLUDE_CONTEXT': '1',
}


class ContextualFilter(logging.Filter):
    """
    Filter that adds contextual information to log records.
    """
    def filter(self, record):
        from flask import has_request_context, has_app_context
        
        # Add correlation ID if available
        if has_app_context():
            from flask import g
            record.correlation_id = getattr(g, 'correlation_id', 'N/A')
        else:
            record.correlation_id = 'N/A'
        
        # Add request information if available
        if has_request_context():
            from flask import g, request
            record.request_id = getattr(g, 'request_id', 'N/A') if has_app_context() else 'N/A'
            record.remote_addr = request.remote_addr
            record.method = request.method
            record.path = request.path
            record.user_agent = request.headers.get('User-Agent', 'N/A')
            record.user_id = getattr(g, 'user_id', 'N/A') if has_app_context() else 'N/A'
        else:
            record.request_id = 'N/A'
            record.remote_addr = 'N/A'
            record.method = 'N/A'
            record.path = 'N/A'
            record.user_agent = 'N/A'
            record.user_id = 'N/A'
            record.method = 'N/A'
            record.path = 'N/A'
            record.user_agent = 'N/A'
            record.user_id = 'N/A'
        
        # Add system information
        record.hostname = socket.gethostname()
        record.process_id = os.getpid()
        
        return True


class CustomJsonFormatter(jsonlogger.JsonFormatter):
    """
    Custom JSON formatter that adds additional fields.
    """
    def add_fields(self, log_record, record, message_dict):
        super(CustomJsonFormatter, self).add_fields(log_record, record, message_dict)
        
        # Add timestamp in ISO format
        log_record['timestamp'] = datetime.datetime.utcnow().isoformat()
        log_record['level'] = record.levelname
        log_record['logger'] = record.name
        
        # Add correlation ID and request information
        log_record['correlation_id'] = getattr(record, 'correlation_id', 'N/A')
        log_record['request_id'] = getattr(record, 'request_id', 'N/A')
        log_record['remote_addr'] = getattr(record, 'remote_addr', 'N/A')
        log_record['method'] = getattr(record, 'method', 'N/A')
        log_record['path'] = getattr(record, 'path', 'N/A')
        log_record['user_agent'] = getattr(record, 'user_agent', 'N/A')
        log_record['user_id'] = getattr(record, 'user_id', 'N/A')
        
        # Add system information
        log_record['hostname'] = getattr(record, 'hostname', socket.gethostname())
        log_record['process_id'] = getattr(record, 'process_id', os.getpid())
        
        # Add exception information if available
        if record.exc_info:
            log_record['exception'] = {
                'type': record.exc_info[0].__name__,
                'message': str(record.exc_info[1]),
                'traceback': traceback.format_exception(*record.exc_info)
            }


class ElasticsearchHandler(logging.Handler):
    """
    Custom handler for sending logs to Elasticsearch.
    This is a simplified implementation - in production, consider using
    a more robust solution like filebeat or logstash.
    """
    def __init__(self, host, index):
        super().__init__()
        self.host = host
        self.index = index
        # Lazy import to avoid dependency issues
        try:
            from elasticsearch import Elasticsearch
            self.es = Elasticsearch([host])
        except ImportError:
            self.es = None
            logging.getLogger(__name__).warning(
                "Elasticsearch package not installed. ELK logging disabled."
            )
    
    def emit(self, record):
        if not self.es:
            return
        
        try:
            log_entry = self.format(record)
            if isinstance(log_entry, str):
                log_entry = json.loads(log_entry)
            
            # Add timestamp for Elasticsearch
            if 'timestamp' not in log_entry:
                log_entry['@timestamp'] = datetime.datetime.utcnow().isoformat()
            
            # Send to Elasticsearch
            self.es.index(index=self.index, body=log_entry)
        except Exception:
            self.handleError(record)


def setup_enhanced_logging(app: Optional[Flask] = None) -> None:
    """
    Set up enhanced logging with JSON formatting, correlation IDs, and ELK integration.
    
    Args:
        app: Optional Flask application instance
    """
    # Create logs directory if it doesn't exist
    log_dir = os.path.dirname(os.getenv('LOG_FILE', DEFAULT_CONFIG['LOG_FILE']))
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    # Get configuration from environment variables or use defaults
    log_level = os.getenv('LOG_LEVEL', DEFAULT_CONFIG['LOG_LEVEL']).upper()
    log_file = os.getenv('LOG_FILE', DEFAULT_CONFIG['LOG_FILE'])
    log_to_file = os.getenv('LOG_TO_FILE', DEFAULT_CONFIG['LOG_TO_FILE']) == '1'
    log_to_console = os.getenv('LOG_TO_CONSOLE', DEFAULT_CONFIG['LOG_TO_CONSOLE']) == '1'
    log_to_elk = os.getenv('LOG_TO_ELK', DEFAULT_CONFIG['LOG_TO_ELK']) == '1'
    log_json_format = os.getenv('LOG_JSON_FORMAT', DEFAULT_CONFIG['LOG_JSON_FORMAT']) == '1'
    log_rotation_size = int(os.getenv('LOG_ROTATION_SIZE', DEFAULT_CONFIG['LOG_ROTATION_SIZE']))
    log_backup_count = int(os.getenv('LOG_BACKUP_COUNT', DEFAULT_CONFIG['LOG_BACKUP_COUNT']))
    elasticsearch_host = os.getenv('ELASTICSEARCH_HOST', DEFAULT_CONFIG['ELASTICSEARCH_HOST'])
    elasticsearch_index = os.getenv('ELASTICSEARCH_INDEX', DEFAULT_CONFIG['ELASTICSEARCH_INDEX'])
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, log_level, logging.INFO))
    
    # Clear existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Create handlers
    handlers = []
    
    # Console handler
    if log_to_console:
        console_handler = logging.StreamHandler(sys.stdout)
        handlers.append(console_handler)
    
    # File handler with rotation
    if log_to_file:
        file_handler = logging.handlers.RotatingFileHandler(
            log_file,
            maxBytes=log_rotation_size,
            backupCount=log_backup_count,
            encoding='utf-8'
        )
        handlers.append(file_handler)
    
    # Elasticsearch handler
    if log_to_elk:
        elk_handler = ElasticsearchHandler(elasticsearch_host, elasticsearch_index)
        handlers.append(elk_handler)
    
    # Configure formatters
    if log_json_format:
        # Use ECS-compliant logging if sending to ELK stack
        if log_to_elk:
            formatter = ecs_logging.StdlibFormatter()
        else:
            formatter = CustomJsonFormatter('%(timestamp)s %(level)s %(name)s %(message)s')
    else:
        formatter = logging.Formatter(
            os.getenv('LOG_FORMAT', DEFAULT_CONFIG['LOG_FORMAT'])
        )
    
    # Add contextual filter and formatter to all handlers
    contextual_filter = ContextualFilter()
    
    for handler in handlers:
        handler.setFormatter(formatter)
        handler.addFilter(contextual_filter)
        root_logger.addHandler(handler)
    
    # Configure Flask app if provided
    if app:
        # Set up request logging middleware
        setup_request_logging_middleware(app)


def setup_request_logging_middleware(app: Flask) -> None:
    """
    Set up middleware for request/response logging.
    
    Args:
        app: Flask application instance
    """
    @app.before_request
    def before_request():
        # Import Flask objects
        from flask import g, request
        
        # Generate request ID and correlation ID
        g.request_id = str(uuid.uuid4())
        
        # Use existing correlation ID from header or generate a new one
        correlation_id_header = os.getenv(
            'LOG_CORRELATION_ID_HEADER',
            DEFAULT_CONFIG['LOG_CORRELATION_ID_HEADER']
        )
        g.correlation_id = request.headers.get(correlation_id_header) or str(uuid.uuid4())
        
        # Store request start time for duration calculation
        g.start_time = datetime.datetime.utcnow()
        
        # Extract user ID if authenticated
        g.user_id = getattr(g, 'user_id', 'anonymous')
        
        # Log request
        app.logger.info(
            f"Request started: {request.method} {request.path}",
            extra={
                'event_type': 'request_started',
                'request': {
                    'method': request.method,
                    'path': request.path,
                    'remote_addr': request.remote_addr,
                    'user_agent': request.headers.get('User-Agent', 'N/A'),
                    'content_length': request.content_length,
                    'content_type': request.content_type,
                    'query_string': request.query_string.decode('utf-8', errors='replace'),
                    'headers': dict(request.headers)
                }
            }
        )
    
    @app.after_request
    def after_request(response):
        # Import Flask objects
        from flask import g, request
        
        # Calculate request duration
        duration = None
        if hasattr(g, 'start_time'):
            duration = (datetime.datetime.utcnow() - g.start_time).total_seconds()
        
        # Log response
        app.logger.info(
            f"Request completed: {request.method} {request.path} {response.status_code}",
            extra={
                'event_type': 'request_completed',
                'request': {
                    'method': request.method,
                    'path': request.path,
                    'remote_addr': request.remote_addr
                },
                'response': {
                    'status_code': response.status_code,
                    'content_length': response.content_length,
                    'content_type': response.content_type,
                    'duration': duration
                }
            }
        )
        
        # Add correlation ID to response headers
        correlation_id_header = os.getenv(
            'LOG_CORRELATION_ID_HEADER', 
            DEFAULT_CONFIG['LOG_CORRELATION_ID_HEADER']
        )
        if hasattr(g, 'correlation_id'):
            response.headers[correlation_id_header] = g.correlation_id
        
        return response
    
    @app.teardown_request
    def teardown_request(exception=None):
        if exception:
            # Import Flask objects
            from flask import request
            
            app.logger.error(
                f"Request failed: {request.method} {request.path}",
                exc_info=exception,
                extra={
                    'event_type': 'request_failed',
                    'request': {
                        'method': request.method,
                        'path': request.path,
                        'remote_addr': request.remote_addr
                    },
                    'error': {
                        'type': type(exception).__name__,
                        'message': str(exception)
                    }
                }
            )


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the given name.
    
    Args:
        name: Logger name
    
    Returns:
        logging.Logger: Configured logger
    """
    return logging.getLogger(name)


def log_with_context(
    logger: logging.Logger,
    level: str,
    message: str,
    context: Optional[Dict[str, Any]] = None,
    exc_info: Optional[Union[bool, Exception]] = None
) -> None:
    """
    Log a message with additional context.
    
    Args:
        logger: Logger instance
        level: Log level (debug, info, warning, error, critical)
        message: Log message
        context: Additional context to include in the log
        exc_info: Exception information
    """
    if context is None:
        context = {}
    
    log_method = getattr(logger, level.lower())
    log_method(message, extra=context, exc_info=exc_info)


if __name__ == "__main__":
    # Example usage
    setup_enhanced_logging()
    logger = get_logger(__name__)
    logger.info("Enhanced logging system initialized")