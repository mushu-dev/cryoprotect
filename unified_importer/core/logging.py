"""
Enhanced logging system for the unified molecular importer.

This module provides structured logging with JSON format support, 
rotating file handlers, and configurable log levels.
"""

import os
import json
import logging
import logging.handlers
import datetime
from typing import Optional, Dict, Any, Union


class StructuredLogRecord(logging.LogRecord):
    """Extended LogRecord with structured data support."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.structured_data = {}


class StructuredLogger(logging.Logger):
    """Logger that supports structured data in log messages."""
    
    def _log(self, level, msg, args, exc_info=None, extra=None, stack_info=False, 
             structured_data=None, **kwargs):
        """
        Enhanced logging with structured data support.
        
        Args:
            structured_data: Dictionary of structured data to include in the log record
        """
        if structured_data is None:
            structured_data = {}
            
        extra_dict = {} if extra is None else extra
        
        if 'structured_data' not in extra_dict:
            extra_dict['structured_data'] = structured_data
        else:
            extra_dict['structured_data'].update(structured_data)
            
        super()._log(level, msg, args, exc_info, extra_dict, stack_info, **kwargs)
    
    def info(self, msg, *args, structured_data=None, **kwargs):
        """Enhanced info logging with structured data."""
        self._log(logging.INFO, msg, args, structured_data=structured_data, **kwargs)
    
    def error(self, msg, *args, structured_data=None, **kwargs):
        """Enhanced error logging with structured data."""
        self._log(logging.ERROR, msg, args, structured_data=structured_data, **kwargs)
    
    def warning(self, msg, *args, structured_data=None, **kwargs):
        """Enhanced warning logging with structured data."""
        self._log(logging.WARNING, msg, args, structured_data=structured_data, **kwargs)
    
    def debug(self, msg, *args, structured_data=None, **kwargs):
        """Enhanced debug logging with structured data."""
        self._log(logging.DEBUG, msg, args, structured_data=structured_data, **kwargs)
    
    def exception(self, msg, *args, structured_data=None, **kwargs):
        """Enhanced exception logging with structured data."""
        exc_info = kwargs.pop('exc_info', True)
        self._log(logging.ERROR, msg, args, exc_info=exc_info, structured_data=structured_data, **kwargs)


class JsonFormatter(logging.Formatter):
    """Formatter that outputs JSON formatted logs."""
    
    def __init__(self, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt)
    
    def format(self, record: Union[logging.LogRecord, StructuredLogRecord]) -> str:
        """Format the log record as JSON."""
        log_data = {
            'timestamp': datetime.datetime.fromtimestamp(record.created).isoformat(),
            'level': record.levelname,
            'logger': record.name,
            'message': record.getMessage(),
            'module': record.module,
            'function': record.funcName,
            'line': record.lineno,
            'thread': record.thread,
            'thread_name': record.threadName,
        }
        
        # Add exception info if available
        if record.exc_info:
            log_data['exception'] = {
                'type': record.exc_info[0].__name__,
                'message': str(record.exc_info[1]),
                'traceback': self.formatException(record.exc_info)
            }
        
        # Add structured data if available
        if hasattr(record, 'structured_data') and record.structured_data:
            log_data['data'] = record.structured_data
            
        return json.dumps(log_data)


def setup_logging(
    logger_name: str,
    log_level: str = "INFO",
    log_file: Optional[str] = None,
    json_format: bool = True,
    max_file_size: int = 10 * 1024 * 1024,  # 10 MB
    backup_count: int = 5
) -> StructuredLogger:
    """
    Setup a logger with console and optional file handlers.
    
    Args:
        logger_name: Name of the logger
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Path to log file (optional)
        json_format: Whether to use JSON formatting
        max_file_size: Maximum size of log file before rotation
        backup_count: Number of backup log files to keep
        
    Returns:
        StructuredLogger instance
    """
    # Register our custom logger class
    logging.setLoggerClass(StructuredLogger)
    
    # Create logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Clear existing handlers
    logger.handlers = []
    
    # Create formatters
    if json_format:
        formatter = JsonFormatter()
    else:
        formatter = logging.Formatter(
            '%(asctime)s [%(levelname)s] %(name)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (if specified)
    if log_file:
        # Ensure log directory exists
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)
            
        file_handler = logging.handlers.RotatingFileHandler(
            log_file,
            maxBytes=max_file_size,
            backupCount=backup_count
        )
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger  # type: ignore


def get_source_logger(source_name: str, config: Dict[str, Any]) -> StructuredLogger:
    """
    Create a logger for a specific data source.
    
    Args:
        source_name: Name of the data source (e.g., 'chembl', 'pubchem')
        config: Configuration dictionary
        
    Returns:
        StructuredLogger instance
    """
    log_level = config.get('log_level', 'INFO')
    log_dir = config.get('log_dir', 'logs')
    json_format = config.get('json_logs', True)
    
    # Ensure log directory exists
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
        
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f"{source_name}_import_{timestamp}.log")
    
    return setup_logging(
        f"cryoprotect.importer.{source_name}",
        log_level=log_level,
        log_file=log_file,
        json_format=json_format
    )