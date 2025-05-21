"""
Comprehensive error handling system for the unified molecular importer.

This module provides advanced error handling capabilities including error
classification, recovery strategies, retry mechanisms, and detailed logging
with context information.
"""

import os
import time
import enum
import json
import logging
import inspect
import traceback
import threading
import asyncio
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Union, Set, Tuple, Callable, Type, TypeVar
from datetime import datetime

# Type variables
T = TypeVar('T')  # Return type for retry functions


class ErrorSeverity(enum.Enum):
    """Severity level of errors."""
    DEBUG = 0
    INFO = 1
    WARNING = 2
    ERROR = 3
    CRITICAL = 4


class ErrorCategory(enum.Enum):
    """Categories of errors for classification."""
    # System errors
    SYSTEM = "system"              # OS, hardware, environment
    NETWORK = "network"            # Connectivity, DNS, TCP/IP
    DATABASE = "database"          # Database connection, query execution
    RESOURCE = "resource"          # Memory, disk, CPU limitations
    
    # API errors
    API = "api"                    # External API access
    API_RATE_LIMIT = "api_rate_limit"  # Rate limiting from external APIs
    API_TIMEOUT = "api_timeout"    # API request timeout
    API_UNAVAILABLE = "api_unavailable"  # API service unavailable
    
    # Data errors
    DATA_FORMAT = "data_format"    # Invalid data format
    DATA_VALIDATION = "data_validation"  # Data validation failures
    DATA_MISSING = "data_missing"  # Required data not found
    DATA_DUPLICATE = "data_duplicate"  # Duplicate data
    
    # Application errors
    CONFIGURATION = "configuration"  # Configuration errors
    PERMISSION = "permission"      # Permission, authorization issues
    AUTHENTICATION = "authentication"  # Authentication issues
    BUSINESS_LOGIC = "business_logic"  # Business logic errors
    
    # Import-specific errors
    MOLECULE_PARSING = "molecule_parsing"  # Molecule parsing failures
    PROPERTY_CALCULATION = "property_calculation"  # Property calculation failures
    TRANSFORMATION = "transformation"  # Transformation failures
    IMPORT_CONSISTENCY = "import_consistency"  # Import consistency issues
    
    # Unknown/uncategorized
    UNKNOWN = "unknown"            # Uncategorized errors


class RecoveryStrategy(enum.Enum):
    """Strategies for recovering from errors."""
    RETRY = "retry"                # Retry the operation
    SKIP = "skip"                  # Skip the failed item and continue
    FALLBACK = "fallback"          # Use fallback/default values
    DELAY = "delay"                # Delay and retry
    ALTERNATE = "alternate"        # Try an alternate approach
    NONE = "none"                  # No recovery possible
    ABORT = "abort"                # Abort the operation


@dataclass
class ErrorContext:
    """
    Context information for an error.
    
    Stores detailed information about the context in which an error occurred,
    including source location, input data, and environment details.
    """
    # Basic information
    timestamp: float = field(default_factory=time.time)
    function_name: str = ""
    file_name: str = ""
    line_number: int = 0
    
    # Input context
    input_data: Optional[Dict[str, Any]] = None
    operation: str = ""
    
    # Environment context
    environment: Dict[str, Any] = field(default_factory=dict)
    
    # Additional context
    additional_context: Dict[str, Any] = field(default_factory=dict)
    
    @classmethod
    def current(cls, input_data: Optional[Dict[str, Any]] = None, operation: str = "", additional_context: Optional[Dict[str, Any]] = None) -> 'ErrorContext':
        """
        Create an ErrorContext from the current execution context.
        
        Args:
            input_data: Input data that led to the error
            operation: Operation being performed
            additional_context: Additional context information
            
        Returns:
            ErrorContext instance
        """
        # Get caller frame
        frame = inspect.currentframe()
        if frame is not None:
            frame = frame.f_back  # Get the caller's frame
        
        if frame is None:
            return cls(
                input_data=input_data,
                operation=operation,
                additional_context=additional_context or {}
            )
        
        # Extract frame information
        file_name = frame.f_code.co_filename
        function_name = frame.f_code.co_name
        line_number = frame.f_lineno
        
        # Create environment context
        environment = {
            "python_version": os.environ.get("PYTHON_VERSION", ""),
            "os_name": os.name,
            "pid": os.getpid(),
            "thread_id": threading.get_ident(),
        }
        
        # Try to get asyncio task name if in asyncio context
        try:
            task = asyncio.current_task()
            if task:
                environment["asyncio_task"] = task.get_name()
        except RuntimeError:
            # Not in asyncio context
            pass
        
        return cls(
            function_name=function_name,
            file_name=file_name,
            line_number=line_number,
            input_data=input_data,
            operation=operation,
            environment=environment,
            additional_context=additional_context or {}
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the error context to a dictionary."""
        return {
            "timestamp": self.timestamp,
            "function_name": self.function_name,
            "file_name": self.file_name,
            "line_number": self.line_number,
            "input_data": self.input_data,
            "operation": self.operation,
            "environment": self.environment,
            "additional_context": self.additional_context
        }
    
    def __str__(self) -> str:
        """String representation of the error context."""
        parts = [
            f"Time: {datetime.fromtimestamp(self.timestamp).isoformat()}",
            f"Location: {self.file_name}:{self.line_number} in {self.function_name}()"
        ]
        
        if self.operation:
            parts.append(f"Operation: {self.operation}")
        
        if self.input_data:
            # Limit input data size for readability
            input_str = str(self.input_data)
            if len(input_str) > 500:
                input_str = input_str[:500] + "..."
            parts.append(f"Input: {input_str}")
        
        # Add environment info
        if self.environment:
            env_parts = []
            for key, value in self.environment.items():
                env_parts.append(f"{key}={value}")
            parts.append(f"Environment: {', '.join(env_parts)}")
        
        # Add additional context
        if self.additional_context:
            context_parts = []
            for key, value in self.additional_context.items():
                value_str = str(value)
                if len(value_str) > 100:
                    value_str = value_str[:100] + "..."
                context_parts.append(f"{key}={value_str}")
            parts.append(f"Context: {', '.join(context_parts)}")
        
        return "\n".join(parts)


@dataclass
class ClassifiedError:
    """
    Classified error with context and recovery strategy.
    
    This class enriches exceptions with classification, context,
    and recovery strategy information.
    """
    # Basic error information
    message: str
    exception: Optional[Exception] = None
    traceback_str: str = ""
    error_id: str = field(default_factory=lambda: f"error-{int(time.time() * 1000)}")
    
    # Classification
    category: ErrorCategory = ErrorCategory.UNKNOWN
    severity: ErrorSeverity = ErrorSeverity.ERROR
    recoverable: bool = False
    
    # Context
    context: Optional[ErrorContext] = None
    
    # Recovery
    strategy: RecoveryStrategy = RecoveryStrategy.NONE
    retry_count: int = 0
    max_retries: int = 0
    retry_delay: float = 0.0
    
    def __post_init__(self):
        """Initialize derived fields after initialization."""
        # Capture traceback if not provided
        if not self.traceback_str and self.exception is not None:
            self.traceback_str = ''.join(traceback.format_exception(
                type(self.exception),
                self.exception,
                self.exception.__traceback__
            ))
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the classified error to a dictionary."""
        result = {
            "error_id": self.error_id,
            "message": self.message,
            "category": self.category.value,
            "severity": self.severity.name,
            "recoverable": self.recoverable,
            "strategy": self.strategy.value,
            "retry_count": self.retry_count,
            "max_retries": self.max_retries,
            "retry_delay": self.retry_delay,
        }
        
        if self.exception:
            result["exception_type"] = type(self.exception).__name__
        
        if self.context:
            result["context"] = self.context.to_dict()
        
        return result
    
    def to_json(self) -> str:
        """Convert the classified error to a JSON string."""
        return json.dumps(self.to_dict(), default=str, indent=2)
    
    def log(self, logger: logging.Logger) -> None:
        """
        Log the error with the appropriate severity level.
        
        Args:
            logger: Logger to use for logging
        """
        # Map severity to log level
        log_level = {
            ErrorSeverity.DEBUG: logging.DEBUG,
            ErrorSeverity.INFO: logging.INFO,
            ErrorSeverity.WARNING: logging.WARNING,
            ErrorSeverity.ERROR: logging.ERROR,
            ErrorSeverity.CRITICAL: logging.CRITICAL
        }.get(self.severity, logging.ERROR)
        
        # Create log message
        log_message = f"{self.category.value.upper()}: {self.message}"
        
        # Add context information
        if self.context:
            log_message += f"\nContext: {str(self.context)}"
        
        # Add recovery information for recoverable errors
        if self.recoverable:
            recovery_info = f"Recovery strategy: {self.strategy.value}"
            if self.strategy == RecoveryStrategy.RETRY:
                recovery_info += f" (attempt {self.retry_count}/{self.max_retries})"
            log_message += f"\n{recovery_info}"
        
        # Log the error
        logger.log(log_level, log_message, exc_info=self.exception)
    
    def can_retry(self) -> bool:
        """
        Check if the error can be retried.
        
        Returns:
            True if the error can be retried, False otherwise
        """
        return (
            self.recoverable and
            self.strategy == RecoveryStrategy.RETRY and
            self.retry_count < self.max_retries
        )
    
    def increment_retry(self) -> 'ClassifiedError':
        """
        Increment the retry count.
        
        Returns:
            Updated ClassifiedError instance
        """
        self.retry_count += 1
        return self
    
    def with_strategy(self, strategy: RecoveryStrategy) -> 'ClassifiedError':
        """
        Update the recovery strategy.
        
        Args:
            strategy: New recovery strategy
            
        Returns:
            Updated ClassifiedError instance
        """
        self.strategy = strategy
        return self
    
    def is_transient(self) -> bool:
        """
        Check if the error is likely transient (temporary).
        
        Returns:
            True if the error is likely transient, False otherwise
        """
        # Categories that are likely to be transient
        transient_categories = {
            ErrorCategory.NETWORK,
            ErrorCategory.API_RATE_LIMIT,
            ErrorCategory.API_TIMEOUT,
            ErrorCategory.API_UNAVAILABLE,
            ErrorCategory.RESOURCE
        }
        
        return self.category in transient_categories


class ErrorClassifier:
    """
    Classifies exceptions into error categories with recovery strategies.
    
    The classifier analyzes exceptions and assigns them to categories,
    determines severity and recoverability, and suggests recovery strategies.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the error classifier.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
        # Default classification rules
        self.rules = self._default_rules()
        
        # Retry configuration
        self.default_max_retries = 3
        self.default_retry_delay = 2.0
        self.retry_backoff_factor = 2.0  # For exponential backoff
        
        # Statistics
        self.stats = {
            "classified_count": 0,
            "categories": {},
            "recoverable_count": 0,
            "unrecoverable_count": 0,
            "retry_success_count": 0,
            "retry_failure_count": 0
        }
    
    def _default_rules(self) -> Dict[Type[Exception], Tuple[ErrorCategory, ErrorSeverity, bool, RecoveryStrategy]]:
        """
        Define default classification rules for common exceptions.
        
        Returns:
            Dictionary mapping exception types to classification tuples
        """
        import socket
        import urllib.error
        import http.client
        import ssl
        import sqlite3
        import json
        
        # Each rule is: (category, severity, recoverable, strategy)
        return {
            # Network errors
            socket.error: (ErrorCategory.NETWORK, ErrorSeverity.ERROR, True, RecoveryStrategy.RETRY),
            socket.timeout: (ErrorCategory.NETWORK, ErrorSeverity.WARNING, True, RecoveryStrategy.RETRY),
            urllib.error.URLError: (ErrorCategory.NETWORK, ErrorSeverity.ERROR, True, RecoveryStrategy.RETRY),
            http.client.HTTPException: (ErrorCategory.API, ErrorSeverity.ERROR, True, RecoveryStrategy.RETRY),
            ssl.SSLError: (ErrorCategory.NETWORK, ErrorSeverity.ERROR, True, RecoveryStrategy.RETRY),
            ConnectionError: (ErrorCategory.NETWORK, ErrorSeverity.ERROR, True, RecoveryStrategy.RETRY),
            TimeoutError: (ErrorCategory.API_TIMEOUT, ErrorSeverity.WARNING, True, RecoveryStrategy.RETRY),
            
            # Database errors
            sqlite3.Error: (ErrorCategory.DATABASE, ErrorSeverity.ERROR, False, RecoveryStrategy.NONE),
            
            # Data errors
            ValueError: (ErrorCategory.DATA_VALIDATION, ErrorSeverity.WARNING, False, RecoveryStrategy.SKIP),
            TypeError: (ErrorCategory.DATA_FORMAT, ErrorSeverity.WARNING, False, RecoveryStrategy.SKIP),
            json.JSONDecodeError: (ErrorCategory.DATA_FORMAT, ErrorSeverity.WARNING, False, RecoveryStrategy.SKIP),
            
            # System errors
            MemoryError: (ErrorCategory.RESOURCE, ErrorSeverity.CRITICAL, False, RecoveryStrategy.ABORT),
            PermissionError: (ErrorCategory.PERMISSION, ErrorSeverity.ERROR, False, RecoveryStrategy.NONE),
            FileNotFoundError: (ErrorCategory.DATA_MISSING, ErrorSeverity.ERROR, False, RecoveryStrategy.NONE),
            
            # Default for unknown exceptions
            Exception: (ErrorCategory.UNKNOWN, ErrorSeverity.ERROR, False, RecoveryStrategy.NONE)
        }
    
    def add_rule(
        self,
        exception_type: Type[Exception],
        category: ErrorCategory,
        severity: ErrorSeverity,
        recoverable: bool,
        strategy: RecoveryStrategy
    ) -> None:
        """
        Add a classification rule for an exception type.
        
        Args:
            exception_type: Type of exception to classify
            category: Error category
            severity: Error severity
            recoverable: Whether the error is recoverable
            strategy: Recovery strategy
        """
        self.rules[exception_type] = (category, severity, recoverable, strategy)
    
    def _get_most_specific_rule(self, exception: Exception) -> Tuple[ErrorCategory, ErrorSeverity, bool, RecoveryStrategy]:
        """
        Get the most specific classification rule for an exception.
        
        This finds the rule for the most specific exception type in the
        exception's class hierarchy.
        
        Args:
            exception: Exception to classify
            
        Returns:
            Tuple of (category, severity, recoverable, strategy)
        """
        # Get exception's class and its base classes
        exception_class = type(exception)
        class_hierarchy = [exception_class]
        base_class = exception_class.__base__
        
        while base_class is not object:
            class_hierarchy.append(base_class)
            base_class = base_class.__base__
        
        # Find the most specific rule
        for exc_type in class_hierarchy:
            if exc_type in self.rules:
                return self.rules[exc_type]
        
        # If no specific rule found, use the default for Exception
        return self.rules.get(Exception, (
            ErrorCategory.UNKNOWN,
            ErrorSeverity.ERROR,
            False,
            RecoveryStrategy.NONE
        ))
    
    def _classify_by_message(self, exception: Exception, category: ErrorCategory) -> ErrorCategory:
        """
        Refine the classification based on exception message.
        
        This looks for keywords in the exception message to further
        refine the classification.
        
        Args:
            exception: Exception to classify
            category: Initial category classification
            
        Returns:
            Refined error category
        """
        message = str(exception).lower()
        
        # Rate limit keywords
        rate_limit_keywords = {"rate limit", "too many requests", "429", "throttl"}
        if category == ErrorCategory.API and any(kw in message for kw in rate_limit_keywords):
            return ErrorCategory.API_RATE_LIMIT
        
        # Timeout keywords
        timeout_keywords = {"timeout", "timed out", "deadline exceeded"}
        if category == ErrorCategory.API and any(kw in message for kw in timeout_keywords):
            return ErrorCategory.API_TIMEOUT
        
        # Service unavailable keywords
        unavailable_keywords = {"unavailable", "down", "maintenance", "503"}
        if category == ErrorCategory.API and any(kw in message for kw in unavailable_keywords):
            return ErrorCategory.API_UNAVAILABLE
        
        # Molecule parsing keywords
        molecule_keywords = {"smiles", "molecule", "atom", "bond", "structure"}
        if category == ErrorCategory.DATA_FORMAT and any(kw in message for kw in molecule_keywords):
            return ErrorCategory.MOLECULE_PARSING
        
        # Property calculation keywords
        property_keywords = {"property", "calculation", "descriptor", "fingerprint"}
        if category == ErrorCategory.DATA_FORMAT and any(kw in message for kw in property_keywords):
            return ErrorCategory.PROPERTY_CALCULATION
        
        return category
    
    def classify(
        self,
        error: Union[Exception, str],
        context: Optional[ErrorContext] = None,
        max_retries: Optional[int] = None,
        retry_delay: Optional[float] = None
    ) -> ClassifiedError:
        """
        Classify an error with context and recovery strategy.
        
        Args:
            error: Exception or error message to classify
            context: Error context
            max_retries: Maximum number of retries (None for default)
            retry_delay: Delay between retries in seconds (None for default)
            
        Returns:
            ClassifiedError instance
        """
        # Convert string errors to exceptions
        if isinstance(error, str):
            exception = Exception(error)
            message = error
        else:
            exception = error
            message = str(error)
        
        # Get the classification rule
        category, severity, recoverable, strategy = self._get_most_specific_rule(exception)
        
        # Refine the category based on the message
        category = self._classify_by_message(exception, category)
        
        # Create context if not provided
        if context is None:
            context = ErrorContext.current()
        
        # Create classified error
        classified_error = ClassifiedError(
            message=message,
            exception=exception,
            category=category,
            severity=severity,
            recoverable=recoverable,
            strategy=strategy,
            context=context,
            max_retries=max_retries or self.default_max_retries,
            retry_delay=retry_delay or self.default_retry_delay
        )
        
        # Update statistics
        self.stats["classified_count"] += 1
        
        if category.value not in self.stats["categories"]:
            self.stats["categories"][category.value] = 0
        self.stats["categories"][category.value] += 1
        
        if recoverable:
            self.stats["recoverable_count"] += 1
        else:
            self.stats["unrecoverable_count"] += 1
        
        return classified_error
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about error classification.
        
        Returns:
            Dictionary with classification statistics
        """
        return self.stats


class RetryManager:
    """
    Manages retrying operations with exponential backoff and jitter.
    
    This class provides retry capabilities for operations that may fail
    transiently, with customizable backoff and jitter.
    """
    
    def __init__(
        self,
        max_retries: int = 3,
        base_delay: float = 1.0,
        max_delay: float = 60.0,
        backoff_factor: float = 2.0,
        jitter: bool = True,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the retry manager.
        
        Args:
            max_retries: Maximum number of retry attempts
            base_delay: Base delay between retries in seconds
            max_delay: Maximum delay between retries in seconds
            backoff_factor: Exponential backoff factor
            jitter: Whether to add jitter to delay times
            logger: Logger instance
        """
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.backoff_factor = backoff_factor
        self.jitter = jitter
        self.logger = logger or logging.getLogger(__name__)
        
        # Statistics
        self.stats = {
            "attempts": 0,
            "successes": 0,
            "failures": 0,
            "retries": 0
        }
    
    def _calculate_delay(self, attempt: int) -> float:
        """
        Calculate the delay for a retry attempt with exponential backoff.
        
        Args:
            attempt: Current attempt number (0-based)
            
        Returns:
            Delay in seconds
        """
        # Calculate exponential backoff
        delay = self.base_delay * (self.backoff_factor ** attempt)
        
        # Cap at max delay
        delay = min(delay, self.max_delay)
        
        # Add jitter if enabled (±25%)
        if self.jitter:
            import random
            jitter_factor = random.uniform(0.75, 1.25)
            delay *= jitter_factor
        
        return delay
    
    async def retry_async(
        self,
        operation: Callable[..., T],
        *args,
        retry_condition: Optional[Callable[[Exception], bool]] = None,
        error_classifier: Optional[ErrorClassifier] = None,
        context_provider: Optional[Callable[[], ErrorContext]] = None,
        **kwargs
    ) -> T:
        """
        Retry an asynchronous operation with exponential backoff.
        
        Args:
            operation: Async function to retry
            *args: Arguments for the operation
            retry_condition: Function that returns True if the exception is retryable
            error_classifier: Error classifier for categorizing exceptions
            context_provider: Function that returns error context
            **kwargs: Keyword arguments for the operation
            
        Returns:
            Result of the operation
            
        Raises:
            Exception: The last exception if all retries fail
        """
        attempt = 0
        last_error = None
        classified_error = None
        
        # Initialize error classifier if needed
        if error_classifier is None and retry_condition is None:
            error_classifier = ErrorClassifier(logger=self.logger)
        
        while True:
            try:
                # Update statistics
                self.stats["attempts"] += 1
                
                # Attempt the operation
                if attempt > 0:
                    self.logger.info(f"Retry attempt {attempt}/{self.max_retries} for {operation.__name__}")
                
                result = await operation(*args, **kwargs)
                
                # Success
                self.stats["successes"] += 1
                
                # Update retry success statistics if this was a retry
                if attempt > 0 and error_classifier is not None:
                    error_classifier.stats["retry_success_count"] += 1
                
                return result
            except Exception as e:
                last_error = e
                
                # Determine if we should retry
                should_retry = False
                
                # If error_classifier is provided, use it to classify the error
                if error_classifier is not None:
                    # Get context if provided
                    context = None
                    if context_provider is not None:
                        context = context_provider()
                    
                    # Classify the error
                    classified_error = error_classifier.classify(
                        e,
                        context=context,
                        max_retries=self.max_retries
                    )
                    
                    # Log the classified error
                    classified_error.log(self.logger)
                    
                    # Check if retryable
                    should_retry = classified_error.can_retry()
                    
                    # Increment retry count if retrying
                    if should_retry:
                        classified_error.increment_retry()
                elif retry_condition is not None:
                    # Use provided retry condition
                    should_retry = retry_condition(e)
                
                # Check if we've reached the maximum retries
                if not should_retry or attempt >= self.max_retries:
                    # Failed after all retries
                    self.stats["failures"] += 1
                    
                    # Update retry failure statistics if this was a retry
                    if attempt > 0 and error_classifier is not None:
                        error_classifier.stats["retry_failure_count"] += 1
                    
                    # Re-raise the last error
                    raise
                
                # Increment attempt count
                attempt += 1
                self.stats["retries"] += 1
                
                # Calculate and apply delay
                delay = self._calculate_delay(attempt - 1)
                
                # Log retry attempt
                self.logger.warning(
                    f"Operation {operation.__name__} failed: {str(e)}. "
                    f"Retrying in {delay:.2f}s (attempt {attempt}/{self.max_retries})"
                )
                
                # Wait before retrying
                await asyncio.sleep(delay)
    
    def retry(
        self,
        operation: Callable[..., T],
        *args,
        retry_condition: Optional[Callable[[Exception], bool]] = None,
        error_classifier: Optional[ErrorClassifier] = None,
        context_provider: Optional[Callable[[], ErrorContext]] = None,
        **kwargs
    ) -> T:
        """
        Retry a synchronous operation with exponential backoff.
        
        Args:
            operation: Function to retry
            *args: Arguments for the operation
            retry_condition: Function that returns True if the exception is retryable
            error_classifier: Error classifier for categorizing exceptions
            context_provider: Function that returns error context
            **kwargs: Keyword arguments for the operation
            
        Returns:
            Result of the operation
            
        Raises:
            Exception: The last exception if all retries fail
        """
        attempt = 0
        last_error = None
        classified_error = None
        
        # Initialize error classifier if needed
        if error_classifier is None and retry_condition is None:
            error_classifier = ErrorClassifier(logger=self.logger)
        
        while True:
            try:
                # Update statistics
                self.stats["attempts"] += 1
                
                # Attempt the operation
                if attempt > 0:
                    self.logger.info(f"Retry attempt {attempt}/{self.max_retries} for {operation.__name__}")
                
                result = operation(*args, **kwargs)
                
                # Success
                self.stats["successes"] += 1
                
                # Update retry success statistics if this was a retry
                if attempt > 0 and error_classifier is not None:
                    error_classifier.stats["retry_success_count"] += 1
                
                return result
            except Exception as e:
                last_error = e
                
                # Determine if we should retry
                should_retry = False
                
                # If error_classifier is provided, use it to classify the error
                if error_classifier is not None:
                    # Get context if provided
                    context = None
                    if context_provider is not None:
                        context = context_provider()
                    
                    # Classify the error
                    classified_error = error_classifier.classify(
                        e,
                        context=context,
                        max_retries=self.max_retries
                    )
                    
                    # Log the classified error
                    classified_error.log(self.logger)
                    
                    # Check if retryable
                    should_retry = classified_error.can_retry()
                    
                    # Increment retry count if retrying
                    if should_retry:
                        classified_error.increment_retry()
                elif retry_condition is not None:
                    # Use provided retry condition
                    should_retry = retry_condition(e)
                
                # Check if we've reached the maximum retries
                if not should_retry or attempt >= self.max_retries:
                    # Failed after all retries
                    self.stats["failures"] += 1
                    
                    # Update retry failure statistics if this was a retry
                    if attempt > 0 and error_classifier is not None:
                        error_classifier.stats["retry_failure_count"] += 1
                    
                    # Re-raise the last error
                    raise
                
                # Increment attempt count
                attempt += 1
                self.stats["retries"] += 1
                
                # Calculate and apply delay
                delay = self._calculate_delay(attempt - 1)
                
                # Log retry attempt
                self.logger.warning(
                    f"Operation {operation.__name__} failed: {str(e)}. "
                    f"Retrying in {delay:.2f}s (attempt {attempt}/{self.max_retries})"
                )
                
                # Wait before retrying
                time.sleep(delay)
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about retry operations.
        
        Returns:
            Dictionary with retry statistics
        """
        stats = self.stats.copy()
        
        # Calculate success rate
        if stats["attempts"] > 0:
            stats["success_rate"] = stats["successes"] / stats["attempts"]
        else:
            stats["success_rate"] = 0.0
        
        # Calculate retry rate (percentage of attempts that were retries)
        if stats["attempts"] > 0:
            stats["retry_rate"] = stats["retries"] / stats["attempts"]
        else:
            stats["retry_rate"] = 0.0
        
        return stats


class ValidationErrorHandler:
    """
    Handles validation errors with customizable strategies.
    
    This class provides a way to handle validation errors with
    different strategies, such as skipping invalid items, using
    defaults, or applying corrections.
    """
    
    def __init__(
        self,
        logger: Optional[logging.Logger] = None,
        error_classifier: Optional[ErrorClassifier] = None
    ):
        """
        Initialize the validation error handler.
        
        Args:
            logger: Logger instance
            error_classifier: Error classifier for categorizing validation errors
        """
        self.logger = logger or logging.getLogger(__name__)
        self.error_classifier = error_classifier or ErrorClassifier(logger=self.logger)
        
        # Default handlers for validation errors
        self.handlers = {}
        
        # Statistics
        self.stats = {
            "validation_errors": 0,
            "handled_errors": 0,
            "strategies": {}
        }
    
    def register_handler(
        self,
        error_type: Union[Type[Exception], str],
        handler: Callable[[Any, Exception, Dict[str, Any]], Tuple[bool, Any]]
    ) -> None:
        """
        Register a handler for a specific validation error type.
        
        Args:
            error_type: Exception type or error message substring
            handler: Function that handles the error and returns (success, result)
        """
        self.handlers[error_type] = handler
    
    def _find_handler(self, error: Exception) -> Optional[Callable[[Any, Exception, Dict[str, Any]], Tuple[bool, Any]]]:
        """
        Find the most appropriate handler for an error.
        
        Args:
            error: Exception to handle
            
        Returns:
            Handler function or None if no handler found
        """
        # Check for exact exception type match
        error_type = type(error)
        if error_type in self.handlers:
            return self.handlers[error_type]
        
        # Check for exception type inheritance
        for handled_type, handler in self.handlers.items():
            if isinstance(handled_type, type) and isinstance(error, handled_type):
                return handler
        
        # Check for string message substring match
        error_message = str(error).lower()
        for handled_type, handler in self.handlers.items():
            if isinstance(handled_type, str) and handled_type.lower() in error_message:
                return handler
        
        # No handler found
        return None
    
    def handle_validation_error(
        self,
        data: Any,
        error: Exception,
        context: Optional[Dict[str, Any]] = None,
        default_strategy: RecoveryStrategy = RecoveryStrategy.SKIP
    ) -> Tuple[bool, Any, Optional[ClassifiedError]]:
        """
        Handle a validation error.
        
        Args:
            data: Data that failed validation
            error: Validation exception
            context: Additional context for error handling
            default_strategy: Default strategy if no handler is found
            
        Returns:
            Tuple of (success, result, classified_error)
        """
        # Update statistics
        self.stats["validation_errors"] += 1
        
        # Create error context
        error_context = ErrorContext.current(
            input_data={"data": data},
            operation="validation",
            additional_context=context
        )
        
        # Classify the error
        classified_error = self.error_classifier.classify(
            error,
            context=error_context
        )
        
        # Log the classified error
        classified_error.log(self.logger)
        
        # Find handler
        handler = self._find_handler(error)
        
        if handler is not None:
            # Apply handler
            try:
                success, result = handler(data, error, context or {})
                
                # Update statistics
                self.stats["handled_errors"] += 1
                
                strategy = "custom_handler"
                if strategy not in self.stats["strategies"]:
                    self.stats["strategies"][strategy] = 0
                self.stats["strategies"][strategy] += 1
                
                return success, result, classified_error
            except Exception as handler_error:
                # Handler failed
                self.logger.error(
                    f"Validation error handler failed: {str(handler_error)}",
                    exc_info=handler_error
                )
        
        # Apply default strategy
        classified_error.with_strategy(default_strategy)
        
        if default_strategy not in self.stats["strategies"]:
            self.stats["strategies"][default_strategy.value] = 0
        self.stats["strategies"][default_strategy.value] += 1
        
        if default_strategy == RecoveryStrategy.SKIP:
            return False, None, classified_error
        elif default_strategy == RecoveryStrategy.FALLBACK:
            # Try to create a minimal valid result
            if isinstance(data, dict):
                result = {k: data.get(k) for k in data if k in context.get("required_fields", [])}
                return True, result, classified_error
            return False, None, classified_error
        else:
            # No recovery
            return False, None, classified_error
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about validation error handling.
        
        Returns:
            Dictionary with validation statistics
        """
        stats = self.stats.copy()
        
        # Calculate handling rate
        if stats["validation_errors"] > 0:
            stats["handling_rate"] = stats["handled_errors"] / stats["validation_errors"]
        else:
            stats["handling_rate"] = 0.0
        
        return stats


# Utility functions

def safe_execute(
    func: Callable[..., T],
    *args,
    default: Optional[T] = None,
    logger: Optional[logging.Logger] = None,
    error_classifier: Optional[ErrorClassifier] = None,
    **kwargs
) -> Tuple[bool, Optional[T], Optional[ClassifiedError]]:
    """
    Execute a function safely, catching and classifying any exceptions.
    
    Args:
        func: Function to execute
        *args: Arguments for the function
        default: Default value to return on error
        logger: Logger instance
        error_classifier: Error classifier
        **kwargs: Keyword arguments for the function
        
    Returns:
        Tuple of (success, result, classified_error)
    """
    # Initialize logger and error classifier if not provided
    logger = logger or logging.getLogger(__name__)
    error_classifier = error_classifier or ErrorClassifier(logger=logger)
    
    try:
        # Execute the function
        result = func(*args, **kwargs)
        return True, result, None
    except Exception as e:
        # Create error context
        context = ErrorContext.current(
            operation=func.__name__,
            additional_context={"args": args, "kwargs": kwargs}
        )
        
        # Classify the error
        classified_error = error_classifier.classify(e, context=context)
        
        # Log the error
        classified_error.log(logger)
        
        # Return failure with default value and classified error
        return False, default, classified_error


async def safe_execute_async(
    func: Callable[..., T],
    *args,
    default: Optional[T] = None,
    logger: Optional[logging.Logger] = None,
    error_classifier: Optional[ErrorClassifier] = None,
    **kwargs
) -> Tuple[bool, Optional[T], Optional[ClassifiedError]]:
    """
    Execute an async function safely, catching and classifying any exceptions.
    
    Args:
        func: Async function to execute
        *args: Arguments for the function
        default: Default value to return on error
        logger: Logger instance
        error_classifier: Error classifier
        **kwargs: Keyword arguments for the function
        
    Returns:
        Tuple of (success, result, classified_error)
    """
    # Initialize logger and error classifier if not provided
    logger = logger or logging.getLogger(__name__)
    error_classifier = error_classifier or ErrorClassifier(logger=logger)
    
    try:
        # Execute the function
        result = await func(*args, **kwargs)
        return True, result, None
    except Exception as e:
        # Create error context
        context = ErrorContext.current(
            operation=func.__name__,
            additional_context={"args": args, "kwargs": kwargs}
        )
        
        # Classify the error
        classified_error = error_classifier.classify(e, context=context)
        
        # Log the error
        classified_error.log(logger)
        
        # Return failure with default value and classified error
        return False, default, classified_error


# Error management class

class ErrorManager:
    """
    Centralized error management system.
    
    This class provides a unified interface for error classification,
    retry management, validation error handling, and other error-related
    functionality.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the error manager.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
        # Components
        self.classifier = ErrorClassifier(logger=self.logger)
        self.retry_manager = RetryManager(logger=self.logger)
        self.validation_handler = ValidationErrorHandler(
            logger=self.logger,
            error_classifier=self.classifier
        )
        
        # Error log
        self.error_log: List[ClassifiedError] = []
        self.max_error_log_size = 1000  # Maximum number of errors to keep in the log
    
    def classify(self, error: Union[Exception, str], context: Optional[ErrorContext] = None) -> ClassifiedError:
        """
        Classify an error.
        
        Args:
            error: Exception or error message
            context: Error context
            
        Returns:
            Classified error
        """
        classified_error = self.classifier.classify(error, context=context)
        
        # Add to error log
        self.add_to_error_log(classified_error)
        
        return classified_error
    
    def retry(self, operation: Callable[..., T], *args, **kwargs) -> T:
        """
        Retry an operation with the configured retry policy.
        
        Args:
            operation: Function to retry
            *args: Arguments for the operation
            **kwargs: Keyword arguments for the operation
            
        Returns:
            Result of the operation
        """
        return self.retry_manager.retry(
            operation,
            *args,
            error_classifier=self.classifier,
            **kwargs
        )
    
    async def retry_async(self, operation: Callable[..., T], *args, **kwargs) -> T:
        """
        Retry an async operation with the configured retry policy.
        
        Args:
            operation: Async function to retry
            *args: Arguments for the operation
            **kwargs: Keyword arguments for the operation
            
        Returns:
            Result of the operation
        """
        return await self.retry_manager.retry_async(
            operation,
            *args,
            error_classifier=self.classifier,
            **kwargs
        )
    
    def handle_validation_error(
        self,
        data: Any,
        error: Exception,
        context: Optional[Dict[str, Any]] = None
    ) -> Tuple[bool, Any]:
        """
        Handle a validation error.
        
        Args:
            data: Data that failed validation
            error: Validation exception
            context: Additional context for error handling
            
        Returns:
            Tuple of (success, result)
        """
        success, result, classified_error = self.validation_handler.handle_validation_error(
            data,
            error,
            context=context
        )
        
        # Add to error log
        if classified_error is not None:
            self.add_to_error_log(classified_error)
        
        return success, result
    
    def add_to_error_log(self, error: ClassifiedError) -> None:
        """
        Add an error to the error log.
        
        Args:
            error: Classified error to add
        """
        self.error_log.append(error)
        
        # Trim log if it exceeds the maximum size
        if len(self.error_log) > self.max_error_log_size:
            self.error_log = self.error_log[-self.max_error_log_size:]
    
    def get_error_log(self, limit: Optional[int] = None) -> List[ClassifiedError]:
        """
        Get the error log.
        
        Args:
            limit: Maximum number of errors to return (newest first)
            
        Returns:
            List of classified errors
        """
        if limit is None:
            return self.error_log.copy()
        
        return self.error_log[-limit:]
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get error management statistics.
        
        Returns:
            Dictionary with error statistics
        """
        return {
            "classifier": self.classifier.get_stats(),
            "retry": self.retry_manager.get_stats(),
            "validation": self.validation_handler.get_stats(),
            "error_log_size": len(self.error_log)
        }
    
    def clear_error_log(self) -> None:
        """Clear the error log."""
        self.error_log.clear()


# Decorator for error handling

def with_error_handling(
    error_manager: Optional[ErrorManager] = None,
    retry: bool = False,
    max_retries: Optional[int] = None,
    log_level: int = logging.ERROR
):
    """
    Decorator for functions to add error handling.
    
    Args:
        error_manager: Error manager to use
        retry: Whether to retry the function on failure
        max_retries: Maximum number of retries
        log_level: Log level for errors
        
    Returns:
        Decorated function
    """
    # Create error manager if not provided
    if error_manager is None:
        logger = logging.getLogger(__name__)
        error_manager = ErrorManager(logger=logger)
    
    def decorator(func):
        """The actual decorator."""
        
        # For async functions
        if asyncio.iscoroutinefunction(func):
            @functools.wraps(func)
            async def async_wrapper(*args, **kwargs):
                """Async wrapper function."""
                if retry:
                    try:
                        return await error_manager.retry_async(func, *args, **kwargs)
                    except Exception as e:
                        # Classify the error
                        context = ErrorContext.current(
                            operation=func.__name__,
                            additional_context={"args": args, "kwargs": kwargs}
                        )
                        classified_error = error_manager.classify(e, context=context)
                        
                        # Log the error
                        error_manager.logger.log(log_level, f"Error in {func.__name__}: {str(e)}", exc_info=e)
                        
                        # Re-raise the exception
                        raise
                else:
                    try:
                        return await func(*args, **kwargs)
                    except Exception as e:
                        # Classify the error
                        context = ErrorContext.current(
                            operation=func.__name__,
                            additional_context={"args": args, "kwargs": kwargs}
                        )
                        classified_error = error_manager.classify(e, context=context)
                        
                        # Log the error
                        error_manager.logger.log(log_level, f"Error in {func.__name__}: {str(e)}", exc_info=e)
                        
                        # Re-raise the exception
                        raise
            
            return async_wrapper
        
        # For synchronous functions
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """Wrapper function."""
            if retry:
                try:
                    return error_manager.retry(func, *args, **kwargs)
                except Exception as e:
                    # Classify the error
                    context = ErrorContext.current(
                        operation=func.__name__,
                        additional_context={"args": args, "kwargs": kwargs}
                    )
                    classified_error = error_manager.classify(e, context=context)
                    
                    # Log the error
                    error_manager.logger.log(log_level, f"Error in {func.__name__}: {str(e)}", exc_info=e)
                    
                    # Re-raise the exception
                    raise
            else:
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    # Classify the error
                    context = ErrorContext.current(
                        operation=func.__name__,
                        additional_context={"args": args, "kwargs": kwargs}
                    )
                    classified_error = error_manager.classify(e, context=context)
                    
                    # Log the error
                    error_manager.logger.log(log_level, f"Error in {func.__name__}: {str(e)}", exc_info=e)
                    
                    # Re-raise the exception
                    raise
        
        return wrapper
    
    return decorator