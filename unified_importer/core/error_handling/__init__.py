"""
Error handling system for the CryoProtect Unified Importer.

This module provides a comprehensive framework for managing errors in a consistent,
recoverable, and informative way throughout the application.

Components:
- Error Classification: Categorizes errors and determines recovery strategies
- Enhanced Retry Mechanism: Implements retry logic with circuit breakers
- Validation Error Handling: Manages validation errors with customizable strategies
"""

# Import from error classification
from .error_classification import (
    ErrorCategory,
    ErrorSeverity,
    RecoveryStrategy,
    ErrorContext,
    ClassifiedError,
    ErrorClassifier
)

# Import from retry enhancement
from .retry_enhancement import (
    RetryConfig,
    CircuitBreakerConfig,
    CircuitState,
    CircuitBreakerState,
    CircuitBreakerRegistry,
    EnhancedRetryManager
)

# Import from validation handling
from .validation_handling import (
    ValidationResult,
    ValidationError,
    ValidationReport,
    ValidationErrorHandler,
    ChemicalValidationRules
)

# Create a unified error management interface
class ErrorManager:
    """Unified interface for all error handling functionality."""
    
    def __init__(self, logger=None):
        """Initialize the error manager.
        
        Args:
            logger: Optional logger instance
        """
        from logging import getLogger
        self.logger = logger or getLogger(__name__)
        self.error_classifier = ErrorClassifier(logger=self.logger)
        self.retry_manager = EnhancedRetryManager(
            error_classifier=self.error_classifier,
            logger=self.logger
        )
        self.validation_handler = ValidationErrorHandler(logger=self.logger)
    
    def classify_error(self, error, context):
        """Classify an error using the error classifier.
        
        Args:
            error: Exception to classify
            context: ErrorContext with information about where the error occurred
            
        Returns:
            ClassifiedError: Classified error with category, severity, and recovery strategy
        """
        return self.error_classifier.classify(error, context)
    
    def register_error_rule(self, error_class, category, severity, recovery_strategy):
        """Register a custom error classification rule.
        
        Args:
            error_class: Exception class to match
            category: ErrorCategory for this error
            severity: ErrorSeverity for this error
            recovery_strategy: RecoveryStrategy for this error
        """
        self.error_classifier.register_rule(
            error_class, category, severity, recovery_strategy
        )
    
    def retry(self, func, component, operation, *args, 
              max_attempts=None, initial_delay=None, max_delay=None,
              backoff_factor=None, jitter_factor=None, 
              use_circuit_breaker=True, **kwargs):
        """Execute a function with retry logic.
        
        Args:
            func: Function to execute
            component: Component name for context
            operation: Operation name for context
            *args: Arguments to pass to the function
            max_attempts: Maximum retry attempts
            initial_delay: Initial delay before retry
            max_delay: Maximum delay between retries
            backoff_factor: Backoff multiplier
            jitter_factor: Jitter factor
            use_circuit_breaker: Whether to use circuit breaker pattern
            **kwargs: Keyword arguments to pass to the function
            
        Returns:
            Result from the function
        """
        return self.retry_manager.retry(
            func, component, operation, *args,
            max_attempts=max_attempts,
            initial_delay=initial_delay,
            max_delay=max_delay,
            backoff_factor=backoff_factor,
            jitter_factor=jitter_factor,
            use_circuit_breaker=use_circuit_breaker,
            **kwargs
        )
    
    def retry_decorator(self, component, operation, 
                        max_attempts=None, initial_delay=None, max_delay=None,
                        backoff_factor=None, jitter_factor=None, 
                        use_circuit_breaker=True):
        """Create a retry decorator.
        
        Args:
            component: Component name for context
            operation: Operation name for context
            max_attempts: Maximum retry attempts
            initial_delay: Initial delay before retry
            max_delay: Maximum delay between retries
            backoff_factor: Backoff multiplier
            jitter_factor: Jitter factor
            use_circuit_breaker: Whether to use circuit breaker pattern
            
        Returns:
            Decorator function
        """
        return self.retry_manager.retry_decorator(
            component, operation,
            max_attempts=max_attempts,
            initial_delay=initial_delay,
            max_delay=max_delay,
            backoff_factor=backoff_factor,
            jitter_factor=jitter_factor,
            use_circuit_breaker=use_circuit_breaker
        )
    
    def validate(self, validation_func, data, component, operation,
                 error_strategy=RecoveryStrategy.SKIP, default_value=None,
                 fallback_func=None, record_errors=True):
        """Handle validation with appropriate error strategy.
        
        Args:
            validation_func: Function to validate data
            data: Data to validate
            component: Component name for context
            operation: Operation name for context
            error_strategy: Strategy for handling validation errors
            default_value: Default value to return on error with SKIP strategy
            fallback_func: Function to call with data on error with FALLBACK strategy
            record_errors: Whether to record errors in the validation report
            
        Returns:
            Validated data or default value
        """
        return self.validation_handler.handle_validation(
            validation_func, data,
            component=component,
            operation=operation,
            error_strategy=error_strategy,
            default_value=default_value,
            fallback_func=fallback_func,
            record_errors=record_errors
        )
    
    def validate_batch(self, validation_func, items, component, operation,
                      error_strategy=RecoveryStrategy.SKIP, default_value=None,
                      fallback_func=None, parallel=False, max_workers=None,
                      continue_on_error=True):
        """Validate a batch of items with error handling.
        
        Args:
            validation_func: Function to validate each item
            items: List of items to validate
            component: Component name for context
            operation: Operation name for context
            error_strategy: Strategy for handling errors
            default_value: Default value for invalid items with SKIP strategy
            fallback_func: Function for invalid items with FALLBACK strategy
            parallel: Whether to process in parallel
            max_workers: Number of parallel workers
            continue_on_error: Whether to continue on error or abort batch
            
        Returns:
            Tuple of (list of validated items, validation report)
        """
        return self.validation_handler.validate_batch(
            validation_func, items,
            component=component,
            operation=operation,
            error_strategy=error_strategy,
            default_value=default_value,
            fallback_func=fallback_func,
            parallel=parallel,
            max_workers=max_workers,
            continue_on_error=continue_on_error
        )
    
    def get_validation_report(self):
        """Get the current validation report.
        
        Returns:
            Current ValidationReport
        """
        return self.validation_handler.get_report()
    
    def reset_validation_report(self):
        """Reset the current validation report."""
        self.validation_handler.reset_report()
    
    def save_validation_report(self, file_path):
        """Save the current validation report to a file.
        
        Args:
            file_path: Path to save the report
        """
        self.validation_handler.save_report(file_path)
    
    def process_bulk(self, process_func, items, component, operation,
                     error_strategy=RecoveryStrategy.SKIP, parallel=False,
                     max_workers=None, error_callback=None):
        """Process items in bulk with error handling.
        
        Args:
            process_func: Function to process each item
            items: List of items to process
            component: Component name for context
            operation: Operation name for context
            error_strategy: Strategy for handling errors
            parallel: Whether to process in parallel
            max_workers: Number of parallel workers
            error_callback: Function to call with (error, item) on error
            
        Returns:
            List of results (None for failed items with SKIP strategy)
        """
        results = [None] * len(items)
        
        if parallel:
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = []
                for i, item in enumerate(items):
                    futures.append(executor.submit(
                        self._process_item, process_func, item, component, operation,
                        error_strategy, error_callback, i
                    ))
                
                for future in futures:
                    # Result is a tuple of (result, index)
                    result, index = future.result()
                    results[index] = result
        else:
            for i, item in enumerate(items):
                result, _ = self._process_item(
                    process_func, item, component, operation,
                    error_strategy, error_callback, i
                )
                results[i] = result
        
        return results
    
    def _process_item(self, process_func, item, component, operation,
                     error_strategy, error_callback, index):
        """Process a single item with error handling.
        
        Args:
            process_func: Function to process the item
            item: Item to process
            component: Component name for context
            operation: Operation name for context
            error_strategy: Strategy for handling errors
            error_callback: Function to call with (error, item) on error
            index: Index of the item in the original list
            
        Returns:
            Tuple of (result, index)
        """
        try:
            result = process_func(item)
            return result, index
        except Exception as e:
            context = ErrorContext(
                component=component,
                operation=operation,
                data={"item": str(item)[:100], "index": index}
            )
            
            classified_error = self.classify_error(e, context)
            
            self.logger.warning(
                f"Error processing item {index} in {component}:{operation}: "
                f"{type(e).__name__}: {str(e)}. "
                f"Using strategy: {error_strategy.name}"
            )
            
            if error_callback:
                error_callback(e, item)
            
            if error_strategy == RecoveryStrategy.SKIP:
                return None, index
            elif error_strategy == RecoveryStrategy.LOG_ONLY:
                return item, index
            else:
                raise
    
    def execute_with_handling(self, func, component, operation, *args,
                             default_value=None, **kwargs):
        """Execute a function with comprehensive error handling.
        
        Args:
            func: Function to execute
            component: Component name for context
            operation: Operation name for context
            *args: Arguments to pass to the function
            default_value: Default value to return on error with SKIP strategy
            **kwargs: Keyword arguments to pass to the function
            
        Returns:
            Result from the function or default value on error
        """
        try:
            return func(*args, **kwargs)
        except Exception as e:
            context = ErrorContext(
                component=component,
                operation=operation,
                data={"args": str(args)[:100], "kwargs": str(kwargs)[:100]}
            )
            
            classified_error = self.classify_error(e, context)
            
            self.logger.warning(
                f"Error in {component}:{operation}: "
                f"{type(e).__name__}: {str(e)}. "
                f"Category: {classified_error.category.name}, "
                f"Strategy: {classified_error.recovery_strategy.name}"
            )
            
            if classified_error.recovery_strategy in {
                RecoveryStrategy.RETRY, 
                RecoveryStrategy.DELAYED_RETRY,
                RecoveryStrategy.CIRCUIT_BREAKER
            }:
                # Try with retry
                try:
                    return self.retry(func, component, operation, *args, **kwargs)
                except Exception as retry_e:
                    self.logger.error(
                        f"Retry failed in {component}:{operation}: {str(retry_e)}"
                    )
                    return default_value
            
            elif classified_error.recovery_strategy == RecoveryStrategy.SKIP:
                return default_value
            
            elif classified_error.recovery_strategy == RecoveryStrategy.LOG_ONLY:
                return None
            
            else:
                # For ABORT and other strategies, re-raise
                raise
    
    def configure_from_dict(self, config_dict):
        """Configure error handling from a dictionary.
        
        Args:
            config_dict: Dictionary with configuration parameters
        """
        self.retry_manager.configure_from_dict(config_dict)
    
    def get_circuit_breaker_status(self):
        """Get the status of all circuit breakers.
        
        Returns:
            Dictionary of circuit breaker states
        """
        return self.retry_manager.get_circuit_breaker_status()
    
    def reset_circuit_breaker(self, service, operation):
        """Reset a circuit breaker to closed state.
        
        Args:
            service: Service identifier
            operation: Operation identifier
        """
        self.retry_manager.reset_circuit_breaker(service, operation)


# For convenience, expose key symbols at the module level
__all__ = [
    # Error Classification
    'ErrorCategory',
    'ErrorSeverity',
    'RecoveryStrategy',
    'ErrorContext',
    'ClassifiedError',
    'ErrorClassifier',
    
    # Retry Enhancement
    'RetryConfig',
    'CircuitBreakerConfig',
    'CircuitState',
    'CircuitBreakerState',
    'CircuitBreakerRegistry',
    'EnhancedRetryManager',
    
    # Validation Handling
    'ValidationResult',
    'ValidationError',
    'ValidationReport',
    'ValidationErrorHandler',
    'ChemicalValidationRules',
    
    # Unified Interface
    'ErrorManager'
]