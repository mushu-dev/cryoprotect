import unittest
import time
from unittest.mock import Mock, patch
import logging
import sys
import io
from contextlib import contextmanager

# Import the error handling system
from unified_importer.core.error_handling import (
    ErrorCategory, ErrorSeverity, RecoveryStrategy, ErrorContext,
    ClassifiedError, ErrorClassifier, RetryManager, 
    ValidationErrorHandler, ErrorManager
)

class TestErrorHandling(unittest.TestCase):
    
    def setUp(self):
        self.logger = Mock(spec=logging.Logger)
    
    def test_error_category_enum(self):
        """Test that the ErrorCategory enum contains all expected categories."""
        self.assertTrue(hasattr(ErrorCategory, "NETWORK"))
        self.assertTrue(hasattr(ErrorCategory, "DATABASE"))
        self.assertTrue(hasattr(ErrorCategory, "VALIDATION"))
        self.assertTrue(hasattr(ErrorCategory, "RESOURCE"))
        self.assertTrue(hasattr(ErrorCategory, "AUTHORIZATION"))
        self.assertTrue(hasattr(ErrorCategory, "CONFIGURATION"))
        self.assertTrue(hasattr(ErrorCategory, "EXTERNAL_SERVICE"))
        self.assertTrue(hasattr(ErrorCategory, "TIMEOUT"))
        self.assertTrue(hasattr(ErrorCategory, "DATA_FORMAT"))
        self.assertTrue(hasattr(ErrorCategory, "UNKNOWN"))
    
    def test_error_severity_enum(self):
        """Test that the ErrorSeverity enum contains all expected severity levels."""
        self.assertTrue(hasattr(ErrorSeverity, "CRITICAL"))
        self.assertTrue(hasattr(ErrorSeverity, "HIGH"))
        self.assertTrue(hasattr(ErrorSeverity, "MEDIUM"))
        self.assertTrue(hasattr(ErrorSeverity, "LOW"))
        self.assertTrue(hasattr(ErrorSeverity, "INFO"))
    
    def test_recovery_strategy_enum(self):
        """Test that the RecoveryStrategy enum contains all expected strategies."""
        self.assertTrue(hasattr(RecoveryStrategy, "RETRY"))
        self.assertTrue(hasattr(RecoveryStrategy, "FALLBACK"))
        self.assertTrue(hasattr(RecoveryStrategy, "SKIP"))
        self.assertTrue(hasattr(RecoveryStrategy, "ABORT"))
        self.assertTrue(hasattr(RecoveryStrategy, "LOG_ONLY"))
        self.assertTrue(hasattr(RecoveryStrategy, "DELAYED_RETRY"))
        self.assertTrue(hasattr(RecoveryStrategy, "CIRCUIT_BREAKER"))
    
    def test_error_context_creation(self):
        """Test creation of ErrorContext with data and metadata."""
        context = ErrorContext(
            component="ChEMBL Importer",
            operation="fetch_molecule",
            data={"molecule_id": "CHEMBL123"},
            metadata={"attempt": 2, "timestamp": time.time()}
        )
        
        self.assertEqual(context.component, "ChEMBL Importer")
        self.assertEqual(context.operation, "fetch_molecule")
        self.assertEqual(context.data["molecule_id"], "CHEMBL123")
        self.assertEqual(context.metadata["attempt"], 2)
        self.assertIn("timestamp", context.metadata)
    
    def test_classified_error_creation(self):
        """Test creation of ClassifiedError with exception and classification."""
        error = ValueError("Invalid molecular structure")
        context = ErrorContext(
            component="PubChem Importer",
            operation="validate_structure"
        )
        
        classified_error = ClassifiedError(
            error=error,
            category=ErrorCategory.VALIDATION,
            severity=ErrorSeverity.MEDIUM,
            recovery_strategy=RecoveryStrategy.SKIP,
            context=context
        )
        
        self.assertIs(classified_error.error, error)
        self.assertEqual(classified_error.category, ErrorCategory.VALIDATION)
        self.assertEqual(classified_error.severity, ErrorSeverity.MEDIUM)
        self.assertEqual(classified_error.recovery_strategy, RecoveryStrategy.SKIP)
        self.assertIs(classified_error.context, context)
        
        # Test string representation
        error_str = str(classified_error)
        self.assertIn("ValueError", error_str)
        self.assertIn("VALIDATION", error_str)
        self.assertIn("MEDIUM", error_str)
        self.assertIn("SKIP", error_str)
        self.assertIn("PubChem Importer", error_str)
    
    def test_error_classifier_network_errors(self):
        """Test that ErrorClassifier correctly classifies network errors."""
        classifier = ErrorClassifier(logger=self.logger)
        
        # Test with ConnectionError
        error = ConnectionError("Failed to connect to API")
        context = ErrorContext(component="ChEMBL API", operation="fetch_data")
        
        classified = classifier.classify(error, context)
        
        self.assertEqual(classified.category, ErrorCategory.NETWORK)
        self.assertEqual(classified.recovery_strategy, RecoveryStrategy.RETRY)
        
        # Test with TimeoutError
        error = TimeoutError("Connection timed out")
        classified = classifier.classify(error, context)
        
        self.assertEqual(classified.category, ErrorCategory.TIMEOUT)
        self.assertEqual(classified.recovery_strategy, RecoveryStrategy.DELAYED_RETRY)
    
    def test_error_classifier_database_errors(self):
        """Test that ErrorClassifier correctly classifies database errors."""
        classifier = ErrorClassifier(logger=self.logger)
        
        # Mock a database error
        class DatabaseError(Exception):
            pass
        
        error = DatabaseError("Database connection failed")
        context = ErrorContext(component="Supabase", operation="insert_molecule")
        
        # Register custom rule for DatabaseError
        classifier.register_rule(
            error_class=DatabaseError,
            category=ErrorCategory.DATABASE,
            severity=ErrorSeverity.HIGH,
            recovery_strategy=RecoveryStrategy.RETRY
        )
        
        classified = classifier.classify(error, context)
        
        self.assertEqual(classified.category, ErrorCategory.DATABASE)
        self.assertEqual(classified.severity, ErrorSeverity.HIGH)
        self.assertEqual(classified.recovery_strategy, RecoveryStrategy.RETRY)
    
    def test_error_classifier_validation_errors(self):
        """Test that ErrorClassifier correctly classifies validation errors."""
        classifier = ErrorClassifier(logger=self.logger)
        
        error = ValueError("Invalid SMILES format")
        context = ErrorContext(
            component="Molecule Validator", 
            operation="validate_smiles",
            data={"smiles": "C1=C(C=CC=C1)C"}
        )
        
        classified = classifier.classify(error, context)
        
        self.assertEqual(classified.category, ErrorCategory.VALIDATION)
        self.assertEqual(classified.recovery_strategy, RecoveryStrategy.SKIP)
    
    def test_error_classifier_unknown_errors(self):
        """Test that ErrorClassifier correctly handles unknown errors."""
        classifier = ErrorClassifier(logger=self.logger)
        
        class CustomUnknownError(Exception):
            pass
        
        error = CustomUnknownError("Something weird happened")
        context = ErrorContext(component="Custom Component", operation="custom_operation")
        
        classified = classifier.classify(error, context)
        
        self.assertEqual(classified.category, ErrorCategory.UNKNOWN)
        self.assertEqual(classified.recovery_strategy, RecoveryStrategy.ABORT)
    
    def test_retry_manager_successful_retry(self):
        """Test RetryManager with a function that succeeds on retry."""
        retry_manager = RetryManager(
            max_attempts=3,
            initial_delay=0.01,
            max_delay=0.1,
            backoff_factor=2,
            logger=self.logger
        )
        
        mock_function = Mock(side_effect=[ConnectionError("Failed"), ValueError("Still failed"), "Success"])
        
        result = retry_manager.execute(
            mock_function,
            component="Test Component",
            operation="test_operation"
        )
        
        self.assertEqual(result, "Success")
        self.assertEqual(mock_function.call_count, 3)
    
    def test_retry_manager_max_attempts_exceeded(self):
        """Test RetryManager when max attempts are exceeded."""
        retry_manager = RetryManager(
            max_attempts=3,
            initial_delay=0.01,
            max_delay=0.1,
            backoff_factor=2,
            logger=self.logger
        )
        
        mock_function = Mock(side_effect=[
            ConnectionError("Failed 1"),
            ConnectionError("Failed 2"),
            ConnectionError("Failed 3"),
            "Success"  # This should not be reached
        ])
        
        with self.assertRaises(Exception) as context:
            retry_manager.execute(
                mock_function,
                component="Test Component",
                operation="test_operation"
            )
        
        self.assertEqual(mock_function.call_count, 3)
        self.assertIn("Max retry attempts (3) exceeded", str(context.exception))
    
    def test_retry_manager_non_retryable_error(self):
        """Test RetryManager with non-retryable errors."""
        retry_manager = RetryManager(
            max_attempts=3,
            initial_delay=0.01,
            max_delay=0.1,
            backoff_factor=2,
            logger=self.logger
        )
        
        # Register a non-retryable error
        retry_manager.register_non_retryable_error(ValueError)
        
        mock_function = Mock(side_effect=[ValueError("Non-retryable error")])
        
        with self.assertRaises(ValueError):
            retry_manager.execute(
                mock_function,
                component="Test Component",
                operation="test_operation"
            )
        
        self.assertEqual(mock_function.call_count, 1)
    
    def test_validation_error_handler_basic(self):
        """Test ValidationErrorHandler basic functionality."""
        handler = ValidationErrorHandler(logger=self.logger)
        
        # Define a simple validation function that raises ValueError for invalid inputs
        def validate_positive(value):
            if value <= 0:
                raise ValueError("Value must be positive")
            return value
        
        # Test with valid input
        result = handler.handle_validation(
            validate_positive, 42,
            error_strategy=RecoveryStrategy.SKIP,
            component="Number Validator",
            operation="validate_positive"
        )
        self.assertEqual(result, 42)
        
        # Test with invalid input and SKIP strategy
        result = handler.handle_validation(
            validate_positive, -1,
            error_strategy=RecoveryStrategy.SKIP,
            component="Number Validator",
            operation="validate_positive",
            default_value=0
        )
        self.assertEqual(result, 0)
    
    def test_validation_error_handler_strategies(self):
        """Test ValidationErrorHandler with different strategies."""
        handler = ValidationErrorHandler(logger=self.logger)
        
        def validate_in_range(value):
            if not (0 <= value <= 100):
                raise ValueError(f"Value {value} must be between 0 and 100")
            return value
        
        # Test FALLBACK strategy
        result = handler.handle_validation(
            validate_in_range, 150,
            error_strategy=RecoveryStrategy.FALLBACK,
            component="Range Validator",
            operation="validate_in_range",
            fallback_func=lambda x: min(max(0, x), 100)  # Clamp value to range
        )
        self.assertEqual(result, 100)
        
        # Test ABORT strategy
        with self.assertRaises(ValueError):
            handler.handle_validation(
                validate_in_range, 150,
                error_strategy=RecoveryStrategy.ABORT,
                component="Range Validator",
                operation="validate_in_range"
            )
        
        # Test LOG_ONLY strategy
        result = handler.handle_validation(
            validate_in_range, 150,
            error_strategy=RecoveryStrategy.LOG_ONLY,
            component="Range Validator",
            operation="validate_in_range"
        )
        self.assertEqual(result, 150)  # Original value returned despite being invalid
    
    def test_error_manager_integration(self):
        """Test ErrorManager as an integrated interface for error handling."""
        error_manager = ErrorManager(logger=self.logger)
        
        # Test classifying an error
        error = ConnectionError("API connection failed")
        context = ErrorContext(
            component="External API",
            operation="fetch_data"
        )
        
        classified = error_manager.classify_error(error, context)
        self.assertEqual(classified.category, ErrorCategory.NETWORK)
        
        # Test retry with a gradually succeeding function
        counter = [0]  # Using list for mutable closure
        
        def flaky_function():
            counter[0] += 1
            if counter[0] < 3:
                raise ConnectionError(f"Failed attempt {counter[0]}")
            return "Success"
        
        result = error_manager.retry(
            flaky_function,
            component="Flaky Component",
            operation="flaky_operation",
            max_attempts=5
        )
        
        self.assertEqual(result, "Success")
        self.assertEqual(counter[0], 3)
        
        # Test validation handling
        def validate_email(email):
            if "@" not in email:
                raise ValueError("Invalid email format")
            return email
        
        # Valid email
        valid_result = error_manager.validate(
            validate_email, "user@example.com",
            component="User Validator",
            operation="validate_email",
            error_strategy=RecoveryStrategy.ABORT
        )
        self.assertEqual(valid_result, "user@example.com")
        
        # Invalid email with SKIP strategy
        invalid_result = error_manager.validate(
            validate_email, "invalid-email",
            component="User Validator",
            operation="validate_email",
            error_strategy=RecoveryStrategy.SKIP,
            default_value="unknown@example.com"
        )
        self.assertEqual(invalid_result, "unknown@example.com")

    def test_circuit_breaker_functionality(self):
        """Test the circuit breaker functionality of the error handling system."""
        error_manager = ErrorManager(logger=self.logger)
        
        # Configure circuit breaker
        error_manager.configure_circuit_breaker(
            failure_threshold=3,
            reset_timeout=0.1  # Short timeout for testing
        )
        
        # Function that always fails
        def always_fails():
            raise ConnectionError("Service unavailable")
        
        # First 3 calls should attempt and fail
        for _ in range(3):
            with self.assertRaises(ConnectionError):
                error_manager.retry(
                    always_fails,
                    component="External Service",
                    operation="api_call",
                    recovery_strategy=RecoveryStrategy.CIRCUIT_BREAKER
                )
        
        # Next call should fail fast with CircuitBreakerOpenError
        with self.assertRaises(Exception) as context:
            error_manager.retry(
                always_fails,
                component="External Service",
                operation="api_call",
                recovery_strategy=RecoveryStrategy.CIRCUIT_BREAKER
            )
        self.assertIn("Circuit breaker open", str(context.exception))
        
        # Wait for reset timeout
        time.sleep(0.15)
        
        # After timeout, should try again (half-open state)
        with self.assertRaises(ConnectionError):
            error_manager.retry(
                always_fails,
                component="External Service",
                operation="api_call",
                recovery_strategy=RecoveryStrategy.CIRCUIT_BREAKER
            )
    
    def test_bulk_error_handling(self):
        """Test handling errors in bulk operations."""
        error_manager = ErrorManager(logger=self.logger)
        
        # List of items to process
        items = [1, 2, "not_a_number", 4, "another_string", 6]
        
        def process_item(item):
            if isinstance(item, int):
                return item * 2
            raise TypeError(f"Cannot process {item}, expected int")
        
        # Process with SKIP strategy for errors
        results = error_manager.process_bulk(
            process_item,
            items,
            component="Bulk Processor",
            operation="process_numbers",
            error_strategy=RecoveryStrategy.SKIP
        )
        
        # Should have results for valid items and None for invalid ones
        self.assertEqual(len(results), 6)
        self.assertEqual(results[0], 2)  # 1 * 2
        self.assertEqual(results[1], 4)  # 2 * 2
        self.assertIsNone(results[2])    # Error, returned None
        self.assertEqual(results[3], 8)  # 4 * 2
        self.assertIsNone(results[4])    # Error, returned None
        self.assertEqual(results[5], 12) # 6 * 2
        
        # Collect errors
        errors = []
        
        # Process with custom error callback
        results = error_manager.process_bulk(
            process_item,
            items,
            component="Bulk Processor",
            operation="process_numbers",
            error_strategy=RecoveryStrategy.SKIP,
            error_callback=lambda err, item: errors.append((err, item))
        )
        
        # Should have collected errors
        self.assertEqual(len(errors), 2)
        self.assertIsInstance(errors[0][0], TypeError)
        self.assertEqual(errors[0][1], "not_a_number")
        self.assertIsInstance(errors[1][0], TypeError)
        self.assertEqual(errors[1][1], "another_string")


if __name__ == "__main__":
    unittest.main()