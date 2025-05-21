"""
Unit tests for the enhanced retry mechanism.
"""

import unittest
import time
import asyncio
import os
import json
import tempfile
from unittest.mock import Mock, patch, MagicMock
import logging

from unified_importer.core.error_handling.error_classification import ErrorCategory, RecoveryStrategy, ErrorContext
from unified_importer.core.error_handling.retry_enhancement import (
    RetryConfig, CircuitBreakerConfig, CircuitBreakerRegistry, 
    CircuitState, EnhancedRetryManager
)

class TestRetryEnhancement(unittest.TestCase):
    """Tests for the enhanced retry mechanism."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.logger = logging.getLogger("test")
        self.logger.setLevel(logging.DEBUG)
        
        # Avoid actual sleep in tests
        self.sleep_patcher = patch('time.sleep')
        self.mock_sleep = self.sleep_patcher.start()
        
        # Avoid actual asyncio.sleep in tests
        self.async_sleep_patcher = patch('asyncio.sleep')
        self.mock_async_sleep = self.async_sleep_patcher.start()
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.sleep_patcher.stop()
        self.async_sleep_patcher.stop()
    
    def test_retry_config_defaults(self):
        """Test RetryConfig default values."""
        config = RetryConfig()
        
        self.assertEqual(config.max_attempts, 3)
        self.assertEqual(config.initial_delay, 0.1)
        self.assertEqual(config.max_delay, 60.0)
        self.assertEqual(config.backoff_factor, 2.0)
        self.assertEqual(config.jitter_factor, 0.1)
        self.assertEqual(config.non_retryable_errors, set())
        self.assertEqual(config.retryable_error_categories, set())
    
    def test_circuit_breaker_config_defaults(self):
        """Test CircuitBreakerConfig default values."""
        config = CircuitBreakerConfig()
        
        self.assertEqual(config.failure_threshold, 5)
        self.assertEqual(config.reset_timeout, 60.0)
        self.assertEqual(config.half_open_max_calls, 1)
        self.assertEqual(config.success_threshold, 1)
        self.assertEqual(config.open_timeout_increment_factor, 1.5)
        self.assertEqual(config.max_open_timeout, 300.0)
    
    def test_circuit_breaker_registry_creation(self):
        """Test creating a circuit breaker registry."""
        registry = CircuitBreakerRegistry()
        self.assertEqual(len(registry._breakers), 0)
    
    def test_circuit_breaker_registry_get_breaker(self):
        """Test getting a circuit breaker from registry."""
        registry = CircuitBreakerRegistry()
        
        # Get a new circuit breaker
        breaker = registry.get_breaker("TestService", "test_operation")
        
        self.assertEqual(breaker.state, CircuitState.CLOSED)
        self.assertEqual(breaker.failure_count, 0)
        self.assertEqual(breaker.service_identifier, "TestService")
        self.assertEqual(breaker.operation_identifier, "test_operation")
        
        # Get the same circuit breaker again
        breaker2 = registry.get_breaker("TestService", "test_operation")
        self.assertIs(breaker, breaker2)
    
    def test_circuit_breaker_should_execute(self):
        """Test circuit breaker execution decisions."""
        registry = CircuitBreakerRegistry()
        breaker = registry.get_breaker("TestService", "test_operation")
        
        # Closed state should allow execution
        self.assertEqual(breaker.state, CircuitState.CLOSED)
        self.assertTrue(breaker.should_execute())
        
        # Open state should prevent execution
        breaker.state = CircuitState.OPEN
        breaker.open_until = time.time() + 10.0
        self.assertFalse(breaker.should_execute())
        
        # Open state should transition to half-open after timeout
        breaker.open_until = time.time() - 1.0
        self.assertTrue(breaker.should_execute())
        self.assertEqual(breaker.state, CircuitState.HALF_OPEN)
        
        # Half-open should limit calls
        breaker.state = CircuitState.HALF_OPEN
        breaker.half_open_calls = 0
        breaker.half_open_max_calls = 2  # Allow 2 calls in half-open
        
        self.assertTrue(breaker.should_execute())  # First call allowed
        self.assertEqual(breaker.half_open_calls, 1)
        
        self.assertTrue(breaker.should_execute())  # Second call allowed
        self.assertEqual(breaker.half_open_calls, 2)
        
        self.assertFalse(breaker.should_execute())  # Third call blocked
    
    def test_circuit_breaker_report_success(self):
        """Test reporting success to circuit breaker."""
        registry = CircuitBreakerRegistry()
        config = CircuitBreakerConfig()
        
        # Test reporting success in closed state
        breaker = registry.get_breaker("TestService", "closed_test")
        breaker.failure_count = 2
        registry.report_success("TestService", "closed_test")
        
        self.assertEqual(breaker.failure_count, 0)
        self.assertEqual(breaker.state, CircuitState.CLOSED)
        
        # Test reporting success in half-open state
        breaker = registry.get_breaker("TestService", "half_open_test")
        breaker.state = CircuitState.HALF_OPEN
        breaker.success_count = 0
        breaker.success_threshold = 2
        
        registry.report_success("TestService", "half_open_test")
        self.assertEqual(breaker.success_count, 1)
        self.assertEqual(breaker.state, CircuitState.HALF_OPEN)
        
        registry.report_success("TestService", "half_open_test")
        self.assertEqual(breaker.success_count, 2)
        self.assertEqual(breaker.state, CircuitState.CLOSED)
        self.assertEqual(breaker.failure_count, 0)
    
    def test_circuit_breaker_report_failure(self):
        """Test reporting failure to circuit breaker."""
        registry = CircuitBreakerRegistry()
        config = CircuitBreakerConfig(
            failure_threshold=3,
            reset_timeout=30.0,
            open_timeout_increment_factor=2.0,
            max_open_timeout=120.0
        )
        
        # Test reporting failure in closed state
        breaker = registry.get_breaker("TestService", "closed_test")
        breaker.failure_count = 0
        
        # First failure
        registry.report_failure("TestService", "closed_test", config)
        self.assertEqual(breaker.failure_count, 1)
        self.assertEqual(breaker.state, CircuitState.CLOSED)
        
        # Second failure
        registry.report_failure("TestService", "closed_test", config)
        self.assertEqual(breaker.failure_count, 2)
        self.assertEqual(breaker.state, CircuitState.CLOSED)
        
        # Third failure should open the circuit
        registry.report_failure("TestService", "closed_test", config)
        self.assertEqual(breaker.failure_count, 3)
        self.assertEqual(breaker.state, CircuitState.OPEN)
        self.assertEqual(breaker.current_timeout, 30.0)
        
        # Test reporting failure in half-open state
        breaker = registry.get_breaker("TestService", "half_open_test")
        breaker.state = CircuitState.HALF_OPEN
        breaker.current_timeout = 30.0
        
        registry.report_failure("TestService", "half_open_test", config)
        self.assertEqual(breaker.state, CircuitState.OPEN)
        self.assertEqual(breaker.current_timeout, 60.0)  # Doubled
        
        # Test timeout capping at max
        breaker.current_timeout = 100.0
        registry.report_failure("TestService", "half_open_test", config)
        self.assertEqual(breaker.current_timeout, 120.0)  # Capped at max
    
    def test_circuit_breaker_state_persistence(self):
        """Test circuit breaker state persistence."""
        # Create a temporary file for state persistence
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            state_file = temp_file.name
        
        try:
            # Create registry with persistence
            registry = CircuitBreakerRegistry(state_file_path=state_file, logger=self.logger)
            
            # Set up circuit breaker state
            breaker = registry.get_breaker("TestService", "persist_test")
            breaker.state = CircuitState.OPEN
            breaker.failure_count = 5
            breaker.open_until = time.time() + 60.0
            breaker.current_timeout = 30.0
            
            # Save state
            registry._save_states()
            
            # Create new registry to load state
            registry2 = CircuitBreakerRegistry(state_file_path=state_file, logger=self.logger)
            
            # Get breaker and check state
            breaker2 = registry2.get_breaker("TestService", "persist_test")
            self.assertEqual(breaker2.state, CircuitState.OPEN)
            self.assertEqual(breaker2.failure_count, 5)
            self.assertAlmostEqual(breaker2.open_until, breaker.open_until, delta=1.0)
            self.assertEqual(breaker2.current_timeout, 30.0)
        finally:
            # Clean up temporary file
            if os.path.exists(state_file):
                os.remove(state_file)
    
    def test_enhanced_retry_manager_creation(self):
        """Test creating an enhanced retry manager."""
        retry_manager = EnhancedRetryManager()
        
        # Check default configurations
        self.assertIsNotNone(retry_manager.default_config)
        self.assertIsNotNone(retry_manager.circuit_breaker_config)
        self.assertIsNotNone(retry_manager.circuit_breaker_registry)
        self.assertIsNotNone(retry_manager.error_classifier)
        self.assertIsNotNone(retry_manager.logger)
        
        # Check default retryable error categories
        self.assertIn(ErrorCategory.NETWORK, retry_manager.default_config.retryable_error_categories)
        self.assertIn(ErrorCategory.TIMEOUT, retry_manager.default_config.retryable_error_categories)
        self.assertIn(ErrorCategory.DATABASE, retry_manager.default_config.retryable_error_categories)
        self.assertIn(ErrorCategory.EXTERNAL_SERVICE, retry_manager.default_config.retryable_error_categories)
        
        # Check default non-retryable errors
        self.assertIn(ValueError, retry_manager.default_config.non_retryable_errors)
        self.assertIn(TypeError, retry_manager.default_config.non_retryable_errors)
    
    def test_execute_success(self):
        """Test successful execution without retries."""
        retry_manager = EnhancedRetryManager(logger=self.logger)
        
        # Mock function that succeeds
        mock_func = Mock(return_value="success")
        
        # Execute function
        result = retry_manager.execute(
            mock_func, "TestComponent", "test_operation", 
            args=("arg1", "arg2"), kwargs={"kwarg1": "value1"}
        )
        
        # Verify results
        self.assertEqual(result, "success")
        mock_func.assert_called_once_with("arg1", "arg2", kwarg1="value1")
        self.mock_sleep.assert_not_called()
    
    def test_execute_with_retries(self):
        """Test execution with retries before success."""
        retry_manager = EnhancedRetryManager(logger=self.logger)
        
        # Mock function that fails twice then succeeds
        mock_func = Mock(side_effect=[
            ConnectionError("Network error"),
            TimeoutError("Timeout error"),
            "success"
        ])
        
        # Execute function
        result = retry_manager.execute(
            mock_func, "TestComponent", "test_operation", 
            config=RetryConfig(max_attempts=5, initial_delay=0.1, backoff_factor=2.0)
        )
        
        # Verify results
        self.assertEqual(result, "success")
        self.assertEqual(mock_func.call_count, 3)
        
        # Verify sleep was called twice with increasing delays
        self.assertEqual(self.mock_sleep.call_count, 2)
        sleep_args = [call_args[0][0] for call_args in self.mock_sleep.call_args_list]
        self.assertAlmostEqual(sleep_args[0], 0.1, delta=0.05)  # Initial delay with jitter
        self.assertAlmostEqual(sleep_args[1], 0.2, delta=0.1)   # Second delay with jitter
    
    def test_execute_max_retries_exceeded(self):
        """Test execution with max retries exceeded."""
        retry_manager = EnhancedRetryManager(logger=self.logger)
        
        # Mock function that always fails
        mock_func = Mock(side_effect=ConnectionError("Network error"))
        
        # Execute function with retries
        with self.assertRaises(Exception) as context:
            retry_manager.execute(
                mock_func, "TestComponent", "test_operation", 
                config=RetryConfig(max_attempts=3, initial_delay=0.1)
            )
        
        # Verify results
        self.assertEqual(mock_func.call_count, 3)
        self.assertIn("Max retry attempts (3) exceeded", str(context.exception))
        self.assertIn("ConnectionError: Network error", str(context.exception))
    
    def test_execute_non_retryable_error(self):
        """Test execution with non-retryable error."""
        retry_manager = EnhancedRetryManager(logger=self.logger)
        
        # Mock function that raises a non-retryable error
        mock_func = Mock(side_effect=ValueError("Invalid value"))
        
        # Execute function
        with self.assertRaises(ValueError) as context:
            retry_manager.execute(
                mock_func, "TestComponent", "test_operation", 
                config=RetryConfig(max_attempts=3, initial_delay=0.1)
            )
        
        # Verify results
        mock_func.assert_called_once()
        self.assertEqual(str(context.exception), "Invalid value")
        self.mock_sleep.assert_not_called()
    
    def test_execute_with_circuit_breaker(self):
        """Test execution with circuit breaker."""
        registry = CircuitBreakerRegistry(logger=self.logger)
        retry_manager = EnhancedRetryManager(
            circuit_breaker_registry=registry,
            circuit_breaker_config=CircuitBreakerConfig(failure_threshold=2),
            logger=self.logger
        )
        
        # Mock function that always fails
        mock_func = Mock(side_effect=ConnectionError("Network error"))
        
        # First execution - should try and fail after retries
        with self.assertRaises(Exception):
            retry_manager.execute(
                mock_func, "TestComponent", "circuit_test", 
                config=RetryConfig(max_attempts=2, initial_delay=0.1)
            )
        
        # Second execution - should also try and fail
        with self.assertRaises(Exception):
            retry_manager.execute(
                mock_func, "TestComponent", "circuit_test", 
                config=RetryConfig(max_attempts=2, initial_delay=0.1)
            )
        
        # Third execution - circuit should be open
        with self.assertRaises(Exception) as context:
            retry_manager.execute(
                mock_func, "TestComponent", "circuit_test", 
                config=RetryConfig(max_attempts=2, initial_delay=0.1)
            )
        
        # Verify circuit breaker opened
        self.assertIn("Circuit breaker open", str(context.exception))
        
        # Verify function was not called on third execution
        self.assertEqual(mock_func.call_count, 4)  # 2 calls x 2 executions
    
    def test_circuit_breaker_half_open_recovery(self):
        """Test circuit breaker half-open state and recovery."""
        registry = CircuitBreakerRegistry(logger=self.logger)
        retry_manager = EnhancedRetryManager(
            circuit_breaker_registry=registry,
            circuit_breaker_config=CircuitBreakerConfig(
                failure_threshold=2,
                reset_timeout=60.0
            ),
            logger=self.logger
        )
        
        # Get circuit breaker directly from registry
        circuit = registry.get_breaker("TestComponent", "recovery_test")
        
        # Manually set circuit to open but expired
        circuit.state = CircuitState.OPEN
        circuit.open_until = time.time() - 1.0  # Already expired
        
        # Mock function that now succeeds
        mock_func = Mock(return_value="success")
        
        # Execute function - should transition to half-open and succeed
        result = retry_manager.execute(
            mock_func, "TestComponent", "recovery_test"
        )
        
        # Verify results
        self.assertEqual(result, "success")
        self.assertEqual(circuit.state, CircuitState.CLOSED)
    
    async def test_execute_async_success(self):
        """Test successful async execution without retries."""
        retry_manager = EnhancedRetryManager(logger=self.logger)
        
        # Mock async function that succeeds
        mock_func = AsyncMock(return_value="async success")
        
        # Execute async function
        result = await retry_manager.execute_async(
            mock_func, "TestComponent", "async_operation", 
            args=("arg1", "arg2"), kwargs={"kwarg1": "value1"}
        )
        
        # Verify results
        self.assertEqual(result, "async success")
        mock_func.assert_called_once_with("arg1", "arg2", kwarg1="value1")
        self.mock_async_sleep.assert_not_called()
    
    async def test_execute_async_with_retries(self):
        """Test async execution with retries before success."""
        retry_manager = EnhancedRetryManager(logger=self.logger)
        
        # Mock async function that fails twice then succeeds
        mock_func = AsyncMock(side_effect=[
            ConnectionError("Network error"),
            TimeoutError("Timeout error"),
            "async success"
        ])
        
        # Execute async function
        result = await retry_manager.execute_async(
            mock_func, "TestComponent", "async_operation", 
            config=RetryConfig(max_attempts=5, initial_delay=0.1, backoff_factor=2.0)
        )
        
        # Verify results
        self.assertEqual(result, "async success")
        self.assertEqual(mock_func.call_count, 3)
        
        # Verify sleep was called twice with increasing delays
        self.assertEqual(self.mock_async_sleep.call_count, 2)
    
    def test_retry_decorator(self):
        """Test retry decorator functionality."""
        retry_manager = EnhancedRetryManager(logger=self.logger)
        
        # Create a decorated function
        counter = [0]
        
        @retry_manager.retry_decorator(
            component="TestComponent", 
            operation="decorated_func",
            max_attempts=3,
            initial_delay=0.1
        )
        def flaky_function(arg1, arg2=None):
            counter[0] += 1
            if counter[0] < 3:
                raise ConnectionError(f"Failed attempt {counter[0]}")
            return f"Success with {arg1} and {arg2}"
        
        # Call the decorated function
        result = flaky_function("test1", arg2="test2")
        
        # Verify results
        self.assertEqual(result, "Success with test1 and test2")
        self.assertEqual(counter[0], 3)
        self.assertEqual(self.mock_sleep.call_count, 2)
    
    def test_configure_from_dict(self):
        """Test configuring retry manager from dictionary."""
        retry_manager = EnhancedRetryManager(logger=self.logger)
        
        # Configure from dictionary
        retry_manager.configure_from_dict({
            'retry': {
                'max_attempts': 5,
                'initial_delay': 0.5,
                'max_delay': 30.0,
                'backoff_factor': 1.5,
                'jitter_factor': 0.2
            },
            'circuit_breaker': {
                'failure_threshold': 3,
                'reset_timeout': 120.0,
                'half_open_max_calls': 2,
                'success_threshold': 2,
                'open_timeout_increment_factor': 2.0,
                'max_open_timeout': 600.0
            }
        })
        
        # Verify retry configuration
        self.assertEqual(retry_manager.default_config.max_attempts, 5)
        self.assertEqual(retry_manager.default_config.initial_delay, 0.5)
        self.assertEqual(retry_manager.default_config.max_delay, 30.0)
        self.assertEqual(retry_manager.default_config.backoff_factor, 1.5)
        self.assertEqual(retry_manager.default_config.jitter_factor, 0.2)
        
        # Verify circuit breaker configuration
        self.assertEqual(retry_manager.circuit_breaker_config.failure_threshold, 3)
        self.assertEqual(retry_manager.circuit_breaker_config.reset_timeout, 120.0)
        self.assertEqual(retry_manager.circuit_breaker_config.half_open_max_calls, 2)
        self.assertEqual(retry_manager.circuit_breaker_config.success_threshold, 2)
        self.assertEqual(retry_manager.circuit_breaker_config.open_timeout_increment_factor, 2.0)
        self.assertEqual(retry_manager.circuit_breaker_config.max_open_timeout, 600.0)
    
    def test_get_circuit_breaker_status(self):
        """Test getting circuit breaker status."""
        registry = CircuitBreakerRegistry(logger=self.logger)
        retry_manager = EnhancedRetryManager(
            circuit_breaker_registry=registry,
            logger=self.logger
        )
        
        # Set up some circuit breakers
        circuit1 = registry.get_breaker("Service1", "operation1")
        circuit1.state = CircuitState.CLOSED
        
        circuit2 = registry.get_breaker("Service2", "operation2")
        circuit2.state = CircuitState.OPEN
        circuit2.failure_count = 5
        circuit2.open_until = time.time() + 60.0
        
        # Get status
        status = retry_manager.get_circuit_breaker_status()
        
        # Verify status
        self.assertEqual(len(status), 2)
        self.assertEqual(status["Service1:operation1"]["state"], "CLOSED")
        self.assertEqual(status["Service2:operation2"]["state"], "OPEN")
        self.assertEqual(status["Service2:operation2"]["failure_count"], 5)
    
    def test_reset_circuit_breaker(self):
        """Test resetting a circuit breaker."""
        registry = CircuitBreakerRegistry(logger=self.logger)
        retry_manager = EnhancedRetryManager(
            circuit_breaker_registry=registry,
            logger=self.logger
        )
        
        # Set up a circuit breaker in open state
        circuit = registry.get_breaker("Service", "reset_test")
        circuit.state = CircuitState.OPEN
        circuit.failure_count = 5
        circuit.current_timeout = 60.0
        
        # Reset the circuit breaker
        retry_manager.reset_circuit_breaker("Service", "reset_test")
        
        # Verify reset
        self.assertEqual(circuit.state, CircuitState.CLOSED)
        self.assertEqual(circuit.failure_count, 0)
        self.assertEqual(circuit.current_timeout, 0.0)


class AsyncMock(MagicMock):
    """Mock class for async functions."""
    
    async def __call__(self, *args, **kwargs):
        return super(AsyncMock, self).__call__(*args, **kwargs)


if __name__ == '__main__':
    unittest.main()