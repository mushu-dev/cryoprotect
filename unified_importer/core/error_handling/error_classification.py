"""
Error classification system for CryoProtect Unified Importer.

This module provides a framework for classifying errors, determining appropriate
recovery strategies, and capturing detailed context for troubleshooting.
"""

import logging
import time
from enum import Enum, auto
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, Type, Set, List

class ErrorCategory(Enum):
    """Categories of errors that can occur in the application."""
    NETWORK = auto()           # Connection issues, DNS failures
    DATABASE = auto()          # Database connectivity or query problems
    VALIDATION = auto()        # Data validation failures
    RESOURCE = auto()          # Resource unavailability (memory, disk space)
    AUTHORIZATION = auto()     # Permission and authentication errors
    CONFIGURATION = auto()     # Configuration errors, missing settings
    EXTERNAL_SERVICE = auto()  # Failures in external services/APIs
    TIMEOUT = auto()           # Operations that exceed time limits
    DATA_FORMAT = auto()       # Issues with data structure/format
    UNKNOWN = auto()           # Unclassified errors

class ErrorSeverity(Enum):
    """Severity levels for errors."""
    CRITICAL = auto()  # System-wide failures requiring immediate attention
    HIGH = auto()      # Severe issues affecting major functionality
    MEDIUM = auto()    # Significant issues affecting some functionality
    LOW = auto()       # Minor issues with limited impact
    INFO = auto()      # Informational issues with negligible impact

class RecoveryStrategy(Enum):
    """Strategies for recovering from errors."""
    RETRY = auto()           # Attempt the operation again immediately
    FALLBACK = auto()        # Use an alternative implementation or data source
    SKIP = auto()            # Skip the problematic item and continue
    ABORT = auto()           # Abort the current operation entirely
    LOG_ONLY = auto()        # Record the error but continue normally
    DELAYED_RETRY = auto()   # Retry after a delay with exponential backoff
    CIRCUIT_BREAKER = auto() # Prevent cascading failures by stopping attempts

@dataclass
class ErrorContext:
    """Context information about an error."""
    component: str  # Component where the error occurred (e.g., "ChEMBL Importer")
    operation: str  # Operation that was being performed (e.g., "fetch_molecule")
    data: Dict[str, Any] = field(default_factory=dict)  # Data related to the error
    metadata: Dict[str, Any] = field(default_factory=dict)  # Additional metadata

    def __post_init__(self):
        """Initialize default metadata if not provided."""
        if 'timestamp' not in self.metadata:
            self.metadata['timestamp'] = time.time()

@dataclass
class ClassifiedError:
    """An exception with classification information."""
    error: Exception  # The original exception
    category: ErrorCategory  # Category of the error
    severity: ErrorSeverity  # Severity level
    recovery_strategy: RecoveryStrategy  # Recommended recovery strategy
    context: ErrorContext  # Context information about the error

    def __str__(self):
        """String representation of the classified error."""
        return (
            f"{type(self.error).__name__}: {str(self.error)} in "
            f"{self.context.component}:{self.context.operation} "
            f"[Category: {self.category.name}, "
            f"Severity: {self.severity.name}, "
            f"Strategy: {self.recovery_strategy.name}]"
        )

class ErrorClassifier:
    """Classifies exceptions and determines recovery strategies."""

    def __init__(self, logger: Optional[logging.Logger] = None):
        """Initialize the error classifier.
        
        Args:
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        self.rules: Dict[Type[Exception], Dict[str, Any]] = {}
        
        # Register default rules
        self._register_default_rules()

    def _register_default_rules(self):
        """Register default classification rules for common exception types."""
        # Network-related errors
        self.register_rule(
            ConnectionError,
            ErrorCategory.NETWORK,
            ErrorSeverity.HIGH,
            RecoveryStrategy.RETRY
        )
        self.register_rule(
            TimeoutError,
            ErrorCategory.TIMEOUT,
            ErrorSeverity.MEDIUM,
            RecoveryStrategy.DELAYED_RETRY
        )
        
        # Data validation errors
        self.register_rule(
            ValueError,
            ErrorCategory.VALIDATION,
            ErrorSeverity.MEDIUM,
            RecoveryStrategy.SKIP
        )
        self.register_rule(
            TypeError,
            ErrorCategory.VALIDATION,
            ErrorSeverity.MEDIUM,
            RecoveryStrategy.SKIP
        )
        
        # Resource errors
        self.register_rule(
            MemoryError,
            ErrorCategory.RESOURCE,
            ErrorSeverity.CRITICAL,
            RecoveryStrategy.ABORT
        )
        self.register_rule(
            OSError,
            ErrorCategory.RESOURCE,
            ErrorSeverity.HIGH,
            RecoveryStrategy.FALLBACK
        )
        
        # Permission errors
        self.register_rule(
            PermissionError,
            ErrorCategory.AUTHORIZATION,
            ErrorSeverity.HIGH,
            RecoveryStrategy.ABORT
        )

    def register_rule(self, error_class: Type[Exception], category: ErrorCategory,
                      severity: ErrorSeverity, recovery_strategy: RecoveryStrategy):
        """Register a classification rule for an exception type.
        
        Args:
            error_class: Exception class to match
            category: ErrorCategory for this error
            severity: ErrorSeverity for this error
            recovery_strategy: RecoveryStrategy for this error
        """
        self.rules[error_class] = {
            'category': category,
            'severity': severity,
            'recovery_strategy': recovery_strategy
        }
        
        self.logger.debug(
            f"Registered error rule: {error_class.__name__} -> "
            f"{category.name}, {severity.name}, {recovery_strategy.name}"
        )

    def classify(self, error: Exception, context: ErrorContext) -> ClassifiedError:
        """Classify an exception based on registered rules.
        
        Args:
            error: Exception to classify
            context: ErrorContext with information about where the error occurred
            
        Returns:
            ClassifiedError: Classified error with category, severity, and recovery strategy
        """
        # Find the most specific matching rule
        error_type = type(error)
        rule = None
        
        # Check for exact match first
        if error_type in self.rules:
            rule = self.rules[error_type]
        else:
            # Look for match in parent classes
            for rule_type, rule_info in self.rules.items():
                if isinstance(error, rule_type):
                    rule = rule_info
                    break
        
        if rule:
            category = rule['category']
            severity = rule['severity']
            recovery_strategy = rule['recovery_strategy']
        else:
            # Default classification for unknown errors
            category = ErrorCategory.UNKNOWN
            severity = ErrorSeverity.HIGH
            recovery_strategy = RecoveryStrategy.ABORT
            
            self.logger.warning(
                f"No classification rule found for {error_type.__name__}. "
                f"Using default classification: {category.name}, {severity.name}, {recovery_strategy.name}"
            )
        
        # Create classified error
        classified_error = ClassifiedError(
            error=error,
            category=category,
            severity=severity,
            recovery_strategy=recovery_strategy,
            context=context
        )
        
        # Log the classification
        self.logger.debug(f"Classified error: {classified_error}")
        
        return classified_error