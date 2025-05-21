"""
Validation error handling system for CryoProtect Unified Importer.

This module provides a framework for handling data validation errors with
customizable strategies and detailed reporting.
"""

import logging
import time
import json
import os
from enum import Enum
from typing import Dict, Any, Optional, List, Callable, TypeVar, Generic, Union, Set, Tuple
from dataclasses import dataclass, field
from datetime import datetime

from .error_classification import ErrorCategory, ErrorSeverity, RecoveryStrategy, ErrorContext, ClassifiedError

# Type variables for generic functions
T = TypeVar('T')
U = TypeVar('U')

@dataclass
class ValidationResult:
    """Result of a validation operation."""
    is_valid: bool
    data: Any
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    context: Dict[str, Any] = field(default_factory=dict)
    
    def __bool__(self):
        """Allow using validation result in boolean context."""
        return self.is_valid

@dataclass
class ValidationError:
    """Detailed information about a validation error."""
    message: str
    data: Any
    field: Optional[str] = None
    code: Optional[str] = None
    severity: ErrorSeverity = ErrorSeverity.MEDIUM
    context: Dict[str, Any] = field(default_factory=dict)
    timestamp: float = field(default_factory=time.time)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert validation error to dictionary for serialization."""
        return {
            'message': self.message,
            'data': str(self.data)[:100],  # Truncate for readability
            'field': self.field,
            'code': self.code,
            'severity': self.severity.name,
            'context': self.context,
            'timestamp': self.timestamp,
            'time': datetime.fromtimestamp(self.timestamp).isoformat()
        }

@dataclass
class ValidationReport:
    """Report of validation errors for a batch of data."""
    valid_count: int = 0
    invalid_count: int = 0
    error_count: int = 0
    warning_count: int = 0
    errors: List[ValidationError] = field(default_factory=list)
    warnings: List[ValidationError] = field(default_factory=list)
    by_field: Dict[str, List[ValidationError]] = field(default_factory=dict)
    by_code: Dict[str, List[ValidationError]] = field(default_factory=dict)
    by_severity: Dict[str, List[ValidationError]] = field(default_factory=dict)
    
    def add_error(self, error: ValidationError) -> None:
        """Add a validation error to the report.
        
        Args:
            error: ValidationError to add
        """
        self.errors.append(error)
        self.error_count += 1
        
        # Index by field
        if error.field:
            if error.field not in self.by_field:
                self.by_field[error.field] = []
            self.by_field[error.field].append(error)
        
        # Index by code
        if error.code:
            if error.code not in self.by_code:
                self.by_code[error.code] = []
            self.by_code[error.code].append(error)
        
        # Index by severity
        severity_name = error.severity.name
        if severity_name not in self.by_severity:
            self.by_severity[severity_name] = []
        self.by_severity[severity_name].append(error)
    
    def add_warning(self, warning: ValidationError) -> None:
        """Add a validation warning to the report.
        
        Args:
            warning: ValidationError to add as warning
        """
        self.warnings.append(warning)
        self.warning_count += 1
        
        # Warnings are also indexed like errors
        if warning.field:
            if warning.field not in self.by_field:
                self.by_field[warning.field] = []
            self.by_field[warning.field].append(warning)
        
        if warning.code:
            if warning.code not in self.by_code:
                self.by_code[warning.code] = []
            self.by_code[warning.code].append(warning)
        
        severity_name = warning.severity.name
        if severity_name not in self.by_severity:
            self.by_severity[severity_name] = []
        self.by_severity[severity_name].append(warning)
    
    def mark_valid(self) -> None:
        """Increment the valid item count."""
        self.valid_count += 1
    
    def mark_invalid(self) -> None:
        """Increment the invalid item count."""
        self.invalid_count += 1
    
    def is_empty(self) -> bool:
        """Check if the report contains no errors or warnings."""
        return self.error_count == 0 and self.warning_count == 0
    
    def get_summary(self) -> Dict[str, Any]:
        """Get a summary of the validation report.
        
        Returns:
            Dictionary with validation summary statistics
        """
        return {
            'total_items': self.valid_count + self.invalid_count,
            'valid_items': self.valid_count,
            'invalid_items': self.invalid_count,
            'error_count': self.error_count,
            'warning_count': self.warning_count,
            'fields_with_errors': list(self.by_field.keys()),
            'error_codes': list(self.by_code.keys()),
            'error_counts_by_severity': {
                severity: len(errors)
                for severity, errors in self.by_severity.items()
            }
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert validation report to dictionary for serialization.
        
        Returns:
            Dictionary representation of the validation report
        """
        return {
            'valid_count': self.valid_count,
            'invalid_count': self.invalid_count,
            'error_count': self.error_count,
            'warning_count': self.warning_count,
            'errors': [e.to_dict() for e in self.errors],
            'warnings': [w.to_dict() for w in self.warnings],
            'by_field': {
                field: [e.to_dict() for e in errors]
                for field, errors in self.by_field.items()
            },
            'by_code': {
                code: [e.to_dict() for e in errors]
                for code, errors in self.by_code.items()
            },
            'by_severity': {
                severity: [e.to_dict() for e in errors]
                for severity, errors in self.by_severity.items()
            },
            'summary': self.get_summary()
        }
    
    def save_to_file(self, file_path: str) -> None:
        """Save the validation report to a JSON file.
        
        Args:
            file_path: Path to save the report
        """
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        
        with open(file_path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
    
    @classmethod
    def load_from_file(cls, file_path: str) -> 'ValidationReport':
        """Load a validation report from a JSON file.
        
        Args:
            file_path: Path to load the report from
            
        Returns:
            ValidationReport loaded from file
        """
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        report = cls(
            valid_count=data['valid_count'],
            invalid_count=data['invalid_count'],
            error_count=data['error_count'],
            warning_count=data['warning_count']
        )
        
        # Reconstruct errors
        for error_dict in data['errors']:
            error = ValidationError(
                message=error_dict['message'],
                data=error_dict['data'],
                field=error_dict['field'],
                code=error_dict['code'],
                severity=ErrorSeverity[error_dict['severity']],
                context=error_dict['context'],
                timestamp=error_dict['timestamp']
            )
            report.errors.append(error)
        
        # Reconstruct warnings
        for warning_dict in data['warnings']:
            warning = ValidationError(
                message=warning_dict['message'],
                data=warning_dict['data'],
                field=warning_dict['field'],
                code=warning_dict['code'],
                severity=ErrorSeverity[warning_dict['severity']],
                context=warning_dict['context'],
                timestamp=warning_dict['timestamp']
            )
            report.warnings.append(warning)
        
        # Indexes are rebuilt by the add_error and add_warning methods
        for error in report.errors:
            if error.field:
                if error.field not in report.by_field:
                    report.by_field[error.field] = []
                report.by_field[error.field].append(error)
            
            if error.code:
                if error.code not in report.by_code:
                    report.by_code[error.code] = []
                report.by_code[error.code].append(error)
            
            severity_name = error.severity.name
            if severity_name not in report.by_severity:
                report.by_severity[severity_name] = []
            report.by_severity[severity_name].append(error)
        
        for warning in report.warnings:
            if warning.field:
                if warning.field not in report.by_field:
                    report.by_field[warning.field] = []
                report.by_field[warning.field].append(warning)
            
            if warning.code:
                if warning.code not in report.by_code:
                    report.by_code[warning.code] = []
                report.by_code[warning.code].append(warning)
            
            severity_name = warning.severity.name
            if severity_name not in report.by_severity:
                report.by_severity[severity_name] = []
            report.by_severity[severity_name].append(warning)
        
        return report
    
    def merge(self, other: 'ValidationReport') -> None:
        """Merge another validation report into this one.
        
        Args:
            other: ValidationReport to merge
        """
        self.valid_count += other.valid_count
        self.invalid_count += other.invalid_count
        self.error_count += other.error_count
        self.warning_count += other.warning_count
        
        # Merge errors
        self.errors.extend(other.errors)
        
        # Merge warnings
        self.warnings.extend(other.warnings)
        
        # Merge indexes
        for field, errors in other.by_field.items():
            if field not in self.by_field:
                self.by_field[field] = []
            self.by_field[field].extend(errors)
        
        for code, errors in other.by_code.items():
            if code not in self.by_code:
                self.by_code[code] = []
            self.by_code[code].extend(errors)
        
        for severity, errors in other.by_severity.items():
            if severity not in self.by_severity:
                self.by_severity[severity] = []
            self.by_severity[severity].extend(errors)


class ValidationErrorHandler:
    """Handles validation errors with configurable strategies."""
    
    def __init__(self, logger=None):
        """Initialize the validation error handler.
        
        Args:
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        self.current_report = ValidationReport()
    
    def handle_validation(
        self, 
        validation_func: Callable[[T], U], 
        data: T,
        component: str,
        operation: str,
        error_strategy: RecoveryStrategy = RecoveryStrategy.SKIP,
        default_value: Optional[U] = None,
        fallback_func: Optional[Callable[[T], U]] = None,
        record_errors: bool = True
    ) -> Union[U, None]:
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
            Validated data, default value, or result of fallback function
        """
        try:
            result = validation_func(data)
            
            if record_errors:
                self.current_report.mark_valid()
                
            return result
            
        except Exception as e:
            context = ErrorContext(
                component=component,
                operation=operation,
                data={"validation_data": str(data)[:100]}
            )
            
            if record_errors:
                validation_error = ValidationError(
                    message=str(e),
                    data=data,
                    code=f"{component}.{operation}",
                    severity=ErrorSeverity.MEDIUM,
                    context={
                        "component": component,
                        "operation": operation,
                        "error_type": type(e).__name__
                    }
                )
                
                self.current_report.add_error(validation_error)
                self.current_report.mark_invalid()
            
            self.logger.warning(
                f"Validation error in {component}:{operation}: "
                f"{type(e).__name__}: {str(e)}. "
                f"Using strategy: {error_strategy.name}"
            )
            
            if error_strategy == RecoveryStrategy.SKIP:
                return default_value
                
            elif error_strategy == RecoveryStrategy.FALLBACK and fallback_func:
                try:
                    return fallback_func(data)
                except Exception as fallback_e:
                    self.logger.error(
                        f"Fallback function failed in {component}:{operation}: {str(fallback_e)}"
                    )
                    return default_value
                    
            elif error_strategy == RecoveryStrategy.LOG_ONLY:
                return data
                
            else:
                raise
    
    def validate_batch(
        self, 
        validation_func: Callable[[T], U], 
        items: List[T],
        component: str,
        operation: str,
        error_strategy: RecoveryStrategy = RecoveryStrategy.SKIP,
        default_value: Optional[U] = None,
        fallback_func: Optional[Callable[[T], U]] = None,
        parallel: bool = False,
        max_workers: Optional[int] = None,
        continue_on_error: bool = True
    ) -> Tuple[List[Optional[U]], ValidationReport]:
        """Validate a batch of items with error handling.
        
        Args:
            validation_func: Function to validate each item
            items: List of items to validate
            component: Component name for context
            operation: Operation name for context
            error_strategy: Strategy for handling validation errors
            default_value: Default value for invalid items with SKIP strategy
            fallback_func: Function for invalid items with FALLBACK strategy
            parallel: Whether to process in parallel
            max_workers: Number of parallel workers
            continue_on_error: Whether to continue on error or abort batch
            
        Returns:
            Tuple of (list of validated items, validation report)
        """
        # Create a fresh report for this batch
        report = ValidationReport()
        results = []
        
        if parallel:
            from concurrent.futures import ThreadPoolExecutor
            
            def validate_item(item):
                try:
                    result = validation_func(item)
                    report.mark_valid()
                    return result, None
                except Exception as e:
                    report.mark_invalid()
                    validation_error = ValidationError(
                        message=str(e),
                        data=item,
                        code=f"{component}.{operation}",
                        severity=ErrorSeverity.MEDIUM,
                        context={
                            "component": component,
                            "operation": operation,
                            "error_type": type(e).__name__
                        }
                    )
                    report.add_error(validation_error)
                    
                    if error_strategy == RecoveryStrategy.SKIP:
                        return default_value, None
                    elif error_strategy == RecoveryStrategy.FALLBACK and fallback_func:
                        try:
                            return fallback_func(item), None
                        except Exception as fallback_e:
                            self.logger.error(
                                f"Fallback function failed in {component}:{operation}: {str(fallback_e)}"
                            )
                            return default_value, None
                    elif error_strategy == RecoveryStrategy.LOG_ONLY:
                        return item, None
                    else:
                        return None, e
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = []
                for item in items:
                    futures.append(executor.submit(validate_item, item))
                
                for future in futures:
                    result, error = future.result()
                    if error and not continue_on_error:
                        raise error
                    results.append(result)
        else:
            for item in items:
                try:
                    result = validation_func(item)
                    report.mark_valid()
                    results.append(result)
                except Exception as e:
                    report.mark_invalid()
                    validation_error = ValidationError(
                        message=str(e),
                        data=item,
                        code=f"{component}.{operation}",
                        severity=ErrorSeverity.MEDIUM,
                        context={
                            "component": component,
                            "operation": operation,
                            "error_type": type(e).__name__
                        }
                    )
                    report.add_error(validation_error)
                    
                    if not continue_on_error:
                        raise
                    
                    if error_strategy == RecoveryStrategy.SKIP:
                        results.append(default_value)
                    elif error_strategy == RecoveryStrategy.FALLBACK and fallback_func:
                        try:
                            results.append(fallback_func(item))
                        except Exception as fallback_e:
                            self.logger.error(
                                f"Fallback function failed in {component}:{operation}: {str(fallback_e)}"
                            )
                            results.append(default_value)
                    elif error_strategy == RecoveryStrategy.LOG_ONLY:
                        results.append(item)
                    else:
                        results.append(None)
        
        return results, report
    
    def reset_report(self) -> None:
        """Reset the current validation report."""
        self.current_report = ValidationReport()
    
    def get_report(self) -> ValidationReport:
        """Get the current validation report.
        
        Returns:
            Current ValidationReport
        """
        return self.current_report
    
    def save_report(self, file_path: str) -> None:
        """Save the current validation report to a file.
        
        Args:
            file_path: Path to save the report
        """
        self.current_report.save_to_file(file_path)


class ChemicalValidationRules:
    """Common validation rules for chemical data."""
    
    @staticmethod
    def validate_smiles(smiles: str) -> str:
        """Validate a SMILES string.
        
        Args:
            smiles: SMILES string to validate
            
        Returns:
            Validated SMILES string
            
        Raises:
            ValueError: If SMILES string is invalid
        """
        if not smiles or not isinstance(smiles, str):
            raise ValueError("SMILES string is required")
        
        # Basic format check
        if not set(smiles).intersection(set('CNOPSFIBr[]()-=#')):
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        return smiles
    
    @staticmethod
    def validate_inchi(inchi: str) -> str:
        """Validate an InChI string.
        
        Args:
            inchi: InChI string to validate
            
        Returns:
            Validated InChI string
            
        Raises:
            ValueError: If InChI string is invalid
        """
        if not inchi or not isinstance(inchi, str):
            raise ValueError("InChI string is required")
        
        # Check prefix
        if not inchi.startswith("InChI="):
            raise ValueError(f"Invalid InChI string (must start with 'InChI='): {inchi}")
        
        return inchi
    
    @staticmethod
    def validate_molecule_name(name: str) -> str:
        """Validate a molecule name.
        
        Args:
            name: Molecule name to validate
            
        Returns:
            Validated molecule name
            
        Raises:
            ValueError: If name is invalid
        """
        if not name or not isinstance(name, str):
            raise ValueError("Molecule name is required")
        
        if len(name) < 2:
            raise ValueError(f"Molecule name too short: {name}")
        
        return name
    
    @staticmethod
    def validate_cid(cid: Union[str, int]) -> str:
        """Validate a PubChem CID.
        
        Args:
            cid: PubChem CID to validate
            
        Returns:
            Validated CID as string
            
        Raises:
            ValueError: If CID is invalid
        """
        if isinstance(cid, int):
            cid = str(cid)
            
        if not cid or not isinstance(cid, str):
            raise ValueError("CID is required")
        
        if not cid.isdigit():
            raise ValueError(f"Invalid CID (must be numeric): {cid}")
        
        return cid
    
    @staticmethod
    def validate_chembl_id(chembl_id: str) -> str:
        """Validate a ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL ID to validate
            
        Returns:
            Validated ChEMBL ID
            
        Raises:
            ValueError: If ChEMBL ID is invalid
        """
        if not chembl_id or not isinstance(chembl_id, str):
            raise ValueError("ChEMBL ID is required")
        
        if not chembl_id.startswith("CHEMBL"):
            raise ValueError(f"Invalid ChEMBL ID (must start with 'CHEMBL'): {chembl_id}")
        
        # Check that the rest is numeric
        if not chembl_id[6:].isdigit():
            raise ValueError(f"Invalid ChEMBL ID (must be CHEMBL followed by numbers): {chembl_id}")
        
        return chembl_id
    
    @staticmethod
    def validate_molecular_formula(formula: str) -> str:
        """Validate a molecular formula.
        
        Args:
            formula: Molecular formula to validate
            
        Returns:
            Validated molecular formula
            
        Raises:
            ValueError: If formula is invalid
        """
        if not formula or not isinstance(formula, str):
            raise ValueError("Molecular formula is required")
        
        # Basic format check: must contain at least one element symbol
        if not any(c.isupper() for c in formula):
            raise ValueError(f"Invalid molecular formula (must contain at least one element symbol): {formula}")
        
        # TODO: Add more sophisticated formula validation if needed
        
        return formula
    
    @staticmethod
    def validate_molecular_weight(weight: Union[float, str]) -> float:
        """Validate a molecular weight.
        
        Args:
            weight: Molecular weight to validate
            
        Returns:
            Validated molecular weight as float
            
        Raises:
            ValueError: If weight is invalid
        """
        if isinstance(weight, str):
            try:
                weight = float(weight)
            except ValueError:
                raise ValueError(f"Invalid molecular weight (must be numeric): {weight}")
        
        if not isinstance(weight, (int, float)):
            raise ValueError(f"Invalid molecular weight (must be numeric): {weight}")
        
        if weight <= 0 or weight > 5000:
            raise ValueError(f"Molecular weight out of realistic range (0-5000): {weight}")
        
        return float(weight)
    
    @staticmethod
    def validate_logp(logp: Union[float, str]) -> float:
        """Validate a logP value.
        
        Args:
            logp: LogP value to validate
            
        Returns:
            Validated logP as float
            
        Raises:
            ValueError: If logP is invalid
        """
        if isinstance(logp, str):
            try:
                logp = float(logp)
            except ValueError:
                raise ValueError(f"Invalid logP (must be numeric): {logp}")
        
        if not isinstance(logp, (int, float)):
            raise ValueError(f"Invalid logP (must be numeric): {logp}")
        
        if logp < -10 or logp > 20:
            raise ValueError(f"LogP out of realistic range (-10 to 20): {logp}")
        
        return float(logp)