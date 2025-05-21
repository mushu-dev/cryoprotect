#!/usr/bin/env python3
"""
Validation Service for CryoProtect Enhanced Experimental Data System.

This module provides services for validating experimental data against
scientific standards, statistical criteria, and business rules.
"""

from typing import Dict, List, Any, Optional, Union, Tuple, Callable
import uuid
from datetime import datetime
import json
import numpy as np
from scipy import stats

from ..models import (
    Experiment,
    ExperimentResult,
    Protocol,
    TissueType,
    ValidationError
)

class RuleViolation:
    """Class representing a validation rule violation."""
    
    def __init__(
        self,
        rule_id: str,
        severity: str,
        message: str,
        entity_id: Optional[str] = None,
        entity_type: Optional[str] = None,
        field_name: Optional[str] = None,
        value: Optional[Any] = None,
        expected_value: Optional[Any] = None,
        suggestion: Optional[str] = None
    ):
        """
        Initialize a rule violation.
        
        Args:
            rule_id: ID of the violated rule
            severity: Severity level ('info', 'warning', 'error', 'critical')
            message: Description of the violation
            entity_id: ID of the entity that violated the rule
            entity_type: Type of entity ('experiment', 'result', etc.)
            field_name: Name of the field that violated the rule
            value: Actual value that violated the rule
            expected_value: Expected value according to the rule
            suggestion: Suggestion for fixing the violation
        """
        self.rule_id = rule_id
        self.severity = severity
        self.message = message
        self.entity_id = entity_id
        self.entity_type = entity_type
        self.field_name = field_name
        self.value = value
        self.expected_value = expected_value
        self.suggestion = suggestion
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "rule_id": self.rule_id,
            "severity": self.severity,
            "message": self.message,
            "entity_id": self.entity_id,
            "entity_type": self.entity_type,
            "field_name": self.field_name,
            "value": self.value,
            "expected_value": self.expected_value,
            "suggestion": self.suggestion
        }

class ValidationRule:
    """Base class for validation rules."""
    
    def __init__(
        self,
        rule_id: str,
        description: str,
        severity: str = "warning",
        entity_type: str = "any"
    ):
        """
        Initialize a validation rule.
        
        Args:
            rule_id: Unique identifier for the rule
            description: Description of what the rule checks
            severity: Default severity level for violations
            entity_type: Type of entity this rule applies to
        """
        self.rule_id = rule_id
        self.description = description
        self.severity = severity
        self.entity_type = entity_type
    
    async def validate(self, entity: Any) -> List[RuleViolation]:
        """
        Validate an entity against this rule.
        
        Args:
            entity: Entity to validate
            
        Returns:
            List of rule violations (empty if validation passes)
        """
        # Base implementation - override in subclasses
        return []
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "rule_id": self.rule_id,
            "description": self.description,
            "severity": self.severity,
            "entity_type": self.entity_type
        }

class RangeRule(ValidationRule):
    """Rule that validates a numeric value is within a specified range."""
    
    def __init__(
        self,
        rule_id: str,
        description: str,
        field_name: str,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None,
        **kwargs
    ):
        """
        Initialize a range rule.
        
        Args:
            rule_id: Unique identifier for the rule
            description: Description of what the rule checks
            field_name: Name of the field to validate
            min_value: Minimum acceptable value (None for no lower bound)
            max_value: Maximum acceptable value (None for no upper bound)
            **kwargs: Additional arguments for ValidationRule
        """
        super().__init__(rule_id, description, **kwargs)
        self.field_name = field_name
        self.min_value = min_value
        self.max_value = max_value
    
    async def validate(self, entity: Any) -> List[RuleViolation]:
        """
        Validate an entity against this range rule.
        
        Args:
            entity: Entity to validate
            
        Returns:
            List of rule violations (empty if validation passes)
        """
        violations = []
        
        # Get field value
        if not hasattr(entity, self.field_name):
            return []
        
        value = getattr(entity, self.field_name)
        
        # Skip if value is None
        if value is None:
            return []
        
        # Validate against range
        if self.min_value is not None and value < self.min_value:
            violations.append(RuleViolation(
                rule_id=self.rule_id,
                severity=self.severity,
                message=f"{self.field_name} is below minimum value",
                entity_id=getattr(entity, 'id', None),
                entity_type=self.entity_type,
                field_name=self.field_name,
                value=value,
                expected_value=f">= {self.min_value}",
                suggestion=f"Increase {self.field_name} to at least {self.min_value}"
            ))
        
        if self.max_value is not None and value > self.max_value:
            violations.append(RuleViolation(
                rule_id=self.rule_id,
                severity=self.severity,
                message=f"{self.field_name} is above maximum value",
                entity_id=getattr(entity, 'id', None),
                entity_type=self.entity_type,
                field_name=self.field_name,
                value=value,
                expected_value=f"<= {self.max_value}",
                suggestion=f"Decrease {self.field_name} to at most {self.max_value}"
            ))
        
        return violations
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "field_name": self.field_name,
            "min_value": self.min_value,
            "max_value": self.max_value
        })
        return result

class RequiredFieldRule(ValidationRule):
    """Rule that validates a required field is present and non-empty."""
    
    def __init__(
        self,
        rule_id: str,
        description: str,
        field_name: str,
        **kwargs
    ):
        """
        Initialize a required field rule.
        
        Args:
            rule_id: Unique identifier for the rule
            description: Description of what the rule checks
            field_name: Name of the field to validate
            **kwargs: Additional arguments for ValidationRule
        """
        super().__init__(rule_id, description, **kwargs)
        self.field_name = field_name
    
    async def validate(self, entity: Any) -> List[RuleViolation]:
        """
        Validate an entity against this required field rule.
        
        Args:
            entity: Entity to validate
            
        Returns:
            List of rule violations (empty if validation passes)
        """
        violations = []
        
        # Get field value
        if not hasattr(entity, self.field_name):
            violations.append(RuleViolation(
                rule_id=self.rule_id,
                severity=self.severity,
                message=f"Required field {self.field_name} is missing",
                entity_id=getattr(entity, 'id', None),
                entity_type=self.entity_type,
                field_name=self.field_name,
                suggestion=f"Add {self.field_name} field"
            ))
            return violations
        
        value = getattr(entity, self.field_name)
        
        # Check if value is empty
        if value is None:
            violations.append(RuleViolation(
                rule_id=self.rule_id,
                severity=self.severity,
                message=f"Required field {self.field_name} is null",
                entity_id=getattr(entity, 'id', None),
                entity_type=self.entity_type,
                field_name=self.field_name,
                suggestion=f"Provide a value for {self.field_name}"
            ))
        elif isinstance(value, str) and not value.strip():
            violations.append(RuleViolation(
                rule_id=self.rule_id,
                severity=self.severity,
                message=f"Required field {self.field_name} is empty",
                entity_id=getattr(entity, 'id', None),
                entity_type=self.entity_type,
                field_name=self.field_name,
                suggestion=f"Provide a non-empty value for {self.field_name}"
            ))
        elif isinstance(value, (list, dict)) and not value:
            violations.append(RuleViolation(
                rule_id=self.rule_id,
                severity=self.severity,
                message=f"Required field {self.field_name} is empty",
                entity_id=getattr(entity, 'id', None),
                entity_type=self.entity_type,
                field_name=self.field_name,
                suggestion=f"Provide non-empty {type(value).__name__} for {self.field_name}"
            ))
        
        return violations
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "field_name": self.field_name
        })
        return result

class FormatRule(ValidationRule):
    """Rule that validates a field matches a specific format."""
    
    def __init__(
        self,
        rule_id: str,
        description: str,
        field_name: str,
        pattern: Optional[str] = None,
        custom_validator: Optional[Callable[[Any], bool]] = None,
        **kwargs
    ):
        """
        Initialize a format rule.
        
        Args:
            rule_id: Unique identifier for the rule
            description: Description of what the rule checks
            field_name: Name of the field to validate
            pattern: Regular expression pattern to match
            custom_validator: Custom validation function
            **kwargs: Additional arguments for ValidationRule
        """
        super().__init__(rule_id, description, **kwargs)
        self.field_name = field_name
        self.pattern = pattern
        self.custom_validator = custom_validator
        
        if not pattern and not custom_validator:
            raise ValueError("Either pattern or custom_validator must be provided")
    
    async def validate(self, entity: Any) -> List[RuleViolation]:
        """
        Validate an entity against this format rule.
        
        Args:
            entity: Entity to validate
            
        Returns:
            List of rule violations (empty if validation passes)
        """
        violations = []
        
        # Get field value
        if not hasattr(entity, self.field_name):
            return []
        
        value = getattr(entity, self.field_name)
        
        # Skip if value is None
        if value is None:
            return []
        
        # Check pattern if provided
        if self.pattern:
            import re
            pattern = re.compile(self.pattern)
            
            if isinstance(value, str) and not pattern.match(value):
                violations.append(RuleViolation(
                    rule_id=self.rule_id,
                    severity=self.severity,
                    message=f"Field {self.field_name} does not match required format",
                    entity_id=getattr(entity, 'id', None),
                    entity_type=self.entity_type,
                    field_name=self.field_name,
                    value=value,
                    expected_value=f"Pattern: {self.pattern}",
                    suggestion=f"Format {self.field_name} according to the required pattern"
                ))
        
        # Use custom validator if provided
        if self.custom_validator and not self.custom_validator(value):
            violations.append(RuleViolation(
                rule_id=self.rule_id,
                severity=self.severity,
                message=f"Field {self.field_name} does not meet format requirements",
                entity_id=getattr(entity, 'id', None),
                entity_type=self.entity_type,
                field_name=self.field_name,
                value=value,
                suggestion=f"Ensure {self.field_name} meets the format requirements: {self.description}"
            ))
        
        return violations
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "field_name": self.field_name,
            "pattern": self.pattern,
            "has_custom_validator": self.custom_validator is not None
        })
        return result

class StatisticalOutlierRule(ValidationRule):
    """Rule that validates a value is not a statistical outlier."""
    
    def __init__(
        self,
        rule_id: str,
        description: str,
        field_name: str,
        z_score_threshold: float = 3.0,
        min_sample_size: int = 5,
        **kwargs
    ):
        """
        Initialize a statistical outlier rule.
        
        Args:
            rule_id: Unique identifier for the rule
            description: Description of what the rule checks
            field_name: Name of the field to validate
            z_score_threshold: Z-score threshold for outlier detection
            min_sample_size: Minimum sample size for outlier detection
            **kwargs: Additional arguments for ValidationRule
        """
        super().__init__(rule_id, description, **kwargs)
        self.field_name = field_name
        self.z_score_threshold = z_score_threshold
        self.min_sample_size = min_sample_size
    
    async def validate(self, entity: Any, reference_data: List[Any]) -> List[RuleViolation]:
        """
        Validate an entity against this statistical outlier rule.
        
        Args:
            entity: Entity to validate
            reference_data: Reference data points for comparison
            
        Returns:
            List of rule violations (empty if validation passes)
        """
        violations = []
        
        # Get field value
        if not hasattr(entity, self.field_name):
            return []
        
        value = getattr(entity, self.field_name)
        
        # Skip if value is None
        if value is None:
            return []
        
        # Check if we have enough reference data
        if len(reference_data) < self.min_sample_size:
            return []
        
        # Extract reference values
        reference_values = []
        for ref in reference_data:
            if hasattr(ref, self.field_name):
                ref_value = getattr(ref, self.field_name)
                if ref_value is not None:
                    reference_values.append(ref_value)
        
        if len(reference_values) < self.min_sample_size:
            return []
        
        # Calculate z-score
        mean = np.mean(reference_values)
        std = np.std(reference_values, ddof=1)
        
        if std == 0:
            return []
        
        z_score = abs((value - mean) / std)
        
        # Check if outlier
        if z_score > self.z_score_threshold:
            violations.append(RuleViolation(
                rule_id=self.rule_id,
                severity=self.severity,
                message=f"Field {self.field_name} is a statistical outlier (z-score: {z_score:.2f})",
                entity_id=getattr(entity, 'id', None),
                entity_type=self.entity_type,
                field_name=self.field_name,
                value=value,
                expected_value=f"Within {self.z_score_threshold} standard deviations of mean ({mean:.2f} ± {self.z_score_threshold * std:.2f})",
                suggestion="Verify the measurement accuracy or document the exceptional condition"
            ))
        
        return violations
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "field_name": self.field_name,
            "z_score_threshold": self.z_score_threshold,
            "min_sample_size": self.min_sample_size
        })
        return result

class ConsistencyRule(ValidationRule):
    """Rule that validates consistency between related fields."""
    
    def __init__(
        self,
        rule_id: str,
        description: str,
        field_names: List[str],
        consistency_check: Callable[[List[Any]], bool],
        error_message: str,
        suggestion: str,
        **kwargs
    ):
        """
        Initialize a consistency rule.
        
        Args:
            rule_id: Unique identifier for the rule
            description: Description of what the rule checks
            field_names: Names of the fields to check for consistency
            consistency_check: Function that takes field values and returns True if consistent
            error_message: Message to display when consistency check fails
            suggestion: Suggestion for fixing inconsistency
            **kwargs: Additional arguments for ValidationRule
        """
        super().__init__(rule_id, description, **kwargs)
        self.field_names = field_names
        self.consistency_check = consistency_check
        self.error_message = error_message
        self.suggestion = suggestion
    
    async def validate(self, entity: Any) -> List[RuleViolation]:
        """
        Validate an entity against this consistency rule.
        
        Args:
            entity: Entity to validate
            
        Returns:
            List of rule violations (empty if validation passes)
        """
        violations = []
        
        # Get field values
        field_values = []
        for field_name in self.field_names:
            if not hasattr(entity, field_name):
                return []
            
            field_values.append(getattr(entity, field_name))
        
        # Check consistency
        if not self.consistency_check(field_values):
            violations.append(RuleViolation(
                rule_id=self.rule_id,
                severity=self.severity,
                message=self.error_message,
                entity_id=getattr(entity, 'id', None),
                entity_type=self.entity_type,
                field_name=", ".join(self.field_names),
                value=field_values,
                suggestion=self.suggestion
            ))
        
        return violations
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "field_names": self.field_names,
            "error_message": self.error_message,
            "suggestion": self.suggestion
        })
        return result

class ValidationService:
    """Service for validating experimental data."""
    
    def __init__(self, db_adapter):
        """
        Initialize validation service.
        
        Args:
            db_adapter: Database adapter for persistence
        """
        self.db_adapter = db_adapter
        self.rules = self._initialize_rules()
    
    def _initialize_rules(self) -> Dict[str, List[ValidationRule]]:
        """
        Initialize validation rules.
        
        Returns:
            Dictionary of rules by entity type
        """
        rules = {
            "experiment": [],
            "experiment_result": [],
            "protocol": [],
            "protocol_step": [],
            "tissue_type": [],
            "time_series": [],
            "time_series_data": []
        }
        
        # Experiment rules
        rules["experiment"].extend([
            RequiredFieldRule(
                rule_id="exp-req-title",
                description="Experiment must have a title",
                field_name="title",
                severity="error",
                entity_type="experiment"
            ),
            RequiredFieldRule(
                rule_id="exp-req-type",
                description="Experiment must have an experiment type",
                field_name="experiment_type_id",
                severity="error",
                entity_type="experiment"
            ),
            RangeRule(
                rule_id="exp-temp-range",
                description="Temperature should be within -200°C to 100°C",
                field_name="temperature",
                min_value=-200,
                max_value=100,
                severity="warning",
                entity_type="experiment"
            ),
            RangeRule(
                rule_id="exp-cooling-rate",
                description="Cooling rate should be negative",
                field_name="cooling_rate",
                max_value=0,
                severity="warning",
                entity_type="experiment"
            ),
            RangeRule(
                rule_id="exp-thawing-rate",
                description="Thawing rate should be positive",
                field_name="thawing_rate",
                min_value=0,
                severity="warning",
                entity_type="experiment"
            )
        ])
        
        # Experiment result rules
        rules["experiment_result"].extend([
            RequiredFieldRule(
                rule_id="result-req-exp",
                description="Result must be associated with an experiment",
                field_name="experiment_id",
                severity="error",
                entity_type="experiment_result"
            ),
            RequiredFieldRule(
                rule_id="result-req-tissue",
                description="Result must specify a tissue type",
                field_name="tissue_type_id",
                severity="error",
                entity_type="experiment_result"
            ),
            RangeRule(
                rule_id="result-viability-range",
                description="Viability percentage should be between 0 and 100",
                field_name="viability_percentage",
                min_value=0,
                max_value=100,
                severity="error",
                entity_type="experiment_result"
            ),
            RangeRule(
                rule_id="result-recovery-range",
                description="Recovery rate should be non-negative",
                field_name="recovery_rate",
                min_value=0,
                severity="error",
                entity_type="experiment_result"
            ),
            ConsistencyRule(
                rule_id="result-molecule-mixture",
                description="Either molecule_id or mixture_id must be set, but not both",
                field_names=["molecule_id", "mixture_id"],
                consistency_check=lambda values: (values[0] is None) != (values[1] is None),
                error_message="Either molecule_id or mixture_id must be set, but not both",
                suggestion="Specify either a molecule or a mixture, not both",
                severity="error",
                entity_type="experiment_result"
            ),
            ConsistencyRule(
                rule_id="result-concentration-unit",
                description="If concentration is provided, concentration_unit is required",
                field_names=["concentration", "concentration_unit"],
                consistency_check=lambda values: values[0] is None or values[1] is not None,
                error_message="Concentration unit is required when concentration is provided",
                suggestion="Specify a concentration unit (e.g., 'M', '%w/v')",
                severity="error",
                entity_type="experiment_result"
            )
        ])
        
        # Protocol rules
        rules["protocol"].extend([
            RequiredFieldRule(
                rule_id="protocol-req-name",
                description="Protocol must have a name",
                field_name="name",
                severity="error",
                entity_type="protocol"
            ),
            RangeRule(
                rule_id="protocol-version",
                description="Protocol version must be positive",
                field_name="version",
                min_value=1,
                severity="error",
                entity_type="protocol"
            )
        ])
        
        # Protocol step rules
        rules["protocol_step"].extend([
            RequiredFieldRule(
                rule_id="step-req-name",
                description="Protocol step must have a name",
                field_name="name",
                severity="error",
                entity_type="protocol_step"
            ),
            FormatRule(
                rule_id="step-type-format",
                description="Step type must be one of: manual, automated, measurement, calculation",
                field_name="step_type",
                custom_validator=lambda v: v in ["manual", "automated", "measurement", "calculation"],
                severity="error",
                entity_type="protocol_step"
            )
        ])
        
        # Time series rules
        rules["time_series"].extend([
            RequiredFieldRule(
                rule_id="timeseries-req-name",
                description="Time series must have a name",
                field_name="name",
                severity="error",
                entity_type="time_series"
            ),
            RequiredFieldRule(
                rule_id="timeseries-req-property",
                description="Time series must specify a property type",
                field_name="property_type",
                severity="error",
                entity_type="time_series"
            ),
            RequiredFieldRule(
                rule_id="timeseries-req-exp",
                description="Time series must be associated with an experiment",
                field_name="experiment_id",
                severity="error",
                entity_type="time_series"
            )
        ])
        
        # Time series data rules
        rules["time_series_data"].extend([
            RequiredFieldRule(
                rule_id="timepoint-req-timestamp",
                description="Time point must have a timestamp",
                field_name="timestamp",
                severity="error",
                entity_type="time_series_data"
            ),
            RequiredFieldRule(
                rule_id="timepoint-req-value",
                description="Time point must have a value",
                field_name="value",
                severity="error",
                entity_type="time_series_data"
            ),
            RequiredFieldRule(
                rule_id="timepoint-req-timeseries",
                description="Time point must be associated with a time series",
                field_name="time_series_id",
                severity="error",
                entity_type="time_series_data"
            )
        ])
        
        return rules
    
    async def register_rule(self, rule: ValidationRule) -> None:
        """
        Register a custom validation rule.
        
        Args:
            rule: Validation rule to register
        """
        entity_type = rule.entity_type
        
        if entity_type == "any":
            # Add to all entity types
            for rules_list in self.rules.values():
                rules_list.append(rule)
        elif entity_type in self.rules:
            self.rules[entity_type].append(rule)
        else:
            # Create new entity type
            self.rules[entity_type] = [rule]
    
    async def validate_experiment(self, experiment: Experiment) -> Dict[str, Any]:
        """
        Validate an experiment against applicable rules.
        
        Args:
            experiment: Experiment to validate
            
        Returns:
            Validation results
        """
        violations = []
        
        # Apply experiment rules
        for rule in self.rules["experiment"]:
            rule_violations = await rule.validate(experiment)
            violations.extend(rule_violations)
        
        # Apply generic rules
        if "any" in self.rules:
            for rule in self.rules["any"]:
                rule_violations = await rule.validate(experiment)
                violations.extend(rule_violations)
        
        # Prepare validation report
        validation_report = {
            "entity_id": experiment.id,
            "entity_type": "experiment",
            "timestamp": datetime.now().isoformat(),
            "is_valid": len([v for v in violations if v.severity in ["error", "critical"]]) == 0,
            "violations": [v.to_dict() for v in violations],
            "warning_count": len([v for v in violations if v.severity == "warning"]),
            "error_count": len([v for v in violations if v.severity == "error"]),
            "critical_count": len([v for v in violations if v.severity == "critical"])
        }
        
        return validation_report
    
    async def validate_experiment_result(
        self, 
        result: ExperimentResult,
        reference_results: Optional[List[ExperimentResult]] = None
    ) -> Dict[str, Any]:
        """
        Validate an experiment result against applicable rules.
        
        Args:
            result: Result to validate
            reference_results: Reference results for statistical validation
            
        Returns:
            Validation results
        """
        violations = []
        
        # Apply experiment result rules
        for rule in self.rules["experiment_result"]:
            rule_violations = await rule.validate(result)
            violations.extend(rule_violations)
        
        # Apply generic rules
        if "any" in self.rules:
            for rule in self.rules["any"]:
                rule_violations = await rule.validate(result)
                violations.extend(rule_violations)
        
        # Apply statistical rules if reference data is provided
        if reference_results:
            for rule in self.rules["experiment_result"]:
                if isinstance(rule, StatisticalOutlierRule):
                    rule_violations = await rule.validate(result, reference_results)
                    violations.extend(rule_violations)
        
        # Prepare validation report
        validation_report = {
            "entity_id": result.id,
            "entity_type": "experiment_result",
            "timestamp": datetime.now().isoformat(),
            "is_valid": len([v for v in violations if v.severity in ["error", "critical"]]) == 0,
            "violations": [v.to_dict() for v in violations],
            "warning_count": len([v for v in violations if v.severity == "warning"]),
            "error_count": len([v for v in violations if v.severity == "error"]),
            "critical_count": len([v for v in violations if v.severity == "critical"])
        }
        
        return validation_report
    
    async def validate_protocol(self, protocol: Protocol) -> Dict[str, Any]:
        """
        Validate a protocol against applicable rules.
        
        Args:
            protocol: Protocol to validate
            
        Returns:
            Validation results
        """
        violations = []
        
        # Apply protocol rules
        for rule in self.rules["protocol"]:
            rule_violations = await rule.validate(protocol)
            violations.extend(rule_violations)
        
        # Apply generic rules
        if "any" in self.rules:
            for rule in self.rules["any"]:
                rule_violations = await rule.validate(protocol)
                violations.extend(rule_violations)
        
        # Validate protocol steps
        step_validations = []
        for step in protocol.steps:
            step_violations = []
            
            # Apply protocol step rules
            for rule in self.rules["protocol_step"]:
                rule_violations = await rule.validate(step)
                step_violations.extend(rule_violations)
            
            # Apply generic rules to step
            if "any" in self.rules:
                for rule in self.rules["any"]:
                    rule_violations = await rule.validate(step)
                    step_violations.extend(rule_violations)
            
            if step_violations:
                step_validations.append({
                    "step_id": step.id,
                    "step_name": step.name,
                    "violations": [v.to_dict() for v in step_violations]
                })
                violations.extend(step_violations)
        
        # Check step prerequisites
        step_ids = {step.id for step in protocol.steps}
        for step in protocol.steps:
            for prereq_id in step.prerequisites:
                if prereq_id not in step_ids:
                    violations.append(RuleViolation(
                        rule_id="protocol-step-prereq",
                        severity="error",
                        message=f"Step references nonexistent prerequisite step",
                        entity_id=protocol.id,
                        entity_type="protocol",
                        field_name="steps.prerequisites",
                        value=prereq_id,
                        suggestion="Remove or correct the prerequisite reference"
                    ))
        
        # Prepare validation report
        validation_report = {
            "entity_id": protocol.id,
            "entity_type": "protocol",
            "timestamp": datetime.now().isoformat(),
            "is_valid": len([v for v in violations if v.severity in ["error", "critical"]]) == 0,
            "violations": [v.to_dict() for v in violations],
            "warning_count": len([v for v in violations if v.severity == "warning"]),
            "error_count": len([v for v in violations if v.severity == "error"]),
            "critical_count": len([v for v in violations if v.severity == "critical"]),
            "step_validations": step_validations
        }
        
        return validation_report
    
    async def validate_batch(self, entities: List[Any]) -> Dict[str, Any]:
        """
        Validate a batch of entities.
        
        Args:
            entities: List of entities to validate
            
        Returns:
            Batch validation results
        """
        validation_results = []
        entity_counts = {}
        violation_counts = {
            "info": 0,
            "warning": 0,
            "error": 0,
            "critical": 0
        }
        
        # Validate each entity
        for entity in entities:
            entity_type = type(entity).__name__.lower()
            
            # Count entity types
            if entity_type not in entity_counts:
                entity_counts[entity_type] = 0
            entity_counts[entity_type] += 1
            
            # Validate entity
            if isinstance(entity, Experiment):
                result = await self.validate_experiment(entity)
            elif isinstance(entity, ExperimentResult):
                result = await self.validate_experiment_result(entity)
            elif isinstance(entity, Protocol):
                result = await self.validate_protocol(entity)
            else:
                # Skip unsupported entity types
                continue
            
            # Count violations
            for violation in result.get("violations", []):
                severity = violation.get("severity", "info")
                if severity in violation_counts:
                    violation_counts[severity] += 1
            
            validation_results.append(result)
        
        # Prepare batch validation report
        batch_report = {
            "timestamp": datetime.now().isoformat(),
            "entity_counts": entity_counts,
            "violation_counts": violation_counts,
            "validation_count": len(validation_results),
            "is_valid": violation_counts["error"] == 0 and violation_counts["critical"] == 0,
            "results": validation_results
        }
        
        return batch_report
    
    async def create_validation_rule(self, rule_data: Dict[str, Any]) -> ValidationRule:
        """
        Create a validation rule from data.
        
        Args:
            rule_data: Rule configuration data
            
        Returns:
            Created validation rule
            
        Raises:
            ValidationError: If rule data is invalid
        """
        rule_type = rule_data.get("rule_type")
        
        if not rule_type:
            raise ValidationError("Rule type is required")
        
        if rule_type == "range":
            rule = RangeRule(
                rule_id=rule_data.get("rule_id", str(uuid.uuid4())),
                description=rule_data.get("description", ""),
                field_name=rule_data.get("field_name"),
                min_value=rule_data.get("min_value"),
                max_value=rule_data.get("max_value"),
                severity=rule_data.get("severity", "warning"),
                entity_type=rule_data.get("entity_type", "any")
            )
        elif rule_type == "required":
            rule = RequiredFieldRule(
                rule_id=rule_data.get("rule_id", str(uuid.uuid4())),
                description=rule_data.get("description", ""),
                field_name=rule_data.get("field_name"),
                severity=rule_data.get("severity", "warning"),
                entity_type=rule_data.get("entity_type", "any")
            )
        elif rule_type == "format":
            rule = FormatRule(
                rule_id=rule_data.get("rule_id", str(uuid.uuid4())),
                description=rule_data.get("description", ""),
                field_name=rule_data.get("field_name"),
                pattern=rule_data.get("pattern"),
                severity=rule_data.get("severity", "warning"),
                entity_type=rule_data.get("entity_type", "any")
            )
        elif rule_type == "statistical":
            rule = StatisticalOutlierRule(
                rule_id=rule_data.get("rule_id", str(uuid.uuid4())),
                description=rule_data.get("description", ""),
                field_name=rule_data.get("field_name"),
                z_score_threshold=rule_data.get("z_score_threshold", 3.0),
                min_sample_size=rule_data.get("min_sample_size", 5),
                severity=rule_data.get("severity", "warning"),
                entity_type=rule_data.get("entity_type", "any")
            )
        else:
            raise ValidationError(f"Unsupported rule type: {rule_type}")
        
        # Register the rule
        await self.register_rule(rule)
        
        # Persist to database if possible
        if hasattr(self.db_adapter, "create") and callable(self.db_adapter.create):
            rule_dict = rule.to_dict()
            rule_dict["rule_type"] = rule_type
            await self.db_adapter.create("validation_rules", rule_dict)
        
        return rule
    
    async def load_rules_from_database(self) -> None:
        """Load validation rules from the database."""
        if not hasattr(self.db_adapter, "list") or not callable(self.db_adapter.list):
            return
        
        # Fetch rules from database
        rules_data, _ = await self.db_adapter.list(
            "validation_rules",
            filters={},
            offset=0,
            limit=1000  # Arbitrary large limit
        )
        
        # Reset rules
        self.rules = {
            "experiment": [],
            "experiment_result": [],
            "protocol": [],
            "protocol_step": [],
            "tissue_type": [],
            "time_series": [],
            "time_series_data": [],
            "any": []
        }
        
        # Create rules from data
        for rule_data in rules_data:
            try:
                await self.create_validation_rule(rule_data)
            except ValidationError as e:
                print(f"Error creating rule from database: {e}")
    
    async def get_validation_rules(
        self,
        entity_type: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        Get validation rules.
        
        Args:
            entity_type: Filter rules by entity type
            
        Returns:
            List of validation rules
        """
        rules = []
        
        if entity_type:
            # Get rules for specific entity type
            if entity_type in self.rules:
                rules.extend([rule.to_dict() for rule in self.rules[entity_type]])
            
            # Include generic rules
            if "any" in self.rules:
                rules.extend([rule.to_dict() for rule in self.rules["any"]])
        else:
            # Get all rules
            for entity_rules in self.rules.values():
                rules.extend([rule.to_dict() for rule in entity_rules])
        
        return rules