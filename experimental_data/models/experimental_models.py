#!/usr/bin/env python3
"""
Experimental data models for the CryoProtect Enhanced Experimental Data System.

This module implements the core models for experiments, protocols, results, and
related entities in the enhanced experimental data system.
"""

from typing import Dict, List, Any, Optional, Union, Set, Tuple
from datetime import datetime
import json
import uuid

from .base_model import (
    BaseModel, Uncertainty, Provenance, ValidationError, 
    DataIntegrityError, JSONDict, UUID
)

class ProtocolStep:
    """Class representing a single step in an experimental protocol."""
    
    def __init__(
        self,
        id: Optional[str] = None,
        name: str = "",
        description: str = "",
        step_type: str = "manual",
        parameters: Optional[Dict[str, Any]] = None,
        duration: Optional[Dict[str, Any]] = None,
        temperature: Optional[Dict[str, Any]] = None,
        equipment: Optional[List[Dict[str, Any]]] = None,
        prerequisites: Optional[List[str]] = None,
        validation_rules: Optional[List[Dict[str, Any]]] = None,
        expected_outcomes: Optional[List[Dict[str, Any]]] = None
    ):
        """
        Initialize a protocol step.
        
        Args:
            id: Unique identifier for the step
            name: Name of the step
            description: Detailed description of the step
            step_type: Type of step (manual, automated, measurement)
            parameters: Parameters for the step
            duration: Duration information (value and unit)
            temperature: Temperature information (value, unit, tolerance)
            equipment: Equipment required for the step
            prerequisites: IDs of prerequisite steps
            validation_rules: Rules to validate step execution
            expected_outcomes: Expected outcomes of the step
        """
        self.id = id or str(uuid.uuid4())
        self.name = name
        self.description = description
        self.step_type = step_type
        self.parameters = parameters or {}
        self.duration = duration or {"value": 0, "unit": "min"}
        self.temperature = temperature or {"value": 22, "unit": "°C", "tolerance": 2}
        self.equipment = equipment or []
        self.prerequisites = prerequisites or []
        self.validation_rules = validation_rules or []
        self.expected_outcomes = expected_outcomes or []
    
    def validate(self) -> None:
        """
        Validate step data.
        
        Raises:
            ValidationError: If validation fails
        """
        if not self.name:
            raise ValidationError("Step name is required")
        
        if self.step_type not in ["manual", "automated", "measurement", "calculation"]:
            raise ValidationError(f"Invalid step type: {self.step_type}")
        
        # Validate duration if present
        if self.duration:
            if "value" not in self.duration or "unit" not in self.duration:
                raise ValidationError("Duration must include value and unit")
            
            if not isinstance(self.duration["value"], (int, float)) or self.duration["value"] < 0:
                raise ValidationError("Duration value must be a non-negative number")
            
            if self.duration["unit"] not in ["s", "min", "h", "d"]:
                raise ValidationError(f"Invalid duration unit: {self.duration['unit']}")
        
        # Validate temperature if present
        if self.temperature:
            if "value" not in self.temperature or "unit" not in self.temperature:
                raise ValidationError("Temperature must include value and unit")
            
            if not isinstance(self.temperature["value"], (int, float)):
                raise ValidationError("Temperature value must be a number")
            
            if self.temperature["unit"] not in ["°C", "K", "°F"]:
                raise ValidationError(f"Invalid temperature unit: {self.temperature['unit']}")
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "step_type": self.step_type,
            "parameters": self.parameters,
            "duration": self.duration,
            "temperature": self.temperature,
            "equipment": self.equipment,
            "prerequisites": self.prerequisites,
            "validation_rules": self.validation_rules,
            "expected_outcomes": self.expected_outcomes
        }
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'ProtocolStep':
        """Create from dictionary."""
        return cls(
            id=data.get("id"),
            name=data.get("name", ""),
            description=data.get("description", ""),
            step_type=data.get("step_type", "manual"),
            parameters=data.get("parameters", {}),
            duration=data.get("duration"),
            temperature=data.get("temperature"),
            equipment=data.get("equipment"),
            prerequisites=data.get("prerequisites"),
            validation_rules=data.get("validation_rules"),
            expected_outcomes=data.get("expected_outcomes")
        )

class Protocol(BaseModel):
    """Class representing an experimental protocol."""
    
    def __init__(
        self,
        name: str,
        description: str = "",
        steps: Optional[List[Union[ProtocolStep, Dict[str, Any]]]] = None,
        version: int = 1,
        parent_id: Optional[UUID] = None,
        template_type: str = "standard",
        target_application: Optional[str] = None,
        author: Optional[str] = None,
        parameters: Optional[Dict[str, Any]] = None,
        **kwargs
    ):
        """
        Initialize a protocol.
        
        Args:
            name: Name of the protocol
            description: Detailed description of the protocol
            steps: List of protocol steps
            version: Version number
            parent_id: ID of parent protocol (for derived protocols)
            template_type: Type of protocol template
            target_application: Target application area
            author: Protocol author
            parameters: Global protocol parameters
            **kwargs: Additional arguments for BaseModel
        """
        super().__init__(**kwargs)
        self.name = name
        self.description = description
        self.version = version
        self.parent_id = parent_id
        self.template_type = template_type
        self.target_application = target_application
        self.author = author
        self.parameters = parameters or {}
        
        # Process steps
        self.steps = []
        if steps:
            for step in steps:
                if isinstance(step, ProtocolStep):
                    self.steps.append(step)
                elif isinstance(step, dict):
                    self.steps.append(ProtocolStep.from_dict(step))
                else:
                    raise ValidationError(f"Invalid step type: {type(step)}")
    
    def validate(self) -> None:
        """
        Validate protocol data.
        
        Raises:
            ValidationError: If validation fails
        """
        super().validate()
        
        if not self.name:
            raise ValidationError("Protocol name is required")
        
        if self.version < 1:
            raise ValidationError("Version must be >= 1")
        
        # Check step prerequisites
        step_ids = {step.id for step in self.steps}
        for step in self.steps:
            for prereq_id in step.prerequisites:
                if prereq_id not in step_ids:
                    raise ValidationError(f"Step {step.id} has missing prerequisite: {prereq_id}")
        
        # Validate each step
        for step in self.steps:
            step.validate()
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "name": self.name,
            "description": self.description,
            "steps": [step.to_dict() for step in self.steps],
            "version": self.version,
            "parent_id": self.parent_id,
            "template_type": self.template_type,
            "target_application": self.target_application,
            "author": self.author,
            "parameters": self.parameters
        })
        return result
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'Protocol':
        """Create from dictionary."""
        steps_data = data.pop("steps", []) if "steps" in data else []
        
        protocol = super().from_dict(data)
        protocol.name = data.get("name", "")
        protocol.description = data.get("description", "")
        protocol.version = data.get("version", 1)
        protocol.parent_id = data.get("parent_id")
        protocol.template_type = data.get("template_type", "standard")
        protocol.target_application = data.get("target_application")
        protocol.author = data.get("author")
        protocol.parameters = data.get("parameters", {})
        
        # Process steps
        protocol.steps = [ProtocolStep.from_dict(step) for step in steps_data]
        
        return protocol
    
    def add_step(self, step: Union[ProtocolStep, Dict[str, Any]]) -> None:
        """
        Add a step to the protocol.
        
        Args:
            step: Protocol step to add
        """
        if isinstance(step, ProtocolStep):
            self.steps.append(step)
        elif isinstance(step, dict):
            self.steps.append(ProtocolStep.from_dict(step))
        else:
            raise ValidationError(f"Invalid step type: {type(step)}")
        
        # Update timestamp
        self.updated_at = datetime.now()
    
    def remove_step(self, step_id: str) -> None:
        """
        Remove a step from the protocol.
        
        Args:
            step_id: ID of the step to remove
        """
        self.steps = [step for step in self.steps if step.id != step_id]
        
        # Update timestamp
        self.updated_at = datetime.now()
    
    def update_step(self, step_id: str, data: Dict[str, Any]) -> None:
        """
        Update a step in the protocol.
        
        Args:
            step_id: ID of the step to update
            data: Dictionary of fields to update
        """
        for i, step in enumerate(self.steps):
            if step.id == step_id:
                for key, value in data.items():
                    if hasattr(step, key):
                        setattr(step, key, value)
                
                # Validate step
                step.validate()
                
                # Update timestamp
                self.updated_at = datetime.now()
                return
        
        raise ValidationError(f"Step not found: {step_id}")
    
    def create_version(self) -> 'Protocol':
        """
        Create a new version of the protocol.
        
        Returns:
            New Protocol instance with incremented version
        """
        new_protocol = Protocol(
            name=self.name,
            description=self.description,
            steps=[step.to_dict() for step in self.steps],
            version=self.version + 1,
            parent_id=self.id,
            template_type=self.template_type,
            target_application=self.target_application,
            author=self.author,
            parameters=self.parameters.copy(),
            metadata=self.metadata.copy(),
            provenance=self.provenance.to_dict()
        )
        
        # Add version change to provenance
        new_protocol.add_provenance_step(
            step_name="version_increment",
            parameters={"previous_version": self.version}
        )
        
        return new_protocol

class ExperimentType(BaseModel):
    """Class representing a type of experiment."""
    
    def __init__(
        self,
        name: str,
        description: str = "",
        protocol_details: Optional[Dict[str, Any]] = None,
        **kwargs
    ):
        """
        Initialize an experiment type.
        
        Args:
            name: Name of the experiment type
            description: Detailed description of the experiment type
            protocol_details: Protocol details specific to this experiment type
            **kwargs: Additional arguments for BaseModel
        """
        super().__init__(**kwargs)
        self.name = name
        self.description = description
        self.protocol_details = protocol_details or {}
    
    def validate(self) -> None:
        """
        Validate experiment type data.
        
        Raises:
            ValidationError: If validation fails
        """
        super().validate()
        
        if not self.name:
            raise ValidationError("Experiment type name is required")
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "name": self.name,
            "description": self.description,
            "protocol_details": self.protocol_details
        })
        return result
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'ExperimentType':
        """Create from dictionary."""
        instance = super().from_dict(data)
        instance.name = data.get("name", "")
        instance.description = data.get("description", "")
        instance.protocol_details = data.get("protocol_details", {})
        return instance

class TissueType(BaseModel):
    """Class representing a type of biological tissue or sample."""
    
    def __init__(
        self,
        name: str,
        description: str = "",
        species: Optional[str] = None,
        taxonomy_id: Optional[int] = None,
        properties: Optional[Dict[str, Any]] = None,
        **kwargs
    ):
        """
        Initialize a tissue type.
        
        Args:
            name: Name of the tissue type
            description: Detailed description of the tissue type
            species: Species of origin
            taxonomy_id: Taxonomy ID (e.g., NCBI Taxonomy ID)
            properties: Additional properties of the tissue
            **kwargs: Additional arguments for BaseModel
        """
        super().__init__(**kwargs)
        self.name = name
        self.description = description
        self.species = species
        self.taxonomy_id = taxonomy_id
        self.properties = properties or {}
    
    def validate(self) -> None:
        """
        Validate tissue type data.
        
        Raises:
            ValidationError: If validation fails
        """
        super().validate()
        
        if not self.name:
            raise ValidationError("Tissue type name is required")
        
        if self.taxonomy_id is not None and not isinstance(self.taxonomy_id, int):
            raise ValidationError("Taxonomy ID must be an integer")
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "name": self.name,
            "description": self.description,
            "species": self.species,
            "taxonomy_id": self.taxonomy_id,
            "properties": self.properties
        })
        return result
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'TissueType':
        """Create from dictionary."""
        instance = super().from_dict(data)
        instance.name = data.get("name", "")
        instance.description = data.get("description", "")
        instance.species = data.get("species")
        instance.taxonomy_id = data.get("taxonomy_id")
        instance.properties = data.get("properties", {})
        return instance

class Experiment(BaseModel):
    """Class representing an experiment."""
    
    def __init__(
        self,
        title: str,
        experiment_type_id: UUID,
        description: str = "",
        protocol_id: Optional[UUID] = None,
        date_performed: Optional[Union[datetime, str]] = None,
        temperature: Optional[float] = None,
        cooling_rate: Optional[float] = None,
        thawing_rate: Optional[float] = None,
        parameters: Optional[Dict[str, Any]] = None,
        version: int = 1,
        **kwargs
    ):
        """
        Initialize an experiment.
        
        Args:
            title: Title of the experiment
            experiment_type_id: ID of the experiment type
            description: Detailed description of the experiment
            protocol_id: ID of the protocol used
            date_performed: Date when the experiment was performed
            temperature: Temperature of the experiment
            cooling_rate: Cooling rate used
            thawing_rate: Thawing rate used
            parameters: Additional parameters for the experiment
            version: Version number
            **kwargs: Additional arguments for BaseModel
        """
        super().__init__(**kwargs)
        self.title = title
        self.experiment_type_id = experiment_type_id
        self.description = description
        self.protocol_id = protocol_id
        
        # Handle date_performed
        if date_performed is None:
            self.date_performed = None
        elif isinstance(date_performed, datetime):
            self.date_performed = date_performed
        elif isinstance(date_performed, str):
            try:
                self.date_performed = datetime.fromisoformat(date_performed)
            except ValueError:
                self.date_performed = None
        else:
            self.date_performed = None
        
        self.temperature = temperature
        self.cooling_rate = cooling_rate
        self.thawing_rate = thawing_rate
        self.parameters = parameters or {}
        self.version = version
    
    def validate(self) -> None:
        """
        Validate experiment data.
        
        Raises:
            ValidationError: If validation fails
        """
        super().validate()
        
        if not self.title:
            raise ValidationError("Experiment title is required")
        
        if not self.experiment_type_id:
            raise ValidationError("Experiment type ID is required")
        
        if self.version < 1:
            raise ValidationError("Version must be >= 1")
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "title": self.title,
            "experiment_type_id": self.experiment_type_id,
            "description": self.description,
            "protocol_id": self.protocol_id,
            "date_performed": self.date_performed.isoformat() if self.date_performed else None,
            "temperature": self.temperature,
            "cooling_rate": self.cooling_rate,
            "thawing_rate": self.thawing_rate,
            "parameters": self.parameters,
            "version": self.version
        })
        return result
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'Experiment':
        """Create from dictionary."""
        instance = super().from_dict(data)
        instance.title = data.get("title", "")
        instance.experiment_type_id = data.get("experiment_type_id")
        instance.description = data.get("description", "")
        instance.protocol_id = data.get("protocol_id")
        
        # Handle date_performed
        date_performed = data.get("date_performed")
        if date_performed is not None:
            if isinstance(date_performed, str):
                try:
                    instance.date_performed = datetime.fromisoformat(date_performed)
                except ValueError:
                    instance.date_performed = None
            else:
                instance.date_performed = None
        else:
            instance.date_performed = None
        
        instance.temperature = data.get("temperature")
        instance.cooling_rate = data.get("cooling_rate")
        instance.thawing_rate = data.get("thawing_rate")
        instance.parameters = data.get("parameters", {})
        instance.version = data.get("version", 1)
        
        return instance

class ExperimentResult(BaseModel):
    """Class representing the result of an experiment."""
    
    def __init__(
        self,
        experiment_id: UUID,
        tissue_type_id: UUID,
        molecule_id: Optional[UUID] = None,
        mixture_id: Optional[UUID] = None,
        concentration: Optional[float] = None,
        concentration_unit: Optional[str] = None,
        viability_percentage: Optional[float] = None,
        recovery_rate: Optional[float] = None,
        functionality_score: Optional[float] = None,
        uncertainty: Optional[Union[Uncertainty, Dict[str, Any], Dict[str, Dict[str, Any]]]] = None,
        result_details: Optional[Dict[str, Any]] = None,
        notes: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize an experiment result.
        
        Args:
            experiment_id: ID of the experiment
            tissue_type_id: ID of the tissue type
            molecule_id: ID of the molecule (mutually exclusive with mixture_id)
            mixture_id: ID of the mixture (mutually exclusive with molecule_id)
            concentration: Concentration used
            concentration_unit: Unit of concentration
            viability_percentage: Viability percentage result
            recovery_rate: Recovery rate result
            functionality_score: Functionality score result
            uncertainty: Uncertainty values for the measurements
            result_details: Additional result details
            notes: Notes about the result
            **kwargs: Additional arguments for BaseModel
        """
        super().__init__(**kwargs)
        self.experiment_id = experiment_id
        self.tissue_type_id = tissue_type_id
        self.molecule_id = molecule_id
        self.mixture_id = mixture_id
        self.concentration = concentration
        self.concentration_unit = concentration_unit
        self.viability_percentage = viability_percentage
        self.recovery_rate = recovery_rate
        self.functionality_score = functionality_score
        self.result_details = result_details or {}
        self.notes = notes
        
        # Process uncertainty data
        self._process_uncertainty(uncertainty)
    
    def _process_uncertainty(self, uncertainty: Optional[Union[Uncertainty, Dict[str, Any], Dict[str, Dict[str, Any]]]]) -> None:
        """
        Process uncertainty data.
        
        Args:
            uncertainty: Uncertainty data in various formats
        """
        self.uncertainty = {}
        
        if uncertainty is None:
            return
        
        if isinstance(uncertainty, Uncertainty):
            # Single uncertainty for all measurements
            self.uncertainty = {
                "viability_percentage": uncertainty,
                "recovery_rate": uncertainty,
                "functionality_score": uncertainty
            }
        elif isinstance(uncertainty, dict):
            # First check if it's a simple uncertainty definition
            if "value" in uncertainty:
                u = Uncertainty.from_dict(uncertainty)
                self.uncertainty = {
                    "viability_percentage": u,
                    "recovery_rate": u,
                    "functionality_score": u
                }
            else:
                # Process per-metric uncertainty values
                for key, value in uncertainty.items():
                    if isinstance(value, dict):
                        self.uncertainty[key] = Uncertainty.from_dict(value)
                    elif isinstance(value, Uncertainty):
                        self.uncertainty[key] = value
    
    def validate(self) -> None:
        """
        Validate experiment result data.
        
        Raises:
            ValidationError: If validation fails
        """
        super().validate()
        
        if not self.experiment_id:
            raise ValidationError("Experiment ID is required")
        
        if not self.tissue_type_id:
            raise ValidationError("Tissue type ID is required")
        
        # Check molecule/mixture mutual exclusivity
        if self.molecule_id and self.mixture_id:
            raise ValidationError("Only one of molecule_id or mixture_id can be set")
        
        if not self.molecule_id and not self.mixture_id:
            raise ValidationError("Either molecule_id or mixture_id must be set")
        
        # Check numeric ranges
        if self.viability_percentage is not None and (self.viability_percentage < 0 or self.viability_percentage > 100):
            raise ValidationError("Viability percentage must be between 0 and 100")
        
        if self.recovery_rate is not None and self.recovery_rate < 0:
            raise ValidationError("Recovery rate cannot be negative")
        
        # Check concentration unit
        if self.concentration is not None and not self.concentration_unit:
            raise ValidationError("Concentration unit is required when concentration is provided")
        
        if self.concentration_unit and self.concentration_unit not in ["M", "mM", "%w/v", "%v/v", "mg/mL", "g/L"]:
            raise ValidationError(f"Invalid concentration unit: {self.concentration_unit}")
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        
        # Convert uncertainty objects to dictionaries
        uncertainty_dict = {}
        for key, value in self.uncertainty.items():
            if isinstance(value, Uncertainty):
                uncertainty_dict[key] = value.to_dict()
            else:
                uncertainty_dict[key] = value
        
        result.update({
            "experiment_id": self.experiment_id,
            "tissue_type_id": self.tissue_type_id,
            "molecule_id": self.molecule_id,
            "mixture_id": self.mixture_id,
            "concentration": self.concentration,
            "concentration_unit": self.concentration_unit,
            "viability_percentage": self.viability_percentage,
            "recovery_rate": self.recovery_rate,
            "functionality_score": self.functionality_score,
            "uncertainty": uncertainty_dict,
            "result_details": self.result_details,
            "notes": self.notes
        })
        return result
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'ExperimentResult':
        """Create from dictionary."""
        # Extract uncertainty data
        uncertainty_data = data.pop("uncertainty", {}) if "uncertainty" in data else {}
        
        instance = super().from_dict(data)
        instance.experiment_id = data.get("experiment_id")
        instance.tissue_type_id = data.get("tissue_type_id")
        instance.molecule_id = data.get("molecule_id")
        instance.mixture_id = data.get("mixture_id")
        instance.concentration = data.get("concentration")
        instance.concentration_unit = data.get("concentration_unit")
        instance.viability_percentage = data.get("viability_percentage")
        instance.recovery_rate = data.get("recovery_rate")
        instance.functionality_score = data.get("functionality_score")
        instance.result_details = data.get("result_details", {})
        instance.notes = data.get("notes")
        
        # Process uncertainty data
        instance._process_uncertainty(uncertainty_data)
        
        return instance

class TimeSeries(BaseModel):
    """Class representing a time series of experimental data."""
    
    def __init__(
        self,
        experiment_id: UUID,
        name: str,
        property_type: str,
        description: str = "",
        unit: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        **kwargs
    ):
        """
        Initialize a time series.
        
        Args:
            experiment_id: ID of the experiment
            name: Name of the time series
            property_type: Type of property measured
            description: Detailed description of the time series
            unit: Unit of measurement
            metadata: Additional metadata
            **kwargs: Additional arguments for BaseModel
        """
        super().__init__(**kwargs)
        self.experiment_id = experiment_id
        self.name = name
        self.property_type = property_type
        self.description = description
        self.unit = unit
        self.metadata = metadata or {}
    
    def validate(self) -> None:
        """
        Validate time series data.
        
        Raises:
            ValidationError: If validation fails
        """
        super().validate()
        
        if not self.experiment_id:
            raise ValidationError("Experiment ID is required")
        
        if not self.name:
            raise ValidationError("Time series name is required")
        
        if not self.property_type:
            raise ValidationError("Property type is required")
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "experiment_id": self.experiment_id,
            "name": self.name,
            "property_type": self.property_type,
            "description": self.description,
            "unit": self.unit
        })
        return result
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'TimeSeries':
        """Create from dictionary."""
        instance = super().from_dict(data)
        instance.experiment_id = data.get("experiment_id")
        instance.name = data.get("name", "")
        instance.property_type = data.get("property_type", "")
        instance.description = data.get("description", "")
        instance.unit = data.get("unit")
        
        # Metadata is already handled by BaseModel
        
        return instance

class TimeSeriesDataPoint(BaseModel):
    """Class representing a single data point in a time series."""
    
    def __init__(
        self,
        time_series_id: UUID,
        timestamp: Union[datetime, str],
        value: float,
        uncertainty: Optional[Union[float, Uncertainty, Dict[str, Any]]] = None,
        **kwargs
    ):
        """
        Initialize a time series data point.
        
        Args:
            time_series_id: ID of the time series
            timestamp: Timestamp of the data point
            value: Measured value
            uncertainty: Uncertainty of the measurement
            **kwargs: Additional arguments for BaseModel
        """
        super().__init__(**kwargs)
        self.time_series_id = time_series_id
        
        # Handle timestamp
        if isinstance(timestamp, datetime):
            self.timestamp = timestamp
        elif isinstance(timestamp, str):
            try:
                self.timestamp = datetime.fromisoformat(timestamp)
            except ValueError:
                raise ValidationError(f"Invalid timestamp format: {timestamp}")
        else:
            raise ValidationError(f"Invalid timestamp type: {type(timestamp)}")
        
        self.value = value
        
        # Process uncertainty
        if uncertainty is None:
            self.uncertainty = None
        elif isinstance(uncertainty, float):
            self.uncertainty = Uncertainty(value=uncertainty)
        elif isinstance(uncertainty, Uncertainty):
            self.uncertainty = uncertainty
        elif isinstance(uncertainty, dict):
            self.uncertainty = Uncertainty.from_dict(uncertainty)
        else:
            raise ValidationError(f"Invalid uncertainty type: {type(uncertainty)}")
    
    def validate(self) -> None:
        """
        Validate time series data point.
        
        Raises:
            ValidationError: If validation fails
        """
        super().validate()
        
        if not self.time_series_id:
            raise ValidationError("Time series ID is required")
        
        if not self.timestamp:
            raise ValidationError("Timestamp is required")
        
        if not isinstance(self.value, (int, float)):
            raise ValidationError("Value must be a number")
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        result = super().to_dict()
        result.update({
            "time_series_id": self.time_series_id,
            "timestamp": self.timestamp.isoformat(),
            "value": self.value,
            "uncertainty": self.uncertainty.to_dict() if self.uncertainty else None
        })
        return result
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'TimeSeriesDataPoint':
        """Create from dictionary."""
        uncertainty_data = data.pop("uncertainty", None) if "uncertainty" in data else None
        timestamp_data = data.pop("timestamp") if "timestamp" in data else None
        value_data = data.pop("value") if "value" in data else None
        
        instance = super().from_dict(data)
        instance.time_series_id = data.get("time_series_id")
        
        # Set timestamp
        if timestamp_data:
            if isinstance(timestamp_data, str):
                try:
                    instance.timestamp = datetime.fromisoformat(timestamp_data)
                except ValueError:
                    raise ValidationError(f"Invalid timestamp format: {timestamp_data}")
            elif isinstance(timestamp_data, datetime):
                instance.timestamp = timestamp_data
            else:
                raise ValidationError(f"Invalid timestamp type: {type(timestamp_data)}")
        else:
            raise ValidationError("Timestamp is required")
        
        # Set value
        instance.value = value_data
        
        # Process uncertainty
        if uncertainty_data is None:
            instance.uncertainty = None
        elif isinstance(uncertainty_data, float):
            instance.uncertainty = Uncertainty(value=uncertainty_data)
        elif isinstance(uncertainty_data, dict):
            instance.uncertainty = Uncertainty.from_dict(uncertainty_data)
        else:
            raise ValidationError(f"Invalid uncertainty type: {type(uncertainty_data)}")
        
        return instance