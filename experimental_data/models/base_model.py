#!/usr/bin/env python3
"""
Base experimental data models for the CryoProtect Enhanced Experimental Data System.

This module provides the foundation classes for all experimental data models,
including common utilities for validation, serialization, and metadata handling.
"""

import uuid
from typing import Dict, List, Any, Optional, Union, TypeVar, Generic, Tuple
from datetime import datetime
import json

# Type definitions
T = TypeVar('T')
UUID = str  # Type alias for UUID strings
JSONDict = Dict[str, Any]  # Type alias for JSON-like dictionaries

class ModelError(Exception):
    """Base exception for all model errors."""
    pass

class ValidationError(ModelError):
    """Exception raised for validation errors."""
    pass

class DataIntegrityError(ModelError):
    """Exception raised for data integrity errors."""
    pass

class Uncertainty:
    """Class representing measurement uncertainty."""
    
    def __init__(
        self,
        value: float,
        type: str = "standard",
        confidence: float = 0.95,
        distribution: str = "normal"
    ):
        """
        Initialize an uncertainty value.
        
        Args:
            value: The uncertainty value (standard deviation or range)
            type: Type of uncertainty ('standard', 'expanded', 'range')
            confidence: Confidence level (typically 0.95 for 95%)
            distribution: Statistical distribution ('normal', 'uniform', 'triangular')
        """
        self.value = value
        self.type = type
        self.confidence = confidence
        self.distribution = distribution
        
        # Validate inputs
        self._validate()
    
    def _validate(self) -> None:
        """Validate uncertainty parameters."""
        if self.value < 0:
            raise ValidationError("Uncertainty value cannot be negative")
        
        if self.type not in ["standard", "expanded", "range"]:
            raise ValidationError(f"Invalid uncertainty type: {self.type}")
        
        if self.confidence <= 0 or self.confidence > 1:
            raise ValidationError(f"Confidence level must be between 0 and 1")
        
        if self.distribution not in ["normal", "uniform", "triangular"]:
            raise ValidationError(f"Invalid distribution: {self.distribution}")
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        return {
            "value": self.value,
            "type": self.type,
            "confidence": self.confidence,
            "distribution": self.distribution
        }
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'Uncertainty':
        """Create from dictionary."""
        return cls(
            value=data.get("value", 0),
            type=data.get("type", "standard"),
            confidence=data.get("confidence", 0.95),
            distribution=data.get("distribution", "normal")
        )
    
    def __repr__(self) -> str:
        return f"Uncertainty(value={self.value}, type={self.type}, confidence={self.confidence})"

class Provenance:
    """
    Class representing data provenance information.
    
    Tracks the origin, processing history, and metadata for experimental data.
    """
    
    def __init__(
        self,
        source: str,
        timestamp: Optional[datetime] = None,
        method: Optional[str] = None,
        operator: Optional[str] = None,
        equipment: Optional[Dict[str, Any]] = None,
        processing_steps: Optional[List[Dict[str, Any]]] = None,
        references: Optional[List[Dict[str, Any]]] = None,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize provenance information.
        
        Args:
            source: Source of the data (e.g., "experiment", "literature")
            timestamp: When the data was collected/created
            method: Method used to collect/create the data
            operator: Person who collected/created the data
            equipment: Equipment used to collect the data
            processing_steps: List of processing steps applied to the data
            references: Literature or data references
            metadata: Additional metadata
        """
        self.source = source
        self.timestamp = timestamp or datetime.now()
        self.method = method
        self.operator = operator
        self.equipment = equipment or {}
        self.processing_steps = processing_steps or []
        self.references = references or []
        self.metadata = metadata or {}
    
    def add_processing_step(
        self,
        step_name: str,
        timestamp: Optional[datetime] = None,
        parameters: Optional[Dict[str, Any]] = None,
        software: Optional[str] = None,
        operator: Optional[str] = None
    ) -> None:
        """
        Add a processing step to the provenance history.
        
        Args:
            step_name: Name of the processing step
            timestamp: When the step was performed
            parameters: Parameters used in the processing step
            software: Software used to perform the step
            operator: Person who performed the step
        """
        step = {
            "step_name": step_name,
            "timestamp": timestamp or datetime.now(),
            "parameters": parameters or {},
            "software": software,
            "operator": operator
        }
        self.processing_steps.append(step)
    
    def add_reference(
        self,
        reference_type: str,
        identifier: str,
        description: Optional[str] = None
    ) -> None:
        """
        Add a reference to the provenance information.
        
        Args:
            reference_type: Type of reference (e.g., "doi", "pmid", "url")
            identifier: Reference identifier
            description: Optional description of the reference
        """
        reference = {
            "type": reference_type,
            "identifier": identifier,
            "description": description
        }
        self.references.append(reference)
    
    def to_dict(self) -> JSONDict:
        """Convert to dictionary for serialization."""
        return {
            "source": self.source,
            "timestamp": self.timestamp.isoformat() if self.timestamp else None,
            "method": self.method,
            "operator": self.operator,
            "equipment": self.equipment,
            "processing_steps": [
                {
                    **step,
                    "timestamp": step["timestamp"].isoformat() if isinstance(step["timestamp"], datetime) else step["timestamp"]
                }
                for step in self.processing_steps
            ],
            "references": self.references,
            "metadata": self.metadata
        }
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'Provenance':
        """Create from dictionary."""
        # Parse timestamp
        timestamp = None
        if data.get("timestamp"):
            try:
                timestamp = datetime.fromisoformat(data["timestamp"])
            except ValueError:
                timestamp = datetime.now()
        
        # Parse processing steps
        processing_steps = []
        for step in data.get("processing_steps", []):
            step_copy = step.copy()
            if isinstance(step_copy.get("timestamp"), str):
                try:
                    step_copy["timestamp"] = datetime.fromisoformat(step_copy["timestamp"])
                except ValueError:
                    step_copy["timestamp"] = datetime.now()
            processing_steps.append(step_copy)
        
        return cls(
            source=data.get("source", "unknown"),
            timestamp=timestamp,
            method=data.get("method"),
            operator=data.get("operator"),
            equipment=data.get("equipment", {}),
            processing_steps=processing_steps,
            references=data.get("references", []),
            metadata=data.get("metadata", {})
        )

class BaseModel:
    """Base class for all experimental data models."""
    
    def __init__(
        self,
        id: Optional[UUID] = None,
        created_by: Optional[UUID] = None,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        metadata: Optional[Dict[str, Any]] = None,
        provenance: Optional[Union[Provenance, Dict[str, Any]]] = None
    ):
        """
        Initialize base model.
        
        Args:
            id: Unique identifier for the model instance
            created_by: User ID who created the instance
            created_at: Creation timestamp
            updated_at: Last update timestamp
            metadata: Additional metadata
            provenance: Provenance information
        """
        self.id = id or str(uuid.uuid4())
        self.created_by = created_by
        self.created_at = created_at or datetime.now()
        self.updated_at = updated_at or datetime.now()
        self.metadata = metadata or {}
        
        # Handle provenance
        if provenance is None:
            self.provenance = Provenance(source="system")
        elif isinstance(provenance, dict):
            self.provenance = Provenance.from_dict(provenance)
        else:
            self.provenance = provenance
    
    def validate(self) -> None:
        """
        Validate model data.
        
        Raises:
            ValidationError: If validation fails
        """
        # Base validation (override in subclasses)
        if not self.id:
            raise ValidationError("ID is required")
    
    def to_dict(self) -> JSONDict:
        """
        Convert model to dictionary for serialization.
        
        Returns:
            Dictionary representation of the model
        """
        return {
            "id": self.id,
            "created_by": self.created_by,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "updated_at": self.updated_at.isoformat() if self.updated_at else None,
            "metadata": self.metadata,
            "provenance": self.provenance.to_dict()
        }
    
    @classmethod
    def from_dict(cls, data: JSONDict) -> 'BaseModel':
        """
        Create model from dictionary.
        
        Args:
            data: Dictionary representation of the model
            
        Returns:
            Model instance
        """
        # Parse timestamps
        created_at = None
        if data.get("created_at"):
            try:
                created_at = datetime.fromisoformat(data["created_at"])
            except ValueError:
                created_at = datetime.now()
        
        updated_at = None
        if data.get("updated_at"):
            try:
                updated_at = datetime.fromisoformat(data["updated_at"])
            except ValueError:
                updated_at = datetime.now()
        
        # Parse provenance
        provenance = None
        if data.get("provenance"):
            if isinstance(data["provenance"], dict):
                provenance = Provenance.from_dict(data["provenance"])
            elif isinstance(data["provenance"], str):
                try:
                    provenance = Provenance.from_dict(json.loads(data["provenance"]))
                except json.JSONDecodeError:
                    provenance = Provenance(source="unknown")
        
        return cls(
            id=data.get("id"),
            created_by=data.get("created_by"),
            created_at=created_at,
            updated_at=updated_at,
            metadata=data.get("metadata", {}),
            provenance=provenance
        )
    
    def update(self, data: Dict[str, Any]) -> None:
        """
        Update model with new data.
        
        Args:
            data: Dictionary of fields to update
            
        Raises:
            ValidationError: If validation fails after update
        """
        # Update fields
        for key, value in data.items():
            if hasattr(self, key):
                setattr(self, key, value)
        
        # Update timestamp
        self.updated_at = datetime.now()
        
        # Validate after update
        self.validate()
    
    def add_provenance_step(
        self,
        step_name: str,
        parameters: Optional[Dict[str, Any]] = None,
        software: Optional[str] = None,
        operator: Optional[str] = None
    ) -> None:
        """
        Add a processing step to the provenance history.
        
        Args:
            step_name: Name of the processing step
            parameters: Parameters used in the processing step
            software: Software used to perform the step
            operator: Person who performed the step
        """
        self.provenance.add_processing_step(
            step_name=step_name,
            timestamp=datetime.now(),
            parameters=parameters,
            software=software,
            operator=operator
        )
        
        # Update timestamp
        self.updated_at = datetime.now()