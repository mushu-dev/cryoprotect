"""
Concentration-dependent scientific models.

This module implements various models for predicting how cryoprotectant properties
change with concentration.
"""

import numpy as np
import math
from typing import Dict, Any, List, Optional, Tuple, Union
import logging

from scientific_models.base import (
    ScientificModel, 
    MoleculePropertyModel,
    ModelValidationError,
    ModelParameterError,
    ModelCalculationError
)

logger = logging.getLogger(__name__)

class ConcentrationModel(MoleculePropertyModel):
    """
    Base class for concentration-dependent models.
    
    This class provides common functionality for modeling how
    cryoprotectant properties change with concentration.
    """
    
    def __init__(self, parameters: Dict[str, Any] = None, name: str = None,
                 description: str = None, property_name: str = None,
                 valid_range: Tuple[float, float] = None):
        """
        Initialize the concentration model.
        
        Args:
            parameters: Dictionary of model parameters
            name: Optional name for the model
            description: Optional description of the model
            property_name: Name of the property being modeled
            valid_range: Tuple of (min_concentration, max_concentration) in M
        """
        super().__init__(parameters, name, description, property_name)
        self.valid_range = valid_range or (0.0, float('inf'))
    
    def validate_parameters(self) -> None:
        """
        Validate that the model parameters are valid.
        
        Raises:
            ModelParameterError: If parameters are invalid
        """
        # Base validation - specific models will extend this
        if self.valid_range[0] < 0:
            raise ModelParameterError("Minimum concentration cannot be negative")
        if self.valid_range[0] >= self.valid_range[1]:
            raise ModelParameterError("Minimum concentration must be less than maximum concentration")
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate property value at the given concentration.
        
        Args:
            inputs: Dictionary containing:
                - concentration: Concentration value in M
                - molecule_id or smiles: Molecule identifier
                
        Returns:
            Dictionary with calculated property and metadata
            
        Raises:
            ModelValidationError: If inputs are invalid
            ModelCalculationError: If calculation fails
        """
        # Validate inputs
        required_keys = ['concentration']
        self.validate_inputs(inputs, required_keys)
        self.validate_molecule_input(inputs)
        
        # Check concentration is within valid range
        concentration = float(inputs['concentration'])
        if concentration < self.valid_range[0] or concentration > self.valid_range[1]:
            raise ModelValidationError(
                f"Concentration {concentration} M is outside valid range "
                f"{self.valid_range[0]}-{self.valid_range[1]} M"
            )
        
        # Calculate property (to be implemented by subclasses)
        value = self._calculate_property(concentration)
        
        return {
            "property_name": self.property_name,
            "concentration": concentration,
            "value": value,
            "units": self._get_property_units(),
            "model_type": self.__class__.__name__,
            "valid_range": self.valid_range
        }
    
    @staticmethod
    def _get_property_units() -> str:
        """Get the units for the property (to be implemented by subclasses)."""
        return ""
    
    def _calculate_property(self, concentration: float) -> float:
        """
        Calculate the property value at the given concentration.
        
        This method should be implemented by subclasses.
        
        Args:
            concentration: Concentration in M
            
        Returns:
            Calculated property value
        """
        raise NotImplementedError("Subclasses must implement _calculate_property")
    
    def to_json(self) -> Dict[str, Any]:
        """
        Convert the model to a JSON-serializable dictionary.
        
        Returns:
            Dictionary representation of the model
        """
        data = super().to_json()
        data["valid_range"] = self.valid_range
        return data
    
    @classmethod
    def from_database(cls, db_record: Dict[str, Any]) -> 'ConcentrationModel':
        """
        Create a model instance from a database record.
        
        Args:
            db_record: Dictionary containing model data from the database
            
        Returns:
            Instance of the appropriate concentration model
        """
        model_type = db_record.get('model_type')
        parameters = db_record.get('parameters', {})
        
        # Convert range from database format to tuple
        valid_range = db_record.get('valid_range')
        if valid_range and isinstance(valid_range, str):
            # Parse range string like "[0,10)" or "(0,10]"
            valid_range = valid_range.strip('()[]')
            min_val, max_val = map(float, valid_range.split(','))
            valid_range = (min_val, max_val)
        else:
            valid_range = (0.0, float('inf'))
        
        # Create the appropriate model type
        if model_type == 'linear':
            return LinearConcentrationModel(
                parameters=parameters,
                name=db_record.get('name'),
                description=db_record.get('description'),
                property_name=db_record.get('property_name'),
                valid_range=valid_range
            )
        elif model_type == 'exponential':
            return ExponentialConcentrationModel(
                parameters=parameters,
                name=db_record.get('name'),
                description=db_record.get('description'),
                property_name=db_record.get('property_name'),
                valid_range=valid_range
            )
        elif model_type == 'sigmoid':
            return SigmoidConcentrationModel(
                parameters=parameters,
                name=db_record.get('name'),
                description=db_record.get('description'),
                property_name=db_record.get('property_name'),
                valid_range=valid_range
            )
        else:
            # Default to linear model
            logger.warning(f"Unknown model type '{model_type}', defaulting to linear")
            return LinearConcentrationModel(
                parameters=parameters,
                name=db_record.get('name'),
                description=db_record.get('description'),
                property_name=db_record.get('property_name'),
                valid_range=valid_range
            )


class LinearConcentrationModel(ConcentrationModel):
    """
    Linear concentration model: value = slope * concentration + intercept
    
    This model assumes a linear relationship between concentration and property value.
    """
    
    def validate_parameters(self) -> None:
        """
        Validate that the model parameters are valid.
        
        Raises:
            ModelParameterError: If parameters are invalid
        """
        super().validate_parameters()
        
        required_params = ['slope', 'intercept']
        missing_params = [param for param in required_params if param not in self.parameters]
        if missing_params:
            raise ModelParameterError(f"Missing required parameters: {', '.join(missing_params)}")
        
        try:
            float(self.parameters['slope'])
            float(self.parameters['intercept'])
        except (ValueError, TypeError):
            raise ModelParameterError("Slope and intercept must be numeric values")
    
    def _calculate_property(self, concentration: float) -> float:
        """
        Calculate the property value using a linear model.
        
        Args:
            concentration: Concentration in M
            
        Returns:
            Calculated property value
        """
        slope = float(self.parameters['slope'])
        intercept = float(self.parameters['intercept'])
        return slope * concentration + intercept


class ExponentialConcentrationModel(ConcentrationModel):
    """
    Exponential concentration model: value = scale * exp(rate * concentration) + offset
    
    This model assumes an exponential relationship between concentration and property value.
    """
    
    def validate_parameters(self) -> None:
        """
        Validate that the model parameters are valid.
        
        Raises:
            ModelParameterError: If parameters are invalid
        """
        super().validate_parameters()
        
        required_params = ['scale', 'rate', 'offset']
        missing_params = [param for param in required_params if param not in self.parameters]
        if missing_params:
            raise ModelParameterError(f"Missing required parameters: {', '.join(missing_params)}")
        
        try:
            float(self.parameters['scale'])
            float(self.parameters['rate'])
            float(self.parameters['offset'])
        except (ValueError, TypeError):
            raise ModelParameterError("Scale, rate, and offset must be numeric values")
    
    def _calculate_property(self, concentration: float) -> float:
        """
        Calculate the property value using an exponential model.
        
        Args:
            concentration: Concentration in M
            
        Returns:
            Calculated property value
        """
        scale = float(self.parameters['scale'])
        rate = float(self.parameters['rate'])
        offset = float(self.parameters['offset'])
        
        try:
            return scale * math.exp(rate * concentration) + offset
        except OverflowError:
            raise ModelCalculationError(f"Exponential calculation overflow at concentration {concentration}")


class SigmoidConcentrationModel(ConcentrationModel):
    """
    Sigmoid (logistic) concentration model: value = bottom + (top - bottom) / (1 + exp((ec50 - concentration) / hillslope))
    
    This model is useful for properties that approach asymptotic values at low and high concentrations.
    """
    
    def validate_parameters(self) -> None:
        """
        Validate that the model parameters are valid.
        
        Raises:
            ModelParameterError: If parameters are invalid
        """
        super().validate_parameters()
        
        required_params = ['bottom', 'top', 'ec50', 'hillslope']
        missing_params = [param for param in required_params if param not in self.parameters]
        if missing_params:
            raise ModelParameterError(f"Missing required parameters: {', '.join(missing_params)}")
        
        try:
            float(self.parameters['bottom'])
            float(self.parameters['top'])
            float(self.parameters['ec50'])
            float(self.parameters['hillslope'])
        except (ValueError, TypeError):
            raise ModelParameterError("Sigmoid parameters must be numeric values")
        
        if float(self.parameters['hillslope']) == 0:
            raise ModelParameterError("Hill slope cannot be zero")
    
    def _calculate_property(self, concentration: float) -> float:
        """
        Calculate the property value using a sigmoid model.
        
        Args:
            concentration: Concentration in M
            
        Returns:
            Calculated property value
        """
        bottom = float(self.parameters['bottom'])
        top = float(self.parameters['top'])
        ec50 = float(self.parameters['ec50'])
        hillslope = float(self.parameters['hillslope'])
        
        try:
            return bottom + (top - bottom) / (1 + math.exp((ec50 - concentration) / hillslope))
        except OverflowError:
            # Handle very large values - return asymptotic value
            if (ec50 - concentration) / hillslope > 0:
                return bottom  # Concentration much less than EC50
            else:
                return top  # Concentration much greater than EC50
        except ZeroDivisionError:
            raise ModelCalculationError("Division by zero in sigmoid calculation")