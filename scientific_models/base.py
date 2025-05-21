"""
Base classes for scientific models.

This module defines the base interfaces and abstract classes for all scientific
models in the CryoProtect system.
"""

from abc import ABC, abstractmethod
import logging
from typing import Dict, Any, Optional, List, Tuple, Union
import json
import numpy as np

logger = logging.getLogger(__name__)

class ModelValidationError(Exception):
    """Exception raised when model validation fails."""
    pass

class ModelParameterError(Exception):
    """Exception raised when model parameters are invalid."""
    pass

class ModelCalculationError(Exception):
    """Exception raised when a model calculation fails."""
    pass

class ScientificModel(ABC):
    """
    Base abstract class for all scientific models.
    
    This class provides common functionality for model validation,
    parameter handling, and error handling.
    """
    
    def __init__(self, parameters: Dict[str, Any] = None, name: str = None, description: str = None):
        """
        Initialize the scientific model.
        
        Args:
            parameters: Dictionary of model parameters
            name: Optional name for the model
            description: Optional description of the model
        """
        self.name = name or self.__class__.__name__
        self.description = description or f"{self.__class__.__name__} scientific model"
        self.parameters = parameters or {}
        self.validate_parameters()
        
    @abstractmethod
    def validate_parameters(self) -> None:
        """
        Validate the model parameters.
        
        This method should be implemented by subclasses to validate
        that all required parameters are present and within valid ranges.
        
        Raises:
            ModelParameterError: If parameters are invalid
        """
        pass
    
    @abstractmethod
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Perform the model calculation.
        
        Args:
            inputs: Dictionary of input values
            
        Returns:
            Dictionary of calculation results
            
        Raises:
            ModelValidationError: If inputs are invalid
            ModelCalculationError: If calculation fails
        """
        pass
    
    def to_json(self) -> Dict[str, Any]:
        """
        Convert the model to a JSON-serializable dictionary.
        
        Returns:
            Dictionary representation of the model
        """
        return {
            "name": self.name,
            "description": self.description,
            "model_type": self.__class__.__name__,
            "parameters": self.parameters
        }
    
    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> 'ScientificModel':
        """
        Create a model instance from JSON data.
        
        Args:
            data: Dictionary containing model data
            
        Returns:
            Instance of the model
            
        Raises:
            ModelParameterError: If parameters are invalid
        """
        return cls(
            parameters=data.get("parameters", {}),
            name=data.get("name"),
            description=data.get("description")
        )
    
    def validate_inputs(self, inputs: Dict[str, Any], required_keys: List[str]) -> None:
        """
        Validate that inputs contain all required keys.
        
        Args:
            inputs: Dictionary of input values
            required_keys: List of required keys
            
        Raises:
            ModelValidationError: If required keys are missing
        """
        missing_keys = [key for key in required_keys if key not in inputs]
        if missing_keys:
            raise ModelValidationError(f"Missing required inputs: {', '.join(missing_keys)}")
    
    def safe_calculate(self, inputs: Dict[str, Any]) -> Tuple[Optional[Dict[str, Any]], Optional[Exception]]:
        """
        Perform the calculation with error handling.
        
        Args:
            inputs: Dictionary of input values
            
        Returns:
            Tuple of (results, error) where error is None if successful
        """
        try:
            results = self.calculate(inputs)
            return results, None
        except Exception as e:
            logger.error(f"Error in {self.name} calculation: {str(e)}")
            return None, e
    
    def __repr__(self) -> str:
        """String representation of the model."""
        return f"{self.__class__.__name__}(name='{self.name}', parameters={self.parameters})"


class MoleculePropertyModel(ScientificModel):
    """
    Base class for models that calculate molecular properties.
    
    This class extends ScientificModel with functionality specific
    to molecular property calculations.
    """
    
    def __init__(self, parameters: Dict[str, Any] = None, name: str = None, 
                 description: str = None, property_name: str = None):
        """
        Initialize the molecular property model.
        
        Args:
            parameters: Dictionary of model parameters
            name: Optional name for the model
            description: Optional description of the model
            property_name: Name of the property being modeled
        """
        super().__init__(parameters, name, description)
        self.property_name = property_name or "unknown_property"
    
    def validate_molecule_input(self, inputs: Dict[str, Any]) -> None:
        """
        Validate molecule input data.
        
        Args:
            inputs: Dictionary containing molecule data
            
        Raises:
            ModelValidationError: If molecule data is invalid
        """
        # Check for either SMILES or molecule_id
        if 'smiles' not in inputs and 'molecule_id' not in inputs:
            raise ModelValidationError("Either 'smiles' or 'molecule_id' must be provided")
    
    def to_json(self) -> Dict[str, Any]:
        """
        Convert the model to a JSON-serializable dictionary.
        
        Returns:
            Dictionary representation of the model
        """
        data = super().to_json()
        data["property_name"] = self.property_name
        return data
    
    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> 'MoleculePropertyModel':
        """
        Create a model instance from JSON data.
        
        Args:
            data: Dictionary containing model data
            
        Returns:
            Instance of the model
        """
        instance = super().from_json(data)
        instance.property_name = data.get("property_name", "unknown_property")
        return instance


class MixtureModel(ScientificModel):
    """
    Base class for models that analyze or optimize mixtures.
    
    This class extends ScientificModel with functionality specific
    to mixture analysis and optimization.
    """
    
    def validate_mixture_input(self, inputs: Dict[str, Any]) -> None:
        """
        Validate mixture input data.
        
        Args:
            inputs: Dictionary containing mixture data
            
        Raises:
            ModelValidationError: If mixture data is invalid
        """
        if 'components' not in inputs:
            raise ModelValidationError("Mixture components must be provided")
        
        components = inputs['components']
        if not isinstance(components, list) or len(components) == 0:
            raise ModelValidationError("Components must be a non-empty list")
        
        for component in components:
            if not isinstance(component, dict):
                raise ModelValidationError("Each component must be a dictionary")
            
            if 'molecule_id' not in component and 'smiles' not in component:
                raise ModelValidationError("Each component must have either 'molecule_id' or 'smiles'")
            
            if 'concentration' not in component:
                raise ModelValidationError("Each component must have a 'concentration'")
            
            try:
                concentration = float(component['concentration'])
                if concentration <= 0:
                    raise ModelValidationError(f"Concentration must be positive: {concentration}")
            except (ValueError, TypeError):
                raise ModelValidationError(f"Invalid concentration value: {component.get('concentration')}")


class ExperimentalDataModel(ScientificModel):
    """
    Base class for models that incorporate experimental data.
    
    This class extends ScientificModel with functionality for
    handling experimental data and validation against experiments.
    """
    
    def __init__(self, parameters: Dict[str, Any] = None, name: str = None,
                 description: str = None, experimental_data: List[Dict[str, Any]] = None):
        """
        Initialize the experimental data model.
        
        Args:
            parameters: Dictionary of model parameters
            name: Optional name for the model
            description: Optional description of the model
            experimental_data: List of experimental data points
        """
        super().__init__(parameters, name, description)
        self.experimental_data = experimental_data or []
    
    def validate_against_experiments(self, prediction: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate a model prediction against experimental data.
        
        Args:
            prediction: Model prediction
            
        Returns:
            Dictionary with validation results (error metrics, etc.)
        """
        if not self.experimental_data:
            return {"validation": "No experimental data available"}
        
        # Implementation will depend on the specific model and data
        # This is a placeholder for subclasses to implement
        return {"validation": "Not implemented"}
    
    def to_json(self) -> Dict[str, Any]:
        """
        Convert the model to a JSON-serializable dictionary.
        
        Returns:
            Dictionary representation of the model
        """
        data = super().to_json()
        data["experimental_data_count"] = len(self.experimental_data)
        return data