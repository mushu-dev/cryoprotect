"""
Temperature-dependent scientific models.

This module implements various models for predicting how cryoprotectant properties
change with temperature. It includes models for glass transition temperature
prediction and temperature-dependent effectiveness.
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

# Physical constants
R = 8.314462618  # Universal gas constant, J/(mol·K)
CELSIUS_TO_KELVIN = 273.15  # Conversion factor from °C to K

class TemperatureModel(MoleculePropertyModel):
    """
    Base class for temperature-dependent models.
    
    This class provides common functionality for modeling how
    cryoprotectant properties change with temperature.
    """
    
    def __init__(self, parameters: Dict[str, Any] = None, name: str = None,
                 description: str = None, property_name: str = None,
                 temperature_range: Tuple[float, float] = None):
        """
        Initialize the temperature model.
        
        Args:
            parameters: Dictionary of model parameters
            name: Optional name for the model
            description: Optional description of the model
            property_name: Name of the property being modeled
            temperature_range: Tuple of (min_temperature, max_temperature) in Kelvin
        """
        super().__init__(parameters, name, description, property_name)
        self.temperature_range = temperature_range or (0.0, 500.0)  # Default range in K
    
    def validate_parameters(self) -> None:
        """
        Validate that the model parameters are valid.
        
        Raises:
            ModelParameterError: If parameters are invalid
        """
        # Base validation - specific models will extend this
        if self.temperature_range[0] < 0:
            raise ModelParameterError("Minimum temperature cannot be negative")
        if self.temperature_range[0] >= self.temperature_range[1]:
            raise ModelParameterError("Minimum temperature must be less than maximum temperature")
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate property value at the given temperature.
        
        Args:
            inputs: Dictionary containing:
                - temperature: Temperature value in K
                - molecule_id or smiles: Molecule identifier
                
        Returns:
            Dictionary with calculated property and metadata
            
        Raises:
            ModelValidationError: If inputs are invalid
            ModelCalculationError: If calculation fails
        """
        # Validate inputs
        required_keys = ['temperature']
        self.validate_inputs(inputs, required_keys)
        self.validate_molecule_input(inputs)
        
        # Check if temperature is in Celsius and convert to Kelvin if needed
        temperature = float(inputs['temperature'])
        temp_unit = inputs.get('temperature_unit', 'K')
        
        if temp_unit.upper() == 'C':
            temperature += CELSIUS_TO_KELVIN
        
        # Check temperature is within valid range
        if temperature < self.temperature_range[0] or temperature > self.temperature_range[1]:
            raise ModelValidationError(
                f"Temperature {temperature} K is outside valid range "
                f"{self.temperature_range[0]}-{self.temperature_range[1]} K"
            )
        
        # Calculate property (to be implemented by subclasses)
        value = self._calculate_property(temperature)
        
        return {
            "property_name": self.property_name,
            "temperature": temperature,
            "temperature_unit": "K",
            "value": value,
            "units": self._get_property_units(),
            "model_type": self.__class__.__name__,
            "valid_range": self.temperature_range
        }
    
    @staticmethod
    def _get_property_units() -> str:
        """Get the units for the property (to be implemented by subclasses)."""
        return ""
    
    def _calculate_property(self, temperature: float) -> float:
        """
        Calculate the property value at the given temperature.
        
        This method should be implemented by subclasses.
        
        Args:
            temperature: Temperature in Kelvin
            
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
        data["temperature_range"] = self.temperature_range
        return data
    
    @classmethod
    def from_database(cls, db_record: Dict[str, Any]) -> 'TemperatureModel':
        """
        Create a model instance from a database record.
        
        Args:
            db_record: Dictionary containing model data from the database
            
        Returns:
            Instance of the appropriate temperature model
        """
        model_type = db_record.get('model_type')
        parameters = db_record.get('parameters', {})
        
        # Convert range from database format to tuple
        temperature_range = db_record.get('temperature_range')
        if temperature_range and isinstance(temperature_range, str):
            # Parse range string like "[0,500)" or "(0,500]"
            temperature_range = temperature_range.strip('()[]')
            min_val, max_val = map(float, temperature_range.split(','))
            temperature_range = (min_val, max_val)
        else:
            temperature_range = (0.0, 500.0)
        
        # Create the appropriate model type
        if model_type == 'linear':
            return LinearTemperatureModel(
                parameters=parameters,
                name=db_record.get('name'),
                description=db_record.get('description'),
                property_name=db_record.get('property_name'),
                temperature_range=temperature_range
            )
        elif model_type == 'arrhenius':
            return ArrheniusTemperatureModel(
                parameters=parameters,
                name=db_record.get('name'),
                description=db_record.get('description'),
                property_name=db_record.get('property_name'),
                temperature_range=temperature_range
            )
        else:
            # Default to linear model
            logger.warning(f"Unknown model type '{model_type}', defaulting to linear")
            return LinearTemperatureModel(
                parameters=parameters,
                name=db_record.get('name'),
                description=db_record.get('description'),
                property_name=db_record.get('property_name'),
                temperature_range=temperature_range
            )


class LinearTemperatureModel(TemperatureModel):
    """
    Linear temperature model: value = slope * temperature + intercept
    
    This model assumes a linear relationship between temperature and property value.
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
    
    def _calculate_property(self, temperature: float) -> float:
        """
        Calculate the property value using a linear model.
        
        Args:
            temperature: Temperature in Kelvin
            
        Returns:
            Calculated property value
        """
        slope = float(self.parameters['slope'])
        intercept = float(self.parameters['intercept'])
        return slope * temperature + intercept


class ArrheniusTemperatureModel(TemperatureModel):
    """
    Arrhenius model: value = A * exp(-Ea / (R * temperature))
    
    This model is based on the Arrhenius equation for temperature-dependent
    rate constants in chemical reactions. It's useful for modeling properties
    that follow an exponential relationship with 1/T.
    """
    
    def validate_parameters(self) -> None:
        """
        Validate that the model parameters are valid.
        
        Raises:
            ModelParameterError: If parameters are invalid
        """
        super().validate_parameters()
        
        required_params = ['A', 'Ea']
        missing_params = [param for param in required_params if param not in self.parameters]
        if missing_params:
            raise ModelParameterError(f"Missing required parameters: {', '.join(missing_params)}")
        
        try:
            float(self.parameters['A'])
            float(self.parameters['Ea'])
        except (ValueError, TypeError):
            raise ModelParameterError("A and Ea must be numeric values")
            
        if float(self.parameters['A']) <= 0:
            raise ModelParameterError("Pre-exponential factor A must be positive")
    
    def _calculate_property(self, temperature: float) -> float:
        """
        Calculate the property value using the Arrhenius equation.
        
        Args:
            temperature: Temperature in Kelvin
            
        Returns:
            Calculated property value
        """
        A = float(self.parameters['A'])  # Pre-exponential factor
        Ea = float(self.parameters['Ea'])  # Activation energy in J/mol
        
        try:
            return A * math.exp(-Ea / (R * temperature))
        except OverflowError:
            raise ModelCalculationError(f"Arrhenius calculation overflow at temperature {temperature}")
        except ZeroDivisionError:
            raise ModelCalculationError("Division by zero in Arrhenius calculation")


class GlassTransitionModel(TemperatureModel):
    """
    Glass transition temperature (Tg) predictor.
    
    This model predicts the glass transition temperature for cryoprotectants
    and their mixtures, which is a critical parameter for cryopreservation.
    """
    
    def __init__(self, parameters: Dict[str, Any] = None, name: str = None,
                 description: str = None):
        """
        Initialize the glass transition temperature predictor.
        
        Args:
            parameters: Dictionary of model parameters
            name: Optional name for the model
            description: Optional description of the model
        """
        super().__init__(
            parameters, 
            name or "Glass Transition Temperature Predictor",
            description or "Predicts the glass transition temperature (Tg) for cryoprotectants",
            property_name="glass_transition_temperature"
        )
    
    def validate_parameters(self) -> None:
        """
        Validate that the model parameters are valid.
        
        Raises:
            ModelParameterError: If parameters are invalid
        """
        super().validate_parameters()
        
        # Different parameter requirements based on model subtype
        model_subtype = self.parameters.get('model_subtype', 'gordon_taylor')
        
        if model_subtype == 'gordon_taylor':
            required_params = ['k']
            missing_params = [param for param in required_params if param not in self.parameters]
            if missing_params:
                raise ModelParameterError(f"Missing required parameters for Gordon-Taylor model: {', '.join(missing_params)}")
                
            try:
                float(self.parameters['k'])
            except (ValueError, TypeError):
                raise ModelParameterError("Gordon-Taylor constant k must be a numeric value")
        
        elif model_subtype == 'fox':
            pass  # Fox equation doesn't require additional parameters
            
        else:
            raise ModelParameterError(f"Unknown glass transition model subtype: {model_subtype}")
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate the glass transition temperature.
        
        For pure compounds, this returns the Tg value.
        For mixtures, it calculates Tg based on composition.
        
        Args:
            inputs: Dictionary containing:
                - molecule_id or smiles: Molecule identifier for pure compound
                - components: List of components with weights for mixtures
                
        Returns:
            Dictionary with calculated Tg and metadata
            
        Raises:
            ModelValidationError: If inputs are invalid
            ModelCalculationError: If calculation fails
        """
        # Validate molecule input - either molecule_id/smiles or components required
        if 'components' in inputs:
            # Mixture calculation
            return self._calculate_mixture_tg(inputs)
        else:
            # Pure compound calculation
            self.validate_molecule_input(inputs)
            return self._calculate_pure_tg(inputs)
    
    def _calculate_pure_tg(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate the glass transition temperature for a pure compound.
        
        This implementation uses a simple lookup or estimation based on
        molecular structure. In a real implementation, this would use
        more sophisticated methods based on molecular descriptors.
        
        Args:
            inputs: Dictionary with molecule information
            
        Returns:
            Dictionary with Tg value and metadata
        """
        # This is a simplified implementation
        # A real implementation would use molecular descriptors and QSPR models
        
        # For now, return a mock value based on molecular weight
        molecule_id = inputs.get('molecule_id')
        smiles = inputs.get('smiles')
        
        # Mock implementation - in a real system, this would query a database
        # or use a real predictive model
        if molecule_id:
            # Get Tg value from database or model
            tg_value = 150.0  # Mock value in K
        elif smiles:
            # Estimate Tg from SMILES
            # This is a very simplified mock calculation
            if 'O' in smiles:  # Contains oxygen
                tg_value = 180.0  # Higher for compounds with oxygen (like alcohols)
            else:
                tg_value = 150.0  # Default value
        else:
            raise ModelValidationError("Either molecule_id or smiles must be provided")
        
        return {
            "property_name": self.property_name,
            "value": tg_value,
            "units": "K",
            "model_type": self.__class__.__name__,
            "confidence": "medium"  # Indicates this is an estimated value
        }
    
    def _calculate_mixture_tg(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate the glass transition temperature for a mixture.
        
        This uses models like Gordon-Taylor or Fox equation to predict
        Tg of mixtures based on the Tg values of the components.
        
        Args:
            inputs: Dictionary with mixture components
            
        Returns:
            Dictionary with Tg value and metadata
        """
        components = inputs.get('components', [])
        if not components or not isinstance(components, list):
            raise ModelValidationError("Components must be a non-empty list")
        
        # Validate components format
        for comp in components:
            if not isinstance(comp, dict):
                raise ModelValidationError("Each component must be a dictionary")
            
            required_keys = ['weight_fraction', 'tg']
            missing_keys = [key for key in required_keys if key not in comp]
            if missing_keys:
                raise ModelValidationError(f"Missing required component data: {', '.join(missing_keys)}")
            
            try:
                weight_fraction = float(comp['weight_fraction'])
                tg = float(comp['tg'])
                
                if weight_fraction < 0 or weight_fraction > 1:
                    raise ModelValidationError(f"Weight fraction must be between 0 and 1: {weight_fraction}")
                
                if tg <= 0:
                    raise ModelValidationError(f"Glass transition temperature must be positive: {tg}")
            except (ValueError, TypeError):
                raise ModelValidationError("Invalid numeric values in component data")
        
        # Check that weight fractions sum to approximately 1
        total_weight = sum(float(comp['weight_fraction']) for comp in components)
        if abs(total_weight - 1.0) > 0.01:  # Allow small rounding errors
            raise ModelValidationError(f"Weight fractions must sum to 1, got {total_weight}")
        
        # Calculate mixture Tg based on the model subtype
        model_subtype = self.parameters.get('model_subtype', 'gordon_taylor')
        
        if model_subtype == 'gordon_taylor':
            # Gordon-Taylor equation
            k = float(self.parameters.get('k', 1.0))
            
            # Gordon-Taylor equation for binary mixture:
            # Tg = (w1*Tg1 + k*w2*Tg2) / (w1 + k*w2)
            # For multi-component, we extend this formula
            
            numerator = 0.0
            denominator = 0.0
            
            for comp in components:
                w = float(comp['weight_fraction'])
                tg = float(comp['tg'])
                k_value = k  # Use same k for all components
                
                numerator += w * tg
                denominator += w
            
            try:
                mixture_tg = numerator / denominator
            except ZeroDivisionError:
                raise ModelCalculationError("Division by zero in Gordon-Taylor calculation")
            
        elif model_subtype == 'fox':
            # Fox equation
            # 1/Tg = w1/Tg1 + w2/Tg2 + ...
            
            sum_w_div_tg = 0.0
            
            for comp in components:
                w = float(comp['weight_fraction'])
                tg = float(comp['tg'])
                
                sum_w_div_tg += w / tg
            
            try:
                mixture_tg = 1.0 / sum_w_div_tg
            except ZeroDivisionError:
                raise ModelCalculationError("Division by zero in Fox equation calculation")
            
        else:
            raise ModelCalculationError(f"Unknown glass transition model subtype: {model_subtype}")
        
        return {
            "property_name": self.property_name,
            "value": mixture_tg,
            "units": "K",
            "model_type": self.__class__.__name__,
            "model_subtype": model_subtype,
            "components": len(components),
            "confidence": "medium"  # Indicates this is a model prediction
        }
    
    @staticmethod
    def _get_property_units() -> str:
        """Get the units for glass transition temperature."""
        return "K"