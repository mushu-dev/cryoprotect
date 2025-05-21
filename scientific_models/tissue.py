"""
Tissue Compatibility Models

This module contains models for predicting compatibility of cryoprotectants
with different tissue types, including penetration and toxicity models.
"""

import math
import numpy as np
from typing import Dict, List, Optional, Any, Union

from scientific_models.base import ScientificModel, ModelValidationError


class TissueCompatibilityModel(ScientificModel):
    """Base class for tissue compatibility models.
    
    This model predicts how well a cryoprotectant or mixture will work
    with a specific tissue type, considering various factors such as
    penetration, toxicity, and overall effectiveness.
    """
    
    def __init__(self, model_params: Dict[str, Any] = None):
        """Initialize the model.
        
        Args:
            model_params: Model parameters.
        """
        super().__init__(model_params or {})
        
        # Default tissue types if not specified
        self.tissue_types = self.model_params.get('tissue_types', [
            'liver', 'kidney', 'heart', 'lung', 'pancreas', 
            'skin', 'cornea', 'nerve', 'bone', 'muscle'
        ])
    
    def validate_inputs(self, inputs: Dict[str, Any]) -> None:
        """Validate model inputs.
        
        Args:
            inputs: Model inputs including molecule and tissue_type.
            
        Raises:
            ModelValidationError: If inputs are invalid.
        """
        if 'molecule_id' not in inputs and 'smiles' not in inputs:
            raise ModelValidationError("Either molecule_id or SMILES string is required")
            
        if 'tissue_type' not in inputs:
            raise ModelValidationError("Tissue type is required")
            
        tissue_type = inputs['tissue_type']
        if tissue_type not in self.tissue_types:
            raise ModelValidationError(
                f"Invalid tissue type: {tissue_type}. "
                f"Must be one of: {', '.join(self.tissue_types)}"
            )
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate tissue compatibility.
        
        This is a base implementation that should be overridden by subclasses.
        
        Args:
            inputs: Model inputs including molecule and tissue_type.
            
        Returns:
            Dict[str, Any]: Model outputs including compatibility score.
        """
        self.validate_inputs(inputs)
        
        # This base implementation returns a random compatibility score
        # Subclasses should implement more sophisticated models
        tissue_type = inputs['tissue_type']
        compatibility = np.random.uniform(0.5, 1.0)
        
        return {
            'tissue_type': tissue_type,
            'compatibility_score': compatibility,
            'confidence': 0.5,  # Low confidence for random score
            'notes': "Base implementation with random score"
        }


class PenetrationModel(TissueCompatibilityModel):
    """Model for predicting penetration of cryoprotectants into tissues.
    
    This model estimates how well a cryoprotectant can penetrate into
    a specific tissue type, which is a critical factor for effectiveness.
    """
    
    def __init__(self, model_params: Dict[str, Any] = None):
        """Initialize the model.
        
        Args:
            model_params: Model parameters.
        """
        super().__init__(model_params or {})
        
        # Parameters related to penetration
        self.molecular_weight_factor = self.model_params.get('molecular_weight_factor', 0.005)
        self.hydrophobicity_factor = self.model_params.get('hydrophobicity_factor', 0.2)
        self.tissue_permeability = {
            'liver': 0.8,
            'kidney': 0.7,
            'heart': 0.6,
            'lung': 0.9,
            'pancreas': 0.7,
            'skin': 0.4,
            'cornea': 0.9,
            'nerve': 0.5,
            'bone': 0.2,
            'muscle': 0.6
        }
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate penetration rate and effectiveness.
        
        Args:
            inputs: Model inputs including molecule properties and tissue_type.
            
        Returns:
            Dict[str, Any]: Model outputs including penetration rate.
        """
        self.validate_inputs(inputs)
        
        tissue_type = inputs['tissue_type']
        
        # Get molecular properties (in a real implementation, these would be
        # retrieved from a database or calculated using RDKit)
        mol_weight = inputs.get('molecular_weight', 200)
        log_p = inputs.get('logP', 0.0)
        
        # Calculate penetration rate
        size_factor = math.exp(-self.molecular_weight_factor * mol_weight)
        hydrophobicity_contribution = 1.0 - abs(log_p - 0.5) * self.hydrophobicity_factor
        tissue_factor = self.tissue_permeability.get(tissue_type, 0.5)
        
        penetration_rate = size_factor * hydrophobicity_contribution * tissue_factor
        
        # Ensure rate is between 0 and 1
        penetration_rate = max(0.0, min(1.0, penetration_rate))
        
        return {
            'tissue_type': tissue_type,
            'penetration_rate': penetration_rate,
            'contributions': {
                'size_factor': size_factor,
                'hydrophobicity_contribution': hydrophobicity_contribution,
                'tissue_factor': tissue_factor
            },
            'confidence': 0.7,
            'notes': "Model based on molecular weight and hydrophobicity"
        }


class TissueToxicityModel(TissueCompatibilityModel):
    """Model for predicting toxicity of cryoprotectants to tissues.
    
    This model estimates how toxic a cryoprotectant is to a specific
    tissue type, which is important for safety and effectiveness.
    """
    
    def __init__(self, model_params: Dict[str, Any] = None):
        """Initialize the model.
        
        Args:
            model_params: Model parameters.
        """
        super().__init__(model_params or {})
        
        # Parameters related to toxicity
        self.concentration_factor = self.model_params.get('concentration_factor', 0.5)
        self.exposure_time_factor = self.model_params.get('exposure_time_factor', 0.1)
        self.tissue_sensitivity = {
            'liver': 0.6,
            'kidney': 0.8,
            'heart': 0.7,
            'lung': 0.6,
            'pancreas': 0.7,
            'skin': 0.4,
            'cornea': 0.5,
            'nerve': 0.9,
            'bone': 0.3,
            'muscle': 0.5
        }
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate toxicity for a tissue.
        
        Args:
            inputs: Model inputs including molecule, concentration, exposure_time,
                and tissue_type.
            
        Returns:
            Dict[str, Any]: Model outputs including toxicity score.
        """
        self.validate_inputs(inputs)
        
        tissue_type = inputs['tissue_type']
        concentration = inputs.get('concentration', 1.0)  # Default to 1.0 M
        exposure_time = inputs.get('exposure_time', 10.0)  # Default to 10 minutes
        
        # Get intrinsic toxicity (in a real implementation, this would be
        # retrieved from a database or calculated using QSAR models)
        intrinsic_toxicity = inputs.get('intrinsic_toxicity', 0.3)
        
        # Calculate toxicity score
        concentration_contribution = math.pow(concentration, self.concentration_factor)
        time_contribution = math.pow(exposure_time / 10.0, self.exposure_time_factor)
        tissue_factor = self.tissue_sensitivity.get(tissue_type, 0.5)
        
        toxicity_score = intrinsic_toxicity * concentration_contribution * time_contribution * tissue_factor
        
        # Ensure score is between 0 and 1
        toxicity_score = max(0.0, min(1.0, toxicity_score))
        
        return {
            'tissue_type': tissue_type,
            'toxicity_score': toxicity_score,
            'safe_concentration': max(0.1, 1.0 - toxicity_score),
            'contributions': {
                'intrinsic_toxicity': intrinsic_toxicity,
                'concentration_contribution': concentration_contribution,
                'time_contribution': time_contribution,
                'tissue_factor': tissue_factor
            },
            'confidence': 0.7,
            'notes': "Model based on concentration, exposure time, and tissue sensitivity"
        }