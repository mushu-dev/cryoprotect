"""
Enhanced RDKit integration module for scientific models.

This module provides a Python interface to the enhanced RDKit service,
allowing scientific models to leverage advanced molecular modeling capabilities.
"""
import os
import json
import logging
import requests
from typing import Dict, List, Optional, Union, Any, Tuple
import numpy as np
from functools import lru_cache

from .base import ScientificModel, ModelValidationError, ModelCalculationError, ModelParameterError

# Configure logging
logger = logging.getLogger(__name__)

# Constants
DEFAULT_RDKIT_SERVICE_URL = os.environ.get(
    'RDKIT_ENHANCED_SERVICE_URL', 'http://localhost:5001'
)

class EnhancedRDKitError(ModelValidationError):
    """Error raised when the enhanced RDKit service encounters an error."""
    pass

class EnhancedRDKitCalculator:
    """Interface to the enhanced RDKit service."""
    
    def __init__(self, service_url: str = None):
        """Initialize the calculator.
        
        Args:
            service_url: URL of the enhanced RDKit service. Defaults to
                the value of the RDKIT_ENHANCED_SERVICE_URL environment
                variable, or http://localhost:5001 if not set.
        """
        self.service_url = service_url or DEFAULT_RDKIT_SERVICE_URL
        self._check_connection()
        
    def _check_connection(self) -> bool:
        """Check if the enhanced RDKit service is available.
        
        Returns:
            bool: True if the service is available, False otherwise.
        
        Raises:
            EnhancedRDKitError: If the service cannot be reached.
        """
        try:
            response = requests.get(f"{self.service_url}/health", timeout=5)
            if response.status_code != 200:
                raise EnhancedRDKitError(
                    f"RDKit service health check failed with status {response.status_code}"
                )
            return True
        except requests.RequestException as e:
            raise EnhancedRDKitError(f"Cannot connect to RDKit service: {str(e)}")
    
    @lru_cache(maxsize=1)
    def get_service_info(self) -> Dict[str, Any]:
        """Get information about the RDKit service.
        
        Returns:
            Dict[str, Any]: Information about the service.
        
        Raises:
            EnhancedRDKitError: If the service request fails.
        """
        try:
            response = requests.get(f"{self.service_url}/api/v1/rdkit/info", timeout=10)
            if response.status_code != 200:
                raise EnhancedRDKitError(
                    f"Failed to get service info: {response.status_code}"
                )
            return response.json()
        except requests.RequestException as e:
            raise EnhancedRDKitError(f"Service info request failed: {str(e)}")
    
    def standardize_molecule(self, smiles: str) -> Dict[str, str]:
        """Standardize a molecule.
        
        Args:
            smiles: SMILES string of the molecule.
            
        Returns:
            Dict[str, str]: Standardized molecule information including SMILES,
            InChI, and InChIKey.
            
        Raises:
            EnhancedRDKitError: If the standardization fails.
        """
        try:
            response = requests.post(
                f"{self.service_url}/api/v1/rdkit/standardize",
                json={"molecule": smiles},
                timeout=30
            )
            if response.status_code != 200:
                raise EnhancedRDKitError(
                    f"Standardization failed with status {response.status_code}: {response.text}"
                )
            return response.json()
        except requests.RequestException as e:
            raise EnhancedRDKitError(f"Standardization request failed: {str(e)}")
    
    def calculate_descriptors(
        self, 
        smiles: str, 
        include_3d: bool = False
    ) -> Dict[str, float]:
        """Calculate molecular descriptors.
        
        Args:
            smiles: SMILES string of the molecule.
            include_3d: Whether to include 3D descriptors.
            
        Returns:
            Dict[str, float]: Calculated descriptors.
            
        Raises:
            EnhancedRDKitError: If the calculation fails.
        """
        try:
            response = requests.post(
                f"{self.service_url}/api/v1/rdkit/descriptors",
                json={"molecule": smiles, "include_3d": include_3d},
                timeout=60 if include_3d else 30
            )
            if response.status_code != 200:
                raise EnhancedRDKitError(
                    f"Descriptor calculation failed with status {response.status_code}: {response.text}"
                )
            result = response.json()
            return result.get('descriptors', {})
        except requests.RequestException as e:
            raise EnhancedRDKitError(f"Descriptor calculation request failed: {str(e)}")
    
    def generate_conformers(
        self, 
        smiles: str, 
        n_conformers: int = 10
    ) -> List[Dict[str, Any]]:
        """Generate conformers for a molecule.
        
        Args:
            smiles: SMILES string of the molecule.
            n_conformers: Number of conformers to generate.
            
        Returns:
            List[Dict[str, Any]]: Generated conformers with their energies.
            
        Raises:
            EnhancedRDKitError: If the conformer generation fails.
        """
        try:
            response = requests.post(
                f"{self.service_url}/api/v1/rdkit/conformers",
                json={"molecule": smiles, "n_conformers": n_conformers},
                timeout=120
            )
            if response.status_code != 200:
                raise EnhancedRDKitError(
                    f"Conformer generation failed with status {response.status_code}: {response.text}"
                )
            result = response.json()
            return result.get('conformers', [])
        except requests.RequestException as e:
            raise EnhancedRDKitError(f"Conformer generation request failed: {str(e)}")
    
    def calculate_3d_similarity(
        self, 
        query_smiles: str, 
        target_smiles: str
    ) -> float:
        """Calculate 3D similarity between two molecules.
        
        Args:
            query_smiles: SMILES string of the query molecule.
            target_smiles: SMILES string of the target molecule.
            
        Returns:
            float: 3D similarity score.
            
        Raises:
            EnhancedRDKitError: If the similarity calculation fails.
        """
        try:
            response = requests.post(
                f"{self.service_url}/api/v1/rdkit/similarity/3d",
                json={"query": query_smiles, "target": target_smiles},
                timeout=120
            )
            if response.status_code != 200:
                raise EnhancedRDKitError(
                    f"3D similarity calculation failed with status {response.status_code}: {response.text}"
                )
            result = response.json()
            return result.get('similarity_score', 0.0)
        except requests.RequestException as e:
            raise EnhancedRDKitError(f"3D similarity calculation request failed: {str(e)}")
    
    def get_pharmacophore(self, smiles: str) -> Dict[str, Any]:
        """Get pharmacophore features for a molecule.
        
        Args:
            smiles: SMILES string of the molecule.
            
        Returns:
            Dict[str, Any]: Pharmacophore features.
            
        Raises:
            EnhancedRDKitError: If the pharmacophore calculation fails.
        """
        try:
            response = requests.post(
                f"{self.service_url}/api/v1/rdkit/pharmacophore",
                json={"molecule": smiles},
                timeout=60
            )
            if response.status_code != 200:
                raise EnhancedRDKitError(
                    f"Pharmacophore calculation failed with status {response.status_code}: {response.text}"
                )
            result = response.json()
            return result.get('pharmacophore_features', [])
        except requests.RequestException as e:
            raise EnhancedRDKitError(f"Pharmacophore calculation request failed: {str(e)}")
    
    def batch_calculate_descriptors(
        self, 
        smiles_list: List[str], 
        include_3d: bool = False
    ) -> List[Dict[str, Any]]:
        """Calculate descriptors for multiple molecules in batch.
        
        Args:
            smiles_list: List of SMILES strings.
            include_3d: Whether to include 3D descriptors.
            
        Returns:
            List[Dict[str, Any]]: Calculated descriptors for each molecule.
            
        Raises:
            EnhancedRDKitError: If the batch calculation fails.
        """
        try:
            response = requests.post(
                f"{self.service_url}/api/v1/rdkit/batch",
                json={
                    "molecules": smiles_list, 
                    "operation": "descriptors",
                    "include_3d": include_3d
                },
                timeout=300 if include_3d else 120
            )
            if response.status_code != 200:
                raise EnhancedRDKitError(
                    f"Batch descriptor calculation failed with status {response.status_code}: {response.text}"
                )
            result = response.json()
            return result.get('results', [])
        except requests.RequestException as e:
            raise EnhancedRDKitError(f"Batch descriptor calculation request failed: {str(e)}")
    
    def minimize_molecule(self, smiles: str) -> Dict[str, Any]:
        """Energy minimize a molecule.
        
        Args:
            smiles: SMILES string of the molecule.
            
        Returns:
            Dict[str, Any]: Minimized molecule with energy information.
            
        Raises:
            EnhancedRDKitError: If the minimization fails.
        """
        try:
            response = requests.post(
                f"{self.service_url}/api/v1/rdkit/minimize",
                json={"molecule": smiles},
                timeout=120
            )
            if response.status_code != 200:
                raise EnhancedRDKitError(
                    f"Molecule minimization failed with status {response.status_code}: {response.text}"
                )
            return response.json()
        except requests.RequestException as e:
            raise EnhancedRDKitError(f"Molecule minimization request failed: {str(e)}")
    
    def predict_property(
        self, 
        property_name: str, 
        smiles: str,
        use_3d: bool = True
    ) -> float:
        """Predict a molecular property using descriptor-based models.
        
        This method leverages the descriptors calculated by the enhanced RDKit
        service to predict properties using pre-trained models.
        
        Args:
            property_name: Name of the property to predict.
                Supported properties depend on the models available in the service.
            smiles: SMILES string of the molecule.
            use_3d: Whether to use 3D descriptors for the prediction.
            
        Returns:
            float: Predicted property value.
            
        Raises:
            EnhancedRDKitError: If the property prediction fails.
            ValueError: If the property is not supported.
        """
        # For now, we'll implement a simple mapping of properties to descriptors
        # In a real implementation, this would use machine learning models
        
        # Calculate descriptors
        descriptors = self.calculate_descriptors(smiles, include_3d=use_3d)
        
        # Define property mappings (placeholder implementation)
        property_mappings = {
            "solubility": lambda d: d.get("MolLogP", 0) * -0.5 + 2.3,
            "permeability": lambda d: 0.2 * d.get("MolLogP", 0) + 0.1 * d.get("NumRotatableBonds", 0),
            "glass_transition": lambda d: 150 - 2 * d.get("NumRotatableBonds", 0) + 0.5 * d.get("NumHDonors", 0),
            "toxicity": lambda d: 0.3 * d.get("MolLogP", 0) + 0.1 * d.get("NumAromaticRings", 0),
            "stability": lambda d: 100 - 2 * d.get("NumRotatableBonds", 0) + 5 * d.get("NumRings", 0)
        }
        
        if property_name not in property_mappings:
            raise ValueError(f"Property '{property_name}' is not supported")
            
        # Apply the property mapping
        return property_mappings[property_name](descriptors)
        
class DescriptorBasedModel(ScientificModel):
    """Scientific model based on molecular descriptors.
    
    This model uses molecular descriptors calculated by the enhanced RDKit service
    to predict properties of interest for cryoprotectants.
    """
    
    def __init__(self, model_params: Dict[str, Any] = None, name: str = None, description: str = None):
        """Initialize the model.
        
        Args:
            model_params: Model parameters.
            name: Optional name for the model.
            description: Optional description of the model.
        """
        super().__init__(model_params or {}, name, description)
        self.calculator = EnhancedRDKitCalculator()
    
    def validate_parameters(self) -> None:
        """Validate model parameters.
        
        Raises:
            ModelParameterError: If parameters are invalid.
        """
        # This base implementation has no special parameters to validate
        pass
        
    def validate_inputs(self, inputs: Dict[str, Any]) -> None:
        """Validate model inputs.
        
        Args:
            inputs: Model inputs.
            
        Raises:
            ModelValidationError: If the inputs are invalid.
        """
        if 'smiles' not in inputs:
            raise ModelValidationError("SMILES string is required")
            
        # Check if the SMILES string is valid by standardizing it
        try:
            self.calculator.standardize_molecule(inputs['smiles'])
        except EnhancedRDKitError as e:
            raise ModelValidationError(f"Invalid SMILES string: {str(e)}")
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate model outputs.
        
        Args:
            inputs: Model inputs.
            
        Returns:
            Dict[str, Any]: Model outputs.
            
        Raises:
            ModelCalculationError: If the calculation fails.
        """
        # This is a base implementation that should be overridden by subclasses
        self.validate_inputs(inputs)
        
        # Calculate descriptors
        try:
            descriptors = self.calculator.calculate_descriptors(
                inputs['smiles'],
                include_3d=self.parameters.get('use_3d', False)
            )
            return {'descriptors': descriptors}
        except EnhancedRDKitError as e:
            raise ModelCalculationError(f"Descriptor calculation failed: {str(e)}")

class MolecularSimilarityModel(DescriptorBasedModel):
    """Model for calculating molecular similarity.
    
    This model uses the enhanced RDKit service to calculate similarity
    between molecules using both 2D and 3D methods.
    """
    
    def validate_inputs(self, inputs: Dict[str, Any]) -> None:
        """Validate model inputs.
        
        Args:
            inputs: Model inputs.
            
        Raises:
            ModelValidationError: If the inputs are invalid.
        """
        if 'query_smiles' not in inputs:
            raise ModelValidationError("Query SMILES string is required")
            
        if 'target_smiles' not in inputs and 'target_smiles_list' not in inputs:
            raise ModelValidationError("Target SMILES string or list is required")
            
        # Check if the SMILES strings are valid
        try:
            self.calculator.standardize_molecule(inputs['query_smiles'])
            
            if 'target_smiles' in inputs:
                self.calculator.standardize_molecule(inputs['target_smiles'])
            elif 'target_smiles_list' in inputs:
                if not isinstance(inputs['target_smiles_list'], list):
                    raise ModelValidationError("target_smiles_list must be a list")
                
                if len(inputs['target_smiles_list']) == 0:
                    raise ModelValidationError("target_smiles_list cannot be empty")
                    
                # Validate only first 3 for efficiency
                for smiles in inputs['target_smiles_list'][:3]:
                    self.calculator.standardize_molecule(smiles)
        except EnhancedRDKitError as e:
            raise ModelValidationError(f"Invalid SMILES string: {str(e)}")
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate similarity between molecules.
        
        Args:
            inputs: Model inputs including query_smiles and target_smiles or target_smiles_list.
            
        Returns:
            Dict[str, Any]: Model outputs including similarity scores.
            
        Raises:
            ModelCalculationError: If the calculation fails.
        """
        self.validate_inputs(inputs)
        
        # Determine similarity mode
        mode = self.parameters.get('similarity_mode', '3d')
        
        try:
            if 'target_smiles' in inputs:
                # Single target comparison
                if mode == '3d':
                    similarity = self.calculator.calculate_3d_similarity(
                        inputs['query_smiles'],
                        inputs['target_smiles']
                    )
                    return {'similarity_score': similarity}
                else:
                    # Not implemented in this example
                    raise ModelValidationError(f"Similarity mode '{mode}' not implemented")
            else:
                # Multiple target comparison
                results = []
                for target_smiles in inputs['target_smiles_list']:
                    if mode == '3d':
                        similarity = self.calculator.calculate_3d_similarity(
                            inputs['query_smiles'],
                            target_smiles
                        )
                    else:
                        # Not implemented in this example
                        raise ModelValidationError(f"Similarity mode '{mode}' not implemented")
                        
                    results.append({
                        'target_smiles': target_smiles,
                        'similarity_score': similarity
                    })
                
                # Sort by similarity score in descending order
                results.sort(key=lambda x: x['similarity_score'], reverse=True)
                
                return {'results': results}
        except EnhancedRDKitError as e:
            raise ModelCalculationError(f"Similarity calculation failed: {str(e)}")

class ConformerAnalysisModel(DescriptorBasedModel):
    """Model for analyzing molecular conformers.
    
    This model generates and analyzes conformers for a molecule to predict
    properties related to flexibility and stability.
    """
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Generate and analyze conformers for a molecule.
        
        Args:
            inputs: Model inputs including SMILES.
            
        Returns:
            Dict[str, Any]: Model outputs including conformer analysis.
            
        Raises:
            ModelCalculationError: If the calculation fails.
        """
        self.validate_inputs(inputs)
        
        n_conformers = self.parameters.get('n_conformers', 10)
        
        try:
            # Generate conformers
            conformers = self.calculator.generate_conformers(
                inputs['smiles'],
                n_conformers=n_conformers
            )
            
            # Calculate statistics
            if conformers:
                energies = [conf['energy'] for conf in conformers]
                energy_min = min(energies)
                energy_max = max(energies)
                energy_range = energy_max - energy_min
                energy_avg = sum(energies) / len(energies)
                
                # Calculate metrics
                flexibility_score = energy_range / max(1.0, len(conformers))
                stability_score = 100.0 / (1.0 + energy_min / 100.0)
                
                return {
                    'n_conformers': len(conformers),
                    'energy_min': energy_min,
                    'energy_max': energy_max,
                    'energy_range': energy_range,
                    'energy_avg': energy_avg,
                    'flexibility_score': flexibility_score,
                    'stability_score': stability_score,
                    'conformer_details': conformers[:3]  # Return only first 3 for brevity
                }
            else:
                raise ModelValidationError("No conformers generated")
        except EnhancedRDKitError as e:
            raise ModelCalculationError(f"Conformer analysis failed: {str(e)}")

class PharmacophoreModel(DescriptorBasedModel):
    """Model for pharmacophore analysis.
    
    This model analyzes the pharmacophore features of a molecule to predict
    interactions with biological targets.
    """
    
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Analyze pharmacophore features of a molecule.
        
        Args:
            inputs: Model inputs including SMILES.
            
        Returns:
            Dict[str, Any]: Model outputs including pharmacophore analysis.
            
        Raises:
            ModelCalculationError: If the calculation fails.
        """
        self.validate_inputs(inputs)
        
        try:
            # Get pharmacophore features
            features = self.calculator.get_pharmacophore(inputs['smiles'])
            
            # Count feature types
            feature_counts = {}
            for feature in features:
                feature_type = feature['type']
                feature_counts[feature_type] = feature_counts.get(feature_type, 0) + 1
            
            # Calculate derived metrics
            total_features = len(features)
            
            # Calculate feature type percentages
            feature_percentages = {
                k: 100.0 * v / total_features if total_features > 0 else 0.0
                for k, v in feature_counts.items()
            }
            
            # Calculate a simple score based on feature counts
            # This is a placeholder for a more sophisticated model
            interaction_potential = sum([
                feature_counts.get('Donor', 0) * 1.5,
                feature_counts.get('Acceptor', 0) * 1.2,
                feature_counts.get('Hydrophobe', 0) * 0.8,
                feature_counts.get('Aromatic', 0) * 0.7,
                feature_counts.get('PosIonizable', 0) * 1.0,
                feature_counts.get('NegIonizable', 0) * 1.0
            ])
            
            return {
                'total_features': total_features,
                'feature_counts': feature_counts,
                'feature_percentages': feature_percentages,
                'interaction_potential': interaction_potential,
                'features': features
            }
        except EnhancedRDKitError as e:
            raise ModelCalculationError(f"Pharmacophore analysis failed: {str(e)}")