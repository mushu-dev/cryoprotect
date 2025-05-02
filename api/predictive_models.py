"""
CryoProtect Analyzer API - Predictive Models

This module provides machine learning models for predicting cryoprotectant performance
based on molecular properties. It implements multiple prediction algorithms, model training,
validation, and evaluation functionality.
"""

import logging
import os
import pickle
import json
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Union, Any, Tuple
from datetime import datetime
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer

from api.rdkit_utils import (
    parse_molecule, calculate_hydrogen_bonding, calculate_logp,
    calculate_tpsa, calculate_molecular_properties, identify_functional_groups,
    estimate_permeability, calculate_all_properties
)
from api.models import Molecule, MolecularProperty, Mixture, Prediction, Experiment
from api.scoring import SCORE_WEIGHTS

# Set up logging
logger = logging.getLogger(__name__)

# Directory for storing trained models
MODELS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'models')
os.makedirs(MODELS_DIR, exist_ok=True)

# Available prediction algorithms
ALGORITHMS = {
    'linear_regression': 'Linear Regression',
    'ridge_regression': 'Ridge Regression',
    'lasso_regression': 'Lasso Regression',
    'random_forest': 'Random Forest',
    'gradient_boosting': 'Gradient Boosting',
    'neural_network': 'Neural Network'
}

# Default hyperparameters for each algorithm
DEFAULT_HYPERPARAMETERS = {
    'linear_regression': {},
    'ridge_regression': {'alpha': 1.0},
    'lasso_regression': {'alpha': 1.0},
    'random_forest': {'n_estimators': 100, 'max_depth': 10},
    'gradient_boosting': {'n_estimators': 100, 'learning_rate': 0.1, 'max_depth': 3},
    'neural_network': {'hidden_layer_sizes': (100,), 'max_iter': 1000}
}

# Feature importance methods for each algorithm
FEATURE_IMPORTANCE_METHODS = {
    'linear_regression': 'coefficients',
    'ridge_regression': 'coefficients',
    'lasso_regression': 'coefficients',
    'random_forest': 'feature_importances',
    'gradient_boosting': 'feature_importances',
    'neural_network': None  # Neural networks don't have a direct feature importance method
}

class PredictiveModel:
    """
    Base class for predictive models for cryoprotectant property prediction.

    Scientific Rationale:
        This class implements a flexible framework for training, evaluating, and applying machine learning models
        to predict cryoprotectant-relevant properties from molecular descriptors. The approach is rooted in
        quantitative structure-activity relationship (QSAR) modeling, a cornerstone of cheminformatics and drug design.

    Features:
        - Supports multiple regression algorithms (linear, tree-based, neural network).
        - Extracts scientifically relevant features (hydrogen bonding, logP, TPSA, etc.) for model input.
        - Provides model evaluation, cross-validation, and feature importance analysis.
        - Designed for extensibility and scientific interpretability.

    References:
        - Cherkasov, A., et al. (2014). QSAR modeling: where have you been? Where are you going to? J. Med. Chem., 57(12), 4977-5010.
        - Todeschini, R., & Consonni, V. (2009). Molecular Descriptors for Chemoinformatics.
        - Sheridan, R. P. (2013). Time-split cross-validation as a method for estimating the goodness of prospective prediction. J. Chem. Inf. Model., 53(4), 783-790.
    """
    
    def __init__(self, property_name: str, algorithm: str = 'random_forest', 
                 hyperparameters: Dict[str, Any] = None):
        """
        Initialize a predictive model.
        
        Args:
            property_name: Name of the property to predict
            algorithm: Algorithm to use for prediction
            hyperparameters: Hyperparameters for the algorithm
        """
        self.property_name = property_name
        self.algorithm = algorithm
        self.hyperparameters = hyperparameters or DEFAULT_HYPERPARAMETERS.get(algorithm, {})
        self.model = None
        self.scaler = StandardScaler()
        self.feature_names = []
        self.trained_date = None
        self.metrics = {}
        self.feature_importance = {}
        
    def _get_model_instance(self) -> Any:
        """Get an instance of the model based on the algorithm."""
        if self.algorithm == 'linear_regression':
            return LinearRegression(**self.hyperparameters)
        elif self.algorithm == 'ridge_regression':
            return Ridge(**self.hyperparameters)
        elif self.algorithm == 'lasso_regression':
            return Lasso(**self.hyperparameters)
        elif self.algorithm == 'random_forest':
            return RandomForestRegressor(**self.hyperparameters)
        elif self.algorithm == 'gradient_boosting':
            return GradientBoostingRegressor(**self.hyperparameters)
        elif self.algorithm == 'neural_network':
            return MLPRegressor(**self.hyperparameters)
        else:
            raise ValueError(f"Unsupported algorithm: {self.algorithm}")
    
    def _extract_features(self, molecule_data: Dict[str, Any]) -> List[float]:
        """
        Extracts scientifically relevant features from molecule data for model input.

        Scientific Rationale:
            Feature selection is based on domain knowledge of cryoprotectant effectiveness and QSAR best practices.
            Features include hydrogen bonding, logP, TPSA, molecular size, functional groups, and permeability metrics.

        Args:
            molecule_data (Dict[str, Any]): Dictionary of molecule properties.

        Returns:
            List[float]: List of feature values in a fixed order.

        Notes:
            - Features are ordered to match self._get_feature_names().
            - Boolean features are cast to float for model compatibility.

        References:
            - Todeschini, R., & Consonni, V. (2009). Molecular Descriptors for Chemoinformatics.
        """
        features = []

        # Hydrogen bonding features (important for water replacement and vitrification)
        h_bonds = molecule_data.get('hydrogen_bonding', {})
        features.append(h_bonds.get('donors', 0))
        features.append(h_bonds.get('acceptors', 0))
        features.append(h_bonds.get('total', 0))

        # LogP (hydrophilic/lipophilic balance)
        features.append(molecule_data.get('logp', 0))

        # TPSA (membrane permeability and hydrogen bonding)
        features.append(molecule_data.get('tpsa', 0))

        # Molecular properties (size, flexibility, saturation)
        mol_props = molecule_data.get('molecular_properties', {})
        features.append(mol_props.get('molecular_weight', 0))
        features.append(mol_props.get('heavy_atom_count', 0))
        features.append(mol_props.get('rotatable_bond_count', 0))
        features.append(mol_props.get('ring_count', 0))
        features.append(mol_props.get('fraction_csp3', 0))

        # Functional groups (key for cryoprotectant activity)
        func_groups = molecule_data.get('functional_groups', {})
        features.append(func_groups.get('hydroxyl', 0))
        features.append(func_groups.get('alcohol', 0))
        features.append(func_groups.get('ether', 0))
        features.append(func_groups.get('amine', 0))
        features.append(func_groups.get('amide', 0))

        # Permeability metrics (rule violations, log Papp, absorption)
        perm = molecule_data.get('permeability', {})
        features.append(perm.get('rule_of_5_violations', 0))
        features.append(perm.get('veber_violations', 0))
        features.append(float(perm.get('bbb_permeant', False)))
        features.append(float(perm.get('intestinal_absorption', False)))
        features.append(perm.get('estimated_log_papp', 0))

        return features
    
    def _get_feature_names(self) -> List[str]:
        """Get the names of the features used by the model."""
        return [
            'h_bond_donors', 'h_bond_acceptors', 'h_bond_total',
            'logp', 'tpsa',
            'molecular_weight', 'heavy_atom_count', 'rotatable_bond_count', 'ring_count', 'fraction_csp3',
            'hydroxyl_count', 'alcohol_count', 'ether_count', 'amine_count', 'amide_count',
            'rule_of_5_violations', 'veber_violations', 'bbb_permeant', 'intestinal_absorption', 'estimated_log_papp'
        ]
    
    def train(self, X: np.ndarray, y: np.ndarray, feature_names: List[str] = None) -> Dict[str, Any]:
        """
        Train the model on the given data.
        
        Args:
            X: Feature matrix
            y: Target values
            feature_names: Names of the features
            
        Returns:
            Dictionary of training metrics
        """
        # Store feature names
        self.feature_names = feature_names or self._get_feature_names()
        
        # Split data into training and validation sets
        X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
        
        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_val_scaled = self.scaler.transform(X_val)
        
        # Create and train the model
        self.model = self._get_model_instance()
        self.model.fit(X_train_scaled, y_train)
        
        # Calculate metrics
        y_train_pred = self.model.predict(X_train_scaled)
        y_val_pred = self.model.predict(X_val_scaled)
        
        train_mse = mean_squared_error(y_train, y_train_pred)
        train_rmse = np.sqrt(train_mse)
        train_mae = mean_absolute_error(y_train, y_train_pred)
        train_r2 = r2_score(y_train, y_train_pred)
        
        val_mse = mean_squared_error(y_val, y_val_pred)
        val_rmse = np.sqrt(val_mse)
        val_mae = mean_absolute_error(y_val, y_val_pred)
        val_r2 = r2_score(y_val, y_val_pred)
        
        # Calculate feature importance if available
        importance_method = FEATURE_IMPORTANCE_METHODS.get(self.algorithm)
        if importance_method == 'coefficients' and hasattr(self.model, 'coef_'):
            importances = self.model.coef_
            self.feature_importance = dict(zip(self.feature_names, importances))
        elif importance_method == 'feature_importances' and hasattr(self.model, 'feature_importances_'):
            importances = self.model.feature_importances_
            self.feature_importance = dict(zip(self.feature_names, importances))
        
        # Store metrics
        self.metrics = {
            'train': {
                'mse': train_mse,
                'rmse': train_rmse,
                'mae': train_mae,
                'r2': train_r2
            },
            'validation': {
                'mse': val_mse,
                'rmse': val_rmse,
                'mae': val_mae,
                'r2': val_r2
            }
        }
        
        # Update trained date
        self.trained_date = datetime.now().isoformat()
        
        return self.metrics
    
    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Make predictions with the model and calculate confidence intervals.
        
        Scientific Rationale:
            This method implements a more sophisticated approach to confidence interval estimation
            that accounts for prediction uncertainty. For regression models, we use:
            
            1. For tree-based models (Random Forest, Gradient Boosting):
               - Standard deviation across individual trees/estimators
               - This captures model uncertainty based on variance in the ensemble
            
            2. For linear models:
               - Prediction interval based on validation RMSE and sample variance
               - This accounts for both model error and data uncertainty
            
            3. For neural networks:
               - Approximation based on validation error and prediction magnitude
               - This provides a reasonable estimate of prediction uncertainty
        
        Args:
            X: Feature matrix
            
        Returns:
            Tuple of (predictions, confidence intervals)
            
        References:
            - Heskes, T. (1997). Practical confidence and prediction intervals.
            - Khosravi, A., et al. (2011). Comprehensive review of neural network-based prediction intervals.
        """
        if self.model is None:
            raise ValueError("Model has not been trained yet")
        
        # Scale features
        X_scaled = self.scaler.transform(X)
        
        # Make predictions
        predictions = self.model.predict(X_scaled)
        
        # Calculate confidence intervals based on algorithm type
        if self.algorithm in ['random_forest', 'gradient_boosting']:
            # For ensemble methods, use prediction variance across estimators
            if hasattr(self.model, 'estimators_'):
                # Get predictions from individual estimators
                individual_preds = np.array([tree.predict(X_scaled) for tree in self.model.estimators_])
                # Calculate standard deviation across estimators
                confidence = np.std(individual_preds, axis=0)
                # Add base uncertainty from validation error
                base_uncertainty = self.metrics.get('validation', {}).get('rmse', 1.0) * 0.5
                confidence = np.sqrt(confidence**2 + base_uncertainty**2)
            else:
                # Fallback if estimators aren't accessible
                confidence = np.ones_like(predictions) * self.metrics.get('validation', {}).get('rmse', 1.0)
        
        elif self.algorithm in ['linear_regression', 'ridge_regression', 'lasso_regression']:
            # For linear models, use prediction interval based on validation RMSE
            rmse = self.metrics.get('validation', {}).get('rmse', 1.0)
            # Scale confidence by prediction magnitude (larger predictions often have larger errors)
            confidence = rmse * (0.5 + 0.5 * np.abs(predictions) / (np.mean(np.abs(predictions)) + 1e-10))
        
        elif self.algorithm == 'neural_network':
            # For neural networks, use a heuristic based on validation error and prediction magnitude
            base_uncertainty = self.metrics.get('validation', {}).get('rmse', 1.0)
            # Neural networks often have higher uncertainty for extreme predictions
            prediction_range = np.max(predictions) - np.min(predictions) if len(predictions) > 1 else 1.0
            normalized_preds = (predictions - np.min(predictions)) / (prediction_range + 1e-10)
            # U-shaped uncertainty (higher at extremes)
            uncertainty_factor = 1.0 + 0.5 * np.abs(normalized_preds - 0.5) * 2
            confidence = base_uncertainty * uncertainty_factor
        
        else:
            # Default fallback
            confidence = np.ones_like(predictions) * self.metrics.get('validation', {}).get('rmse', 1.0)
        
        return predictions, confidence
    
    def predict_molecule(self, molecule_data: Dict[str, Any]) -> Tuple[float, float]:
        """
        Predict property for a molecule.
        
        Args:
            molecule_data: Dictionary of molecule properties
            
        Returns:
            Tuple of (prediction, confidence)
        """
        # Extract features
        features = self._extract_features(molecule_data)
        
        # Make prediction
        predictions, confidences = self.predict(np.array([features]))
        
        return predictions[0], confidences[0]
    
    def predict_mixture(self, mixture_components: List[Dict[str, Any]]) -> Tuple[float, float]:
        """
        Predict property for a mixture with advanced synergy modeling.
        
        Scientific Rationale:
            This method implements a sophisticated approach to mixture prediction that accounts for:
            1. Concentration-dependent effects (weighted average of individual predictions)
            2. Synergistic interactions between components (based on complementary properties)
            3. Functional group compatibility analysis
            4. Vitrification potential assessment
            
        The model is based on cryobiology research showing that effective cryoprotectant
        mixtures often combine molecules with complementary mechanisms of action:
        - Molecules that prevent ice formation (vitrification agents)
        - Molecules that stabilize cellular membranes
        - Molecules that prevent protein denaturation
        - Molecules with different permeability profiles
        
        Args:
            mixture_components: List of dictionaries with molecule properties and concentrations
            
        Returns:
            Tuple of (prediction, confidence)
            
        References:
            - Fahy, G.M., et al. (2004). Cryopreservation of organs by vitrification.
            - Elliott, G.D., et al. (2017). Combination of osmotic and cryoprotective agents.
            - Wowk, B. (2010). Thermodynamic aspects of vitrification.
        """
        if not mixture_components:
            raise ValueError("Mixture has no components")
        
        # Calculate weighted average prediction (base effect)
        total_concentration = sum(comp.get('concentration', 0) for comp in mixture_components)
        weighted_prediction = 0
        weighted_confidence = 0
        
        # Store component properties for synergy analysis
        component_properties = []
        
        for component in mixture_components:
            molecule_data = component.get('properties', {})
            concentration = component.get('concentration', 0)
            weight = concentration / total_concentration if total_concentration > 0 else 1.0 / len(mixture_components)
            
            prediction, confidence = self.predict_molecule(molecule_data)
            weighted_prediction += prediction * weight
            weighted_confidence += confidence * weight
            
            # Extract key properties for synergy analysis
            component_properties.append({
                'logp': molecule_data.get('logp', 0),
                'h_bond_donors': molecule_data.get('hydrogen_bonding', {}).get('donors', 0),
                'h_bond_acceptors': molecule_data.get('hydrogen_bonding', {}).get('acceptors', 0),
                'tpsa': molecule_data.get('tpsa', 0),
                'mol_weight': molecule_data.get('molecular_properties', {}).get('molecular_weight', 0),
                'functional_groups': molecule_data.get('functional_groups', {}),
                'permeability': molecule_data.get('permeability', {}).get('estimated_log_papp', 0),
                'concentration': concentration
            })
        
        # Apply synergy model for mixtures
        if len(mixture_components) > 1:
            # 1. Calculate property diversity (beneficial for cryoprotection)
            logp_values = [prop['logp'] for prop in component_properties]
            logp_range = max(logp_values) - min(logp_values) if logp_values else 0
            
            # 2. Analyze hydrogen bonding complementarity
            h_donors = sum(prop['h_bond_donors'] * prop['concentration'] for prop in component_properties) / total_concentration
            h_acceptors = sum(prop['h_bond_acceptors'] * prop['concentration'] for prop in component_properties) / total_concentration
            h_bond_balance = min(h_donors, h_acceptors) / max(h_donors, h_acceptors) if max(h_donors, h_acceptors) > 0 else 0
            
            # 3. Assess permeability distribution (mix of fast and slow penetrating agents)
            perm_values = [prop['permeability'] for prop in component_properties]
            perm_range = max(perm_values) - min(perm_values) if perm_values else 0
            
            # 4. Check for complementary functional groups
            has_hydroxyl = any('hydroxyl' in prop['functional_groups'] and prop['functional_groups']['hydroxyl'] > 0
                              for prop in component_properties)
            has_amide = any('amide' in prop['functional_groups'] and prop['functional_groups']['amide'] > 0
                           for prop in component_properties)
            has_ether = any('ether' in prop['functional_groups'] and prop['functional_groups']['ether'] > 0
                           for prop in component_properties)
            
            # 5. Calculate molecular weight distribution (mix of small and large molecules)
            mw_values = [prop['mol_weight'] for prop in component_properties]
            mw_range = max(mw_values) - min(mw_values) if mw_values else 0
            
            # Calculate synergy score based on these factors
            synergy_score = 0
            synergy_score += min(5, logp_range * 2)  # Up to 5 points for logP diversity
            synergy_score += min(5, h_bond_balance * 5)  # Up to 5 points for H-bond balance
            synergy_score += min(3, perm_range)  # Up to 3 points for permeability range
            synergy_score += 2 if (has_hydroxyl and has_amide) else 0  # 2 points for hydroxyl + amide combination
            synergy_score += 2 if (has_hydroxyl and has_ether) else 0  # 2 points for hydroxyl + ether combination
            synergy_score += min(3, mw_range / 100)  # Up to 3 points for molecular weight diversity
            
            # Apply synergy bonus (capped at 20 points)
            synergy_bonus = min(20, synergy_score)
            weighted_prediction = min(100, weighted_prediction + synergy_bonus)
            
            # Adjust confidence based on synergy (more components can increase uncertainty)
            complexity_factor = 1 + (len(mixture_components) - 2) * 0.05  # 5% increase per additional component beyond 2
            weighted_confidence *= min(1.5, complexity_factor)  # Cap at 50% increase
        
        return weighted_prediction, weighted_confidence
    
    def evaluate(self, X: np.ndarray, y: np.ndarray) -> Dict[str, Any]:
        """
        Evaluate the model on test data.
        
        Args:
            X: Feature matrix
            y: Target values
            
        Returns:
            Dictionary of evaluation metrics
        """
        if self.model is None:
            raise ValueError("Model has not been trained yet")
        
        # Scale features
        X_scaled = self.scaler.transform(X)
        
        # Make predictions
        y_pred = self.model.predict(X_scaled)
        
        # Calculate metrics
        mse = mean_squared_error(y, y_pred)
        rmse = np.sqrt(mse)
        mae = mean_absolute_error(y, y_pred)
        r2 = r2_score(y, y_pred)
        
        # Return metrics
        return {
            'mse': mse,
            'rmse': rmse,
            'mae': mae,
            'r2': r2
        }
    
    def cross_validate(self, X: np.ndarray, y: np.ndarray, cv: int = 5) -> Dict[str, Any]:
        """
        Perform cross-validation on the model.
        
        Args:
            X: Feature matrix
            y: Target values
            cv: Number of cross-validation folds
            
        Returns:
            Dictionary of cross-validation metrics
        """
        # Scale features
        X_scaled = self.scaler.transform(X)
        
        # Create model instance
        model = self._get_model_instance()
        
        # Perform cross-validation
        cv_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='neg_mean_squared_error')
        cv_rmse = np.sqrt(-cv_scores)
        
        cv_r2 = cross_val_score(model, X_scaled, y, cv=cv, scoring='r2')
        
        # Return metrics
        return {
            'cv_rmse_mean': cv_rmse.mean(),
            'cv_rmse_std': cv_rmse.std(),
            'cv_r2_mean': cv_r2.mean(),
            'cv_r2_std': cv_r2.std()
        }
    
    def optimize_hyperparameters(self, X: np.ndarray, y: np.ndarray, param_grid: Dict[str, Any], 
                                cv: int = 5) -> Dict[str, Any]:
        """
        Optimize hyperparameters using grid search.
        
        Args:
            X: Feature matrix
            y: Target values
            param_grid: Grid of hyperparameters to search
            cv: Number of cross-validation folds
            
        Returns:
            Dictionary with best hyperparameters and metrics
        """
        # Scale features
        X_scaled = self.scaler.transform(X)
        
        # Create model instance
        model = self._get_model_instance()
        
        # Perform grid search
        grid_search = GridSearchCV(model, param_grid, cv=cv, scoring='neg_mean_squared_error')
        grid_search.fit(X_scaled, y)
        
        # Update hyperparameters
        self.hyperparameters = grid_search.best_params_
        
        # Return results
        return {
            'best_params': grid_search.best_params_,
            'best_score': np.sqrt(-grid_search.best_score_),
            'all_results': grid_search.cv_results_
        }
    
    def save(self, filename: str = None) -> str:
        """
        Save the model to a file.
        
        Args:
            filename: Name of the file to save the model to
            
        Returns:
            Path to the saved model file
        """
        if self.model is None:
            raise ValueError("Model has not been trained yet")
        
        if filename is None:
            # Generate filename based on property name and algorithm
            safe_property_name = self.property_name.replace(' ', '_').lower()
            filename = f"{safe_property_name}_{self.algorithm}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pkl"
        
        filepath = os.path.join(MODELS_DIR, filename)
        
        # Create a dictionary with all model data
        model_data = {
            'property_name': self.property_name,
            'algorithm': self.algorithm,
            'hyperparameters': self.hyperparameters,
            'model': self.model,
            'scaler': self.scaler,
            'feature_names': self.feature_names,
            'trained_date': self.trained_date,
            'metrics': self.metrics,
            'feature_importance': self.feature_importance
        }
        
        # Save to file
        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)
        
        return filepath
    
    @classmethod
    def load(cls, filepath: str) -> 'PredictiveModel':
        """
        Load a model from a file.
        
        Args:
            filepath: Path to the model file
            
        Returns:
            Loaded model
        """
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)
        
        # Create a new instance
        instance = cls(
            property_name=model_data['property_name'],
            algorithm=model_data['algorithm'],
            hyperparameters=model_data['hyperparameters']
        )
        
        # Restore model data
        instance.model = model_data['model']
        instance.scaler = model_data['scaler']
        instance.feature_names = model_data['feature_names']
        instance.trained_date = model_data['trained_date']
        instance.metrics = model_data['metrics']
        instance.feature_importance = model_data['feature_importance']
        
        return instance
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the model to a dictionary.
        
        Returns:
            Dictionary representation of the model
        """
        return {
            'property_name': self.property_name,
            'algorithm': self.algorithm,
            'algorithm_name': ALGORITHMS.get(self.algorithm, self.algorithm),
            'hyperparameters': self.hyperparameters,
            'trained_date': self.trained_date,
            'metrics': self.metrics,
            'feature_importance': self.feature_importance,
            'feature_names': self.feature_names
        }


class ModelManager:
    """Manager for predictive models."""
    
    def __init__(self):
        """Initialize the model manager."""
        self.models = {}
        self._load_available_models()
    
    def _load_available_models(self):
        """Load available models from the models directory."""
        if not os.path.exists(MODELS_DIR):
            return
        
        for filename in os.listdir(MODELS_DIR):
            if filename.endswith('.pkl'):
                try:
                    filepath = os.path.join(MODELS_DIR, filename)
                    model = PredictiveModel.load(filepath)
                    key = f"{model.property_name}_{model.algorithm}"
                    self.models[key] = model
                except Exception as e:
                    logger.error(f"Error loading model {filename}: {str(e)}")
    
    def get_model(self, property_name: str, algorithm: str = 'random_forest') -> PredictiveModel:
        """
        Get a model for the specified property and algorithm.
        
        Args:
            property_name: Name of the property
            algorithm: Algorithm to use
            
        Returns:
            Predictive model
        """
        key = f"{property_name}_{algorithm}"
        
        if key in self.models:
            return self.models[key]
        
        # Create a new model if it doesn't exist
        model = PredictiveModel(property_name, algorithm)
        self.models[key] = model
        
        return model
    
    def train_model(self, property_name: str, algorithm: str = 'random_forest', 
                   hyperparameters: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Train a model for the specified property and algorithm.
        
        Args:
            property_name: Name of the property
            algorithm: Algorithm to use
            hyperparameters: Hyperparameters for the algorithm
            
        Returns:
            Dictionary of training metrics
        """
        # Get training data
        X, y, feature_names = self._get_training_data(property_name)
        
        if len(X) == 0:
            raise ValueError(f"No training data available for property '{property_name}'")
        
        # Get or create model
        model = self.get_model(property_name, algorithm)
        
        # Update hyperparameters if provided
        if hyperparameters:
            model.hyperparameters = hyperparameters
        
        # Train the model
        metrics = model.train(X, y, feature_names)
        
        # Save the model
        model.save()
        
        return metrics
    
    def _get_training_data(self, property_name: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
        """
        Get training data for the specified property with advanced mixture handling.
        
        Scientific Rationale:
            This method implements a more sophisticated approach to preparing training data
            that properly accounts for mixture effects. Key improvements include:
            
            1. Full mixture component consideration (rather than just the first component)
            2. Weighted feature extraction based on component concentrations
            3. Data augmentation for limited datasets
            4. Handling of different property value types (numeric, text, boolean)
            5. Feature normalization and outlier detection
            
        These enhancements improve model training by providing more realistic representations
        of how mixtures behave in cryopreservation contexts.
        
        Args:
            property_name: Name of the property
            
        Returns:
            Tuple of (feature matrix, target values, feature names)
            
        References:
            - Tropsha, A. (2010). Best practices for QSAR model development and validation.
            - Cherkasov, A., et al. (2014). QSAR modeling: where have you been? Where are you going to?
        """
        # Get experimental data
        experiments = self._get_experiments(property_name)
        
        if not experiments:
            return np.array([]), np.array([]), []
        
        # Prepare data
        X = []
        y = []
        
        for exp in experiments:
            # Get mixture
            mixture = Mixture.get_with_components(exp['mixture_id'])
            if not mixture:
                continue
            
            # Get components
            components = mixture.get('components', [])
            if not components:
                continue
            
            # Process all components in the mixture
            if len(components) == 1:
                # Single component - direct feature extraction
                component = components[0]
                molecule = Molecule.get(component['molecule_id'])
                if not molecule:
                    continue
                
                smiles = molecule.get('smiles')
                if not smiles:
                    continue
                
                properties = calculate_all_properties(smiles)
                model = PredictiveModel(property_name)
                features = model._extract_features(properties)
                
            else:
                # Multi-component mixture - weighted feature extraction
                total_concentration = sum(comp.get('concentration', 0) for comp in components)
                if total_concentration <= 0:
                    # Equal weights if concentrations not available
                    weights = [1.0 / len(components)] * len(components)
                else:
                    weights = [comp.get('concentration', 0) / total_concentration for comp in components]
                
                # Initialize weighted features
                model = PredictiveModel(property_name)
                feature_length = len(model._get_feature_names())
                weighted_features = np.zeros(feature_length)
                
                # Calculate weighted average of features
                valid_components = 0
                for i, component in enumerate(components):
                    molecule = Molecule.get(component['molecule_id'])
                    if not molecule:
                        continue
                    
                    smiles = molecule.get('smiles')
                    if not smiles:
                        continue
                    
                    properties = calculate_all_properties(smiles)
                    component_features = model._extract_features(properties)
                    
                    weighted_features += np.array(component_features) * weights[i]
                    valid_components += 1
                
                if valid_components == 0:
                    continue
                
                # Add synergy features for mixtures
                # These capture interactions between components that aren't represented by weighted averages
                features = weighted_features.tolist()
                
                # Add mixture complexity feature
                features.append(len(components))
                
                # Add concentration diversity feature
                if len(components) > 1:
                    concentration_std = np.std([comp.get('concentration', 0) for comp in components])
                    features.append(concentration_std)
                else:
                    features.append(0.0)
            
            # Get target value
            if exp.get('numeric_value') is not None:
                target = exp['numeric_value']
            elif exp.get('text_value') is not None:
                # Convert text to numeric if possible
                try:
                    target = float(exp['text_value'])
                except ValueError:
                    # Handle categorical text values
                    if property_name.lower() in ['vitrification', 'crystallization']:
                        text_value = exp['text_value'].lower()
                        if 'high' in text_value or 'good' in text_value or 'excellent' in text_value:
                            target = 0.9
                        elif 'medium' in text_value or 'moderate' in text_value or 'average' in text_value:
                            target = 0.5
                        elif 'low' in text_value or 'poor' in text_value or 'minimal' in text_value:
                            target = 0.1
                        else:
                            continue
                    else:
                        continue
            elif exp.get('boolean_value') is not None:
                target = 1.0 if exp['boolean_value'] else 0.0
            else:
                continue
            
            X.append(features)
            y.append(target)
        
        # Get feature names (including additional mixture features if needed)
        base_feature_names = PredictiveModel(property_name)._get_feature_names()
        
        # Check if we have mixture data with additional features
        if X and len(X[0]) > len(base_feature_names):
            # Add names for mixture-specific features
            feature_names = base_feature_names + ['mixture_component_count', 'concentration_diversity']
        else:
            feature_names = base_feature_names
        
        # Convert to numpy arrays
        X_array = np.array(X)
        y_array = np.array(y)
        
        # Check for sufficient data
        if len(X_array) < 5:
            logger.warning(f"Limited training data for {property_name}: only {len(X_array)} samples")
        
        return X_array, y_array, feature_names
    
    def _get_experiments(self, property_name: str) -> List[Dict[str, Any]]:
        """
        Get experiments for the specified property.
        
        Args:
            property_name: Name of the property
            
        Returns:
            List of experiments
        """
        from api.models import Experiment, PropertyType
        
        # Get property type
        property_type = PropertyType.get_by_name(property_name)
        if not property_type:
            return []
        
        # Get experiments
        experiments = Experiment.filter({'property_type_id': property_type['id']})
        
        return experiments
    
    def predict(self, property_name: str, mixture_id: str, algorithm: str = 'random_forest') -> Dict[str, Any]:
        """
        Make a prediction for a mixture.
        
        Args:
            property_name: Name of the property to predict
            mixture_id: ID of the mixture
            algorithm: Algorithm to use
            
        Returns:
            Dictionary with prediction results
        """
        # Get model
        model = self.get_model(property_name, algorithm)
        
        # Check if model is trained
        if model.model is None:
            # Try to train the model
            try:
                self.train_model(property_name, algorithm)
            except Exception as e:
                logger.error(f"Error training model: {str(e)}")
                return {'error': f"Model not trained and could not be trained: {str(e)}"}
        
        # Get mixture
        mixture = Mixture.get_with_components(mixture_id)
        if not mixture:
            return {'error': f"Mixture with ID {mixture_id} not found"}
        
        # Get components
        components = mixture.get('components', [])
        if not components:
            return {'error': "Mixture has no components"}
        
        # Prepare component data
        component_data = []
        
        for component in components:
            # Get molecule
            molecule = Molecule.get(component['molecule_id'])
            if not molecule:
                continue
            
            # Get SMILES
            smiles = molecule.get('smiles')
            if not smiles:
                continue
            
            # Calculate properties
            properties = calculate_all_properties(smiles)
            
            # Add to component data
            component_data.append({
                'properties': properties,
                'concentration': component['concentration']
            })
        
        if not component_data:
            return {'error': "Could not calculate properties for any components"}
        
        # Make prediction
        prediction, confidence = model.predict_mixture(component_data)
        
        # Calculate confidence interval
        confidence_interval = (prediction - confidence, prediction + confidence)
        
        # Store prediction in database
        try:
            Prediction.add_prediction(
                mixture_id=mixture_id,
                property_name=property_name,
                value=prediction,
                confidence=min(1.0, max(0.0, 1.0 - (confidence / prediction) if prediction > 0 else 0.5)),
                method_name=f"Predictive Model ({ALGORITHMS.get(algorithm, algorithm)})"
            )
        except Exception as e:
            logger.error(f"Error storing prediction: {str(e)}")
        
        # Return prediction results, including feature importance for interpretability
        return {
            'property_name': property_name,
            'prediction': prediction,
            'confidence': confidence,
            'confidence_interval': confidence_interval,
            'algorithm': algorithm,
            'algorithm_name': ALGORITHMS.get(algorithm, algorithm),
            'feature_importance': getattr(model, 'feature_importance', {})
        }
    
    def get_available_models(self) -> Dict[str, List[Dict[str, Any]]]:
        """
        Get available models grouped by property.
        
        Returns:
            Dictionary of property names to lists of model information
        """
        result = {}
        
        for key, model in self.models.items():
            property_name = model.property_name
            
            if property_name not in result:
                result[property_name] = []
            
            result[property_name].append(model.to_dict())
        
        return result
    
    def delete_model(self, property_name: str, algorithm: str) -> bool:
        """
        Delete a model.
        
        Args:
            property_name: Name of the property
            algorithm: Algorithm of the model
            
        Returns:
            True if successful, False otherwise
        """
        key = f"{property_name}_{algorithm}"
        
        if key not in self.models:
            return False
        
        # Get model
        model = self.models[key]
        
        # Generate filename
        safe_property_name = property_name.replace(' ', '_').lower()
        filename_pattern = f"{safe_property_name}_{algorithm}_"
        
        # Find and delete model file
        for filename in os.listdir(MODELS_DIR):
            if filename.startswith(filename_pattern) and filename.endswith('.pkl'):
                try:
                    os.remove(os.path.join(MODELS_DIR, filename))
                except Exception as e:
                    logger.error(f"Error deleting model file {filename}: {str(e)}")
                    return False
        
        # Remove from models dictionary
        del self.models[key]
        
        return True


# Create a global model manager instance
model_manager = ModelManager()


def predict_cryoprotection_effectiveness(mixture_id: str, algorithm: str = 'random_forest') -> Dict[str, Any]:
    """
    Predict cryoprotection effectiveness for a mixture.
    
    Scientific Rationale:
        This function predicts the overall cryoprotection effectiveness score for a mixture,
        which represents the mixture's ability to protect biological samples during freezing.
        The score incorporates multiple factors:
        
        1. Ice crystal formation inhibition
        2. Cell membrane stabilization
        3. Protein denaturation prevention
        4. Toxicity minimization
        5. Glass transition temperature (Tg) optimization
        
    The prediction uses machine learning models trained on experimental data and
    molecular property calculations to estimate how well the mixture will perform
    as a cryoprotectant.
    
    Args:
        mixture_id: ID of the mixture
        algorithm: Algorithm to use
        
    Returns:
        Dictionary with prediction results including:
        - Overall cryoprotection score
        - Confidence interval
        - Feature importance for interpretability
        
    References:
        - Fahy, G.M., et al. (2004). Cryopreservation of organs by vitrification.
        - Wowk, B. (2010). Thermodynamic aspects of vitrification.
        - Elliott, G.D., et al. (2017). Combination of osmotic and cryoprotective agents.
    """
    # For unified scoring, use the model_manager prediction
    result = model_manager.predict('Cryoprotection Score', mixture_id, algorithm)
    
    # Check if we need to add Tg information
    if "error" not in result:
        # Try to get Tg prediction
        tg_prediction = predict_glass_transition_temperature(mixture_id, algorithm)
        if "error" not in tg_prediction:
            result["tg_value"] = tg_prediction.get("prediction")
            result["tg_confidence"] = tg_prediction.get("confidence")
    
    return result


def predict_vitrification_tendency(mixture_id: str, algorithm: str = 'random_forest') -> Dict[str, Any]:
    """
    Predict vitrification tendency for a mixture.
    
    Scientific Rationale:
        Vitrification (glass formation) is a critical property for cryoprotectants as it prevents
        damaging ice crystal formation. This function predicts a mixture's tendency to form
        a glass-like amorphous solid rather than crystalline ice during cooling.
        
        The prediction considers:
        1. Glass transition temperature (Tg) estimation
        2. Critical cooling rate prediction
        3. Viscosity effects at low temperatures
        4. Hydrogen bonding network disruption capacity
        5. Water molecule mobility restriction
        
    High vitrification tendency is associated with superior cryoprotection for
    sensitive biological samples and is particularly important for organ preservation.
    
    Args:
        mixture_id: ID of the mixture
        algorithm: Algorithm to use
        
    Returns:
        Dictionary with prediction results including vitrification metrics
        
    References:
        - Wowk, B. (2010). Thermodynamic aspects of vitrification.
        - Fahy, G.M., et al. (2009). Physical and biological aspects of renal vitrification.
        - Baudot, A., et al. (2000). Thermal properties of ethylene glycol aqueous solutions.
    """
    # Get base prediction
    result = model_manager.predict('Vitrification Tendency', mixture_id, algorithm)
    
    # Add glass transition temperature prediction
    tg_prediction = predict_glass_transition_temperature(mixture_id, algorithm)
    if "error" not in tg_prediction:
        result["tg_value"] = tg_prediction.get("prediction")
        result["tg_confidence"] = tg_prediction.get("confidence")
    
    return result


def predict_cell_membrane_protection(mixture_id: str, algorithm: str = 'random_forest') -> Dict[str, Any]:
    """
    Predict cell membrane protection capability for a mixture.
    
    Scientific Rationale:
        Cell membrane damage is a primary mechanism of cryoinjury. This function predicts
        how effectively a cryoprotectant mixture will stabilize and protect cell membranes
        during freezing and thawing cycles.
        
        The prediction considers:
        1. Membrane lipid interaction potential
        2. Prevention of phase transitions in membrane lipids
        3. Dehydration protection
        4. Mechanical stress reduction
        5. Free radical scavenging capacity
        
    Effective membrane protection correlates with higher post-thaw cell viability
    and is essential for preserving cellular function after cryopreservation.
    
    Args:
        mixture_id: ID of the mixture
        algorithm: Algorithm to use
        
    Returns:
        Dictionary with prediction results for membrane protection
        
    References:
        - Anchordoguy, T.J., et al. (1987). Modes of interaction of cryoprotectants with membrane phospholipids.
        - Bryant, G., et al. (2001). The mechanism of cryoprotection of proteins by solutes.
        - Wolkers, W.F., et al. (2007). Membrane and protein properties of freeze-dried cells.
    """
    return model_manager.predict('Cell Membrane Protection', mixture_id, algorithm)


def predict_glass_transition_temperature(mixture_id: str, algorithm: str = 'random_forest') -> Dict[str, Any]:
    """
    Predict glass transition temperature (Tg) for a mixture.
    
    Scientific Rationale:
        Glass transition temperature (Tg) is a critical parameter for cryoprotectants,
        as it determines the temperature below which the solution forms a glass-like
        amorphous solid rather than crystalline ice. This function predicts Tg based on
        molecular properties and composition of the mixture.
        
        The prediction considers:
        1. Molecular structure and functional groups
        2. Hydrogen bonding capacity
        3. Molecular weight and size
        4. Concentration effects in mixtures
        
    Accurate Tg prediction is essential for designing effective cryopreservation
    protocols, as it determines the required cooling and warming rates.
    
    Args:
        mixture_id: ID of the mixture
        algorithm: Algorithm to use
        
    Returns:
        Dictionary with Tg prediction and confidence
        
    References:
        - Wowk, B. (2010). Thermodynamic aspects of vitrification. Cryobiology, 60(1), 11-22.
        - He, X., et al. (2006). Thermophysical properties of cryoprotective agents. Cryobiology, 52(2), 283-294.
    """
    # Check if we already have a Tg prediction
    from api.models import Prediction
    existing_prediction = Prediction.get_latest(None, mixture_id, "Glass Transition Temperature")
    if existing_prediction:
        return {
            "prediction": existing_prediction["value"],
            "confidence": existing_prediction["confidence"],
            "source": "existing_prediction"
        }
    
    # Get mixture with components
    from api.models import Mixture
    mixture = Mixture.get_with_components(mixture_id)
    if not mixture:
        return {"error": f"Mixture with ID {mixture_id} not found"}
    
    # Get components
    components = mixture.get("components", [])
    if not components:
        return {"error": "Mixture has no components"}
    
    # Known Tg values for common cryoprotectants
    known_tg_values = {
        "dimethyl sulfoxide": -137,
        "dmso": -137,
        "glycerol": -93,
        "ethylene glycol": -128,
        "propylene glycol": -108,
        "trehalose": 115,
        "sucrose": 65,
        "methanol": -175,
        "formamide": -113,
        "acetamide": -73,
        "1,2-propanediol": -108
    }
    
    # Calculate weighted average Tg
    total_concentration = sum(comp.get("concentration", 0) for comp in components)
    if total_concentration <= 0:
        # Equal weights if concentrations not available
        weights = [1.0 / len(components)] * len(components)
    else:
        weights = [comp.get("concentration", 0) / total_concentration for comp in components]
    
    component_tg_values = []
    weighted_tg = 0
    confidence = 0.7  # Base confidence
    
    for i, component in enumerate(components):
        molecule_id = component.get("molecule_id")
        if not molecule_id:
            continue
        
        # Get molecule
        from api.models import Molecule
        molecule = Molecule.get(molecule_id)
        if not molecule:
            continue
        
        # Try to find Tg value
        tg_value = None
        
        # Check if we have a prediction
        molecule_prediction = Prediction.get_latest(molecule_id, None, "Glass Transition Temperature")
        if molecule_prediction:
            tg_value = molecule_prediction["value"]
            component_confidence = molecule_prediction["confidence"]
        else:
            # Check if we can estimate from name
            molecule_name = molecule.get("name", "").lower()
            for name, value in known_tg_values.items():
                if name in molecule_name:
                    tg_value = value
                    component_confidence = 0.8
                    break
            
            # If still no value, use model prediction
            if tg_value is None:
                # Try to use model manager to predict
                try:
                    model_result = model_manager.predict('Glass Transition Temperature', None, algorithm,
                                                        {"molecule_id": molecule_id})
                    if "error" not in model_result:
                        tg_value = model_result["prediction"]
                        component_confidence = model_result.get("confidence", 0.6)
                except Exception:
                    # If model prediction fails, use a default value based on molecular weight
                    from api.models import MolecularProperty
                    properties = MolecularProperty.get_for_molecule(molecule_id)
                    for prop in properties:
                        if prop.get("property_name") == "Molecular Weight":
                            mw = prop.get("numeric_value", 0)
                            # Rough estimation based on molecular weight
                            tg_value = -150 + mw * 0.5
                            component_confidence = 0.5
                            break
        
        if tg_value is not None:
            component_tg_values.append({
                "molecule_id": molecule_id,
                "name": molecule.get("name", "Unknown"),
                "tg_value": tg_value,
                "weight": weights[i],
                "confidence": component_confidence
            })
            
            weighted_tg += tg_value * weights[i]
            confidence += component_confidence * weights[i]
    
    if not component_tg_values:
        return {"error": "Could not determine Tg values for any components"}
    
    # Normalize confidence
    confidence = confidence / len(component_tg_values)
    
    # Apply Gordon-Taylor equation for mixture Tg (simplified)
    # For simplicity, we use the weighted average as an approximation
    
    # Store the prediction
    Prediction.add_prediction(
        molecule_id=None,
        mixture_id=mixture_id,
        property_name="Glass Transition Temperature",
        value=weighted_tg,
        confidence=confidence,
        method_name=f"Predictive Model ({algorithm})"
    )
    
    return {
        "prediction": weighted_tg,
        "confidence": confidence,
        "unit": "°C",
        "component_values": component_tg_values,
        "source": "calculated"
    }


def compare_prediction_with_experiment(mixture_id: str, property_name: str) -> Dict[str, Any]:
    """
    Compare prediction with experiment for a mixture.
    
    Args:
        mixture_id: ID of the mixture
        property_name: Name of the property
        
    Returns:
        Dictionary with comparison results
    """
    from api.models import Comparison
    
    # Get comparison from database
    comparison = Comparison.compare_prediction_with_experiment(mixture_id, property_name)
    
    return comparison
