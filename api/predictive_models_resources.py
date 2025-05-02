"""
CryoProtect Analyzer API - Predictive Models Resources

This module provides API resources for the predictive models functionality, including:
- Training and managing machine learning models
- Making predictions for mixtures and custom properties
- Evaluating model performance
- Batch training of multiple models

All endpoints follow RESTful design principles and return standardized
response formats. Authentication is required for operations that modify
the model repository.
"""

import logging
from flask import request, current_app
from flask_restful import Resource, marshal_with, reqparse, abort, fields
from marshmallow import ValidationError, Schema, fields as ma_fields, validate
from typing import Dict, List, Any, Tuple, Optional, Union

# Import API documentation utilities
try:
    from flask_apispec import use_kwargs, marshal_with as apispec_marshal_with, doc
    from flask_apispec.views import MethodResource
except ImportError:
    # Create dummy decorators if flask-apispec is not installed
    def use_kwargs(*args, **kwargs): return lambda f: f
    def apispec_marshal_with(*args, **kwargs): return lambda f: f
    def doc(*args, **kwargs): return lambda f: f
    MethodResource = Resource

from api.models import prediction_fields
from api.utils import token_required, handle_supabase_error, get_user_id
from api.predictive_models import model_manager, ALGORITHMS, predict_cryoprotection_effectiveness

# Set up logging
logger = logging.getLogger(__name__)

# Custom fields for model information
model_info_fields = {
    'property_name': fields.String,
    'algorithm': fields.String,
    'algorithm_name': fields.String,
    'hyperparameters': fields.Raw,
    'trained_date': fields.String,
    'metrics': fields.Raw,
    'feature_importance': fields.Raw
}

# Custom fields for prediction results
prediction_result_fields = {
    'property_name': fields.String,
    'prediction': fields.Float,
    'confidence': fields.Float,
    'confidence_interval': fields.Raw,
    'algorithm': fields.String,
    'algorithm_name': fields.String,
    'feature_importance': fields.Raw
}

# Validation schemas
class ModelTrainingSchema(Schema):
    """Schema for model training requests."""
    property_name = ma_fields.String(required=True)
    algorithm = ma_fields.String(validate=validate.OneOf(list(ALGORITHMS.keys())), default='random_forest')
    hyperparameters = ma_fields.Dict()

class BatchModelTrainingSchema(Schema):
    """Schema for batch model training requests."""
    models = ma_fields.List(ma_fields.Dict(keys=ma_fields.String(), values=ma_fields.Field()), required=True)

class ModelEvaluationSchema(Schema):
    """Schema for model evaluation requests."""
    property_name = ma_fields.String(required=True)
    algorithm = ma_fields.String(validate=validate.OneOf(list(ALGORITHMS.keys())), default='random_forest')
    test_size = ma_fields.Float(validate=validate.Range(min=0.1, max=0.5), default=0.2)
    cross_validation = ma_fields.Boolean(default=True)
    cv_folds = ma_fields.Integer(validate=validate.Range(min=2, max=10), default=5)

class PredictionRequestSchema(Schema):
    """Schema for prediction requests."""
    algorithm = ma_fields.String(validate=validate.OneOf(list(ALGORITHMS.keys())), default='random_forest')


class PredictiveModelListResource(MethodResource):
    """Resource for listing and training predictive models.
    
    Scientific Rationale:
        This endpoint provides access to the available predictive models and allows
        for training new models. It supports various machine learning algorithms
        and customizable hyperparameters for optimizing model performance.
        
        The model training process incorporates:
        1. Feature selection based on physicochemical properties relevant to cryoprotection
        2. Molecular descriptor calculation using RDKit and custom algorithms
        3. Data preprocessing including scaling and outlier detection
        4. Hyperparameter optimization for model performance
        5. Cross-validation for robust evaluation
        
        These models enable quantitative structure-activity relationship (QSAR) analysis
        for cryoprotectants, facilitating rational design and optimization of
        cryopreservation protocols.
        
    References:
        - Cherkasov, A., et al. (2014). QSAR modeling: where have you been? Where are you going to?
        - Todeschini, R., & Consonni, V. (2009). Molecular Descriptors for Chemoinformatics.
    """
    
    @doc(description='Get a list of available predictive models',
         tags=['Predictive Models'])
    def get(self) -> Tuple[Dict[str, Any], int]:
        """Get a list of available predictive models.
        
        Retrieves information about all available predictive models, including
        their properties, algorithms, training status, and performance metrics.
        
        Returns:
            Tuple[Dict[str, Any], int]: List of available models and HTTP status code
            
        Raises:
            500: If server error occurs
            
        Response:
            List of model information including:
            - Property name
            - Algorithm
            - Training status
            - Performance metrics (if trained)
            - Feature importance (if available)
        """
        try:
            models = model_manager.get_available_models()
            return models, 200
        except Exception as e:
            current_app.logger.error(f"Error fetching models: {str(e)}")
            abort(500, message=f"Error fetching models: {str(e)}")
    
    @doc(description='Train a new predictive model',
         tags=['Predictive Models'],
         security=[{'Bearer': []}],
         request_body={
             'property_name': {'description': 'Property to predict', 'type': 'string', 'required': True},
             'algorithm': {'description': 'Machine learning algorithm to use', 'type': 'string', 'default': 'random_forest', 'enum': list(ALGORITHMS.keys())},
             'hyperparameters': {'description': 'Algorithm-specific hyperparameters', 'type': 'object'}
         })
    @token_required
    def post(self) -> Tuple[Dict[str, Any], int]:
        """Train a new predictive model.
        
        Trains a machine learning model to predict the specified property using
        the selected algorithm and hyperparameters. The model is trained on
        available data in the system and stored for future predictions.
        
        Returns:
            Tuple[Dict[str, Any], int]: Trained model information and HTTP status code
            
        Raises:
            400: If validation error occurs or training fails
            500: If server error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            property_name (str): Property to predict (required)
            algorithm (str): Machine learning algorithm to use (default: 'random_forest')
                Supported algorithms: 'random_forest', 'gradient_boosting', 'svm', 'neural_network'
            hyperparameters (dict, optional): Algorithm-specific hyperparameters
                
        Response:
            Trained model information including:
            - Property name
            - Algorithm
            - Training metrics
            - Feature importance
        """
        try:
            # Validate request data
            schema = ModelTrainingSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Extract parameters
            property_name = data['property_name']
            algorithm = data.get('algorithm', 'random_forest')
            hyperparameters = data.get('hyperparameters')
            
            # Train model
            try:
                metrics = model_manager.train_model(property_name, algorithm, hyperparameters)
                
                # Return model info
                model = model_manager.get_model(property_name, algorithm)
                return model.to_dict(), 201
            except ValueError as e:
                abort(400, message=str(e))
            except Exception as e:
                current_app.logger.error(f"Error training model: {str(e)}")
                abort(500, message=f"Error training model: {str(e)}")
                
        except Exception as e:
            current_app.logger.error(f"Error processing request: {str(e)}")
            abort(500, message=f"Error processing request: {str(e)}")


class PredictiveModelResource(MethodResource):
    """Resource for managing a specific predictive model.
    
    Scientific Rationale:
        This endpoint provides access to information about a specific predictive model
        and allows for its deletion. It enables detailed inspection of model characteristics
        and performance metrics, which is essential for:
        
        1. Model validation and quality assessment
        2. Scientific reproducibility and transparency
        3. Informed decision-making in cryoprotectant design
        4. Understanding structure-property relationships
        
        The detailed model information includes training metrics, hyperparameters,
        and feature importance, providing insights into the molecular properties
        that drive cryoprotectant effectiveness.
    """
    
    @doc(description='Get information about a specific model',
         tags=['Predictive Models'],
         params={
             'property_name': {'description': 'Property the model predicts', 'type': 'string', 'required': True},
             'algorithm': {'description': 'Algorithm used by the model', 'type': 'string', 'required': True}
         })
    def get(self, property_name: str, algorithm: str) -> Tuple[Dict[str, Any], int]:
        """Get information about a specific model.
        
        Retrieves detailed information about a specific predictive model,
        including its properties, algorithm, training status, and performance metrics.
        
        Args:
            property_name (str): Property the model predicts
            algorithm (str): Algorithm used by the model
            
        Returns:
            Tuple[Dict[str, Any], int]: Model information and HTTP status code
            
        Raises:
            404: If model not found or not trained
            500: If server error occurs
            
        Response:
            Detailed model information including:
            - Property name
            - Algorithm
            - Training date
            - Hyperparameters
            - Performance metrics
            - Feature importance
        """
        try:
            model = model_manager.get_model(property_name, algorithm)
            
            if model.model is None:
                abort(404, message=f"Model for property '{property_name}' with algorithm '{algorithm}' not found or not trained")
            
            return model.to_dict(), 200
        except Exception as e:
            current_app.logger.error(f"Error fetching model: {str(e)}")
            abort(500, message=f"Error fetching model: {str(e)}")
    
    @doc(description='Delete a specific model',
         tags=['Predictive Models'],
         params={
             'property_name': {'description': 'Property the model predicts', 'type': 'string', 'required': True},
             'algorithm': {'description': 'Algorithm used by the model', 'type': 'string', 'required': True}
         },
         security=[{'Bearer': []}])
    @token_required
    def delete(self, property_name: str, algorithm: str) -> Tuple[Dict[str, Any], int]:
        """Delete a specific model.
        
        Permanently removes a predictive model from the system. This operation
        cannot be undone, and the model will need to be retrained if needed again.
        
        Args:
            property_name (str): Property the model predicts
            algorithm (str): Algorithm used by the model
            
        Returns:
            Tuple[Dict[str, Any], int]: Success message and HTTP status code
            
        Raises:
            404: If model not found
            500: If server error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Response:
            Success message confirming the model was deleted
        """
        try:
            success = model_manager.delete_model(property_name, algorithm)
            
            if not success:
                abort(404, message=f"Model for property '{property_name}' with algorithm '{algorithm}' not found")
            
            return {'message': f"Model for property '{property_name}' with algorithm '{algorithm}' deleted successfully"}, 200
        except Exception as e:
            current_app.logger.error(f"Error deleting model: {str(e)}")
            abort(500, message=f"Error deleting model: {str(e)}")


class MixturePredictionResource(MethodResource):
    """Resource for making predictions for a mixture.
    
    Scientific Rationale:
        This endpoint predicts the cryoprotection effectiveness of a mixture using
        trained machine learning models. The prediction integrates multiple factors:
        
        1. Ice crystal formation inhibition
        2. Cell membrane stabilization
        3. Protein denaturation prevention
        4. Toxicity minimization
        5. Glass transition temperature (Tg) optimization
        
        The models are trained on experimental data and molecular property calculations
        to estimate how well the mixture will perform as a cryoprotectant. The prediction
        includes confidence intervals and feature importance to provide scientific
        interpretability.
        
    References:
        - Fahy, G.M., et al. (2004). Cryopreservation of organs by vitrification.
        - Wowk, B. (2010). Thermodynamic aspects of vitrification.
        - Elliott, G.D., et al. (2017). Combination of osmotic and cryoprotective agents.
    """
    
    @doc(description='Get a prediction for a mixture',
         tags=['Predictive Models'],
         params={
             'mixture_id': {'description': 'ID of the mixture to predict for', 'type': 'string', 'required': True},
             'algorithm': {'description': 'Algorithm to use for prediction', 'type': 'string', 'default': 'random_forest', 'enum': list(ALGORITHMS.keys())}
         })
    @apispec_marshal_with(prediction_result_fields, code=200, description='Prediction results')
    @marshal_with(prediction_result_fields)
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """Get a prediction for a mixture.
        
        Predicts the cryoprotection effectiveness of a mixture using the specified
        machine learning algorithm. Returns the predicted value along with confidence
        information and feature importance.
        
        Args:
            mixture_id (str): ID of the mixture to predict for
            
        Returns:
            Tuple[Dict[str, Any], int]: Prediction results and HTTP status code
            
        Raises:
            400: If validation error occurs or prediction fails
            500: If server error occurs
            
        Query Parameters:
            algorithm (str, optional): Algorithm to use for prediction (default: 'random_forest')
                Supported algorithms: 'random_forest', 'gradient_boosting', 'svm', 'neural_network'
                
        Response:
            Prediction results including:
            - Predicted value
            - Confidence score
            - Confidence interval (if available)
            - Feature importance
            - Algorithm information
        """
        try:
            # Validate query parameters
            schema = PredictionRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Extract parameters
            algorithm = args.get('algorithm', 'random_forest')
            
            # Make prediction
            result = predict_cryoprotection_effectiveness(mixture_id, algorithm)
            
            if 'error' in result:
                abort(400, message=result['error'])
            
            return result, 200
        except Exception as e:
            current_app.logger.error(f"Error making prediction: {str(e)}")
            abort(500, message=f"Error making prediction: {str(e)}")


class CustomPropertyPredictionResource(MethodResource):
    """Resource for making predictions for custom properties.
    
    Scientific Rationale:
        This endpoint predicts custom properties for a mixture using trained
        machine learning models. It allows for flexible prediction of various
        properties beyond the standard cryoprotection effectiveness.
        
        Custom properties may include:
        - Glass transition temperature (Tg)
        - Toxicity metrics
        - Permeability coefficients
        - Osmotic stress parameters
        - Protein stabilization effects
        - Specific cell type compatibility
        
        The prediction models integrate molecular descriptors, experimental data,
        and structure-property relationships to estimate how a mixture will
        perform for the specified property.
    """
    
    @doc(description='Get a prediction for a custom property',
         tags=['Predictive Models'],
         params={
             'mixture_id': {'description': 'ID of the mixture to predict for', 'type': 'string', 'required': True},
             'property_name': {'description': 'Property to predict', 'type': 'string', 'required': True},
             'algorithm': {'description': 'Algorithm to use for prediction', 'type': 'string', 'default': 'random_forest', 'enum': list(ALGORITHMS.keys())}
         })
    @apispec_marshal_with(prediction_result_fields, code=200, description='Prediction results')
    @marshal_with(prediction_result_fields)
    def get(self, mixture_id: str, property_name: str) -> Tuple[Dict[str, Any], int]:
        """Get a prediction for a custom property.
        
        Predicts a custom property for a mixture using the specified machine
        learning algorithm. Returns the predicted value along with confidence
        information and feature importance.
        
        Args:
            mixture_id (str): ID of the mixture to predict for
            property_name (str): Property to predict
            
        Returns:
            Tuple[Dict[str, Any], int]: Prediction results and HTTP status code
            
        Raises:
            400: If validation error occurs or prediction fails
            500: If server error occurs
            
        Query Parameters:
            algorithm (str, optional): Algorithm to use for prediction (default: 'random_forest')
                Supported algorithms: 'random_forest', 'gradient_boosting', 'svm', 'neural_network'
                
        Response:
            Prediction results including:
            - Predicted value
            - Confidence score
            - Confidence interval (if available)
            - Feature importance
            - Algorithm information
        """
        try:
            # Validate query parameters
            schema = PredictionRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Extract parameters
            algorithm = args.get('algorithm', 'random_forest')
            
            # Make prediction
            result = model_manager.predict(property_name, mixture_id, algorithm)
            
            if 'error' in result:
                abort(400, message=result['error'])
            
            return result, 200
        except Exception as e:
            current_app.logger.error(f"Error making prediction: {str(e)}")
            abort(500, message=f"Error making prediction: {str(e)}")


class BatchModelTrainingResource(MethodResource):
    """Resource for training multiple predictive models in batch.
    
    This endpoint allows for training multiple predictive models in a single
    request. It is useful for efficiently training models for different properties
    or comparing multiple algorithms for the same property.
    """
    
    @doc(description='Train multiple predictive models in batch',
         tags=['Predictive Models'],
         security=[{'Bearer': []}],
         request_body={
             'models': {
                 'description': 'List of models to train',
                 'type': 'array',
                 'items': {
                     'type': 'object',
                     'properties': {
                         'property_name': {'description': 'Property to predict', 'type': 'string'},
                         'algorithm': {'description': 'Machine learning algorithm to use', 'type': 'string', 'default': 'random_forest', 'enum': list(ALGORITHMS.keys())},
                         'hyperparameters': {'description': 'Algorithm-specific hyperparameters', 'type': 'object'}
                     },
                     'required': ['property_name']
                 },
                 'required': True
             }
         })
    @token_required
    def post(self) -> Tuple[Dict[str, Any], int]:
        """
        Train multiple predictive models in batch.
        
        Trains multiple machine learning models in a single request, allowing for
        efficient training of models for different properties or comparing multiple
        algorithms for the same property.
        
        Returns:
            Tuple[Dict[str, Any], int]: Batch training results and HTTP status code
            
        Raises:
            400: If validation error occurs or training fails
            500: If server error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            models (List[Dict]): List of models to train (required)
                Each model specification includes:
                - property_name (str): Property to predict (required)
                - algorithm (str): Machine learning algorithm to use (default: 'random_forest')
                - hyperparameters (dict, optional): Algorithm-specific hyperparameters
                
        Response:
            Batch training results including:
            - Overall status
            - Individual results for each model
            - Success/error information for each model
        """
        try:
            # Validate request data
            schema = BatchModelTrainingSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Extract models to train
            models_to_train = data['models']
            
            if not models_to_train:
                abort(400, message="No models specified for training")
            
            # Train each model
            results = []
            for model_spec in models_to_train:
                try:
                    # Extract parameters
                    property_name = model_spec.get('property_name')
                    algorithm = model_spec.get('algorithm', 'random_forest')
                    hyperparameters = model_spec.get('hyperparameters')
                    
                    if not property_name:
                        results.append({
                            'property_name': None,
                            'algorithm': algorithm,
                            'status': 'ERROR',
                            'message': 'Property name is required'
                        })
                        continue
                    
                    # Train model
                    try:
                        metrics = model_manager.train_model(property_name, algorithm, hyperparameters)
                        
                        # Get model info
                        model = model_manager.get_model(property_name, algorithm)
                        
                        results.append({
                            'property_name': property_name,
                            'algorithm': algorithm,
                            'status': 'SUCCESS',
                            'metrics': metrics,
                            'model_info': model.to_dict() if model else None
                        })
                    except ValueError as e:
                        results.append({
                            'property_name': property_name,
                            'algorithm': algorithm,
                            'status': 'ERROR',
                            'message': str(e)
                        })
                    except Exception as e:
                        current_app.logger.error(f"Error training model {property_name}/{algorithm}: {str(e)}")
                        results.append({
                            'property_name': property_name,
                            'algorithm': algorithm,
                            'status': 'ERROR',
                            'message': f"Error training model: {str(e)}"
                        })
                except Exception as e:
                    current_app.logger.error(f"Error processing model spec: {str(e)}")
                    results.append({
                        'status': 'ERROR',
                        'message': f"Error processing model specification: {str(e)}"
                    })
            
            # Determine overall status
            overall_status = "SUCCESS"
            if any(r.get('status') == 'ERROR' for r in results):
                overall_status = "COMPLETED_WITH_ERRORS"
            
            return {
                'status': overall_status,
                'results': results
            }, 200
                
        except Exception as e:
            current_app.logger.error(f"Error processing batch training request: {str(e)}")
            abort(500, message=f"Error processing batch training request: {str(e)}")


class ModelEvaluationResource(MethodResource):
    """Resource for evaluating predictive models.
    
    Scientific Rationale:
        Model evaluation is a critical step in the development of reliable predictive models
        for cryoprotectant properties. This endpoint implements rigorous statistical evaluation
        methods including:
        
        1. Train/test splitting with stratification where appropriate
        2. K-fold cross-validation for robust performance estimation
        3. Multiple evaluation metrics (R², RMSE, MAE) for comprehensive assessment
        4. Learning curve analysis to detect overfitting/underfitting
        5. Feature importance analysis for scientific interpretability
        
    These evaluation techniques ensure that models are properly validated before being
    used for cryoprotectant design and optimization, following best practices in
    cheminformatics and QSAR modeling.
    """
    
    @doc(description='Evaluate a predictive model using test data or cross-validation',
         tags=['Predictive Models'],
         security=[{'Bearer': []}],
         request_body={
             'property_name': {'description': 'Property the model predicts', 'type': 'string', 'required': True},
             'algorithm': {'description': 'Algorithm used by the model', 'type': 'string', 'default': 'random_forest', 'enum': list(ALGORITHMS.keys())},
             'test_size': {'description': 'Fraction of data to use for testing', 'type': 'number', 'default': 0.2, 'minimum': 0.1, 'maximum': 0.5},
             'cross_validation': {'description': 'Whether to use cross-validation', 'type': 'boolean', 'default': True},
             'cv_folds': {'description': 'Number of folds for cross-validation', 'type': 'integer', 'default': 5, 'minimum': 2, 'maximum': 10}
         })
    @token_required
    def post(self) -> Tuple[Dict[str, Any], int]:
        """
        Evaluate a predictive model using test data or cross-validation.
        
        Scientific Rationale:
            This method implements a comprehensive evaluation protocol for cryoprotectant
            property prediction models. The evaluation process includes:
            
            1. Data splitting using stratification where appropriate
            2. K-fold cross-validation for robust performance estimation
            3. Multiple evaluation metrics for comprehensive assessment:
               - R² score (coefficient of determination) for explained variance
               - RMSE (root mean squared error) for prediction accuracy
               - MAE (mean absolute error) for prediction precision
            4. Learning curve analysis to detect overfitting/underfitting
            5. Feature importance analysis for scientific interpretability
            
        These evaluation techniques ensure that models are properly validated before
        being used for cryoprotectant design and optimization.
        
        Returns:
            Tuple[Dict[str, Any], int]: Evaluation metrics and HTTP status code
            
        Raises:
            400: If validation error occurs
            404: If model not found or not trained
            500: If server error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            property_name (str): Property the model predicts (required)
            algorithm (str): Algorithm used by the model (default: 'random_forest')
            test_size (float): Fraction of data to use for testing (default: 0.2)
            cross_validation (bool): Whether to use cross-validation (default: True)
            cv_folds (int): Number of folds for cross-validation (default: 5)
                
        Response:
            Evaluation metrics including:
            - R² score
            - Mean absolute error
            - Mean squared error
            - Root mean squared error
            - Cross-validation metrics (if enabled)
            - Feature importance rankings
            - Learning curve data (if available)
            
        References:
            - Tropsha, A. (2010). Best practices for QSAR model development and validation.
            - Sheridan, R. P. (2013). Time-split cross-validation as a method for estimating the goodness of prospective prediction.
            - Cherkasov, A., et al. (2014). QSAR modeling: where have you been? Where are you going to?
        """
        try:
            # Validate request data
            schema = ModelEvaluationSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Extract parameters
            property_name = data['property_name']
            algorithm = data.get('algorithm', 'random_forest')
            test_size = data.get('test_size', 0.2)
            cross_validation = data.get('cross_validation', True)
            cv_folds = data.get('cv_folds', 5)
            
            # Check if model exists
            model = model_manager.get_model(property_name, algorithm)
            if model.model is None:
                abort(404, message=f"Model for property '{property_name}' with algorithm '{algorithm}' not found or not trained")
            
            # Evaluate model
            try:
                # Get training data
                X, y, feature_names = model_manager._get_training_data(property_name)
                
                if len(X) == 0:
                    abort(400, message=f"No training data available for property '{property_name}'")
                
                # Split data for evaluation
                from sklearn.model_selection import train_test_split, cross_val_score, learning_curve
                from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
                import numpy as np
                
                # Scale features using the model's scaler
                X_scaled = model.scaler.transform(X)
                
                # Split data
                X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=test_size, random_state=42)
                
                # Evaluate on test set
                y_pred = model.model.predict(X_test)
                
                # Calculate metrics
                r2 = r2_score(y_test, y_pred)
                mse = mean_squared_error(y_test, y_pred)
                rmse = np.sqrt(mse)
                mae = mean_absolute_error(y_test, y_pred)
                
                # Prepare results
                metrics = {
                    'property_name': property_name,
                    'algorithm': algorithm,
                    'algorithm_name': ALGORITHMS.get(algorithm, algorithm),
                    'test_size': test_size,
                    'cross_validation': cross_validation,
                    'cv_folds': cv_folds if cross_validation else None,
                    'data_points': len(X),
                    'metrics': {
                        'r2_score': round(r2, 4),
                        'mean_absolute_error': round(mae, 4),
                        'mean_squared_error': round(mse, 4),
                        'root_mean_squared_error': round(rmse, 4)
                    }
                }
                
                # Add feature importance if available
                if hasattr(model, 'feature_importance') and model.feature_importance:
                    # Sort features by importance
                    sorted_features = sorted(
                        model.feature_importance.items(),
                        key=lambda x: abs(x[1]),
                        reverse=True
                    )
                    metrics['feature_importance'] = {
                        'ranked_features': [
                            {'feature': feature, 'importance': round(float(importance), 4)}
                            for feature, importance in sorted_features
                        ]
                    }
                
                # Add cross-validation metrics if requested
                if cross_validation:
                    # Perform cross-validation
                    cv_r2 = cross_val_score(model.model, X_scaled, y, cv=cv_folds, scoring='r2')
                    cv_neg_mse = cross_val_score(model.model, X_scaled, y, cv=cv_folds, scoring='neg_mean_squared_error')
                    cv_rmse = np.sqrt(-cv_neg_mse)
                    cv_neg_mae = cross_val_score(model.model, X_scaled, y, cv=cv_folds, scoring='neg_mean_absolute_error')
                    cv_mae = -cv_neg_mae
                    
                    metrics['cv_metrics'] = {
                        'r2_score_cv': [round(score, 4) for score in cv_r2.tolist()],
                        'mean_absolute_error_cv': [round(score, 4) for score in cv_mae.tolist()],
                        'root_mean_squared_error_cv': [round(score, 4) for score in cv_rmse.tolist()],
                        'r2_score_mean': round(cv_r2.mean(), 4),
                        'r2_score_std': round(cv_r2.std(), 4),
                        'mae_mean': round(cv_mae.mean(), 4),
                        'mae_std': round(cv_mae.std(), 4),
                        'rmse_mean': round(cv_rmse.mean(), 4),
                        'rmse_std': round(cv_rmse.std(), 4)
                    }
                
                # Generate learning curve data if dataset is large enough
                if len(X) >= 20:
                    try:
                        # Calculate learning curve
                        train_sizes, train_scores, test_scores = learning_curve(
                            model.model, X_scaled, y,
                            train_sizes=np.linspace(0.2, 1.0, 5),
                            cv=min(3, cv_folds),  # Use fewer folds for efficiency
                            scoring='r2',
                            random_state=42
                        )
                        
                        # Format learning curve data
                        metrics['learning_curve'] = {
                            'train_sizes': train_sizes.tolist(),
                            'train_scores_mean': np.mean(train_scores, axis=1).tolist(),
                            'train_scores_std': np.std(train_scores, axis=1).tolist(),
                            'test_scores_mean': np.mean(test_scores, axis=1).tolist(),
                            'test_scores_std': np.std(test_scores, axis=1).tolist()
                        }
                        
                        # Add overfitting/underfitting assessment
                        train_score_final = metrics['learning_curve']['train_scores_mean'][-1]
                        test_score_final = metrics['learning_curve']['test_scores_mean'][-1]
                        gap = train_score_final - test_score_final
                        
                        if gap > 0.3:
                            metrics['model_diagnosis'] = "Potential overfitting detected (large gap between training and validation performance)"
                        elif test_score_final < 0.5 and train_score_final < 0.7:
                            metrics['model_diagnosis'] = "Potential underfitting detected (poor performance on both training and validation data)"
                        else:
                            metrics['model_diagnosis'] = "Model appears well-balanced (good generalization)"
                            
                    except Exception as e:
                        current_app.logger.warning(f"Could not generate learning curve: {str(e)}")
                
                return metrics, 200
            except Exception as e:
                current_app.logger.error(f"Error evaluating model: {str(e)}")
                abort(500, message=f"Error evaluating model: {str(e)}")
                
        except Exception as e:
            current_app.logger.error(f"Error processing evaluation request: {str(e)}")
            abort(500, message=f"Error processing evaluation request: {str(e)}")


class VitrificationPredictionResource(MethodResource):
    """Resource for predicting vitrification tendency of a mixture.
    
    Scientific Rationale:
        Vitrification (glass formation) is a critical property for cryoprotectants as it prevents
        damaging ice crystal formation. This endpoint predicts a mixture's tendency to form
        a glass-like amorphous solid rather than crystalline ice during cooling.
        
        The prediction considers:
        1. Glass transition temperature (Tg) estimation
        2. Critical cooling rate prediction
        3. Viscosity effects at low temperatures
        4. Hydrogen bonding network disruption capacity
        5. Water molecule mobility restriction
        
        High vitrification tendency is associated with superior cryoprotection for
        sensitive biological samples and is particularly important for organ preservation.
        
    References:
        - Wowk, B. (2010). Thermodynamic aspects of vitrification.
        - Fahy, G.M., et al. (2009). Physical and biological aspects of renal vitrification.
        - Baudot, A., et al. (2000). Thermal properties of ethylene glycol aqueous solutions.
    """
    
    @doc(description='Predict vitrification tendency for a mixture',
         tags=['Predictive Models'],
         params={
             'mixture_id': {'description': 'ID of the mixture to predict for', 'type': 'string', 'required': True},
             'algorithm': {'description': 'Algorithm to use for prediction', 'type': 'string', 'default': 'random_forest', 'enum': list(ALGORITHMS.keys())}
         })
    @apispec_marshal_with(prediction_result_fields, code=200, description='Prediction results')
    @marshal_with(prediction_result_fields)
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """Get a vitrification tendency prediction for a mixture.
        
        Predicts the vitrification tendency of a mixture using the specified
        machine learning algorithm. Returns the predicted value along with confidence
        information and feature importance.
        
        Args:
            mixture_id (str): ID of the mixture to predict for
            
        Returns:
            Tuple[Dict[str, Any], int]: Prediction results and HTTP status code
            
        Raises:
            400: If validation error occurs or prediction fails
            500: If server error occurs
            
        Query Parameters:
            algorithm (str, optional): Algorithm to use for prediction (default: 'random_forest')
                
        Response:
            Prediction results including:
            - Predicted vitrification tendency
            - Confidence score
            - Confidence interval
            - Feature importance
            - Algorithm information
        """
        try:
            # Validate query parameters
            schema = PredictionRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Extract parameters
            algorithm = args.get('algorithm', 'random_forest')
            
            # Import the prediction function from predictive_models
            from api.predictive_models import predict_vitrification_tendency
            
            # Make prediction
            result = predict_vitrification_tendency(mixture_id, algorithm)
            
            if 'error' in result:
                abort(400, message=result['error'])
            
            return result, 200
        except Exception as e:
            current_app.logger.error(f"Error making vitrification prediction: {str(e)}")
            abort(500, message=f"Error making vitrification prediction: {str(e)}")


class MembraneProtectionPredictionResource(MethodResource):
    """Resource for predicting cell membrane protection capability of a mixture.
    
    Scientific Rationale:
        Cell membrane damage is a primary mechanism of cryoinjury. This endpoint predicts
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
        
    References:
        - Anchordoguy, T.J., et al. (1987). Modes of interaction of cryoprotectants with membrane phospholipids.
        - Bryant, G., et al. (2001). The mechanism of cryoprotection of proteins by solutes.
        - Wolkers, W.F., et al. (2007). Membrane and protein properties of freeze-dried cells.
    """
    
    @doc(description='Predict cell membrane protection capability for a mixture',
         tags=['Predictive Models'],
         params={
             'mixture_id': {'description': 'ID of the mixture to predict for', 'type': 'string', 'required': True},
             'algorithm': {'description': 'Algorithm to use for prediction', 'type': 'string', 'default': 'random_forest', 'enum': list(ALGORITHMS.keys())}
         })
    @apispec_marshal_with(prediction_result_fields, code=200, description='Prediction results')
    @marshal_with(prediction_result_fields)
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """Get a cell membrane protection prediction for a mixture.
        
        Predicts how effectively a mixture will protect cell membranes during
        freezing and thawing cycles. Returns the predicted value along with
        confidence information and feature importance.
        
        Args:
            mixture_id (str): ID of the mixture to predict for
            
        Returns:
            Tuple[Dict[str, Any], int]: Prediction results and HTTP status code
            
        Raises:
            400: If validation error occurs or prediction fails
            500: If server error occurs
            
        Query Parameters:
            algorithm (str, optional): Algorithm to use for prediction (default: 'random_forest')
                
        Response:
            Prediction results including:
            - Predicted membrane protection capability
            - Confidence score
            - Confidence interval
            - Feature importance
            - Algorithm information
        """
        try:
            # Validate query parameters
            schema = PredictionRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Extract parameters
            algorithm = args.get('algorithm', 'random_forest')
            
            # Import the prediction function from predictive_models
            from api.predictive_models import predict_cell_membrane_protection
            
            # Make prediction
            result = predict_cell_membrane_protection(mixture_id, algorithm)
            
            if 'error' in result:
                abort(400, message=result['error'])
            
            return result, 200
        except Exception as e:
            current_app.logger.error(f"Error making membrane protection prediction: {str(e)}")
            abort(500, message=f"Error making membrane protection prediction: {str(e)}")


class ProteinStabilizationPredictionResource(MethodResource):
    """Resource for predicting protein stabilization capability of a mixture.
    
    Scientific Rationale:
        Protein denaturation is a major cause of cellular damage during freezing.
        This endpoint predicts how effectively a cryoprotectant mixture will
        stabilize proteins and prevent denaturation during freezing and thawing.
        
        The prediction considers:
        1. Preferential hydration effects
        2. Preferential exclusion from protein surfaces
        3. Hydrogen bonding with protein residues
        4. Reduction of ice-protein interactions
        5. Prevention of protein aggregation
        
        Effective protein stabilization is critical for preserving enzymatic
        activity and cellular function after cryopreservation.
        
    References:
        - Timasheff, S.N. (1998). Control of protein stability and reactions by weakly interacting cosolvents.
        - Arakawa, T., et al. (2001). Factors affecting short-term and long-term stabilities of proteins.
        - Carpenter, J.F., et al. (1990). Stabilization of phosphofructokinase during freezing.
    """
    
    @doc(description='Predict protein stabilization capability for a mixture',
         tags=['Predictive Models'],
         params={
             'mixture_id': {'description': 'ID of the mixture to predict for', 'type': 'string', 'required': True},
             'algorithm': {'description': 'Algorithm to use for prediction', 'type': 'string', 'default': 'random_forest', 'enum': list(ALGORITHMS.keys())}
         })
    @apispec_marshal_with(prediction_result_fields, code=200, description='Prediction results')
    @marshal_with(prediction_result_fields)
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """Get a protein stabilization prediction for a mixture.
        
        Predicts how effectively a mixture will stabilize proteins during
        freezing and thawing cycles. Returns the predicted value along with
        confidence information and feature importance.
        
        Args:
            mixture_id (str): ID of the mixture to predict for
            
        Returns:
            Tuple[Dict[str, Any], int]: Prediction results and HTTP status code
            
        Raises:
            400: If validation error occurs or prediction fails
            500: If server error occurs
            
        Query Parameters:
            algorithm (str, optional): Algorithm to use for prediction (default: 'random_forest')
                
        Response:
            Prediction results including:
            - Predicted protein stabilization capability
            - Confidence score
            - Confidence interval
            - Feature importance
            - Algorithm information
        """
        try:
            # Validate query parameters
            schema = PredictionRequestSchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Extract parameters
            algorithm = args.get('algorithm', 'random_forest')
            
            # For now, we'll use the general prediction function
            # In a future implementation, this would call a specialized function
            result = model_manager.predict('Protein Stabilization', mixture_id, algorithm)
            
            if 'error' in result:
                abort(400, message=result['error'])
            
            return result, 200
        except Exception as e:
            current_app.logger.error(f"Error making protein stabilization prediction: {str(e)}")
            abort(500, message=f"Error making protein stabilization prediction: {str(e)}")


def register_resources(api):
    """
    Register predictive model resources with the API.
    
    Scientific Rationale:
        This function registers all predictive model endpoints with the provided
        Flask-RESTful API instance. The API design follows RESTful principles and
        provides a comprehensive set of endpoints for:
        
        1. Model management (listing, training, evaluation, deletion)
        2. Property prediction for mixtures (general and property-specific)
        3. Batch operations for efficient model training
        4. Specialized endpoints for cryoprotection-specific properties:
           - Vitrification tendency
           - Cell membrane protection
           - Protein stabilization
        
    The API structure is designed to support both research and practical applications
    in cryoprotectant development and optimization.
    
    Args:
        api: Flask-RESTful API instance
    """
    # Core model management endpoints
    api.add_resource(PredictiveModelListResource, '/api/v1/predictive-models')
    api.add_resource(PredictiveModelResource, '/api/v1/predictive-models/<string:property_name>/<string:algorithm>')
    api.add_resource(BatchModelTrainingResource, '/api/v1/predictive-models/train-batch')
    api.add_resource(ModelEvaluationResource, '/api/v1/predictive-models/evaluate')
    
    # Prediction endpoints
    api.add_resource(MixturePredictionResource, '/api/v1/mixtures/<string:mixture_id>/predict')
    api.add_resource(CustomPropertyPredictionResource, '/api/v1/mixtures/<string:mixture_id>/predict/<string:property_name>')
    
    # Specialized cryoprotection prediction endpoints
    api.add_resource(VitrificationPredictionResource, '/api/v1/mixtures/<string:mixture_id>/predict/vitrification')
    api.add_resource(MembraneProtectionPredictionResource, '/api/v1/mixtures/<string:mixture_id>/predict/membrane-protection')
    api.add_resource(ProteinStabilizationPredictionResource, '/api/v1/mixtures/<string:mixture_id>/predict/protein-stabilization')