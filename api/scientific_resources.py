"""
CryoProtect Analyzer API - Scientific Resources

This module provides API resources for the scientific models in the CryoProtect system,
including concentration-dependent models, temperature-dependent models, and mixture optimization.
"""

from flask import request, g, current_app
from flask_restful import Resource, abort
from typing import Dict, List, Optional, Union, Any, Tuple
import uuid
from datetime import datetime
import json

from api.api_standards import (
    create_standard_response,
    create_error_response, 
    create_success_response
)
from api.api_decorators import (
    standardize_response,
    validate_request_schema,
    rate_limit
)
from api.utils import get_supabase_client, token_required
from config import active_config

# Import scientific models
from scientific_models import (
    ConcentrationModel,
    LinearConcentrationModel,
    ExponentialConcentrationModel,
    SigmoidConcentrationModel,
    TemperatureModel,
    LinearTemperatureModel,
    ArrheniusTemperatureModel,
    GlassTransitionModel,
    MixtureOptimizationModel,
    SynergyPredictionModel,
    ComponentInteractionModel,
    ModelValidationError,
    ModelParameterError,
    ModelCalculationError
)

class ConcentrationModelResource(Resource):
    """Resource for concentration-dependent models."""
    
    @token_required
    def get(self, molecule_id: str = None) -> Tuple[Dict[str, Any], int]:
        """
        Get concentration models for a molecule.
        
        Args:
            molecule_id: UUID of the molecule (optional)
            
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            supabase = get_supabase_client()
            
            # Get models for a specific molecule or all models
            if molecule_id:
                # Query models for the specific molecule
                result = (
                    supabase.table('concentration_models')
                    .select('*')
                    .eq('molecule_id', molecule_id)
                    .execute()
                )
                
                if not hasattr(result, 'data') or not result.data:
                    return create_error_response(
                        error=f"No concentration models found for molecule {molecule_id}",
                        status_code=404,
                        context="Concentration models"
                    )
                
                # Get molecule details
                molecule_result = (
                    supabase.table('molecules')
                    .select('name, smiles')
                    .eq('id', molecule_id)
                    .single()
                    .execute()
                )
                
                molecule_name = (
                    molecule_result.data.get('name', "Unknown")
                    if hasattr(molecule_result, 'data') and molecule_result.data
                    else "Unknown"
                )
                
                # Return models with molecule information
                return create_success_response(
                    data={
                        "molecule_id": molecule_id,
                        "molecule_name": molecule_name,
                        "models": result.data
                    },
                    message=f"Retrieved {len(result.data)} concentration models for molecule {molecule_id}",
                    status_code=200
                )
            else:
                # Query all models with pagination
                page = request.args.get('page', 1, type=int)
                per_page = request.args.get('per_page', 20, type=int)
                
                # Validate and constrain pagination values
                page = max(1, page)
                per_page = max(1, min(100, per_page))
                
                # Query with pagination
                result = (
                    supabase.table('concentration_models')
                    .select('*')
                    .order('created_at', desc=True)
                    .range((page - 1) * per_page, page * per_page - 1)
                    .execute()
                )
                
                # Count total
                count_result = (
                    supabase.table('concentration_models')
                    .select('id', count='exact')
                    .execute()
                )
                
                total_count = count_result.count if hasattr(count_result, 'count') else 0
                
                # Calculate pagination values
                total_pages = (total_count + per_page - 1) // per_page if per_page > 0 else 0
                
                # Prepare model data
                models_data = result.data if hasattr(result, 'data') else []
                
                # Return all models with pagination
                return create_success_response(
                    data={
                        "models": models_data,
                        "count": len(models_data),
                        "total": total_count
                    },
                    message=f"Retrieved {len(models_data)} concentration models",
                    status_code=200,
                    pagination={
                        "page": page,
                        "per_page": per_page,
                        "total_items": total_count,
                        "total_pages": total_pages
                    }
                )
        except Exception as e:
            current_app.logger.error(f"Error retrieving concentration models: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Concentration models"
            )
    
    @token_required
    @validate_request_schema({
        'type': 'object',
        'properties': {
            'molecule_id': {'type': 'string'},
            'model_type': {'type': 'string'},
            'parameters': {'type': 'object'},
            'valid_range': {
                'type': 'array',
                'items': {'type': 'number'},
                'minItems': 2,
                'maxItems': 2
            },
            'name': {'type': 'string'},
            'description': {'type': 'string'}
        },
        'required': ['molecule_id', 'model_type', 'parameters']
    })
    def post(self) -> Tuple[Dict[str, Any], int]:
        """
        Create a new concentration model.
        
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            data = request.get_json()
            
            # Extract model data
            molecule_id = data['molecule_id']
            model_type = data['model_type']
            parameters = data['parameters']
            valid_range = data.get('valid_range', [0, 10])
            name = data.get('name', f"{model_type.capitalize()} concentration model")
            description = data.get('description', f"Concentration-dependent model using {model_type} approach")
            
            # Validate molecule exists
            supabase = get_supabase_client()
            molecule_result = (
                supabase.table('molecules')
                .select('id')
                .eq('id', molecule_id)
                .single()
                .execute()
            )
            
            if not hasattr(molecule_result, 'data') or not molecule_result.data:
                return create_error_response(
                    error=f"Molecule with ID {molecule_id} not found",
                    status_code=404,
                    context="Concentration models"
                )
            
            # Validate model parameters
            try:
                # Create model instance to validate parameters
                if model_type == 'linear':
                    model = LinearConcentrationModel(
                        parameters=parameters,
                        name=name,
                        description=description,
                        valid_range=(valid_range[0], valid_range[1])
                    )
                elif model_type == 'exponential':
                    model = ExponentialConcentrationModel(
                        parameters=parameters,
                        name=name,
                        description=description,
                        valid_range=(valid_range[0], valid_range[1])
                    )
                elif model_type == 'sigmoid':
                    model = SigmoidConcentrationModel(
                        parameters=parameters,
                        name=name,
                        description=description,
                        valid_range=(valid_range[0], valid_range[1])
                    )
                else:
                    return create_error_response(
                        error=f"Invalid model type: {model_type}",
                        status_code=400,
                        context="Concentration models"
                    )
                
                # Validate parameters
                model.validate_parameters()
                
            except (ModelValidationError, ModelParameterError) as e:
                return create_error_response(
                    error=str(e),
                    status_code=400,
                    context="Model validation"
                )
            
            # Create model record in database
            model_data = {
                'molecule_id': molecule_id,
                'model_type': model_type,
                'parameters': parameters,
                'valid_range': f"[{valid_range[0]},{valid_range[1]})",
                'name': name,
                'description': description,
                'created_by': g.user['email'] if hasattr(g, 'user') and g.user else None,
                'created_at': datetime.now().isoformat()
            }
            
            result = (
                supabase.table('concentration_models')
                .insert(model_data)
                .execute()
            )
            
            if not hasattr(result, 'data') or not result.data:
                return create_error_response(
                    error="Failed to create concentration model",
                    status_code=500,
                    context="Concentration models"
                )
            
            created_model = result.data[0] if hasattr(result, 'data') and result.data else None
            
            return create_success_response(
                data=created_model,
                message=f"Created new {model_type} concentration model for molecule {molecule_id}",
                status_code=201
            )
            
        except Exception as e:
            current_app.logger.error(f"Error creating concentration model: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Concentration models"
            )
    
    @token_required
    def delete(self, model_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Delete a concentration model.
        
        Args:
            model_id: UUID of the model to delete
            
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            supabase = get_supabase_client()
            
            # Check if model exists
            model_result = (
                supabase.table('concentration_models')
                .select('*')
                .eq('id', model_id)
                .single()
                .execute()
            )
            
            if not hasattr(model_result, 'data') or not model_result.data:
                return create_error_response(
                    error=f"Concentration model with ID {model_id} not found",
                    status_code=404,
                    context="Concentration models"
                )
            
            # Delete the model
            result = (
                supabase.table('concentration_models')
                .delete()
                .eq('id', model_id)
                .execute()
            )
            
            if not hasattr(result, 'data') or not result.data:
                return create_error_response(
                    error=f"Failed to delete concentration model {model_id}",
                    status_code=500,
                    context="Concentration models"
                )
            
            return create_success_response(
                data=None,
                message=f"Deleted concentration model {model_id}",
                status_code=200
            )
            
        except Exception as e:
            current_app.logger.error(f"Error deleting concentration model: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Concentration models"
            )
    
    @token_required
    @validate_request_schema({
        'type': 'object',
        'properties': {
            'concentration': {'type': 'number', 'minimum': 0},
            'temperature': {'type': 'number', 'minimum': 0},
            'temperature_unit': {'type': 'string', 'enum': ['K', 'C']}
        },
        'required': ['concentration']
    })
    def post(self, model_id: str, action: str = 'calculate') -> Tuple[Dict[str, Any], int]:
        """
        Calculate property value using a concentration model.
        
        Args:
            model_id: UUID of the model to use
            action: Action to perform (default: calculate)
            
        Returns:
            Tuple of (response_data, status_code)
        """
        if action != 'calculate':
            return create_error_response(
                error=f"Invalid action: {action}",
                status_code=400,
                context="Concentration models"
            )
        
        try:
            data = request.get_json()
            concentration = data['concentration']
            
            # Get model from database
            supabase = get_supabase_client()
            model_result = (
                supabase.table('concentration_models')
                .select('*')
                .eq('id', model_id)
                .single()
                .execute()
            )
            
            if not hasattr(model_result, 'data') or not model_result.data:
                return create_error_response(
                    error=f"Concentration model with ID {model_id} not found",
                    status_code=404,
                    context="Concentration models"
                )
            
            model_data = model_result.data
            
            # Create model instance
            try:
                # Parse valid range
                valid_range_str = model_data.get('valid_range')
                valid_range = (0, 10)  # Default
                
                if valid_range_str:
                    try:
                        # Remove brackets and split
                        range_values = valid_range_str.strip('()[]').split(',')
                        valid_range = (float(range_values[0]), float(range_values[1]))
                    except (ValueError, IndexError):
                        current_app.logger.warning(f"Invalid valid_range format: {valid_range_str}")
                
                # Create model based on type
                model_type = model_data.get('model_type')
                parameters = model_data.get('parameters', {})
                name = model_data.get('name')
                description = model_data.get('description')
                
                if model_type == 'linear':
                    model = LinearConcentrationModel(
                        parameters=parameters,
                        name=name,
                        description=description,
                        valid_range=valid_range
                    )
                elif model_type == 'exponential':
                    model = ExponentialConcentrationModel(
                        parameters=parameters,
                        name=name,
                        description=description,
                        valid_range=valid_range
                    )
                elif model_type == 'sigmoid':
                    model = SigmoidConcentrationModel(
                        parameters=parameters,
                        name=name,
                        description=description,
                        valid_range=valid_range
                    )
                else:
                    return create_error_response(
                        error=f"Unsupported model type: {model_type}",
                        status_code=400,
                        context="Concentration models"
                    )
                
                # Calculate property value
                result = model.calculate({
                    'concentration': concentration,
                    'molecule_id': model_data.get('molecule_id')
                })
                
                return create_success_response(
                    data=result,
                    message=f"Calculated property value for concentration {concentration}",
                    status_code=200
                )
                
            except (ModelValidationError, ModelParameterError, ModelCalculationError) as e:
                return create_error_response(
                    error=str(e),
                    status_code=400,
                    context="Model calculation"
                )
            
        except Exception as e:
            current_app.logger.error(f"Error calculating with concentration model: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Concentration models"
            )


class TemperatureModelResource(Resource):
    """Resource for temperature-dependent models."""
    
    @token_required
    def get(self, molecule_id: str = None) -> Tuple[Dict[str, Any], int]:
        """
        Get temperature models for a molecule.
        
        Args:
            molecule_id: UUID of the molecule (optional)
            
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            supabase = get_supabase_client()
            
            # Get models for a specific molecule or all models
            if molecule_id:
                # Query models for the specific molecule
                result = (
                    supabase.table('temperature_models')
                    .select('*')
                    .eq('molecule_id', molecule_id)
                    .execute()
                )
                
                if not hasattr(result, 'data') or not result.data:
                    return create_error_response(
                        error=f"No temperature models found for molecule {molecule_id}",
                        status_code=404,
                        context="Temperature models"
                    )
                
                # Get molecule details
                molecule_result = (
                    supabase.table('molecules')
                    .select('name, smiles')
                    .eq('id', molecule_id)
                    .single()
                    .execute()
                )
                
                molecule_name = (
                    molecule_result.data.get('name', "Unknown")
                    if hasattr(molecule_result, 'data') and molecule_result.data
                    else "Unknown"
                )
                
                # Return models with molecule information
                return create_success_response(
                    data={
                        "molecule_id": molecule_id,
                        "molecule_name": molecule_name,
                        "models": result.data
                    },
                    message=f"Retrieved {len(result.data)} temperature models for molecule {molecule_id}",
                    status_code=200
                )
            else:
                # Query all models with pagination
                page = request.args.get('page', 1, type=int)
                per_page = request.args.get('per_page', 20, type=int)
                
                # Validate and constrain pagination values
                page = max(1, page)
                per_page = max(1, min(100, per_page))
                
                # Query with pagination
                result = (
                    supabase.table('temperature_models')
                    .select('*')
                    .order('created_at', desc=True)
                    .range((page - 1) * per_page, page * per_page - 1)
                    .execute()
                )
                
                # Count total
                count_result = (
                    supabase.table('temperature_models')
                    .select('id', count='exact')
                    .execute()
                )
                
                total_count = count_result.count if hasattr(count_result, 'count') else 0
                
                # Calculate pagination values
                total_pages = (total_count + per_page - 1) // per_page if per_page > 0 else 0
                
                # Prepare model data
                models_data = result.data if hasattr(result, 'data') else []
                
                # Return all models with pagination
                return create_success_response(
                    data={
                        "models": models_data,
                        "count": len(models_data),
                        "total": total_count
                    },
                    message=f"Retrieved {len(models_data)} temperature models",
                    status_code=200,
                    pagination={
                        "page": page,
                        "per_page": per_page,
                        "total_items": total_count,
                        "total_pages": total_pages
                    }
                )
        except Exception as e:
            current_app.logger.error(f"Error retrieving temperature models: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Temperature models"
            )


class GlassTransitionResource(Resource):
    """Resource for glass transition temperature predictions."""
    
    @token_required
    @validate_request_schema({
        'type': 'object',
        'properties': {
            'molecule_id': {'type': 'string'},
            'smiles': {'type': 'string'},
            'components': {
                'type': 'array',
                'items': {
                    'type': 'object',
                    'properties': {
                        'weight_fraction': {'type': 'number', 'minimum': 0, 'maximum': 1},
                        'tg': {'type': 'number', 'minimum': 0}
                    },
                    'required': ['weight_fraction', 'tg']
                }
            }
        }
    })
    def post(self) -> Tuple[Dict[str, Any], int]:
        """
        Predict glass transition temperature.
        
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            data = request.get_json()
            
            # Create model
            model = GlassTransitionModel(
                parameters={
                    'model_subtype': 'gordon_taylor',
                    'k': 1.0  # Default Gordon-Taylor parameter
                }
            )
            
            # Calculate Tg
            result = model.calculate(data)
            
            return create_success_response(
                data=result,
                message="Glass transition temperature prediction successful",
                status_code=200
            )
            
        except (ModelValidationError, ModelParameterError, ModelCalculationError) as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Glass transition prediction"
            )
        except Exception as e:
            current_app.logger.error(f"Error predicting glass transition temperature: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Glass transition prediction"
            )


class MixtureOptimizerResource(Resource):
    """Resource for mixture optimization."""
    
    @token_required
    @validate_request_schema({
        'type': 'object',
        'properties': {
            'candidate_molecules': {
                'type': 'array',
                'items': {'type': 'object'}
            },
            'initial_mixture': {'type': 'object'},
            'constraints': {'type': 'object'},
            'algorithm': {'type': 'string'},
            'optimization_parameters': {'type': 'object'}
        },
        'required': ['candidate_molecules']
    })
    def post(self) -> Tuple[Dict[str, Any], int]:
        """
        Optimize a cryoprotectant mixture.
        
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            data = request.get_json()
            
            # Extract optimization parameters
            candidate_molecules = data['candidate_molecules']
            initial_mixture = data.get('initial_mixture')
            constraints = data.get('constraints', {})
            algorithm = data.get('algorithm', 'genetic')
            optimization_parameters = data.get('optimization_parameters', {})
            
            # Create optimizer
            parameters = {
                'algorithm': algorithm,
                'constraints': constraints,
                'objective_function': self._dummy_objective_function  # Simplified for demo
            }
            
            # Add algorithm-specific parameters
            if algorithm == 'genetic':
                parameters['ga_params'] = {
                    'population_size': optimization_parameters.get('population_size', 50),
                    'generations': optimization_parameters.get('generations', 100),
                    'mutation_rate': optimization_parameters.get('mutation_rate', 0.1),
                    'crossover_rate': optimization_parameters.get('crossover_rate', 0.8)
                }
            elif algorithm == 'grid_search':
                parameters['grid_params'] = {
                    'grid_points': optimization_parameters.get('grid_points', 5)
                }
            
            optimizer = MixtureOptimizationModel(parameters=parameters)
            
            # Run optimization
            result = optimizer.calculate({
                'candidate_molecules': candidate_molecules,
                'initial_mixture': initial_mixture,
                'constraints': constraints
            })
            
            # Store optimization result in database
            try:
                supabase = get_supabase_client()
                
                optimization_data = {
                    'name': data.get('name', f"Optimization_{datetime.now().strftime('%Y%m%d_%H%M%S')}"),
                    'components': result['optimized_mixture']['components'],
                    'optimization_parameters': {
                        'algorithm': algorithm,
                        'constraints': constraints,
                        'parameters': optimization_parameters
                    },
                    'effectiveness_score': result['fitness'],
                    'created_by': g.user['email'] if hasattr(g, 'user') and g.user else None,
                    'created_at': datetime.now().isoformat(),
                    'notes': data.get('notes', '')
                }
                
                db_result = (
                    supabase.table('mixture_optimizations')
                    .insert(optimization_data)
                    .execute()
                )
                
                if hasattr(db_result, 'data') and db_result.data:
                    result['optimization_id'] = db_result.data[0]['id']
            except Exception as db_error:
                current_app.logger.error(f"Error storing optimization result: {str(db_error)}", exc_info=True)
                # Continue even if database storage fails
            
            return create_success_response(
                data=result,
                message=f"Mixture optimization completed using {algorithm} algorithm",
                status_code=200
            )
            
        except (ModelValidationError, ModelParameterError, ModelCalculationError) as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Mixture optimization"
            )
        except Exception as e:
            current_app.logger.error(f"Error optimizing mixture: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Mixture optimization"
            )
    
    @staticmethod
    def _dummy_objective_function(mixture: Dict[str, Any]) -> float:
        """
        Dummy objective function for demonstration purposes.
        
        In a real implementation, this would evaluate cryoprotectant effectiveness.
        
        Args:
            mixture: Mixture to evaluate
            
        Returns:
            Effectiveness score
        """
        # Simple scoring based on number of components and concentration
        components = mixture.get('components', [])
        n_components = len(components)
        
        # Basic score: more components is better, up to a point
        base_score = min(n_components, 3) / 3.0
        
        # Adjust based on total concentration
        total_concentration = sum(c.get('concentration', 0) for c in components)
        if total_concentration < 1.0:
            # Too dilute
            concentration_factor = total_concentration
        elif total_concentration > 10.0:
            # Too concentrated
            concentration_factor = max(0, 1.0 - (total_concentration - 10.0) / 10.0)
        else:
            # Ideal range: 1-10 M
            concentration_factor = 1.0
        
        score = base_score * concentration_factor
        
        # Add some randomness for demonstration
        import random
        score += random.uniform(-0.1, 0.1)
        
        # Constrain to [0, 1] range
        score = max(0, min(1, score))
        
        return score


class SynergyPredictorResource(Resource):
    """Resource for synergy predictions between cryoprotectant components."""
    
    @token_required
    @validate_request_schema({
        'type': 'object',
        'properties': {
            'components': {
                'type': 'array',
                'items': {'type': 'object'}
            },
            'synergy_mode': {'type': 'string'},
            'parameters': {'type': 'object'}
        },
        'required': ['components']
    })
    def post(self) -> Tuple[Dict[str, Any], int]:
        """
        Predict synergistic effects between components.
        
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            data = request.get_json()
            
            # Extract parameters
            components = data['components']
            synergy_mode = data.get('synergy_mode', 'property_based')
            parameters = data.get('parameters', {})
            
            # Create predictor
            predictor_params = {
                'synergy_mode': synergy_mode
            }
            
            # Add mode-specific parameters
            if synergy_mode == 'property_based':
                predictor_params['property_weights'] = parameters.get('property_weights', {
                    'molecular_weight': 0.3,
                    'logP': 0.5,
                    'hydrogen_bond_donors': 0.2
                })
            elif synergy_mode == 'mechanism_based':
                predictor_params['mechanisms'] = parameters.get('mechanisms', [
                    'vitrification',
                    'cell_permeation',
                    'antioxidant',
                    'membrane_stabilization'
                ])
            elif synergy_mode == 'hybrid':
                predictor_params['property_weight'] = parameters.get('property_weight', 0.3)
                predictor_params['mechanism_weight'] = parameters.get('mechanism_weight', 0.3)
                predictor_params['empirical_weight'] = parameters.get('empirical_weight', 0.4)
            
            predictor = SynergyPredictionModel(parameters=predictor_params)
            
            # Calculate synergy
            result = predictor.calculate({
                'components': components
            })
            
            return create_success_response(
                data=result,
                message=f"Synergy prediction completed using {synergy_mode} mode",
                status_code=200
            )
            
        except (ModelValidationError, ModelParameterError, ModelCalculationError) as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Synergy prediction"
            )
        except Exception as e:
            current_app.logger.error(f"Error predicting synergy: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Synergy prediction"
            )


def register_scientific_resources(api):
    """
    Register scientific resources with the API.
    
    Args:
        api: Flask-RESTful API instance
    """
    # Concentration model resources
    api.add_resource(
        ConcentrationModelResource,
        '/api/v1/scientific/concentration-models',
        '/api/v1/scientific/concentration-models/<string:molecule_id>',
        '/api/v1/scientific/concentration-models/<string:model_id>/<string:action>',
        endpoint='concentration_models'
    )
    
    # Temperature model resources
    api.add_resource(
        TemperatureModelResource,
        '/api/v1/scientific/temperature-models',
        '/api/v1/scientific/temperature-models/<string:molecule_id>',
        endpoint='temperature_models'
    )
    
    # Glass transition resource
    api.add_resource(
        GlassTransitionResource,
        '/api/v1/scientific/glass-transition',
        endpoint='glass_transition'
    )
    
    # Mixture optimizer resource
    api.add_resource(
        MixtureOptimizerResource,
        '/api/v1/scientific/mixture-optimizer',
        endpoint='mixture_optimizer'
    )
    
    # Synergy predictor resource
    api.add_resource(
        SynergyPredictorResource,
        '/api/v1/scientific/synergy-predictor',
        endpoint='synergy_predictor'
    )


def register_scientific_docs(docs):
    """
    Register scientific resources for documentation.
    
    Args:
        docs: FlaskApiSpec instance
    """
    from api.api_docs import register_resource
    
    register_resource(docs, ConcentrationModelResource, 'concentration_models')
    register_resource(docs, TemperatureModelResource, 'temperature_models')
    register_resource(docs, GlassTransitionResource, 'glass_transition')
    register_resource(docs, MixtureOptimizerResource, 'mixture_optimizer')
    register_resource(docs, SynergyPredictorResource, 'synergy_predictor')

# Aliases for backward compatibility and validation
register_resources = register_scientific_resources