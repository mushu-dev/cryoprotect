"""
CryoProtect Analyzer API - Mixture Analysis Resources

This module provides API endpoints for mixture analysis capabilities, including:
- Predicting mixture properties
- Analyzing component compatibility
- Detecting synergistic effects
- Optimizing mixture compositions
- Providing recommendations for improvement

All endpoints follow RESTful design principles and return standardized
response formats. Authentication is handled via JWT tokens for protected
endpoints.
"""

import logging
from flask import request, jsonify
from flask_restful import Resource, marshal_with, abort, fields
from marshmallow import ValidationError
from typing import Dict, List, Any, Tuple, Optional

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

from api.utils import (
    _handle_json_serialization, token_required,
    handle_error, handle_supabase_error, get_user_id
)
from api.mixture_analysis import (
    MixtureProperty, MixtureCompatibility, MixtureSynergy,
    MixtureOptimization, MixtureRecommendation
)
from api.models import Mixture, FlexibleDateTime

# Set up logging
logger = logging.getLogger(__name__)

# Define field schemas for marshal_with decorators
mixture_properties_fields = {
    'mixture_id': fields.String,
    'mixture_name': fields.String,
    'properties': fields.Raw,
    'raw_properties': fields.Raw
}

mixture_compatibility_fields = {
    'mixture_id': fields.String,
    'mixture_name': fields.String,
    'compatibility': fields.Raw
}

mixture_synergy_fields = {
    'mixture_id': fields.String,
    'mixture_name': fields.String,
    'synergy': fields.Raw
}

mixture_optimization_fields = {
    'mixture_id': fields.String,
    'mixture_name': fields.String,
    'original_composition': fields.Raw,
    'optimized_composition': fields.Raw,
    'target_property': fields.String,
    'target_value': fields.Float,
    'achieved_value': fields.Float,
    'improvement': fields.Float
}

mixture_step_optimization_fields = {
    'mixture_id': fields.String,
    'mixture_name': fields.String,
    'steps': fields.Raw,
    'target_property': fields.String,
    'target_value': fields.Float,
    'final_value': fields.Float
}

mixture_analysis_fields = {
    'mixture_id': fields.String,
    'mixture_name': fields.String,
    'analysis': fields.Raw,
    'recommendations': fields.Raw
}

mixture_component_recommendation_fields = {
    'mixture_id': fields.String,
    'mixture_name': fields.String,
    'target_property': fields.String,
    'target_value': fields.Float,
    'recommended_components': fields.Raw
}

class MixturePropertiesResource(MethodResource):
    """Resource for predicting mixture properties.
    
    This endpoint calculates and returns predicted physicochemical properties
    for a mixture based on its components. Properties include viscosity,
    density, freezing point depression, and other relevant cryoprotection
    characteristics.
    """
    
    @doc(description='Get predicted properties for a mixture',
         tags=['Mixture Analysis'],
         params={'mixture_id': {'description': 'ID of the mixture to analyze', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}])
    @apispec_marshal_with(mixture_properties_fields, code=200, description='Predicted mixture properties')
    @token_required
    @marshal_with(mixture_properties_fields)
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Get predicted properties for a mixture.
        
        Calculates various physicochemical properties of a mixture based on its
        components and their concentrations. Returns both the predicted properties
        and the raw calculated values.
        
        Args:
            mixture_id (str): ID of the mixture to analyze
            
        Returns:
            Tuple[Dict[str, Any], int]: JSON response with predicted properties and HTTP status code
            
        Raises:
            404: If mixture not found or has no components
            500: If calculation error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
        """
        try:
            # Get mixture with components
            mixture = Mixture.get_with_components(mixture_id)
            if not mixture or not mixture.get("components"):
                error_response, error_status = handle_error(
                    "Mixture not found or has no components",
                    context="Fetching mixture properties",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
                
            components = mixture["components"]
            
            # Calculate properties
            properties = MixtureProperty.predict_mixture_properties(components)
            raw_properties = MixtureProperty.calculate_raw_properties(components)
            
            return {
                "mixture_id": mixture_id,
                "mixture_name": mixture.get("name", ""),
                "properties": properties,
                "raw_properties": raw_properties
            }
            
        except Exception as e:
            # Standard error handling
            error_response, error_status = handle_error(
                e,
                context="Predicting mixture properties",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MixtureCompatibilityResource(MethodResource):
    """Resource for analyzing mixture compatibility.
    
    This endpoint analyzes the chemical and physical compatibility between
    components in a mixture. It identifies potential interactions, conflicts,
    or stability issues that might affect the mixture's performance as a
    cryoprotectant.
    """
    
    @doc(description='Analyze compatibility between mixture components',
         tags=['Mixture Analysis'],
         params={'mixture_id': {'description': 'ID of the mixture to analyze', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}])
    @apispec_marshal_with(mixture_compatibility_fields, code=200, description='Compatibility analysis results')
    @token_required
    @marshal_with(mixture_compatibility_fields)
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Analyze compatibility between mixture components.
        
        Evaluates chemical and physical compatibility between all components in the
        mixture, identifying potential adverse interactions or stability issues.
        
        Args:
            mixture_id (str): ID of the mixture to analyze
            
        Returns:
            Tuple[Dict[str, Any], int]: JSON response with compatibility analysis and HTTP status code
            
        Raises:
            404: If mixture not found or has no components
            500: If analysis error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
        """
        try:
            # Get mixture with components
            mixture = Mixture.get_with_components(mixture_id)
            if not mixture or not mixture.get("components"):
                error_response, error_status = handle_error(
                    "Mixture not found or has no components",
                    context="Analyzing mixture compatibility",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
                
            components = mixture["components"]
            
            # Analyze compatibility
            compatibility = MixtureCompatibility.analyze_compatibility(components)
            
            return {
                "mixture_id": mixture_id,
                "mixture_name": mixture.get("name", ""),
                "compatibility": compatibility
            }
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Analyzing mixture compatibility",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MixtureSynergyResource(MethodResource):
    """Resource for analyzing mixture synergy.
    
    This endpoint evaluates potential synergistic or antagonistic effects
    between components in a mixture. It identifies combinations that may
    enhance or diminish the overall cryoprotective effectiveness beyond
    what would be expected from individual components.
    """
    
    @doc(description='Analyze synergistic or antagonistic effects in a mixture',
         tags=['Mixture Analysis'],
         params={'mixture_id': {'description': 'ID of the mixture to analyze', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}])
    @apispec_marshal_with(mixture_synergy_fields, code=200, description='Synergy analysis results')
    @token_required
    @marshal_with(mixture_synergy_fields)
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Analyze synergistic or antagonistic effects in a mixture.
        
        Evaluates how components interact to produce effects that are greater than
        (synergistic) or less than (antagonistic) the sum of their individual effects.
        
        Args:
            mixture_id (str): ID of the mixture to analyze
            
        Returns:
            Tuple[Dict[str, Any], int]: JSON response with synergy analysis and HTTP status code
            
        Raises:
            404: If mixture not found or analysis fails
            500: If server error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
        """
        try:
            # Analyze synergy
            synergy = MixtureSynergy.analyze_synergy(mixture_id)
            
            if "error" in synergy:
                error_response, error_status = handle_error(
                    synergy["error"],
                    context="Analyzing mixture synergy",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
                
            # Get mixture name
            mixture = Mixture.get(mixture_id)
            mixture_name = mixture.get("name", "") if mixture else ""
            
            return {
                "mixture_id": mixture_id,
                "mixture_name": mixture_name,
                "synergy": synergy
            }
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Analyzing mixture synergy",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MixtureOptimizationResource(MethodResource):
    """Resource for optimizing mixture compositions.
    
    This endpoint optimizes the composition of a mixture to achieve a target
    property value or maximize overall effectiveness. It adjusts component
    concentrations while respecting specified constraints.
    """
    
    @doc(description='Optimize the composition of a mixture',
         tags=['Mixture Analysis'],
         params={'mixture_id': {'description': 'ID of the mixture to optimize', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}],
         request_body={
             'target_property': {'description': 'Property to optimize for', 'type': 'string', 'default': 'Cryoprotection Score'},
             'target_value': {'description': 'Target value for the property', 'type': 'number'},
             'constraints': {'description': 'Constraints for optimization', 'type': 'object'}
         })
    @apispec_marshal_with(mixture_optimization_fields, code=200, description='Optimized mixture composition')
    @token_required
    @marshal_with(mixture_optimization_fields)
    def post(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Optimize the composition of a mixture.
        
        Adjusts component concentrations to achieve a target property value or
        maximize overall effectiveness while respecting specified constraints.
        
        Args:
            mixture_id (str): ID of the mixture to optimize
            
        Returns:
            Tuple[Dict[str, Any], int]: JSON response with optimized composition and HTTP status code
            
        Raises:
            400: If validation error occurs
            404: If mixture not found
            500: If optimization error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            target_property (str, optional): Property to optimize for. Default: "Cryoprotection Score"
            target_value (float, optional): Target value for the property
            constraints (dict, optional): Constraints for optimization, such as:
                - min_concentration: Minimum concentration for components
                - max_concentration: Maximum concentration for components
                - fixed_components: List of components that should not be modified
        """
        try:
            # Get request data
            data = request.get_json() or {}
            
            # Extract parameters
            target_property = data.get("target_property", "Cryoprotection Score")
            target_value = data.get("target_value")
            constraints = data.get("constraints", {})
            
            # Optimize composition
            result = MixtureOptimization.optimize_composition(
                mixture_id, target_property, target_value, constraints
            )
            
            if "error" in result:
                error_response, error_status = handle_error(
                    result["error"],
                    context="Optimizing mixture composition",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
                
            return result
            
        except ValidationError as err:
            error_response, error_status = handle_error(
                err,
                context="Validating optimization parameters",
                log_level='error',
                return_response=True,
                status_code=400
            )
            abort(error_status, **error_response)
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Optimizing mixture composition",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MixtureStepOptimizationResource(MethodResource):
    """Resource for step-by-step mixture optimization.
    
    This endpoint performs a step-by-step optimization of a mixture's composition,
    showing the progression of changes and their effects on target properties.
    It provides transparency into the optimization process and allows for
    better understanding of how each change affects the mixture.
    """
    
    @doc(description='Optimize the composition of a mixture step by step',
         tags=['Mixture Analysis'],
         params={'mixture_id': {'description': 'ID of the mixture to optimize', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}],
         request_body={
             'target_property': {'description': 'Property to optimize for', 'type': 'string', 'default': 'Cryoprotection Score'},
             'target_value': {'description': 'Target value for the property', 'type': 'number'},
             'constraints': {'description': 'Constraints for optimization', 'type': 'object'}
         })
    @apispec_marshal_with(mixture_step_optimization_fields, code=200, description='Step-by-step optimization results')
    @token_required
    @marshal_with(mixture_step_optimization_fields)
    def post(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Optimize the composition of a mixture step by step.
        
        Performs incremental optimization of the mixture composition, showing
        each step of the process and how it affects the target property.
        
        Args:
            mixture_id (str): ID of the mixture to optimize
            
        Returns:
            Tuple[Dict[str, Any], int]: JSON response with step-by-step optimization results and HTTP status code
            
        Raises:
            400: If validation error occurs
            404: If mixture not found
            500: If optimization error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            target_property (str, optional): Property to optimize for. Default: "Cryoprotection Score"
            target_value (float, optional): Target value for the property
            constraints (dict, optional): Constraints for optimization, such as:
                - min_concentration: Minimum concentration for components
                - max_concentration: Maximum concentration for components
                - fixed_components: List of components that should not be modified
        """
        try:
            # Get request data
            data = request.get_json() or {}
            
            # Extract parameters
            target_property = data.get("target_property", "Cryoprotection Score")
            target_value = data.get("target_value")
            constraints = data.get("constraints", {})
            
            # Optimize composition
            result = MixtureOptimization.optimize_step_by_step(
                mixture_id, target_property, target_value, constraints
            )
            
            if "error" in result:
                error_response, error_status = handle_error(
                    result["error"],
                    context="Step-by-step mixture optimization",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
                
            return result
            
        except ValidationError as err:
            error_response, error_status = handle_error(
                err,
                context="Validating step optimization parameters",
                log_level='error',
                return_response=True,
                status_code=400
            )
            abort(error_status, **error_response)
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Step-by-step mixture optimization",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MixtureAnalysisResource(MethodResource):
    """Resource for analyzing mixtures and providing recommendations.
    
    This endpoint performs a comprehensive analysis of a mixture and provides
    recommendations for improving its cryoprotective effectiveness. It evaluates
    various aspects including component interactions, property balance, and
    potential areas for enhancement.
    """
    
    @doc(description='Analyze a mixture and provide recommendations for improvement',
         tags=['Mixture Analysis'],
         params={'mixture_id': {'description': 'ID of the mixture to analyze', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}])
    @apispec_marshal_with(mixture_analysis_fields, code=200, description='Mixture analysis and recommendations')
    @token_required
    @marshal_with(mixture_analysis_fields)
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Analyze a mixture and provide recommendations for improvement.
        
        Performs a comprehensive analysis of the mixture's properties, component
        interactions, and effectiveness as a cryoprotectant. Provides specific
        recommendations for improving its performance.
        
        Args:
            mixture_id (str): ID of the mixture to analyze
            
        Returns:
            Tuple[Dict[str, Any], int]: JSON response with analysis and recommendations and HTTP status code
            
        Raises:
            404: If mixture not found
            500: If analysis error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
        """
        try:
            # Analyze mixture
            result = MixtureRecommendation.analyze_mixture(mixture_id)
            
            if "error" in result:
                error_response, error_status = handle_error(
                    result["error"],
                    context="Analyzing mixture",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
                
            return result
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Analyzing mixture",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MixtureComponentRecommendationResource(MethodResource):
    """Resource for recommending new components for a mixture.
    
    This endpoint suggests new components that could be added to a mixture
    to improve its properties or effectiveness. It analyzes the current
    composition and identifies molecules that would complement or enhance
    the mixture's cryoprotective capabilities.
    """
    
    @doc(description='Recommend new components to add to a mixture',
         tags=['Mixture Analysis'],
         params={'mixture_id': {'description': 'ID of the mixture to improve', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}],
         request_body={
             'target_property': {'description': 'Property to improve', 'type': 'string'},
             'target_value': {'description': 'Target value for the property', 'type': 'number'},
             'count': {'description': 'Number of recommendations to return', 'type': 'integer', 'default': 3}
         })
    @apispec_marshal_with(mixture_component_recommendation_fields, code=200, description='Component recommendations')
    @token_required
    @marshal_with(mixture_component_recommendation_fields)
    def post(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Recommend new components to add to a mixture.
        
        Analyzes the current mixture and suggests new components that would
        improve its properties or effectiveness as a cryoprotectant.
        
        Args:
            mixture_id (str): ID of the mixture to improve
            
        Returns:
            Tuple[Dict[str, Any], int]: JSON response with component recommendations and HTTP status code
            
        Raises:
            400: If validation error occurs
            404: If mixture not found
            500: If recommendation error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            target_property (str, optional): Property to improve
            target_value (float, optional): Target value for the property
            count (int, optional): Number of recommendations to return. Default: 3
        """
        try:
            # Get request data
            data = request.get_json() or {}
            
            # Extract parameters
            target_property = data.get("target_property")
            target_value = data.get("target_value")
            count = data.get("count", 3)
            
            # Get recommendations
            result = MixtureRecommendation.recommend_components(
                mixture_id, target_property, target_value, count
            )
            
            if "error" in result:
                error_response, error_status = handle_error(
                    result["error"],
                    context="Recommending mixture components",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
                
            return result
            
        except ValidationError as err:
            error_response, error_status = handle_error(
                err,
                context="Validating recommendation parameters",
                log_level='error',
                return_response=True,
                status_code=400
            )
            abort(error_status, **error_response)
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Recommending mixture components",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


def register_resources(api):
    """
    Register mixture analysis resources with the API.
    
    Args:
        api: Flask-RESTful API instance
    
    This function registers all mixture analysis endpoints with the provided
    Flask-RESTful API instance, making them available for use.
    """
    api.add_resource(MixturePropertiesResource, '/api/v1/mixtures/<string:mixture_id>/properties')
    api.add_resource(MixtureCompatibilityResource, '/api/v1/mixtures/<string:mixture_id>/compatibility')
    api.add_resource(MixtureSynergyResource, '/api/v1/mixtures/<string:mixture_id>/synergy')
    api.add_resource(MixtureOptimizationResource, '/api/v1/mixtures/<string:mixture_id>/optimize')
    api.add_resource(MixtureStepOptimizationResource, '/api/v1/mixtures/<string:mixture_id>/optimize-step')
    api.add_resource(MixtureAnalysisResource, '/api/v1/mixtures/<string:mixture_id>/analyze')
    api.add_resource(MixtureComponentRecommendationResource, '/api/v1/mixtures/<string:mixture_id>/recommend-components')