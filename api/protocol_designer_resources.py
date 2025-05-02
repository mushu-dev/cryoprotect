"""
CryoProtect Analyzer API - Protocol Designer Resources

This module provides API endpoints for the concentration gradient protocol designer.
"""

import logging
from flask import request
from flask_restful import Resource, marshal_with, abort
from marshmallow import ValidationError

from api.utils import token_required, handle_error
from api.protocol_designer import ProtocolDesigner
from api.models import Mixture, protocol_fields

# Set up logging
logger = logging.getLogger(__name__)

class ProtocolDesignerResource(Resource):
    """Resource for designing concentration gradient protocols."""
    
    @token_required
    @marshal_with(protocol_fields)
    def post(self, mixture_id):
        """
        Design a concentration gradient protocol for a mixture.
        
        Args:
            mixture_id: ID of the mixture
            
        Returns:
            JSON response with protocol details
        """
        try:
            # Check if mixture exists
            mixture = Mixture.get(mixture_id)
            if not mixture:
                error_response, error_status = handle_error(
                    f"Mixture with ID {mixture_id} not found",
                    context="ProtocolDesignerResource.post",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Get request data
            data = request.get_json() or {}
            
            # Extract parameters
            target_concentration = data.get("target_concentration", 10.0)
            sample_type = data.get("sample_type", "cell_line")
            starting_temperature = data.get("starting_temperature", 4.0)
            target_temperature = data.get("target_temperature")
            step_count = data.get("step_count")
            custom_sensitivity = data.get("custom_sensitivity")
            
            # Validate parameters
            if target_concentration <= 0:
                error_response, error_status = handle_error(
                    "Target concentration must be greater than 0",
                    context="ProtocolDesignerResource.post",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Design protocol
            protocol = ProtocolDesigner.design_concentration_gradient(
                mixture_id=mixture_id,
                target_concentration=target_concentration,
                sample_type=sample_type,
                starting_temperature=starting_temperature,
                target_temperature=target_temperature,
                step_count=step_count,
                custom_sensitivity=custom_sensitivity
            )
            
            if "error" in protocol:
                error_response, error_status = handle_error(
                    protocol["error"],
                    context="ProtocolDesignerResource.post",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
                
            return protocol, 200
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="ProtocolDesignerResource.post",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)
class SaveProtocolResource(Resource):
    """Resource for saving protocols to the database."""
    
    @token_required
    @marshal_with(protocol_fields)
    def post(self, mixture_id):
        """
        Save a protocol to the database.
        
        Args:
            mixture_id: ID of the mixture
            
        Returns:
            JSON response with saved protocol ID
        """
        try:
            # Check if mixture exists
            mixture = Mixture.get(mixture_id)
            if not mixture:
                error_response, error_status = handle_error(
                    f"Mixture with ID {mixture_id} not found",
                    context="SaveProtocolResource.post",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Get request data
            data = request.get_json() or {}
            
            # Extract protocol metadata
            name = data.get("name")
            if not name:
                error_response, error_status = handle_error(
                    "Protocol name is required",
                    context="SaveProtocolResource.post",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
                
            description = data.get("description", "")
            
            # Extract protocol parameters
            target_concentration = data.get("target_concentration", 10.0)
            sample_type = data.get("sample_type", "cell_line")
            starting_temperature = data.get("starting_temperature", 4.0)
            target_temperature = data.get("target_temperature")
            step_count = data.get("step_count")
            custom_sensitivity = data.get("custom_sensitivity")
            
            # Generate protocol if not provided
            if "protocol" in data:
                protocol_data = data["protocol"]
            else:
                # Generate protocol using production implementation
                protocol_data = ProtocolDesigner.design_concentration_gradient(
                    mixture_id=mixture_id,
                    target_concentration=target_concentration,
                    sample_type=sample_type,
                    starting_temperature=starting_temperature,
                    target_temperature=target_temperature,
                    step_count=step_count,
                    custom_sensitivity=custom_sensitivity
                )
                
                if "error" in protocol_data:
                    error_response, error_status = handle_error(
                        protocol_data["error"],
                        context="SaveProtocolResource.post",
                        log_level='error',
                        return_response=True,
                        status_code=400
                    )
                    abort(error_status, **error_response)
            
            # Save protocol to database
            saved_protocol = ProtocolDesigner.save_protocol(
                mixture_id=mixture_id,
                name=name,
                description=description,
                protocol_data=protocol_data
            )
            
            if "error" in saved_protocol:
                error_response, error_status = handle_error(
                    saved_protocol["error"],
                    context="SaveProtocolResource.post",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
                
            return saved_protocol, 201

        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="SaveProtocolResource.post",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class ProtocolResource(Resource):
    """Resource for retrieving saved protocols."""
    
    @token_required
    @marshal_with(protocol_fields)
    def get(self, protocol_id):
        """
        Get a saved protocol.
        
        Args:
            protocol_id: ID of the protocol
            
        Returns:
            JSON response with protocol details
        """
        try:
            # Get protocol
            protocol = ProtocolDesigner.get_saved_protocol(protocol_id)
            
            if "error" in protocol:
                error_response, error_status = handle_error(
                    protocol["error"],
                    context="ProtocolResource.get",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
                
            return protocol, 200
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="ProtocolResource.get",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)

class MixtureProtocolsResource(Resource):
    """Resource for listing protocols for a mixture."""
    
    @token_required
    @marshal_with({
        'mixture_id': protocol_fields['mixture_id'],
        'mixture_name': protocol_fields['mixture_name'],
        'protocols': protocol_fields
    })
    def get(self, mixture_id):
        """
        List all protocols for a mixture.
        
        Args:
            mixture_id: ID of the mixture
            
        Returns:
            JSON response with list of protocols
        """
        try:
            # Check if mixture exists
            mixture = Mixture.get(mixture_id)
            if not mixture:
                error_response, error_status = handle_error(
                    f"Mixture with ID {mixture_id} not found",
                    context="MixtureProtocolsResource.get",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Get protocols
            protocols = ProtocolDesigner.list_protocols_for_mixture(mixture_id)
            
            # Format response
            return {
                "mixture_id": mixture_id,
                "mixture_name": mixture.get("name", ""),
                "protocols": protocols
            }, 200
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="MixtureProtocolsResource.get",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)

class SampleSensitivityProfilesResource(Resource):
    """Resource for retrieving sample sensitivity profiles."""
    
    @token_required
    @marshal_with({
        'profiles': protocol_fields['custom_sensitivity']
    })
    def get(self):
        """
        Get available sample sensitivity profiles.
        
        Returns:
            JSON response with sensitivity profiles
        """
        try:
            # Get profiles
            profiles = ProtocolDesigner.get_sample_sensitivity_profiles()
            
            return {
                "profiles": profiles
            }, 200
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="SampleSensitivityProfilesResource.get",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)

class ProtocolComparisonResource(Resource):
    """Resource for comparing multiple protocols."""
    
    @token_required
    @marshal_with({
        'protocols': protocol_fields,
        'comparison': protocol_fields['steps']
    })
    def post(self):
        """
        Compare multiple protocols.
        
        Returns:
            JSON response with comparison results
        """
        try:
            # Get request data
            data = request.get_json() or {}
            
            # Extract parameters
            protocol_ids = data.get("protocol_ids", [])
            
            if not protocol_ids:
                error_response, error_status = handle_error(
                    "No protocol IDs provided",
                    context="ProtocolComparisonResource.post",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Compare protocols
            comparison = ProtocolDesigner.compare_protocols(protocol_ids)
            
            if "error" in comparison:
                error_response, error_status = handle_error(
                    comparison["error"],
                    context="ProtocolComparisonResource.post",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
                
            return comparison, 200
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="ProtocolComparisonResource.post",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)

def register_resources(api):
    """
    Register protocol designer resources with the API.
    
    Args:
        api: Flask-RESTful API instance
    """
    api.add_resource(ProtocolDesignerResource, '/api/v1/protocols/design/<string:mixture_id>')
    api.add_resource(SaveProtocolResource, '/api/v1/protocols/save/<string:mixture_id>')
    api.add_resource(ProtocolResource, '/api/v1/protocols/<string:protocol_id>')
    api.add_resource(MixtureProtocolsResource, '/api/v1/mixtures/<string:mixture_id>/protocols')
    api.add_resource(SampleSensitivityProfilesResource, '/api/v1/protocols/sensitivity-profiles')
    api.add_resource(ProtocolComparisonResource, '/api/v1/protocols/compare')