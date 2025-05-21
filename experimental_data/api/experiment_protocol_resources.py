#!/usr/bin/env python3
"""
Experiment-Protocol Linkage Resources for CryoProtect Enhanced Experimental Data System.

This module provides API resources for linking experiments with protocols.
"""

from flask import request, g, current_app
from flask_restful import Resource, abort
from typing import Dict, List, Any, Optional, Union, Tuple
import uuid
from datetime import datetime
import json

from ..services import (
    ExperimentService,
    create_database_adapter,
)
from ..adapters.protocol_adapter import ProtocolAdapter
from ..models import (
    Experiment,
    ValidationError
)

# Create database adapter and services
db_adapter = create_database_adapter()
experiment_service = ExperimentService(db_adapter)
protocol_adapter = ProtocolAdapter()

# API response helpers
def create_success_response(data=None, message=None, status_code=200, **kwargs):
    """Create a standardized success response."""
    response = {
        "status": "success",
        "message": message,
        "data": data
    }
    
    # Add any additional fields
    for key, value in kwargs.items():
        response[key] = value
    
    return response, status_code

def create_error_response(error, status_code=400, context=None):
    """Create a standardized error response."""
    response = {
        "status": "error",
        "message": str(error),
        "context": context
    }
    
    return response, status_code

class ExperimentProtocolLinkResource(Resource):
    """Resource for linking experiments with protocols."""
    
    async def post(self):
        """
        Link an experiment with a protocol.
        
        Request body:
        - experiment_id: ID of the experiment
        - protocol_id: ID of the protocol
        
        Returns:
            JSON response with link details
        """
        try:
            # Parse request body
            data = request.get_json()
            
            # Validate required fields
            experiment_id = data.get("experiment_id")
            protocol_id = data.get("protocol_id")
            
            if not experiment_id:
                return create_error_response(
                    error="Experiment ID is required",
                    status_code=400,
                    context="Experiment-Protocol Link"
                )
            
            if not protocol_id:
                return create_error_response(
                    error="Protocol ID is required",
                    status_code=400,
                    context="Experiment-Protocol Link"
                )
            
            # Get experiment to verify it exists
            experiment = await experiment_service.get_experiment(experiment_id)
            
            if not experiment:
                return create_error_response(
                    error=f"Experiment not found: {experiment_id}",
                    status_code=404,
                    context="Experiment-Protocol Link"
                )
            
            # Create the link
            link = await protocol_adapter.link_protocol_to_experiment(protocol_id, experiment_id)
            
            # Update experiment to set protocol_id
            await experiment_service.update_experiment(
                experiment_id, 
                {"protocol_id": protocol_id}
            )
            
            return create_success_response(
                data=link,
                message="Experiment linked to protocol successfully",
                status_code=201
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Experiment-Protocol Link"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error linking experiment to protocol: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment-Protocol Link"
            )
    
    async def get(self):
        """
        Get linked experiments or protocols.
        
        Query parameters:
        - experiment_id: ID of the experiment to get linked protocols for
        - protocol_id: ID of the protocol to get linked experiments for
        
        Returns:
            JSON response with links
        """
        try:
            # Parse query parameters
            experiment_id = request.args.get("experiment_id")
            protocol_id = request.args.get("protocol_id")
            
            if experiment_id:
                # Get protocols linked to experiment
                links = await protocol_adapter.get_protocols_for_experiment(experiment_id)
                
                return create_success_response(
                    data=links,
                    message=f"Retrieved protocols linked to experiment {experiment_id}",
                    status_code=200
                )
            
            elif protocol_id:
                # Get experiments linked to protocol
                links = await protocol_adapter.get_experiments_for_protocol(protocol_id)
                
                return create_success_response(
                    data=links,
                    message=f"Retrieved experiments linked to protocol {protocol_id}",
                    status_code=200
                )
            
            else:
                return create_error_response(
                    error="Either experiment_id or protocol_id is required",
                    status_code=400,
                    context="Experiment-Protocol Link"
                )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Experiment-Protocol Link"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error getting links: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment-Protocol Link"
            )

class ExperimentProtocolLinkDetailResource(Resource):
    """Resource for managing individual experiment-protocol links."""
    
    async def get(self, link_id):
        """
        Get details of a specific experiment-protocol link.
        
        Args:
            link_id: ID of the link
        
        Returns:
            JSON response with link details
        """
        try:
            # Get link details from Convex API
            async with aiohttp.ClientSession() as session:
                headers = await get_auth_headers()
                url = f"{protocol_adapter.convex_api_url}api/protocol-experiment-links/{link_id}"
                
                async with session.get(url, headers=headers) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        return create_error_response(
                            error=f"Failed to get link: {error_text}",
                            status_code=response.status,
                            context="Experiment-Protocol Link"
                        )
                    
                    link = await response.json()
            
            return create_success_response(
                data=link,
                message=f"Retrieved link {link_id}",
                status_code=200
            )
        
        except Exception as e:
            current_app.logger.error(f"Error getting link: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment-Protocol Link"
            )
    
    async def delete(self, link_id):
        """
        Delete an experiment-protocol link.
        
        Args:
            link_id: ID of the link
        
        Returns:
            JSON response indicating success or failure
        """
        try:
            # Delete link from Convex API
            async with aiohttp.ClientSession() as session:
                headers = await get_auth_headers()
                url = f"{protocol_adapter.convex_api_url}api/protocol-experiment-links/{link_id}"
                
                async with session.delete(url, headers=headers) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        return create_error_response(
                            error=f"Failed to delete link: {error_text}",
                            status_code=response.status,
                            context="Experiment-Protocol Link"
                        )
            
            return create_success_response(
                message=f"Deleted link {link_id}",
                status_code=200
            )
        
        except Exception as e:
            current_app.logger.error(f"Error deleting link: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment-Protocol Link"
            )

class ExperimentWithProtocolResource(Resource):
    """Resource for retrieving experiments with their linked protocols."""
    
    async def get(self, experiment_id):
        """
        Get an experiment with its linked protocol.
        
        Args:
            experiment_id: ID of the experiment
        
        Returns:
            JSON response with experiment and protocol
        """
        try:
            # Get experiment
            experiment = await experiment_service.get_experiment(experiment_id)
            
            if not experiment:
                return create_error_response(
                    error=f"Experiment not found: {experiment_id}",
                    status_code=404,
                    context="Experiment with Protocol"
                )
            
            # Convert experiment to dict
            experiment_dict = experiment.to_dict()
            
            # Check if experiment has a linked protocol
            if experiment.protocol_id:
                try:
                    # Get protocol
                    protocol = await protocol_adapter.get_protocol(experiment.protocol_id)
                    
                    # Add protocol to response
                    experiment_dict["protocol"] = protocol.to_dict()
                except ValidationError:
                    # Protocol not found, but that's okay
                    experiment_dict["protocol"] = None
            else:
                experiment_dict["protocol"] = None
            
            return create_success_response(
                data=experiment_dict,
                message=f"Retrieved experiment {experiment_id} with protocol",
                status_code=200
            )
        
        except Exception as e:
            current_app.logger.error(f"Error getting experiment with protocol: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment with Protocol"
            )

class ProtocolWithExperimentsResource(Resource):
    """Resource for retrieving protocols with their linked experiments."""
    
    async def get(self, protocol_id):
        """
        Get a protocol with its linked experiments.
        
        Args:
            protocol_id: ID of the protocol
        
        Returns:
            JSON response with protocol and experiments
        """
        try:
            # Get protocol
            protocol = await protocol_adapter.get_protocol(protocol_id)
            
            # Convert to dict
            protocol_dict = protocol.to_dict()
            
            # Get linked experiments
            linked_experiments = await protocol_adapter.get_experiments_for_protocol(protocol_id)
            
            # Get full experiment details for each linked experiment
            experiments = []
            for link in linked_experiments:
                experiment_id = link.get("experiment_id")
                if experiment_id:
                    experiment = await experiment_service.get_experiment(experiment_id)
                    if experiment:
                        experiments.append(experiment.to_dict())
            
            # Add experiments to response
            protocol_dict["experiments"] = experiments
            
            return create_success_response(
                data=protocol_dict,
                message=f"Retrieved protocol {protocol_id} with experiments",
                status_code=200
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Protocol with Experiments"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error getting protocol with experiments: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol with Experiments"
            )

def register_experiment_protocol_resources(api):
    """
    Register experiment-protocol resources with the API.
    
    Args:
        api: Flask-RESTful API instance
    """
    api.add_resource(
        ExperimentProtocolLinkResource,
        '/api/v1/experiment-protocol-links',
        endpoint='experiment_protocol_links'
    )
    api.add_resource(
        ExperimentProtocolLinkDetailResource,
        '/api/v1/experiment-protocol-links/<string:link_id>',
        endpoint='experiment_protocol_link'
    )
    api.add_resource(
        ExperimentWithProtocolResource,
        '/api/v1/experiments/<string:experiment_id>/with-protocol',
        endpoint='experiment_with_protocol'
    )
    api.add_resource(
        ProtocolWithExperimentsResource,
        '/api/v1/protocols/<string:protocol_id>/with-experiments',
        endpoint='protocol_with_experiments'
    )