#!/usr/bin/env python3
"""
Protocol Resources for CryoProtect Enhanced Experimental Data System.

This module provides API resources for protocols and protocol steps.
"""

from flask import request, g, current_app
from flask_restful import Resource, abort
from typing import Dict, List, Any, Optional, Union, Tuple
import uuid
from datetime import datetime
import json

from ..services import (
    ProtocolService, 
    ValidationService,
    create_database_adapter
)
from ..models import (
    Protocol,
    ProtocolStep,
    ValidationError
)

# Create database adapter
db_adapter = create_database_adapter()

# Create services
protocol_service = ProtocolService(db_adapter)
validation_service = ValidationService(db_adapter)

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

class ProtocolListResource(Resource):
    """Resource for listing and creating protocols."""
    
    async def get(self):
        """
        Get a list of protocols.
        
        Query parameters:
        - page: Page number (default: 1)
        - per_page: Items per page (default: 20)
        - sort_by: Field to sort by (default: created_at)
        - sort_order: Sort order (asc or desc, default: desc)
        - Any other parameter is used as a filter (e.g., template_type=standard)
        
        Returns:
            JSON response with protocols
        """
        try:
            # Parse query parameters
            page = request.args.get('page', 1, type=int)
            per_page = request.args.get('per_page', 20, type=int)
            sort_by = request.args.get('sort_by', 'created_at')
            sort_order = request.args.get('sort_order', 'desc')
            
            # Build filters from remaining query parameters
            filters = {}
            for key, value in request.args.items():
                if key not in ['page', 'per_page', 'sort_by', 'sort_order']:
                    filters[key] = value
            
            # Get protocols
            protocols, total_count = await protocol_service.list_protocols(
                filters=filters,
                page=page,
                page_size=per_page,
                sort_by=sort_by,
                sort_order=sort_order
            )
            
            # Calculate pagination values
            total_pages = (total_count + per_page - 1) // per_page if per_page > 0 else 0
            
            # Convert protocols to dictionaries
            protocol_dicts = [proto.to_dict() for proto in protocols]
            
            return create_success_response(
                data=protocol_dicts,
                message=f"Retrieved {len(protocols)} protocols",
                status_code=200,
                pagination={
                    "page": page,
                    "per_page": per_page,
                    "total_items": total_count,
                    "total_pages": total_pages
                }
            )
        
        except Exception as e:
            current_app.logger.error(f"Error retrieving protocols: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol list"
            )
    
    async def post(self):
        """
        Create a new protocol.
        
        Request body:
        - JSON object with protocol data
        
        Returns:
            JSON response with created protocol
        """
        try:
            # Parse request body
            data = request.get_json()
            
            # Add current user ID if available
            if hasattr(g, 'user') and g.user and 'id' in g.user:
                data['created_by'] = g.user['id']
            
            # Create protocol
            protocol = await protocol_service.create_protocol(data)
            
            # Validate protocol
            validation_result = await validation_service.validate_protocol(protocol)
            
            return create_success_response(
                data=protocol.to_dict(),
                message="Protocol created successfully",
                status_code=201,
                validation=validation_result
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Protocol validation"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error creating protocol: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol creation"
            )

class ProtocolResource(Resource):
    """Resource for individual protocols."""
    
    async def get(self, protocol_id):
        """
        Get a protocol by ID.
        
        Args:
            protocol_id: ID of the protocol
        
        Returns:
            JSON response with protocol
        """
        try:
            # Get protocol
            protocol = await protocol_service.get_protocol(protocol_id)
            
            if not protocol:
                return create_error_response(
                    error=f"Protocol not found: {protocol_id}",
                    status_code=404,
                    context="Protocol retrieval"
                )
            
            return create_success_response(
                data=protocol.to_dict(),
                message=f"Retrieved protocol {protocol_id}",
                status_code=200
            )
        
        except Exception as e:
            current_app.logger.error(f"Error retrieving protocol: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol retrieval"
            )
    
    async def put(self, protocol_id):
        """
        Update a protocol.
        
        Args:
            protocol_id: ID of the protocol to update
        
        Request body:
        - JSON object with updated protocol data
        
        Returns:
            JSON response with updated protocol
        """
        try:
            # Parse request body
            data = request.get_json()
            
            # Update protocol
            protocol = await protocol_service.update_protocol(protocol_id, data)
            
            if not protocol:
                return create_error_response(
                    error=f"Protocol not found: {protocol_id}",
                    status_code=404,
                    context="Protocol update"
                )
            
            # Validate protocol
            validation_result = await validation_service.validate_protocol(protocol)
            
            return create_success_response(
                data=protocol.to_dict(),
                message=f"Updated protocol {protocol_id}",
                status_code=200,
                validation=validation_result
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Protocol validation"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error updating protocol: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol update"
            )
    
    async def delete(self, protocol_id):
        """
        Delete a protocol.
        
        Args:
            protocol_id: ID of the protocol to delete
        
        Returns:
            JSON response indicating success or failure
        """
        try:
            # Delete protocol
            result = await protocol_service.delete_protocol(protocol_id)
            
            if not result:
                return create_error_response(
                    error=f"Protocol not found: {protocol_id}",
                    status_code=404,
                    context="Protocol deletion"
                )
            
            return create_success_response(
                message=f"Deleted protocol {protocol_id}",
                status_code=200
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Protocol deletion"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error deleting protocol: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol deletion"
            )

class ProtocolVersionResource(Resource):
    """Resource for protocol versioning operations."""
    
    async def post(self, protocol_id):
        """
        Create a new version of a protocol.
        
        Args:
            protocol_id: ID of the protocol to version
        
        Request body (optional):
        - JSON object with changes to apply to the new version
        
        Returns:
            JSON response with the new protocol version
        """
        try:
            # Parse request body
            data = request.get_json() if request.is_json else None
            
            # Create new version
            new_protocol = await protocol_service.create_protocol_version(protocol_id, data)
            
            # Validate new protocol
            validation_result = await validation_service.validate_protocol(new_protocol)
            
            return create_success_response(
                data=new_protocol.to_dict(),
                message=f"Created new version of protocol {protocol_id}",
                status_code=201,
                validation=validation_result
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Protocol versioning"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error creating protocol version: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol versioning"
            )
    
    async def get(self, protocol_id):
        """
        Get version history for a protocol.
        
        Args:
            protocol_id: ID of the protocol
        
        Returns:
            JSON response with version history
        """
        try:
            # Get version history
            history = await protocol_service.get_protocol_version_history(protocol_id)
            
            if not history:
                return create_error_response(
                    error=f"Protocol not found: {protocol_id}",
                    status_code=404,
                    context="Protocol version history"
                )
            
            return create_success_response(
                data=history,
                message=f"Retrieved version history for protocol {protocol_id}",
                status_code=200
            )
        
        except Exception as e:
            current_app.logger.error(f"Error retrieving protocol version history: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol version history"
            )

class ProtocolCompareResource(Resource):
    """Resource for comparing protocols."""
    
    async def get(self, protocol_id1, protocol_id2):
        """
        Compare two protocols.
        
        Args:
            protocol_id1: ID of the first protocol
            protocol_id2: ID of the second protocol
        
        Returns:
            JSON response with comparison results
        """
        try:
            # Compare protocols
            comparison = await protocol_service.compare_protocols(protocol_id1, protocol_id2)
            
            return create_success_response(
                data=comparison,
                message=f"Compared protocols {protocol_id1} and {protocol_id2}",
                status_code=200
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Protocol comparison"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error comparing protocols: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol comparison"
            )

class ProtocolExportResource(Resource):
    """Resource for exporting protocols."""
    
    async def get(self, protocol_id, format='json'):
        """
        Export a protocol.
        
        Args:
            protocol_id: ID of the protocol
            format: Export format (json, yaml, or human)
        
        Returns:
            Exported protocol
        """
        try:
            # Export protocol
            export_data = await protocol_service.export_protocol(protocol_id, format)
            
            # Set response content type based on format
            if format == 'json':
                content_type = 'application/json'
            elif format == 'yaml':
                content_type = 'application/yaml'
            else:
                content_type = 'text/plain'
            
            return export_data, 200, {'Content-Type': content_type}
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Protocol export"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error exporting protocol: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol export"
            )

class ProtocolImportResource(Resource):
    """Resource for importing protocols."""
    
    async def post(self):
        """
        Import a protocol.
        
        Request body:
        - JSON object with protocol data
        
        Returns:
            JSON response with imported protocol
        """
        try:
            # Parse request body
            data = request.get_json()
            
            # Import protocol
            protocol = await protocol_service.import_protocol(data)
            
            # Validate protocol
            validation_result = await validation_service.validate_protocol(protocol)
            
            return create_success_response(
                data=protocol.to_dict(),
                message="Protocol imported successfully",
                status_code=201,
                validation=validation_result
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Protocol import"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error importing protocol: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol import"
            )

class ProtocolValidateResource(Resource):
    """Resource for validating protocols."""
    
    async def post(self, protocol_id):
        """
        Validate a protocol.
        
        Args:
            protocol_id: ID of the protocol to validate
        
        Returns:
            JSON response with validation results
        """
        try:
            # Validate protocol
            validation_result = await protocol_service.validate_protocol(protocol_id)
            
            return create_success_response(
                data=validation_result,
                message=f"Validated protocol {protocol_id}",
                status_code=200
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Protocol validation"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error validating protocol: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Protocol validation"
            )

def register_protocol_resources(api):
    """
    Register protocol resources with the API.
    
    Args:
        api: Flask-RESTful API instance
    """
    # Register protocol resources
    api.add_resource(
        ProtocolListResource,
        '/api/v1/protocols',
        endpoint='protocols'
    )
    api.add_resource(
        ProtocolResource,
        '/api/v1/protocols/<string:protocol_id>',
        endpoint='protocol'
    )
    api.add_resource(
        ProtocolVersionResource,
        '/api/v1/protocols/<string:protocol_id>/versions',
        endpoint='protocol_versions'
    )
    api.add_resource(
        ProtocolCompareResource,
        '/api/v1/protocols/compare/<string:protocol_id1>/<string:protocol_id2>',
        endpoint='protocol_compare'
    )
    api.add_resource(
        ProtocolExportResource,
        '/api/v1/protocols/<string:protocol_id>/export',
        '/api/v1/protocols/<string:protocol_id>/export/<string:format>',
        endpoint='protocol_export'
    )
    api.add_resource(
        ProtocolImportResource,
        '/api/v1/protocols/import',
        endpoint='protocol_import'
    )
    api.add_resource(
        ProtocolValidateResource,
        '/api/v1/protocols/<string:protocol_id>/validate',
        endpoint='protocol_validate'
    )