"""
CryoProtect Analyzer API - Batch Operations

This module provides API endpoints for batch operations on molecules, mixtures, and experiments.
Supported operations include property calculation, mixture optimization, predictive scoring, and export.
Endpoints accept lists of IDs and return structured results for each item, with error handling for partial failures.

Scientific basis and logic are documented in code comments.
"""

import logging
from typing import Dict, List, Any, Tuple, Optional, Union

from flask import request, current_app
from flask_restful import Resource
from api.utils import token_required, _handle_json_serialization

# Import core logic modules
from api.mixture_analysis import MixtureProperty, MixtureOptimization
from api.predictive_models import ModelManager
from api.scoring import calculate_molecule_score, calculate_mixture_score
from api.export_resources import get_data_for_export, generate_csv, generate_excel, generate_json, generate_pdf

# Set up logging
logger = logging.getLogger(__name__)

class BatchOperationService:
    """
    Service class for batch operations on molecules, mixtures, and experiments.
    
    This class contains the business logic for batch operations, separated from
    the Flask-specific code to improve testability.
    """
    
    @staticmethod
    def validate_request(operation: str, entity_type: str, ids: List[str]) -> Optional[Dict[str, Any]]:
        """
        Validate the batch operation request parameters.
        
        Args:
            operation: The operation to perform
            entity_type: The type of entity to operate on
            ids: List of entity IDs
            
        Returns:
            Error response dict if validation fails, None if validation passes
        """
        if not operation or not entity_type or not isinstance(ids, list):
            return {
                "status": "ERROR",
                "results": [],
                "errors": ["Missing or invalid parameters: operation, entity_type, ids"]
            }
        return None
    
    @staticmethod
    def process_property_calculation(entity_type: str, item_id: str) -> Dict[str, Any]:
        """
        Process property calculation for a single entity.
        
        Args:
            entity_type: Type of entity (molecule, mixture)
            item_id: ID of the entity
            
        Returns:
            Result of the property calculation
            
        Raises:
            ValueError: If entity_type is not supported
        """
        if entity_type == "molecule":
            return calculate_molecule_score(item_id)
        elif entity_type == "mixture":
            return calculate_mixture_score(item_id)
        else:
            raise ValueError(f"Unsupported entity_type for property_calculation: {entity_type}")
    
    @staticmethod
    def process_mixture_optimization(entity_type: str, item_id: str) -> Dict[str, Any]:
        """
        Process mixture optimization for a single entity.
        
        Args:
            entity_type: Type of entity (must be mixture)
            item_id: ID of the entity
            
        Returns:
            Result of the mixture optimization
            
        Raises:
            ValueError: If entity_type is not mixture
        """
        if entity_type == "mixture":
            return MixtureOptimization.optimize_composition(item_id)
        else:
            raise ValueError(f"Unsupported entity_type for mixture_optimization: {entity_type}")
    
    @staticmethod
    def process_predictive_scoring(entity_type: str, item_id: str,
                                  property_name: str = 'Cryoprotection Score',
                                  algorithm: str = 'random_forest') -> Dict[str, Any]:
        """
        Process predictive scoring for a single entity.
        
        Args:
            entity_type: Type of entity (molecule, mixture, experiment)
            item_id: ID of the entity
            property_name: Name of the property to predict
            algorithm: Algorithm to use for prediction
            
        Returns:
            Result of the predictive scoring
        """
        model_manager = ModelManager()
        return model_manager.predict(property_name, item_id, algorithm)
    
    @staticmethod
    def process_export(entity_type: str, item_id: str,
                      export_format: str = 'csv') -> Dict[str, Any]:
        """
        Process export for a single entity.
        
        Args:
            entity_type: Type of entity (molecule, mixture, experiment)
            item_id: ID of the entity
            export_format: Format to export (csv, json, excel, pdf)
            
        Returns:
            Result of the export operation
            
        Raises:
            ValueError: If export_format is not supported
        """
        # Ensure entity_type has 's' suffix for the API
        entity_type_plural = entity_type + 's'
        
        # Get data for export
        data = get_data_for_export(entity_type_plural, item_id)
        
        # Generate export in the specified format
        if export_format == 'csv':
            export_result = generate_csv(data, f"{entity_type}_{item_id or 'all'}")
            return export_result.getvalue()
        elif export_format == 'json':
            export_result = generate_json(data, f"{entity_type}_{item_id or 'all'}")
            return export_result.getvalue()
        elif export_format == 'excel':
            export_result = generate_excel(data, f"{entity_type}_{item_id or 'all'}")
            return export_result.getvalue()
        elif export_format == 'pdf':
            export_result = generate_pdf(data, f"{entity_type}_{item_id or 'all'}", entity_type_plural)
            return export_result.getvalue()
        else:
            raise ValueError(f"Unsupported export format: {export_format}")
    
    @staticmethod
    def process_batch_operation(operation: str, entity_type: str, ids: List[str],
                               additional_params: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Process a batch operation on multiple entities.
        
        Args:
            operation: The operation to perform
            entity_type: The type of entity to operate on
            ids: List of entity IDs
            additional_params: Additional parameters for the operation
            
        Returns:
            Dictionary with operation results and status
        """
        if additional_params is None:
            additional_params = {}
            
        results = []
        errors = []
        
        # Validate request parameters
        validation_error = BatchOperationService.validate_request(operation, entity_type, ids)
        if validation_error:
            return validation_error
        
        # Process each item in the batch
        for item_id in ids:
            try:
                if operation == "property_calculation":
                    result = BatchOperationService.process_property_calculation(entity_type, item_id)
                    results.append({"id": item_id, "result": result})
                    
                elif operation == "mixture_optimization":
                    result = BatchOperationService.process_mixture_optimization(entity_type, item_id)
                    results.append({"id": item_id, "result": result})
                    
                elif operation == "predictive_scoring":
                    property_name = additional_params.get('property_name', 'Cryoprotection Score')
                    algorithm = additional_params.get('algorithm', 'random_forest')
                    result = BatchOperationService.process_predictive_scoring(
                        entity_type, item_id, property_name, algorithm
                    )
                    results.append({"id": item_id, "result": result})
                    
                elif operation == "export":
                    export_format = additional_params.get('format', 'csv')
                    export_result = BatchOperationService.process_export(
                        entity_type, item_id, export_format
                    )
                    results.append({"id": item_id, "export": export_result})
                    
                else:
                    raise ValueError(f"Unsupported operation: {operation}")
                    
            except Exception as e:
                logger.error(f"Batch operation error for {item_id}: {str(e)}", exc_info=True)
                errors.append({"id": item_id, "error": str(e)})
        
        # Determine overall status
        status = "SUCCESS"
        if errors and results:
            status = "COMPLETED_WITH_WARNINGS"
        elif errors and not results:
            status = "ERROR"
        
        return {
            "status": status,
            "results": results,
            "errors": errors
        }


class BatchOperationResource(Resource):
    """
    Resource for batch operations on molecules, mixtures, and experiments.
    
    This class handles HTTP requests and responses, delegating business logic
    to the BatchOperationService class.
    """

    @token_required
    def post(self):
        """
        Perform a batch operation.

        Request JSON:
        {
            "operation": "property_calculation" | "mixture_optimization" | "predictive_scoring" | "export",
            "entity_type": "molecule" | "mixture" | "experiment",
            "ids": [ ... ],  # List of IDs
            "format": "csv" | "json" | "excel" | "pdf",  # Optional, for export operation
            "property_name": "...",  # Optional, for predictive_scoring operation
            "algorithm": "..."  # Optional, for predictive_scoring operation
        }

        Returns:
            {
                "status": "SUCCESS" | "COMPLETED_WITH_WARNINGS" | "ERROR",
                "operation_id": "...",  # Unique ID for the batch operation
                "results": [ ... ],  # List of results per item
                "errors": [ ... ]    # List of errors per item (if any)
            }
        """
        try:
            # Parse request data
            data = request.get_json(force=True)
            operation = data.get("operation")
            entity_type = data.get("entity_type")
            ids = data.get("ids", [])
            
            # Extract additional parameters
            additional_params = {
                'format': data.get('format', 'csv'),
                'property_name': data.get('property_name', 'Cryoprotection Score'),
                'algorithm': data.get('algorithm', 'random_forest')
            }
            
            # Process the batch operation using the service
            result = BatchOperationService.process_batch_operation(
                operation, entity_type, ids, additional_params
            )
            
            # Generate a unique operation ID
            import uuid
            operation_id = str(uuid.uuid4())
            
            # Store the operation result for later retrieval
            from api.utils import get_supabase_client, get_user_id
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Store operation metadata
            operation_data = {
                "id": operation_id,
                "operation": operation,
                "entity_type": entity_type,
                "status": result["status"],
                "created_by": user_id,
                "created_at": datetime.now().isoformat(),
                "completed_at": datetime.now().isoformat(),
                "result_count": len(result.get("results", [])),
                "error_count": len(result.get("errors", [])),
                "is_cancelled": False
            }
            
            # Store in database
            try:
                supabase.table("batch_operations").insert(operation_data).execute()
            except Exception as e:
                logger.warning(f"Failed to store batch operation metadata: {str(e)}")
            
            # Add operation ID to result
            result["operation_id"] = operation_id
            
            # Return the serialized result
            return _handle_json_serialization(result)
            
        except Exception as e:
            logger.error(f"Unexpected error in batch operation: {str(e)}", exc_info=True)
            return _handle_json_serialization({
                "status": "ERROR",
                "results": [],
                "errors": [f"Unexpected error: {str(e)}"]
            })

class BatchOperationStatusResource(Resource):
    """
    Resource for checking the status of a batch operation.
    
    This class handles HTTP requests and responses for retrieving
    the status and results of a previously submitted batch operation.
    """
    
    @token_required
    def get(self, operation_id):
        """
        Get the status and results of a batch operation.
        
        Args:
            operation_id: ID of the batch operation
            
        Returns:
            JSON response with operation status and results
        """
        try:
            # Get operation from database
            from api.utils import get_supabase_client, get_user_id
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Query the operation
            response = supabase.table("batch_operations").select("*").eq("id", operation_id).execute()
            
            if not response.data:
                return _handle_json_serialization({
                    "status": "ERROR",
                    "errors": [f"Batch operation with ID {operation_id} not found"]
                }), 404
            
            operation = response.data[0]
            
            # Check if user has access to this operation
            if operation.get("created_by") != user_id:
                # Check if user has admin role
                from api.rbac import UserRoleManager
                user_roles = UserRoleManager.get_user_roles(user_id)
                has_admin = any(role.get("name") == "admin" for role in user_roles)
                
                if not has_admin:
                    return _handle_json_serialization({
                        "status": "ERROR",
                        "errors": ["You do not have permission to access this batch operation"]
                    }), 403
            
            # Get operation results
            results = []
            errors = []
            
            # In a real implementation, results would be stored in a separate table
            # For this example, we'll return the operation metadata
            
            return _handle_json_serialization({
                "operation_id": operation_id,
                "status": operation.get("status"),
                "operation": operation.get("operation"),
                "entity_type": operation.get("entity_type"),
                "created_at": operation.get("created_at"),
                "completed_at": operation.get("completed_at"),
                "result_count": operation.get("result_count"),
                "error_count": operation.get("error_count"),
                "is_cancelled": operation.get("is_cancelled", False)
            })
            
        except Exception as e:
            logger.error(f"Error retrieving batch operation status: {str(e)}", exc_info=True)
            return _handle_json_serialization({
                "status": "ERROR",
                "errors": [f"Error retrieving batch operation status: {str(e)}"]
            }), 500

class BatchOperationCancellationResource(Resource):
    """
    Resource for cancelling a batch operation.
    
    This class handles HTTP requests and responses for cancelling
    a previously submitted batch operation.
    """
    
    @token_required
    def delete(self, operation_id):
        """
        Cancel a batch operation.
        
        Args:
            operation_id: ID of the batch operation to cancel
            
        Returns:
            JSON response with cancellation status
        """
        try:
            # Get operation from database
            from api.utils import get_supabase_client, get_user_id
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Query the operation
            response = supabase.table("batch_operations").select("*").eq("id", operation_id).execute()
            
            if not response.data:
                return _handle_json_serialization({
                    "status": "ERROR",
                    "errors": [f"Batch operation with ID {operation_id} not found"]
                }), 404
            
            operation = response.data[0]
            
            # Check if user has access to this operation
            if operation.get("created_by") != user_id:
                # Check if user has admin role
                from api.rbac import UserRoleManager
                user_roles = UserRoleManager.get_user_roles(user_id)
                has_admin = any(role.get("name") == "admin" for role in user_roles)
                
                if not has_admin:
                    return _handle_json_serialization({
                        "status": "ERROR",
                        "errors": ["You do not have permission to cancel this batch operation"]
                    }), 403
            
            # Check if operation is already completed
            if operation.get("status") in ["SUCCESS", "ERROR", "COMPLETED_WITH_WARNINGS"]:
                return _handle_json_serialization({
                    "status": "ERROR",
                    "errors": ["Cannot cancel a completed batch operation"]
                }), 400
            
            # Mark operation as cancelled
            update_response = supabase.table("batch_operations").update({
                "is_cancelled": True,
                "status": "CANCELLED"
            }).eq("id", operation_id).execute()
            
            if update_response.error:
                return _handle_json_serialization({
                    "status": "ERROR",
                    "errors": [f"Failed to cancel batch operation: {update_response.error.message}"]
                }), 500
            
            return _handle_json_serialization({
                "status": "SUCCESS",
                "message": f"Batch operation {operation_id} has been cancelled"
            })
            
        except Exception as e:
            logger.error(f"Error cancelling batch operation: {str(e)}", exc_info=True)
            return _handle_json_serialization({
                "status": "ERROR",
                "errors": [f"Error cancelling batch operation: {str(e)}"]
            }), 500

# Scientific/logic notes:
# - Property calculation uses existing mixture/molecule property logic.
# - Mixture optimization uses the same optimization logic as single-item endpoints, but in batch.
# - Predictive scoring leverages the predictive models for each entity type.
# - Export returns protocol files/links as in the export API.
# - Partial failures are reported in the 'errors' field, with successful results in 'results'.
# - Batch operations can be tracked and cancelled using their operation_id.

def register_resources(api):
    """Register batch resources with the API."""
    api.add_resource(BatchOperationResource, '/api/v1/batch')
    api.add_resource(BatchOperationStatusResource, '/api/v1/batch/<string:operation_id>')
    api.add_resource(BatchOperationCancellationResource, '/api/v1/batch/<string:operation_id>')