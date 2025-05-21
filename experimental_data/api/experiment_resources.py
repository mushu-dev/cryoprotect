#!/usr/bin/env python3
"""
Experiment Resources for CryoProtect Enhanced Experimental Data System.

This module provides API resources for experiments and their results.
"""

from flask import request, g, current_app
from flask_restful import Resource, abort
from typing import Dict, List, Any, Optional, Union, Tuple
import uuid
from datetime import datetime
import json

from ..services import (
    ExperimentService, 
    ValidationService,
    create_database_adapter
)
from ..models import (
    Experiment,
    ExperimentResult,
    ValidationError
)

# Create database adapter
db_adapter = create_database_adapter()

# Create services
experiment_service = ExperimentService(db_adapter)
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

class ExperimentListResource(Resource):
    """Resource for listing and creating experiments."""
    
    async def get(self):
        """
        Get a list of experiments.
        
        Query parameters:
        - page: Page number (default: 1)
        - per_page: Items per page (default: 20)
        - sort_by: Field to sort by (default: created_at)
        - sort_order: Sort order (asc or desc, default: desc)
        - Any other parameter is used as a filter (e.g., experiment_type_id=123)
        
        Returns:
            JSON response with experiments
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
            
            # Get experiments
            experiments, total_count = await experiment_service.list_experiments(
                filters=filters,
                page=page,
                page_size=per_page,
                sort_by=sort_by,
                sort_order=sort_order
            )
            
            # Calculate pagination values
            total_pages = (total_count + per_page - 1) // per_page if per_page > 0 else 0
            
            # Convert experiments to dictionaries
            experiment_dicts = [exp.to_dict() for exp in experiments]
            
            return create_success_response(
                data=experiment_dicts,
                message=f"Retrieved {len(experiments)} experiments",
                status_code=200,
                pagination={
                    "page": page,
                    "per_page": per_page,
                    "total_items": total_count,
                    "total_pages": total_pages
                }
            )
        
        except Exception as e:
            current_app.logger.error(f"Error retrieving experiments: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment list"
            )
    
    async def post(self):
        """
        Create a new experiment.
        
        Request body:
        - JSON object with experiment data
        
        Returns:
            JSON response with created experiment
        """
        try:
            # Parse request body
            data = request.get_json()
            
            # Add current user ID if available
            if hasattr(g, 'user') and g.user and 'id' in g.user:
                data['created_by'] = g.user['id']
            
            # Create experiment
            experiment = await experiment_service.create_experiment(data)
            
            # Validate experiment
            validation_result = await validation_service.validate_experiment(experiment)
            
            return create_success_response(
                data=experiment.to_dict(),
                message="Experiment created successfully",
                status_code=201,
                validation=validation_result
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Experiment validation"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error creating experiment: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment creation"
            )

class ExperimentResource(Resource):
    """Resource for individual experiments."""
    
    async def get(self, experiment_id):
        """
        Get an experiment by ID.
        
        Args:
            experiment_id: ID of the experiment
        
        Returns:
            JSON response with experiment
        """
        try:
            # Get experiment
            experiment = await experiment_service.get_experiment(experiment_id)
            
            if not experiment:
                return create_error_response(
                    error=f"Experiment not found: {experiment_id}",
                    status_code=404,
                    context="Experiment retrieval"
                )
            
            return create_success_response(
                data=experiment.to_dict(),
                message=f"Retrieved experiment {experiment_id}",
                status_code=200
            )
        
        except Exception as e:
            current_app.logger.error(f"Error retrieving experiment: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment retrieval"
            )
    
    async def put(self, experiment_id):
        """
        Update an experiment.
        
        Args:
            experiment_id: ID of the experiment to update
        
        Request body:
        - JSON object with updated experiment data
        
        Returns:
            JSON response with updated experiment
        """
        try:
            # Parse request body
            data = request.get_json()
            
            # Update experiment
            experiment = await experiment_service.update_experiment(experiment_id, data)
            
            if not experiment:
                return create_error_response(
                    error=f"Experiment not found: {experiment_id}",
                    status_code=404,
                    context="Experiment update"
                )
            
            # Validate experiment
            validation_result = await validation_service.validate_experiment(experiment)
            
            return create_success_response(
                data=experiment.to_dict(),
                message=f"Updated experiment {experiment_id}",
                status_code=200,
                validation=validation_result
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Experiment validation"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error updating experiment: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment update"
            )
    
    async def delete(self, experiment_id):
        """
        Delete an experiment.
        
        Args:
            experiment_id: ID of the experiment to delete
        
        Returns:
            JSON response indicating success or failure
        """
        try:
            # Delete experiment
            result = await experiment_service.delete_experiment(experiment_id)
            
            if not result:
                return create_error_response(
                    error=f"Experiment not found: {experiment_id}",
                    status_code=404,
                    context="Experiment deletion"
                )
            
            return create_success_response(
                message=f"Deleted experiment {experiment_id}",
                status_code=200
            )
        
        except Exception as e:
            current_app.logger.error(f"Error deleting experiment: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment deletion"
            )

class ExperimentResultListResource(Resource):
    """Resource for listing and creating experiment results."""
    
    async def get(self, experiment_id):
        """
        Get results for an experiment.
        
        Args:
            experiment_id: ID of the experiment
        
        Returns:
            JSON response with experiment results
        """
        try:
            # Get experiment results
            results = await experiment_service.get_experiment_results(experiment_id)
            
            # Convert results to dictionaries
            result_dicts = [result.to_dict() for result in results]
            
            return create_success_response(
                data=result_dicts,
                message=f"Retrieved {len(results)} results for experiment {experiment_id}",
                status_code=200
            )
        
        except Exception as e:
            current_app.logger.error(f"Error retrieving experiment results: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment results"
            )
    
    async def post(self, experiment_id):
        """
        Create a new experiment result.
        
        Args:
            experiment_id: ID of the experiment
        
        Request body:
        - JSON object with result data
        
        Returns:
            JSON response with created result
        """
        try:
            # Parse request body
            data = request.get_json()
            
            # Ensure experiment_id is set
            data['experiment_id'] = experiment_id
            
            # Add current user ID if available
            if hasattr(g, 'user') and g.user and 'id' in g.user:
                data['created_by'] = g.user['id']
            
            # Create result
            result = await experiment_service.create_experiment_result(data)
            
            # Validate result
            validation_result = await validation_service.validate_experiment_result(result)
            
            return create_success_response(
                data=result.to_dict(),
                message="Experiment result created successfully",
                status_code=201,
                validation=validation_result
            )
        
        except ValidationError as e:
            return create_error_response(
                error=str(e),
                status_code=400,
                context="Experiment result validation"
            )
        
        except Exception as e:
            current_app.logger.error(f"Error creating experiment result: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment result creation"
            )

class ExperimentAnalysisResource(Resource):
    """Resource for analyzing experiment results."""
    
    async def get(self, experiment_id, analysis_type='basic'):
        """
        Get analysis of experiment results.
        
        Args:
            experiment_id: ID of the experiment
            analysis_type: Type of analysis to perform (basic or detailed)
        
        Returns:
            JSON response with analysis results
        """
        try:
            # Check if analysis type is supported
            if analysis_type not in ['basic', 'detailed']:
                return create_error_response(
                    error=f"Unsupported analysis type: {analysis_type}",
                    status_code=400,
                    context="Experiment analysis"
                )
            
            # Perform analysis
            analysis = await experiment_service.analyze_experiment_results(
                experiment_id=experiment_id,
                analysis_type=analysis_type
            )
            
            return create_success_response(
                data=analysis,
                message=f"Performed {analysis_type} analysis for experiment {experiment_id}",
                status_code=200
            )
        
        except Exception as e:
            current_app.logger.error(f"Error analyzing experiment: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Experiment analysis"
            )

def register_experiment_resources(api):
    """
    Register experiment resources with the API.
    
    Args:
        api: Flask-RESTful API instance
    """
    # Register experiment resources
    api.add_resource(
        ExperimentListResource,
        '/api/v1/experiments',
        endpoint='experiments'
    )
    api.add_resource(
        ExperimentResource,
        '/api/v1/experiments/<string:experiment_id>',
        endpoint='experiment'
    )
    api.add_resource(
        ExperimentResultListResource,
        '/api/v1/experiments/<string:experiment_id>/results',
        endpoint='experiment_results'
    )
    api.add_resource(
        ExperimentAnalysisResource,
        '/api/v1/experiments/<string:experiment_id>/analysis',
        '/api/v1/experiments/<string:experiment_id>/analysis/<string:analysis_type>',
        endpoint='experiment_analysis'
    )