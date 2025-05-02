"""
CryoProtect Analyzer API - Scoring Resources

This module contains API resources for the cryoprotection effectiveness scoring system.
It provides endpoints for:
- Calculating cryoprotection scores for molecules and mixtures
- Storing scores in the database
- Batch scoring of multiple entities

The scoring system evaluates how effective a molecule or mixture is likely to be
as a cryoprotectant based on various physicochemical properties and structural features.
All endpoints follow RESTful design principles and return standardized response formats.
"""

from flask import request, g, current_app
from flask_restful import Resource, marshal_with, abort, fields
from marshmallow import Schema, fields as ma_fields, ValidationError, validate
import logging
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

from api.utils import token_required, get_supabase_client, get_user_id, handle_error
from api.models import Molecule, Mixture
from api.scoring import (
    calculate_molecule_score, calculate_mixture_score,
    store_molecule_score, store_mixture_score
)

# Set up logging
logger = logging.getLogger(__name__)

# Define response fields
score_component_fields = {
    'hydrogen_bonding': fields.Integer,
    'logp': fields.Integer,
    'molecular_size': fields.Integer,
    'tpsa': fields.Integer,
    'functional_groups': fields.Integer,
    'permeability': fields.Integer
}

molecule_score_fields = {
    'overall_score': fields.Integer,
    'component_scores': fields.Nested(score_component_fields),
    'properties': fields.Raw
}

component_score_fields = {
    'molecule_id': fields.String,
    'name': fields.String,
    'concentration': fields.Float,
    'concentration_unit': fields.String,
    'score': fields.Integer
}

mixture_score_fields = {
    'mixture_id': fields.String,
    'name': fields.String,
    'overall_score': fields.Integer,
    'component_scores': fields.List(fields.Nested(component_score_fields))
}

# Define request schemas
class MoleculeScoreInputSchema(Schema):
    """Schema for molecule score input."""
    molecule_data = ma_fields.String(required=True)
    input_format = ma_fields.String(default='smiles')
    store_result = ma_fields.Boolean(default=False)

class MoleculeIdScoreSchema(Schema):
    """Schema for scoring a molecule by ID."""
    store_result = ma_fields.Boolean(default=True)

class MixtureScoreSchema(Schema):
    """Schema for mixture score input."""
    store_result = ma_fields.Boolean(default=True)

class BatchScoringSchema(Schema):
    """Schema for batch scoring input."""
    entity_type = ma_fields.String(required=True, validate=validate.OneOf(['molecule', 'mixture']))
    ids = ma_fields.List(ma_fields.String(), required=True)
    store_results = ma_fields.Boolean(default=True)


class MoleculeScoreResource(MethodResource):
    """Resource for calculating cryoprotection scores for molecules.
    
    This endpoint calculates a cryoprotection effectiveness score for a molecule
    based on its physicochemical properties and structural features. The score
    indicates how likely the molecule is to be effective as a cryoprotectant.
    """
    
    @doc(description='Calculate cryoprotection score for a molecule',
         tags=['Scoring'],
         request_body={
             'molecule_data': {'description': 'Molecule data (SMILES, InChI, etc.)', 'type': 'string', 'required': True},
             'input_format': {'description': 'Format of the input data', 'type': 'string', 'default': 'smiles', 'enum': ['smiles', 'inchi', 'mol']},
             'store_result': {'description': 'Whether to store the result in the database', 'type': 'boolean', 'default': False}
         })
    @apispec_marshal_with(molecule_score_fields, code=200, description='Cryoprotection score results')
    @marshal_with(molecule_score_fields)
    def post(self) -> Tuple[Dict[str, Any], int]:
        """Calculate cryoprotection score for a molecule.
        
        Evaluates a molecule's potential effectiveness as a cryoprotectant based on
        its physicochemical properties and structural features. Returns an overall
        score and component scores for different aspects of cryoprotection.
        
        Returns:
            Tuple[Dict[str, Any], int]: Score results and HTTP status code
            
        Raises:
            400: If request data is invalid or molecule cannot be processed
            500: If server error occurs
            
        Request Body:
            molecule_data (str): Molecule data in the specified format (required)
            input_format (str): Format of the input data (default: 'smiles')
                Supported formats: 'smiles', 'inchi', 'mol'
            store_result (bool): Whether to store the result in the database (default: False)
                Note: Storing requires authentication
                
        Response:
            Score information including:
            - Overall score (0-100)
            - Component scores for different aspects of cryoprotection
            - Properties used in the calculation
        """
        try:
            # Validate request data
            schema = MoleculeScoreInputSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating molecule score input",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Calculate score
            score_data = calculate_molecule_score(
                data["molecule_data"],
                data.get("input_format", "smiles")
            )
            
            if "error" in score_data:
                error_response, error_status = handle_error(
                    score_data["error"],
                    context="Calculating molecule score",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Store result if requested and user is authenticated
            if data.get("store_result", False):
                try:
                    # Create molecule in database if it doesn't exist
                    supabase = get_supabase_client()
                    user_id = get_user_id()
                    
                    if user_id:
                        # This would require additional implementation to store
                        # an arbitrary molecule that's not already in the database
                        pass
                except Exception as e:
                    handle_error(
                        e,
                        context="Storing molecule score",
                        log_level='warning'
                    )
            
            return score_data, 200
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Calculating molecule score",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MoleculeIdScoreResource(MethodResource):
    """Resource for calculating cryoprotection scores for molecules by ID.
    
    This endpoint calculates a cryoprotection effectiveness score for a molecule
    that already exists in the database. It retrieves the molecule by ID, calculates
    the score, and optionally stores the result in the database.
    """
    
    @doc(description='Calculate cryoprotection score for a molecule by ID',
         tags=['Scoring'],
         params={'molecule_id': {'description': 'ID of the molecule to score', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}],
         request_body={
             'store_result': {'description': 'Whether to store the result in the database', 'type': 'boolean', 'default': True}
         })
    @apispec_marshal_with(molecule_score_fields, code=200, description='Cryoprotection score results')
    @token_required
    @marshal_with(molecule_score_fields)
    def post(self, molecule_id: str) -> Tuple[Dict[str, Any], int]:
        """Calculate cryoprotection score for a molecule by ID.
        
        Retrieves a molecule from the database by ID, calculates its cryoprotection
        effectiveness score, and optionally stores the result in the database.
        
        Args:
            molecule_id (str): ID of the molecule to score
            
        Returns:
            Tuple[Dict[str, Any], int]: Score results and HTTP status code
            
        Raises:
            400: If validation error occurs or molecule cannot be processed
            404: If molecule not found
            500: If server error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            store_result (bool): Whether to store the result in the database (default: True)
                
        Response:
            Score information including:
            - Overall score (0-100)
            - Component scores for different aspects of cryoprotection
            - Properties used in the calculation
        """
        try:
            # Validate request data
            schema = MoleculeIdScoreSchema()
            try:
                data = schema.load(request.json or {})
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating molecule score input",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Get molecule from database
            molecule = Molecule.get(molecule_id)
            if not molecule:
                error_response, error_status = handle_error(
                    f"Molecule with ID {molecule_id} not found",
                    context="Calculating molecule score by ID",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Get SMILES
            smiles = molecule.get("smiles")
            if not smiles:
                error_response, error_status = handle_error(
                    "Molecule does not have SMILES data",
                    context=f"Calculating score for molecule {molecule_id}",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Calculate score
            score_data = calculate_molecule_score(smiles)
            
            if "error" in score_data:
                error_response, error_status = handle_error(
                    score_data["error"],
                    context=f"Calculating score for molecule {molecule_id}",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Store result if requested
            if data.get("store_result", True):
                success = store_molecule_score(molecule_id, score_data)
                if not success:
                    handle_error(
                        f"Could not store score for molecule {molecule_id}",
                        context="Storing molecule score",
                        log_level='warning'
                    )
            
            return score_data, 200
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context=f"Calculating score for molecule {molecule_id}",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MixtureScoreResource(MethodResource):
    """Resource for calculating cryoprotection scores for mixtures.
    
    This endpoint calculates a cryoprotection effectiveness score for a mixture
    based on its components and their properties. It evaluates how effective the
    mixture is likely to be as a cryoprotectant and provides component-level scores.
    """
    
    @doc(description='Calculate cryoprotection score for a mixture',
         tags=['Scoring'],
         params={'mixture_id': {'description': 'ID of the mixture to score', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}],
         request_body={
             'store_result': {'description': 'Whether to store the result in the database', 'type': 'boolean', 'default': True}
         })
    @apispec_marshal_with(mixture_score_fields, code=200, description='Mixture cryoprotection score results')
    @token_required
    @marshal_with(mixture_score_fields)
    def post(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """Calculate cryoprotection score for a mixture.
        
        Evaluates a mixture's potential effectiveness as a cryoprotectant based on
        its components and their properties. Returns an overall score and individual
        scores for each component.
        
        Args:
            mixture_id (str): ID of the mixture to score
            
        Returns:
            Tuple[Dict[str, Any], int]: Score results and HTTP status code
            
        Raises:
            400: If validation error occurs or mixture cannot be processed
            404: If mixture not found
            500: If server error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            store_result (bool): Whether to store the result in the database (default: True)
                
        Response:
            Score information including:
            - Overall score (0-100)
            - Individual scores for each component
            - Mixture name and ID
        """
        try:
            # Validate request data
            schema = MixtureScoreSchema()
            try:
                data = schema.load(request.json or {})
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating mixture score input",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Calculate score
            score_data = calculate_mixture_score(mixture_id)
            
            if "error" in score_data:
                error_response, error_status = handle_error(
                    score_data["error"],
                    context=f"Calculating score for mixture {mixture_id}",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Store result if requested
            if data.get("store_result", True):
                success = store_mixture_score(mixture_id, score_data)
                if not success:
                    handle_error(
                        f"Could not store score for mixture {mixture_id}",
                        context="Storing mixture score",
                        log_level='warning'
                    )
            
            return score_data, 200
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context=f"Calculating score for mixture {mixture_id}",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class BatchScoringResource(MethodResource):
    """Resource for batch scoring of molecules or mixtures.
    
    This endpoint calculates cryoprotection effectiveness scores for multiple
    molecules or mixtures in a single request. It is useful for efficiently
    processing and comparing scores for multiple entities.
    """
    
    @doc(description='Calculate cryoprotection scores for multiple molecules or mixtures',
         tags=['Scoring'],
         security=[{'Bearer': []}],
         request_body={
             'entity_type': {'description': 'Type of entities to score', 'type': 'string', 'enum': ['molecule', 'mixture'], 'required': True},
             'ids': {'description': 'List of entity IDs to score', 'type': 'array', 'items': {'type': 'string'}, 'required': True},
             'store_results': {'description': 'Whether to store the results in the database', 'type': 'boolean', 'default': True}
         })
    @token_required
    def post(self) -> Tuple[Dict[str, Any], int]:
        """
        Calculate cryoprotection scores for multiple molecules or mixtures.
        
        Processes multiple entities in a single request, calculating cryoprotection
        effectiveness scores for each one. Returns results for all entities along
        with any errors that occurred during processing.
        
        Returns:
            Tuple[Dict[str, Any], int]: Batch scoring results and HTTP status code
            
        Raises:
            400: If validation error occurs or no valid IDs provided
            500: If server error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            entity_type (str): Type of entities to score (required)
                Supported types: 'molecule', 'mixture'
            ids (List[str]): List of entity IDs to score (required)
            store_results (bool): Whether to store the results in the database (default: True)
                
        Response:
            Batch scoring results including:
            - Overall status ('SUCCESS', 'COMPLETED_WITH_WARNINGS', 'ERROR')
            - Individual results for each entity
            - Errors for failed entities
        """
        try:
            # Validate request data
            schema = BatchScoringSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating batch scoring input",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Extract parameters
            entity_type = data["entity_type"]
            ids = data["ids"]
            store_results = data.get("store_results", True)
            
            if not ids:
                error_response, error_status = handle_error(
                    "No IDs provided",
                    context="Batch scoring",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Process each entity
            results = []
            errors = []
            
            for entity_id in ids:
                try:
                    # Calculate score based on entity type
                    if entity_type == "molecule":
                        # Get molecule
                        molecule = Molecule.get(entity_id)
                        if not molecule:
                            errors.append({
                                "id": entity_id,
                                "error": f"Molecule with ID {entity_id} not found"
                            })
                            continue
                        
                        # Get SMILES
                        smiles = molecule.get("smiles")
                        if not smiles:
                            errors.append({
                                "id": entity_id,
                                "error": "Molecule does not have SMILES data"
                            })
                            continue
                        
                        # Calculate score
                        score_data = calculate_molecule_score(smiles)
                        
                        if "error" in score_data:
                            errors.append({
                                "id": entity_id,
                                "error": score_data["error"]
                            })
                            continue
                        
                        # Store result if requested
                        if store_results:
                            success = store_molecule_score(entity_id, score_data)
                            if not success:
                                logger.warning(f"Could not store score for molecule {entity_id}")
                        
                        # Add to results
                        results.append({
                            "id": entity_id,
                            "name": molecule.get("name", ""),
                            "type": "molecule",
                            "score": score_data
                        })
                        
                    elif entity_type == "mixture":
                        # Calculate score
                        score_data = calculate_mixture_score(entity_id)
                        
                        if "error" in score_data:
                            errors.append({
                                "id": entity_id,
                                "error": score_data["error"]
                            })
                            continue
                        
                        # Store result if requested
                        if store_results:
                            success = store_mixture_score(entity_id, score_data)
                            if not success:
                                logger.warning(f"Could not store score for mixture {entity_id}")
                        
                        # Add to results
                        results.append({
                            "id": entity_id,
                            "name": score_data.get("mixture_name", ""),
                            "type": "mixture",
                            "score": score_data
                        })
                        
                    else:
                        errors.append({
                            "id": entity_id,
                            "error": f"Unsupported entity type: {entity_type}"
                        })
                        
                except Exception as e:
                    errors.append({
                        "id": entity_id,
                        "error": str(e)
                    })
            
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
            }, 200
            
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Batch scoring",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


def register_resources(api):
    """
    Register scoring resources with the API.
    
    Args:
        api: Flask-RESTful API instance
        
    This function registers all scoring endpoints with the provided
    Flask-RESTful API instance, making them available for use.
    """
    api.add_resource(MoleculeScoreResource, '/api/v1/scoring/molecules')
    api.add_resource(MoleculeIdScoreResource, '/api/v1/molecules/<string:molecule_id>/score')
    api.add_resource(MixtureScoreResource, '/api/v1/mixtures/<string:mixture_id>/score')
    api.add_resource(BatchScoringResource, '/api/v1/scoring/batch')