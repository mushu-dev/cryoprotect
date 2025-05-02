"""
CryoProtect Analyzer API - Resources

This module contains the API resource classes that handle HTTP requests.
It implements standardized patterns for error handling, authentication,
response formatting, and request validation.

All endpoints follow RESTful design principles and return standardized
response formats. Authentication is handled via JWT tokens for protected
endpoints.
"""

from flask import request, g, current_app
from flask_restful import Resource, marshal_with, reqparse, abort, fields
from marshmallow import ValidationError
from functools import wraps
from datetime import datetime
from typing import Dict, List, Optional, Union, Any, Tuple

from api.schemas import PropertyComparisonRequestSchema
from api.pagination_utils import parse_pagination_params, paginate_query_response, get_total_count

# Import API documentation utilities
try:
    from flask_apispec import use_kwargs, marshal_with as apispec_marshal_with, doc
    from flask_apispec.views import MethodResource
    from api.docs import auth_required
except ImportError:
    # Create dummy decorators if flask-apispec is not installed
    def use_kwargs(*args, **kwargs): return lambda f: f
    def apispec_marshal_with(*args, **kwargs): return lambda f: f
    def doc(*args, **kwargs): return lambda f: f
    def auth_required(*args, **kwargs): return lambda f: f
    MethodResource = Resource

from api.models import (
    molecule_fields, mixture_fields, prediction_fields, experiment_fields, comparison_fields,
    MoleculeSchema, MixtureSchema, PredictionSchema, ExperimentSchema, ComparisonQuerySchema
)
from api.utils import (
    get_supabase_client, token_required, handle_supabase_error, get_user_id,
    handle_error
)
from api.comparisons import compare_entities

import json
import uuid


class MoleculeListResource(MethodResource):
    """Resource for listing and creating molecules.
    
    Provides endpoints for retrieving a paginated list of molecules and
    importing molecules from PubChem.
    """
    
    # Schema for molecule import
    molecule_parser = reqparse.RequestParser()
    molecule_parser.add_argument('cid', type=int, required=True,
                               help='PubChem Compound ID is required')
    
    # Schema for request validation
    molecule_import_schema = MoleculeSchema(only=["cid"])
    
    @doc(description='Get a paginated list of molecules with their properties and optional filtering.',
         tags=['Molecules'],
         params={
             'limit': {
                 'description': 'Maximum number of records to return (max 500)',
                 'type': 'integer',
                 'default': 100,
                 'minimum': 1,
                 'maximum': 500,
                 'example': 50
             },
             'offset': {
                 'description': 'Number of records to skip for pagination',
                 'type': 'integer',
                 'default': 0,
                 'minimum': 0,
                 'example': 100
             }
         },
         responses={
             200: {
                 'description': 'List of molecules with their properties',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'array',
                             'items': {
                                 'type': 'object',
                                 'properties': {
                                     'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the molecule'},
                                     'cid': {'type': 'integer', 'description': 'PubChem Compound ID'},
                                     'name': {'type': 'string', 'description': 'Name of the molecule'},
                                     'molecular_formula': {'type': 'string', 'description': 'Chemical formula'},
                                     'smiles': {'type': 'string', 'description': 'SMILES notation representing molecular structure'},
                                     'pubchem_link': {'type': 'string', 'description': 'Link to PubChem page'},
                                     'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                     'updated_at': {'type': 'string', 'format': 'date-time', 'description': 'Last update timestamp'},
                                     'properties': {'type': 'object', 'description': 'Molecular properties'}
                                 }
                             }
                         },
                         'example': [
                             {
                                 'id': '123e4567-e89b-12d3-a456-426614174000',
                                 'cid': 2244,
                                 'name': 'Aspirin',
                                 'molecular_formula': 'C9H8O4',
                                 'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                                 'pubchem_link': 'https://pubchem.ncbi.nlm.nih.gov/compound/2244',
                                 'created_at': '2025-01-15T10:30:00Z',
                                 'updated_at': '2025-01-15T10:30:00Z',
                                 'properties': {
                                     'molecular_weight': 180.16,
                                     'logp': 1.43,
                                     'hydrogen_bond_donors': 1,
                                     'hydrogen_bond_acceptors': 4,
                                     'rotatable_bonds': 3
                                 }
                             }
                         ]
                     }
                 }
             },
             400: {
                 'description': 'Invalid pagination parameters',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Invalid pagination parameters',
                             'message': 'Limit must be a positive integer not exceeding 500',
                             'status': 400
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error fetching molecules list',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(molecule_fields, code=200, description='List of molecules')
    @marshal_with(molecule_fields)
    def get(self) -> Tuple[List[Dict[str, Any]], int]:
        """
        Get a paginated list of molecules with their properties.

        This endpoint retrieves a list of molecules from the database with pagination support.
        Each molecule includes its basic information (name, formula, SMILES) and associated
        molecular properties.

        Query Parameters:
            limit (int, optional): Maximum number of records to return. Default is 100, maximum is 500.
                Example: ?limit=50
            offset (int, optional): Number of records to skip for pagination. Default is 0.
                Example: ?offset=100

        Returns:
            Tuple[List[Dict[str, Any]], int]: List of molecule objects with pagination and HTTP status code
            
        Raises:
            400: If pagination parameters are invalid
            500: If database error occurs
            
        Response Schema:
            The response is a list of objects following the MoleculeSchema structure with the following fields:
            - id (UUID): Unique identifier for the molecule
            - cid (int): PubChem Compound ID
            - name (str): Name of the molecule
            - molecular_formula (str): Chemical formula
            - smiles (str): SMILES notation representing molecular structure
            - pubchem_link (str): Link to PubChem page
            - created_at (datetime): Creation timestamp
            - updated_at (datetime): Last update timestamp
            - properties (object): Molecular properties including:
                - molecular_weight (float): Molecular weight in g/mol
                - logp (float): Calculated LogP value
                - hydrogen_bond_donors (int): Number of hydrogen bond donors
                - hydrogen_bond_acceptors (int): Number of hydrogen bond acceptors
                - rotatable_bonds (int): Number of rotatable bonds
                - tpsa (float): Topological polar surface area
                - functional_groups (object): Detected functional groups
        
        Example Request:
            GET /api/molecules?limit=10&offset=0
            
        Example Response:
            [
                {
                    "id": "123e4567-e89b-12d3-a456-426614174000",
                    "cid": 2244,
                    "name": "Aspirin",
                    "molecular_formula": "C9H8O4",
                    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                    "pubchem_link": "https://pubchem.ncbi.nlm.nih.gov/compound/2244",
                    "created_at": "2025-01-15T10:30:00Z",
                    "updated_at": "2025-01-15T10:30:00Z",
                    "properties": {
                        "molecular_weight": 180.16,
                        "logp": 1.43,
                        "hydrogen_bond_donors": 1,
                        "hydrogen_bond_acceptors": 4,
                        "rotatable_bonds": 3,
                        "tpsa": 63.6,
                        "functional_groups": {
                            "carboxylic_acid": 1,
                            "ester": 1,
                            "aromatic": 1
                        }
                    }
                },
                ...
            ]
        
        Notes for Developers:
            - The endpoint uses the 'molecules_with_properties' view which joins the molecules table
              with their properties for efficient retrieval.
            - Properties are dynamically loaded from the molecular_properties table.
            - For large datasets, consider using the pagination parameters to improve performance.
            - The maximum limit is capped at 500 to prevent excessive server load.
        
        Performance Optimization:
            - Added pagination (limit, offset) to avoid loading all molecules at once.
            - This reduces memory usage and improves response time for large datasets.
            - The endpoint uses a database view that pre-joins molecules with their properties.
        """
        try:
            supabase = get_supabase_client()
            
            # Parse pagination parameters using the utility function
            try:
                limit, offset = parse_pagination_params(default_limit=100, max_limit=500)
            except ValueError as e:
                error_response, error_status = handle_error(
                    str(e),
                    context="Invalid pagination parameters",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Get total count for pagination metadata
            total_count = get_total_count(supabase, "molecules_with_properties")
            
            # Execute the query with pagination
            response = supabase.table("molecules_with_properties").select("*").range(offset, offset + limit - 1).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context="Fetching molecules list",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            # Format the response with pagination metadata
            resource_endpoint = "/api/molecules"
            query_params = {k: v for k, v in request.args.items() if k not in ['limit', 'offset', 'page', 'per_page']}
            
            paginated_response = paginate_query_response(
                response.data,
                total_count,
                limit,
                offset,
                resource_endpoint,
                query_params
            )
            
            return paginated_response, 200
        except Exception as e:
            if isinstance(e, ValueError):
                status_code = 400
                context = "Invalid pagination parameters"
            else:
                status_code = 500
                context = "Fetching molecules list"
                
            error_response, error_status = handle_error(
                e,
                context=context,
                log_level='error',
                return_response=True,
                status_code=status_code
            )
            abort(error_status, **error_response)
    
    @doc(description='Import a molecule from PubChem.',
         tags=['Molecules'],
         params={'cid': {'description': 'PubChem Compound ID', 'type': 'integer', 'required': True}},
         security=[{'bearerAuth': []}],
         requestBody={
             'description': 'PubChem Compound ID to import',
             'required': True,
             'content': {
                 'application/json': {
                     'schema': {
                         'type': 'object',
                         'properties': {
                             'cid': {
                                 'type': 'integer',
                                 'description': 'PubChem Compound ID',
                                 'example': 2244
                             }
                         },
                         'required': ['cid']
                     },
                     'example': {
                         'cid': 2244
                     }
                 }
             }
         },
         responses={
             201: {
                 'description': 'Molecule successfully imported',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the molecule'},
                                 'cid': {'type': 'integer', 'description': 'PubChem Compound ID'},
                                 'name': {'type': 'string', 'description': 'Name of the molecule'},
                                 'molecular_formula': {'type': 'string', 'description': 'Chemical formula'},
                                 'smiles': {'type': 'string', 'description': 'SMILES notation representing molecular structure'},
                                 'pubchem_link': {'type': 'string', 'description': 'Link to PubChem page'},
                                 'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                 'updated_at': {'type': 'string', 'format': 'date-time', 'description': 'Last update timestamp'},
                                 'properties': {'type': 'object', 'description': 'Molecular properties'}
                             }
                         },
                         'example': {
                             'id': '123e4567-e89b-12d3-a456-426614174000',
                             'cid': 2244,
                             'name': 'Aspirin',
                             'molecular_formula': 'C9H8O4',
                             'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                             'pubchem_link': 'https://pubchem.ncbi.nlm.nih.gov/compound/2244',
                             'created_at': '2025-01-15T10:30:00Z',
                             'updated_at': '2025-01-15T10:30:00Z',
                             'properties': {
                                 'molecular_weight': 180.16,
                                 'logp': 1.43,
                                 'hydrogen_bond_donors': 1,
                                 'hydrogen_bond_acceptors': 4,
                                 'rotatable_bonds': 3,
                                 'tpsa': 63.6,
                                 'functional_groups': {
                                     'carboxylic_acid': 1,
                                     'ester': 1,
                                     'aromatic': 1
                                 }
                             }
                         }
                     }
                 }
             },
             400: {
                 'description': 'Invalid request parameters',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Invalid request parameters',
                             'message': 'PubChem Compound ID (cid) is required',
                             'status': 400
                         }
                     }
                 }
             },
             404: {
                 'description': 'Molecule not found in PubChem',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Not found',
                             'message': 'Molecule with CID 999999999 not found in PubChem',
                             'status': 404
                         }
                     }
                 }
             },
             401: {
                 'description': 'Authentication required',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Unauthorized',
                             'message': 'Authentication token is missing or invalid',
                             'status': 401
                         }
                     }
                 }
             },
             403: {
                 'description': 'Insufficient permissions',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Forbidden',
                             'message': 'You do not have permission to import molecules',
                             'status': 403
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error importing molecule from PubChem',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @use_kwargs(molecule_import_schema, location='json')
    @apispec_marshal_with(molecule_fields, code=201, description='Imported molecule')
    @auth_required('bearerAuth')
    @token_required
    @marshal_with(molecule_fields)
    def post(self, user=None, **kwargs) -> Tuple[Dict[str, Any], int]:
        """Import a molecule from PubChem.
        
        Args:
            user: User object from token_required decorator
            **kwargs: Validated fields from the molecule_import_schema
            
        Returns:
            Tuple[Dict[str, Any], int]: Imported molecule object and HTTP status code
            
        Raises:
            400: If request data is invalid
            401: If authentication token is missing or invalid
            403: If user doesn't have permission to import molecules
            404: If molecule not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Authorization:
            Any authenticated user can import molecules. No special roles required.
            
        Request Body:
            cid (int): PubChem Compound ID of the molecule to import
            
        Performance Optimization:
            - Uses database function to efficiently import molecule data
            - Retrieves complete molecule data with properties in a single query
        """
        try:
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Get CID from validated schema
            cid = kwargs.get('cid')
            if not cid:
                error_response, error_status = handle_error(
                    "PubChem Compound ID (cid) is required",
                    context="Importing molecule from PubChem",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Call the database function to import the molecule
            response = supabase.rpc(
                "import_molecule_from_pubchem",
                {"p_cid": cid, "p_user_id": user_id}
            ).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Importing molecule from PubChem (CID: {args['cid']})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            molecule_id = response.data
            
            # Get the imported molecule with properties
            response = supabase.table("molecules_with_properties").select("*").eq("id", molecule_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Retrieving imported molecule (ID: {molecule_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            if not response.data:
                error_response, error_status = handle_error(
                    f"Molecule with ID {molecule_id} not found",
                    context="Retrieving imported molecule",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            return response.data[0], 201
        except Exception as e:
            if isinstance(e, ValueError):
                status_code = 400
                context = "Invalid request parameters"
            else:
                status_code = 500
                context = f"Importing molecule from PubChem (CID: {args['cid']})"
                
            error_response, error_status = handle_error(
                e,
                context=context,
                log_level='error',
                return_response=True,
                status_code=status_code
            )
            abort(error_status, **error_response)


class MoleculeResource(MethodResource):
    """Resource for retrieving, updating, and deleting a specific molecule.
    
    Provides endpoints for retrieving, updating, and deleting a molecule by ID.
    """
    
    # Schema for molecule update
    molecule_update_parser = reqparse.RequestParser()
    molecule_update_parser.add_argument('name', type=str, help='Name of the molecule')
    molecule_update_parser.add_argument('formula', type=str, help='Chemical formula')
    molecule_update_parser.add_argument('smiles', type=str, help='SMILES notation')
    
    # Schema for request validation
    molecule_update_schema = MoleculeSchema(partial=True, only=["name", "molecular_formula", "smiles"])
    
    @doc(description='Get a molecule with its properties.',
         tags=['Molecules'],
         params={'molecule_id': {'description': 'ID of the molecule to retrieve', 'type': 'string', 'required': True}})
    @apispec_marshal_with(molecule_fields, code=200, description='Molecule object')
    @marshal_with(molecule_fields)
    def get(self, molecule_id: str) -> Tuple[Dict[str, Any], int]:
        """Get a molecule with its properties.
        
        Args:
            molecule_id (str): ID of the molecule to retrieve
            
        Returns:
            Tuple[Dict[str, Any], int]: Molecule object if found and HTTP status code
            
        Raises:
            404: If molecule not found
            500: If database error occurs
            
        Response Schema:
            The response follows the MoleculeSchema structure with the following fields:
            - id: UUID of the molecule
            - cid: PubChem Compound ID (integer)
            - name: Name of the molecule (string)
            - molecular_formula: Chemical formula (string)
            - smiles: SMILES notation (string)
            - pubchem_link: Link to PubChem page (string)
            - created_at: Creation timestamp (datetime)
            - updated_at: Last update timestamp (datetime)
            - properties: Molecular properties (object)
            
        Performance Optimization:
            - Retrieves molecule with all properties in a single query
            - Uses standardized error handling for consistent responses
        """
        try:
            supabase = get_supabase_client()
            response = supabase.table("molecules_with_properties").select("*").eq("id", molecule_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Fetching molecule (ID: {molecule_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            if not response.data:
                error_response, error_status = handle_error(
                    f"Molecule with ID {molecule_id} not found",
                    context="Fetching molecule",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            return response.data[0], 200
        except Exception as e:
            from werkzeug.exceptions import HTTPException
            if isinstance(e, HTTPException):
                raise
                
            error_response, error_status = handle_error(
                e,
                context=f"Fetching molecule (ID: {molecule_id})",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)
    
    @doc(description='Update a molecule\'s information.',
         tags=['Molecules'],
         params={'molecule_id': {'description': 'ID of the molecule to update', 'type': 'string', 'required': True}},
         security=[{'bearerAuth': []}],
         requestBody={
             'description': 'Updated molecule information',
             'required': True,
             'content': {
                 'application/json': {
                     'schema': {
                         'type': 'object',
                         'properties': {
                             'name': {
                                 'type': 'string',
                                 'description': 'Name of the molecule',
                                 'example': 'Acetylsalicylic acid'
                             },
                             'molecular_formula': {
                                 'type': 'string',
                                 'description': 'Chemical formula of the molecule',
                                 'example': 'C9H8O4'
                             },
                             'smiles': {
                                 'type': 'string',
                                 'description': 'SMILES notation representing molecular structure',
                                 'example': 'CC(=O)OC1=CC=CC=C1C(=O)O'
                             }
                         }
                     },
                     'examples': {
                         'Full update': {
                             'value': {
                                 'name': 'Acetylsalicylic acid',
                                 'molecular_formula': 'C9H8O4',
                                 'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'
                             }
                         },
                         'Partial update': {
                             'value': {
                                 'name': 'Acetylsalicylic acid'
                             }
                         }
                     }
                 }
             }
         },
         responses={
             200: {
                 'description': 'Molecule successfully updated',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the molecule'},
                                 'cid': {'type': 'integer', 'description': 'PubChem Compound ID'},
                                 'name': {'type': 'string', 'description': 'Name of the molecule'},
                                 'molecular_formula': {'type': 'string', 'description': 'Chemical formula'},
                                 'smiles': {'type': 'string', 'description': 'SMILES notation representing molecular structure'},
                                 'pubchem_link': {'type': 'string', 'description': 'Link to PubChem page'},
                                 'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                 'updated_at': {'type': 'string', 'format': 'date-time', 'description': 'Last update timestamp'},
                                 'properties': {'type': 'object', 'description': 'Molecular properties'}
                             }
                         },
                         'example': {
                             'id': '123e4567-e89b-12d3-a456-426614174000',
                             'cid': 2244,
                             'name': 'Acetylsalicylic acid',
                             'molecular_formula': 'C9H8O4',
                             'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                             'pubchem_link': 'https://pubchem.ncbi.nlm.nih.gov/compound/2244',
                             'created_at': '2025-01-15T10:30:00Z',
                             'updated_at': '2025-01-15T11:45:00Z',
                             'properties': {
                                 'molecular_weight': 180.16,
                                 'logp': 1.43,
                                 'hydrogen_bond_donors': 1,
                                 'hydrogen_bond_acceptors': 4,
                                 'rotatable_bonds': 3,
                                 'tpsa': 63.6,
                                 'functional_groups': {
                                     'carboxylic_acid': 1,
                                     'ester': 1,
                                     'aromatic': 1
                                 }
                             }
                         }
                     }
                 }
             },
             400: {
                 'description': 'Invalid request parameters',
                 'content': {
                     'application/json': {
                         'examples': {
                             'No update data': {
                                 'value': {
                                     'error': 'Invalid request parameters',
                                     'message': 'No update data provided',
                                     'status': 400
                                 }
                             },
                             'Invalid data': {
                                 'value': {
                                     'error': 'Invalid request parameters',
                                     'message': 'Invalid SMILES notation',
                                     'status': 400
                                 }
                             }
                         }
                     }
                 }
             },
             404: {
                 'description': 'Molecule not found',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Not found',
                             'message': 'Molecule with ID 123e4567-e89b-12d3-a456-426614174000 not found',
                             'status': 404
                         }
                     }
                 }
             },
             401: {
                 'description': 'Authentication required',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Unauthorized',
                             'message': 'Authentication token is missing or invalid',
                             'status': 401
                         }
                     }
                 }
             },
             403: {
                 'description': 'Insufficient permissions',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Forbidden',
                             'message': 'You do not have permission to update this molecule',
                             'status': 403
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error updating molecule',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @use_kwargs(molecule_update_schema, location='json')
    @apispec_marshal_with(molecule_fields, code=200, description='Updated molecule')
    @auth_required('bearerAuth')
    @token_required
    @marshal_with(molecule_fields)
    def put(self, molecule_id: str, user=None, **kwargs) -> Tuple[Dict[str, Any], int]:
        """Update a molecule's information.
        
        Args:
            molecule_id (str): ID of the molecule to update
            user: User object from token_required decorator
            **kwargs: Validated fields from the molecule_update_schema
            
        Returns:
            Tuple[Dict[str, Any], int]: Updated molecule object and HTTP status code
            
        Raises:
            400: If no update data provided
            401: If authentication token is missing or invalid
            403: If user doesn't have permission to update this molecule
            404: If molecule not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Authorization:
            Any authenticated user can update molecules. No special roles required.
            In a production environment, you might want to restrict updates to:
            - The user who created the molecule
            - Users with specific roles (admin, curator)
            
        Request Body:
            name (str, optional): Name of the molecule
            molecular_formula (str, optional): Chemical formula
            smiles (str, optional): SMILES notation
            
        Performance Optimization:
            - Only updates provided fields
            - Verifies molecule exists before attempting update
            - Retrieves updated molecule with properties in a single query
        """
        try:
            # Verify molecule exists
            supabase = get_supabase_client()
            check_response = supabase.table("molecules").select("id").eq("id", molecule_id).execute()
            
            error_message, status_code = handle_supabase_error(check_response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Verifying molecule (ID: {molecule_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            if not check_response.data:
                error_response, error_status = handle_error(
                    f"Molecule with ID {molecule_id} not found",
                    context="Updating molecule",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Prepare update data from validated schema
            update_data = {}
            if 'name' in kwargs:
                update_data['name'] = kwargs['name']
            if 'molecular_formula' in kwargs:
                update_data['formula'] = kwargs['molecular_formula']
            if 'smiles' in kwargs:
                update_data['smiles'] = kwargs['smiles']
            
            if not update_data:
                error_response, error_status = handle_error(
                    "No update data provided",
                    context="Updating molecule",
                    log_level='warning',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Add updated_at timestamp
            update_data['updated_at'] = datetime.now().isoformat()
            
            # Update the molecule
            response = supabase.table("molecules").update(update_data).eq("id", molecule_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Updating molecule (ID: {molecule_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            # Get the updated molecule with properties
            response = supabase.table("molecules_with_properties").select("*").eq("id", molecule_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Retrieving updated molecule (ID: {molecule_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            current_app.logger.info(f"Molecule {molecule_id} updated successfully")
            return response.data[0], 200
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context=f"Updating molecule (ID: {molecule_id})",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)
    
    @doc(description='Delete a molecule.',
         tags=['Molecules'],
         params={'molecule_id': {'description': 'ID of the molecule to delete', 'type': 'string', 'required': True}},
         security=[{'bearerAuth': []}],
         responses={
             200: {
                 'description': 'Molecule successfully deleted',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'message': {
                                     'type': 'string',
                                     'description': 'Success message'
                                 }
                             }
                         },
                         'example': {
                             'message': 'Molecule with ID 123e4567-e89b-12d3-a456-426614174000 deleted successfully'
                         }
                     }
                 }
             },
             404: {
                 'description': 'Molecule not found',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Not found',
                             'message': 'Molecule with ID 123e4567-e89b-12d3-a456-426614174000 not found',
                             'status': 404
                         }
                     }
                 }
             },
             401: {
                 'description': 'Authentication required',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Unauthorized',
                             'message': 'Authentication token is missing or invalid',
                             'status': 401
                         }
                     }
                 }
             },
             403: {
                 'description': 'Insufficient permissions',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Forbidden',
                             'message': 'You do not have permission to delete this molecule',
                             'status': 403
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error deleting molecule',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(fields.Raw, code=200, description='Success message')
    @auth_required('bearerAuth')
    @token_required
    def delete(self, molecule_id: str, user=None) -> Tuple[Dict[str, str], int]:
        """Delete a molecule.
        
        Args:
            molecule_id (str): ID of the molecule to delete
            user: User object from token_required decorator
            
        Returns:
            Tuple[Dict[str, str], int]: Success message and HTTP status code
            
        Raises:
            401: If authentication token is missing or invalid
            403: If user doesn't have permission to delete this molecule
            404: If molecule not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Authorization:
            Any authenticated user can delete molecules. No special roles required.
            In a production environment, you might want to restrict deletion to:
            - The user who created the molecule
            - Users with specific roles (admin, curator)
            
        Security Considerations:
            - Deletion is permanent and cannot be undone
            - Consider implementing soft delete in production environments
            - Ensure proper authorization checks before deletion
            
        Performance Optimization:
            - Verifies molecule exists before attempting deletion
            - Uses standardized error handling for consistent responses
        """
        try:
            supabase = get_supabase_client()
            
            # Check if molecule exists
            check_response = supabase.table("molecules").select("id").eq("id", molecule_id).execute()
            
            error_message, status_code = handle_supabase_error(check_response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Verifying molecule (ID: {molecule_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            if not check_response.data:
                error_response, error_status = handle_error(
                    f"Molecule with ID {molecule_id} not found",
                    context="Deleting molecule",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Delete the molecule
            response = supabase.table("molecules").delete().eq("id", molecule_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Deleting molecule (ID: {molecule_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            current_app.logger.info(f"Molecule {molecule_id} deleted successfully")
            return {"message": f"Molecule with ID {molecule_id} deleted successfully"}, 200
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context=f"Deleting molecule (ID: {molecule_id})",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MixtureListResource(MethodResource):
    """Resource for listing and creating mixtures.
    
    Provides endpoints for retrieving a paginated list of mixtures and
    creating new mixtures with components.
    """
    
    @doc(description='Get a paginated list of mixtures with their components and optional filtering.',
         tags=['Mixtures'],
         params={
             'limit': {
                 'description': 'Maximum number of records to return (max 500)',
                 'type': 'integer',
                 'default': 100,
                 'minimum': 1,
                 'maximum': 500,
                 'example': 50
             },
             'offset': {
                 'description': 'Number of records to skip for pagination',
                 'type': 'integer',
                 'default': 0,
                 'minimum': 0,
                 'example': 100
             },
             'name': {
                 'description': 'Filter mixtures by name (case-insensitive partial match)',
                 'type': 'string',
                 'example': 'cryo'
             },
             'created_by': {
                 'description': 'Filter mixtures by creator user ID',
                 'type': 'string',
                 'format': 'uuid',
                 'example': '123e4567-e89b-12d3-a456-426614174000'
             },
             'sort': {
                 'description': 'Field to sort by (prefix with - for descending order)',
                 'type': 'string',
                 'enum': ['name', '-name', 'created_at', '-created_at'],
                 'default': '-created_at',
                 'example': 'name'
             }
         },
         security=[{'bearerAuth': []}],
         responses={
             200: {
                 'description': 'List of mixtures with their components',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'array',
                             'items': {
                                 'type': 'object',
                                 'properties': {
                                     'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the mixture'},
                                     'name': {'type': 'string', 'description': 'Name of the mixture'},
                                     'description': {'type': 'string', 'description': 'Description of the mixture'},
                                     'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                     'updated_at': {'type': 'string', 'format': 'date-time', 'description': 'Last update timestamp'},
                                     'created_by': {'type': 'string', 'format': 'uuid', 'description': 'ID of the user who created the mixture'},
                                     'components': {
                                         'type': 'array',
                                         'description': 'List of components in the mixture',
                                         'items': {
                                             'type': 'object',
                                             'properties': {
                                                 'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the component'},
                                                 'mixture_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the mixture this component belongs to'},
                                                 'molecule_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the molecule in this component'},
                                                 'concentration': {'type': 'number', 'format': 'float', 'description': 'Concentration value'},
                                                 'concentration_unit': {'type': 'string', 'description': 'Unit of concentration (e.g., mg/mL, %w/v)'},
                                                 'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                                 'molecule': {
                                                     'type': 'object',
                                                     'description': 'Details of the molecule in this component',
                                                     'properties': {
                                                         'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the molecule'},
                                                         'name': {'type': 'string', 'description': 'Name of the molecule'},
                                                         'molecular_formula': {'type': 'string', 'description': 'Chemical formula'},
                                                         'smiles': {'type': 'string', 'description': 'SMILES notation representing molecular structure'}
                                                     }
                                                 }
                                             }
                                         }
                                     }
                                 }
                             }
                         },
                         'example': [
                             {
                                 'id': '123e4567-e89b-12d3-a456-426614174000',
                                 'name': 'Cryoprotective Solution A',
                                 'description': 'Standard cryoprotective solution for cell preservation',
                                 'created_at': '2025-01-15T10:30:00Z',
                                 'updated_at': '2025-01-15T10:30:00Z',
                                 'created_by': '123e4567-e89b-12d3-a456-426614174001',
                                 'components': [
                                     {
                                         'id': '123e4567-e89b-12d3-a456-426614174002',
                                         'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                         'molecule_id': '123e4567-e89b-12d3-a456-426614174003',
                                         'concentration': 10.5,
                                         'concentration_unit': 'mg/mL',
                                         'created_at': '2025-01-15T10:30:00Z',
                                         'molecule': {
                                             'id': '123e4567-e89b-12d3-a456-426614174003',
                                             'name': 'Glycerol',
                                             'molecular_formula': 'C3H8O3',
                                             'smiles': 'C(C(CO)O)O'
                                         }
                                     },
                                     {
                                         'id': '123e4567-e89b-12d3-a456-426614174004',
                                         'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                         'molecule_id': '123e4567-e89b-12d3-a456-426614174005',
                                         'concentration': 5.0,
                                         'concentration_unit': '%v/v',
                                         'created_at': '2025-01-15T10:30:00Z',
                                         'molecule': {
                                             'id': '123e4567-e89b-12d3-a456-426614174005',
                                             'name': 'Dimethyl sulfoxide',
                                             'molecular_formula': 'C2H6OS',
                                             'smiles': 'CS(=O)C'
                                         }
                                     }
                                 ]
                             }
                         ]
                     }
                 }
             },
             400: {
                 'description': 'Invalid pagination parameters',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Invalid pagination parameters',
                             'message': 'Limit must be a positive integer not exceeding 500',
                             'status': 400
                         }
                     }
                 }
             },
             401: {
                 'description': 'Authentication required',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Unauthorized',
                             'message': 'Authentication token is missing or invalid',
                             'status': 401
                         }
                     }
                 }
             },
             403: {
                 'description': 'Insufficient permissions',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Forbidden',
                             'message': 'You do not have permission to access these mixtures',
                             'status': 403
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error fetching mixtures list',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(mixture_fields, code=200, description='List of mixtures')
    @marshal_with(mixture_fields)
    @auth_required('bearerAuth')
    @token_required
    def get(self) -> Tuple[List[Dict[str, Any]], int]:
        """
        Get a paginated list of mixtures with their components and optional filtering.

        This endpoint retrieves a list of mixtures from the database with pagination support.
        Each mixture includes its basic information (name, description) and associated
        components with their molecule details.

        Query Parameters:
            limit (int, optional): Maximum number of records to return. Default is 100, maximum is 500.
                Example: ?limit=50
            offset (int, optional): Number of records to skip for pagination. Default is 0.
                Example: ?offset=100
            name (str, optional): Filter mixtures by name (case-insensitive partial match).
                Example: ?name=cryo
            created_by (str, optional): Filter mixtures by creator user ID.
                Example: ?created_by=123e4567-e89b-12d3-a456-426614174000
            sort (str, optional): Field to sort by (prefix with - for descending order).
                Allowed values: name, -name, created_at, -created_at
                Default: -created_at
                Example: ?sort=name

        Returns:
            Tuple[List[Dict[str, Any]], int]: List of mixture objects with pagination and HTTP status code
            
        Raises:
            400: If pagination parameters are invalid
            401: If authentication token is missing or invalid
            403: If user doesn't have permission to access mixtures
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Authorization:
            Any authenticated user can view mixtures. Users can see all public mixtures
            and private mixtures they have created or have been granted access to.
            
        Response Schema:
            The response is a list of objects following the MixtureSchema structure with the following fields:
            - id (UUID): Unique identifier for the mixture
            - name (str): Name of the mixture
            - description (str): Description of the mixture
            - created_at (datetime): Creation timestamp
            - updated_at (datetime): Last update timestamp
            - created_by (UUID): ID of the user who created the mixture
            - components (list): List of components in the mixture, each with:
                - id (UUID): Unique identifier for the component
                - mixture_id (UUID): ID of the mixture this component belongs to
                - molecule_id (UUID): ID of the molecule in this component
                - concentration (float): Concentration value
                - concentration_unit (str): Unit of concentration (e.g., mg/mL, %w/v)
                - created_at (datetime): Creation timestamp
                - molecule (object): Details of the molecule in this component, including:
                    - id (UUID): Unique identifier for the molecule
                    - name (str): Name of the molecule
                    - molecular_formula (str): Chemical formula
                    - smiles (str): SMILES notation representing molecular structure
        
        Example Request:
            GET /api/mixtures?limit=10&offset=0&name=cryo&sort=name
            
        Example Response:
            [
                {
                    "id": "123e4567-e89b-12d3-a456-426614174000",
                    "name": "Cryoprotective Solution A",
                    "description": "Standard cryoprotective solution for cell preservation",
                    "created_at": "2025-01-15T10:30:00Z",
                    "updated_at": "2025-01-15T10:30:00Z",
                    "created_by": "123e4567-e89b-12d3-a456-426614174001",
                    "components": [
                        {
                            "id": "123e4567-e89b-12d3-a456-426614174002",
                            "mixture_id": "123e4567-e89b-12d3-a456-426614174000",
                            "molecule_id": "123e4567-e89b-12d3-a456-426614174003",
                            "concentration": 10.5,
                            "concentration_unit": "mg/mL",
                            "created_at": "2025-01-15T10:30:00Z",
                            "molecule": {
                                "id": "123e4567-e89b-12d3-a456-426614174003",
                                "name": "Glycerol",
                                "molecular_formula": "C3H8O3",
                                "smiles": "C(C(CO)O)O"
                            }
                        },
                        {
                            "id": "123e4567-e89b-12d3-a456-426614174004",
                            "mixture_id": "123e4567-e89b-12d3-a456-426614174000",
                            "molecule_id": "123e4567-e89b-12d3-a456-426614174005",
                            "concentration": 5.0,
                            "concentration_unit": "%v/v",
                            "created_at": "2025-01-15T10:30:00Z",
                            "molecule": {
                                "id": "123e4567-e89b-12d3-a456-426614174005",
                                "name": "Dimethyl sulfoxide",
                                "molecular_formula": "C2H6OS",
                                "smiles": "CS(=O)C"
                            }
                        }
                    ]
                },
                ...
            ]
        
        Notes for Developers:
            - The endpoint uses the 'mixtures_with_components' view which joins the mixtures table
              with their components and molecule details for efficient retrieval.
            - Components are dynamically loaded from the mixture_components table.
            - For large datasets, consider using the pagination parameters to improve performance.
            - The maximum limit is capped at 500 to prevent excessive server load.
            - Filtering by name uses a case-insensitive partial match (ILIKE %name% in SQL).
            - The default sort order is by creation date (newest first).
            - The response structure follows the MixtureSchema with nested components.

        Performance Optimization:
            - Added pagination (limit, offset) to avoid loading all mixtures at once.
            - This reduces memory usage and improves response time for large datasets.
            - The endpoint uses a database view that pre-joins mixtures with their components.
        """
        try:
            supabase = get_supabase_client()
            
            # Parse pagination parameters using the utility function
            try:
                limit, offset = parse_pagination_params(default_limit=100, max_limit=500)
            except ValueError as e:
                error_response, error_status = handle_error(
                    str(e),
                    context="Invalid pagination parameters",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Build the query
            query = supabase.table("mixtures_with_components").select("*")
            
            # Apply filters if provided
            name_filter = request.args.get("name")
            if name_filter:
                # Case-insensitive partial match using ILIKE
                query = query.ilike("name", f"%{name_filter}%")
            
            created_by_filter = request.args.get("created_by")
            if created_by_filter:
                query = query.eq("created_by", created_by_filter)
            
            # Apply sorting
            sort_param = request.args.get("sort", "-created_at")
            allowed_sort_fields = ["name", "-name", "created_at", "-created_at"]
            
            if sort_param not in allowed_sort_fields:
                sort_param = "-created_at"  # Default to newest first if invalid
            
            # Handle descending order (prefixed with -)
            if sort_param.startswith("-"):
                sort_field = sort_param[1:]  # Remove the - prefix
                sort_order = "desc"
            else:
                sort_field = sort_param
                sort_order = "asc"
            
            query = query.order(sort_field, sort_order)
            
            # Apply pagination
            query = query.range(offset, offset + limit - 1)
            
            # Execute the query
            # Get total count for pagination metadata
            # Build a count query with the same filters
            count_query = supabase.table("mixtures_with_components").select("*", count="exact")
            
            # Apply the same filters to the count query
            if name_filter:
                count_query = count_query.ilike("name", f"%{name_filter}%")
            
            if created_by_filter:
                count_query = count_query.eq("created_by", created_by_filter)
                
            count_response = count_query.execute()
            total_count = count_response.count if hasattr(count_response, 'count') else 0
            
            # Execute the main query
            response = query.execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context="Fetching mixtures list",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            # Format the response with pagination metadata
            resource_endpoint = "/api/mixtures"
            query_params = {
                k: v for k, v in request.args.items()
                if k not in ['limit', 'offset', 'page', 'per_page']
            }
            
            paginated_response = paginate_query_response(
                response.data,
                total_count,
                limit,
                offset,
                resource_endpoint,
                query_params
            )
            
            return paginated_response, 200
        except Exception as e:
            if isinstance(e, ValueError):
                status_code = 400
                context = "Invalid query parameters"
            else:
                status_code = 500
                context = "Fetching mixtures list"
                
            error_response, error_status = handle_error(
                e,
                context=context,
                log_level='error',
                return_response=True,
                status_code=status_code
            )
            abort(error_status, **error_response)
    
    @doc(description='Create a new mixture with components.',
         tags=['Mixtures'],
         requestBody={
             'description': 'Mixture data with components',
             'required': True,
             'content': {
                 'application/json': {
                     'schema': {
                         'type': 'object',
                         'properties': {
                             'name': {
                                 'type': 'string',
                                 'description': 'Name of the mixture',
                                 'example': 'Cryoprotective Solution A'
                             },
                             'description': {
                                 'type': 'string',
                                 'description': 'Description of the mixture',
                                 'example': 'Standard cryoprotective solution for cell preservation'
                             },
                             'components': {
                                 'type': 'array',
                                 'description': 'List of components in the mixture',
                                 'items': {
                                     'type': 'object',
                                     'properties': {
                                         'molecule_id': {
                                             'type': 'string',
                                             'format': 'uuid',
                                             'description': 'ID of the molecule in this component',
                                             'example': '123e4567-e89b-12d3-a456-426614174003'
                                         },
                                         'concentration': {
                                             'type': 'number',
                                             'format': 'float',
                                             'description': 'Concentration value',
                                             'example': 10.5
                                         },
                                         'concentration_unit': {
                                             'type': 'string',
                                             'description': 'Unit of concentration (e.g., mg/mL, %w/v)',
                                             'example': 'mg/mL'
                                         }
                                     },
                                     'required': ['molecule_id', 'concentration', 'concentration_unit']
                                 }
                             }
                         },
                         'required': ['name', 'components']
                     },
                     'example': {
                         'name': 'Cryoprotective Solution A',
                         'description': 'Standard cryoprotective solution for cell preservation',
                         'components': [
                             {
                                 'molecule_id': '123e4567-e89b-12d3-a456-426614174003',
                                 'concentration': 10.5,
                                 'concentration_unit': 'mg/mL'
                             },
                             {
                                 'molecule_id': '123e4567-e89b-12d3-a456-426614174005',
                                 'concentration': 5.0,
                                 'concentration_unit': '%v/v'
                             }
                         ]
                     }
                 }
             }
         },
         responses={
             201: {
                 'description': 'Mixture successfully created',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the mixture'},
                                 'name': {'type': 'string', 'description': 'Name of the mixture'},
                                 'description': {'type': 'string', 'description': 'Description of the mixture'},
                                 'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                 'updated_at': {'type': 'string', 'format': 'date-time', 'description': 'Last update timestamp'},
                                 'created_by': {'type': 'string', 'format': 'uuid', 'description': 'ID of the user who created the mixture'},
                                 'components': {
                                     'type': 'array',
                                     'description': 'List of components in the mixture',
                                     'items': {
                                         'type': 'object',
                                         'properties': {
                                             'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the component'},
                                             'mixture_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the mixture this component belongs to'},
                                             'molecule_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the molecule in this component'},
                                             'concentration': {'type': 'number', 'format': 'float', 'description': 'Concentration value'},
                                             'concentration_unit': {'type': 'string', 'description': 'Unit of concentration (e.g., mg/mL, %w/v)'},
                                             'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                             'molecule': {
                                                 'type': 'object',
                                                 'description': 'Details of the molecule in this component',
                                                 'properties': {
                                                     'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the molecule'},
                                                     'name': {'type': 'string', 'description': 'Name of the molecule'},
                                                     'molecular_formula': {'type': 'string', 'description': 'Chemical formula'},
                                                     'smiles': {'type': 'string', 'description': 'SMILES notation representing molecular structure'}
                                                 }
                                             }
                                         }
                                     }
                                 }
                             }
                         },
                         'example': {
                             'id': '123e4567-e89b-12d3-a456-426614174000',
                             'name': 'Cryoprotective Solution A',
                             'description': 'Standard cryoprotective solution for cell preservation',
                             'created_at': '2025-01-15T10:30:00Z',
                             'updated_at': '2025-01-15T10:30:00Z',
                             'created_by': '123e4567-e89b-12d3-a456-426614174001',
                             'components': [
                                 {
                                     'id': '123e4567-e89b-12d3-a456-426614174002',
                                     'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                     'molecule_id': '123e4567-e89b-12d3-a456-426614174003',
                                     'concentration': 10.5,
                                     'concentration_unit': 'mg/mL',
                                     'created_at': '2025-01-15T10:30:00Z',
                                     'molecule': {
                                         'id': '123e4567-e89b-12d3-a456-426614174003',
                                         'name': 'Glycerol',
                                         'molecular_formula': 'C3H8O3',
                                         'smiles': 'C(C(CO)O)O'
                                     }
                                 },
                                 {
                                     'id': '123e4567-e89b-12d3-a456-426614174004',
                                     'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                     'molecule_id': '123e4567-e89b-12d3-a456-426614174005',
                                     'concentration': 5.0,
                                     'concentration_unit': '%v/v',
                                     'created_at': '2025-01-15T10:30:00Z',
                                     'molecule': {
                                         'id': '123e4567-e89b-12d3-a456-426614174005',
                                         'name': 'Dimethyl sulfoxide',
                                         'molecular_formula': 'C2H6OS',
                                         'smiles': 'CS(=O)C'
                                     }
                                 }
                             ]
                         }
                     }
                 }
             },
             400: {
                 'description': 'Invalid request parameters',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Validation error',
                             'message': 'Name is required. Components must have at least one item.',
                             'status': 400
                         }
                     }
                 }
             },
             401: {
                 'description': 'Authentication required',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Unauthorized',
                             'message': 'Authentication token is missing or invalid',
                             'status': 401
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error creating mixture',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(mixture_fields, code=201, description='Created mixture')
    @auth_required('bearerAuth')
    @token_required
    @marshal_with(mixture_fields)
    def post(self, user=None) -> Tuple[Dict[str, Any], int]:
        """Create a new mixture with components.
        
        Args:
            user: User object from token_required decorator
            
        Returns:
            Tuple[Dict[str, Any], int]: Created mixture object and HTTP status code
            
        Raises:
            400: If request data is invalid
            401: If authentication is missing or invalid
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            JSON object with mixture data following the MixtureSchema:
            - name (str): Name of the mixture (required)
            - description (str): Description of the mixture (optional)
            - components (list): List of components (required)
              - molecule_id (str): ID of the molecule (required)
              - concentration (float): Concentration value (required)
              - concentration_unit (str): Unit of concentration (required)
              
        Example Request:
            POST /api/mixtures
            Authorization: Bearer <token>
            Content-Type: application/json
            
            {
                "name": "Cryoprotective Solution A",
                "description": "Standard cryoprotective solution for cell preservation",
                "components": [
                    {
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174003",
                        "concentration": 10.5,
                        "concentration_unit": "mg/mL"
                    },
                    {
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174005",
                        "concentration": 5.0,
                        "concentration_unit": "%v/v"
                    }
                ]
            }
            
        Example Response:
            Status: 201 Created
            Content-Type: application/json
            
            {
                "id": "123e4567-e89b-12d3-a456-426614174000",
                "name": "Cryoprotective Solution A",
                "description": "Standard cryoprotective solution for cell preservation",
                "created_at": "2025-01-15T10:30:00Z",
                "updated_at": "2025-01-15T10:30:00Z",
                "created_by": "123e4567-e89b-12d3-a456-426614174001",
                "components": [
                    {
                        "id": "123e4567-e89b-12d3-a456-426614174002",
                        "mixture_id": "123e4567-e89b-12d3-a456-426614174000",
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174003",
                        "concentration": 10.5,
                        "concentration_unit": "mg/mL",
                        "created_at": "2025-01-15T10:30:00Z",
                        "molecule": {
                            "id": "123e4567-e89b-12d3-a456-426614174003",
                            "name": "Glycerol",
                            "molecular_formula": "C3H8O3",
                            "smiles": "C(C(CO)O)O"
                        }
                    },
                    {
                        "id": "123e4567-e89b-12d3-a456-426614174004",
                        "mixture_id": "123e4567-e89b-12d3-a456-426614174000",
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174005",
                        "concentration": 5.0,
                        "concentration_unit": "%v/v",
                        "created_at": "2025-01-15T10:30:00Z",
                        "molecule": {
                            "id": "123e4567-e89b-12d3-a456-426614174005",
                            "name": "Dimethyl sulfoxide",
                            "molecular_formula": "C2H6OS",
                            "smiles": "CS(=O)C"
                        }
                    }
                ]
            }
            
        Notes for Developers:
            - The endpoint creates a new mixture record and its associated components in a transaction-like pattern
            - If component creation fails, the mixture record is automatically rolled back
            - The response includes the complete mixture data with nested component and molecule details
            - Molecule IDs must exist in the database before they can be used in components
            - Concentration units should follow standard notation (e.g., mg/mL, %w/v, %v/v, mM)
            - The created_by field is automatically set to the authenticated user's ID
            
        Performance Optimization:
            - Uses Marshmallow schema for efficient validation
            - Performs database operations in a transaction-like pattern
            - Rolls back mixture creation if component addition fails
        """
        try:
            # Validate request data
            schema = MixtureSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating mixture data",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Create mixture
            response = supabase.table("mixtures").insert({
                "name": data["name"],
                "description": data.get("description", ""),
                "created_by": user_id
            }).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context="Creating mixture",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            mixture_id = response.data[0]["id"]
            
            # Add components
            component_inserts = [{
                "mixture_id": mixture_id,
                "molecule_id": comp["molecule_id"],
                "concentration": comp["concentration"],
                "concentration_unit": comp["concentration_unit"],
                "created_by": user_id
            } for comp in data["components"]]
            
            response = supabase.table("mixture_components").insert(component_inserts).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                # Rollback mixture creation
                supabase.table("mixtures").delete().eq("id", mixture_id).execute()
                error_response, error_status = handle_error(
                    error_message,
                    context="Adding mixture components",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            # Get the created mixture with components
            response = supabase.table("mixtures_with_components").select("*").eq("id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Retrieving created mixture (ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            if not response.data:
                error_response, error_status = handle_error(
                    f"Mixture with ID {mixture_id} not found",
                    context="Retrieving created mixture",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            return response.data[0], 201
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Creating mixture",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class MixtureResource(MethodResource):
    """Resource for retrieving, updating, and deleting a specific mixture.
    
    Provides endpoints for retrieving, updating, and deleting a mixture by ID.
    """
    
    @doc(description='Get a mixture with its components.',
         tags=['Mixtures'],
         params={'mixture_id': {'description': 'ID of the mixture to retrieve', 'type': 'string', 'required': True}},
         security=[{'bearerAuth': []}],
         responses={
             200: {
                 'description': 'Mixture with its components',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the mixture'},
                                 'name': {'type': 'string', 'description': 'Name of the mixture'},
                                 'description': {'type': 'string', 'description': 'Description of the mixture'},
                                 'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                 'updated_at': {'type': 'string', 'format': 'date-time', 'description': 'Last update timestamp'},
                                 'created_by': {'type': 'string', 'format': 'uuid', 'description': 'ID of the user who created the mixture'},
                                 'components': {
                                     'type': 'array',
                                     'description': 'List of components in the mixture',
                                     'items': {
                                         'type': 'object',
                                         'properties': {
                                             'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the component'},
                                             'mixture_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the mixture this component belongs to'},
                                             'molecule_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the molecule in this component'},
                                             'concentration': {'type': 'number', 'format': 'float', 'description': 'Concentration value'},
                                             'concentration_unit': {'type': 'string', 'description': 'Unit of concentration (e.g., mg/mL, %w/v)'},
                                             'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                             'molecule': {
                                                 'type': 'object',
                                                 'description': 'Details of the molecule in this component',
                                                 'properties': {
                                                     'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the molecule'},
                                                     'name': {'type': 'string', 'description': 'Name of the molecule'},
                                                     'molecular_formula': {'type': 'string', 'description': 'Chemical formula'},
                                                     'smiles': {'type': 'string', 'description': 'SMILES notation representing molecular structure'}
                                                 }
                                             }
                                         }
                                     }
                                 }
                             }
                         },
                         'example': {
                             'id': '123e4567-e89b-12d3-a456-426614174000',
                             'name': 'Cryoprotective Solution A',
                             'description': 'Standard cryoprotective solution for cell preservation',
                             'created_at': '2025-01-15T10:30:00Z',
                             'updated_at': '2025-01-15T10:30:00Z',
                             'created_by': '123e4567-e89b-12d3-a456-426614174001',
                             'components': [
                                 {
                                     'id': '123e4567-e89b-12d3-a456-426614174002',
                                     'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                     'molecule_id': '123e4567-e89b-12d3-a456-426614174003',
                                     'concentration': 10.5,
                                     'concentration_unit': 'mg/mL',
                                     'created_at': '2025-01-15T10:30:00Z',
                                     'molecule': {
                                         'id': '123e4567-e89b-12d3-a456-426614174003',
                                         'name': 'Glycerol',
                                         'molecular_formula': 'C3H8O3',
                                         'smiles': 'C(C(CO)O)O'
                                     }
                                 },
                                 {
                                     'id': '123e4567-e89b-12d3-a456-426614174004',
                                     'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                     'molecule_id': '123e4567-e89b-12d3-a456-426614174005',
                                     'concentration': 5.0,
                                     'concentration_unit': '%v/v',
                                     'created_at': '2025-01-15T10:30:00Z',
                                     'molecule': {
                                         'id': '123e4567-e89b-12d3-a456-426614174005',
                                         'name': 'Dimethyl sulfoxide',
                                         'molecular_formula': 'C2H6OS',
                                         'smiles': 'CS(=O)C'
                                     }
                                 }
                             ]
                         }
                     }
                 }
             },
             401: {
                 'description': 'Authentication required',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Unauthorized',
                             'message': 'Authentication token is missing or invalid',
                             'status': 401
                         }
                     }
                 }
             },
             403: {
                 'description': 'Insufficient permissions',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Forbidden',
                             'message': 'You do not have permission to access this mixture',
                             'status': 403
                         }
                     }
                 }
             },
             404: {
                 'description': 'Mixture not found',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Not found',
                             'message': 'Mixture with ID 123e4567-e89b-12d3-a456-426614174000 not found',
                             'status': 404
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error fetching mixture',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(mixture_fields, code=200, description='Mixture object')
    @marshal_with(mixture_fields)
    @auth_required('bearerAuth')
    @token_required
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """Get a mixture with its components.
        
        Args:
            mixture_id (str): ID of the mixture to retrieve
            
        Returns:
            Tuple[Dict[str, Any], int]: Mixture object if found and HTTP status code
            
        Raises:
            401: If authentication token is missing or invalid
            403: If user doesn't have permission to access this mixture
            404: If mixture not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Authorization:
            Any authenticated user can view mixtures they have created or have been granted access to.
            Administrators can view all mixtures.
        """
        try:
            supabase = get_supabase_client()
            response = supabase.table("mixtures_with_components").select("*").eq("id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Fetching mixture (ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            if not response.data:
                error_response, error_status = handle_error(
                    f"Mixture with ID {mixture_id} not found",
                    context="Fetching mixture",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            return response.data[0], 200
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context=f"Fetching mixture (ID: {mixture_id})",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)
    
    @doc(description='Update a mixture\'s information and components.',
         tags=['Mixtures'],
         params={'mixture_id': {'description': 'ID of the mixture to update', 'type': 'string', 'required': True}},
         requestBody={
             'description': 'Updated mixture data with components',
             'required': True,
             'content': {
                 'application/json': {
                     'schema': {
                         'type': 'object',
                         'properties': {
                             'name': {
                                 'type': 'string',
                                 'description': 'Name of the mixture',
                                 'example': 'Improved Cryoprotective Solution A'
                             },
                             'description': {
                                 'type': 'string',
                                 'description': 'Description of the mixture',
                                 'example': 'Enhanced cryoprotective solution for improved cell preservation'
                             },
                             'components': {
                                 'type': 'array',
                                 'description': 'List of components in the mixture',
                                 'items': {
                                     'type': 'object',
                                     'properties': {
                                         'molecule_id': {
                                             'type': 'string',
                                             'format': 'uuid',
                                             'description': 'ID of the molecule in this component',
                                             'example': '123e4567-e89b-12d3-a456-426614174003'
                                         },
                                         'concentration': {
                                             'type': 'number',
                                             'format': 'float',
                                             'description': 'Concentration value',
                                             'example': 12.0
                                         },
                                         'concentration_unit': {
                                             'type': 'string',
                                             'description': 'Unit of concentration (e.g., mg/mL, %w/v)',
                                             'example': 'mg/mL'
                                         }
                                     },
                                     'required': ['molecule_id', 'concentration', 'concentration_unit']
                                 }
                             }
                         },
                         'required': ['name', 'components']
                     },
                     'example': {
                         'name': 'Improved Cryoprotective Solution A',
                         'description': 'Enhanced cryoprotective solution for improved cell preservation',
                         'components': [
                             {
                                 'molecule_id': '123e4567-e89b-12d3-a456-426614174003',
                                 'concentration': 12.0,
                                 'concentration_unit': 'mg/mL'
                             },
                             {
                                 'molecule_id': '123e4567-e89b-12d3-a456-426614174005',
                                 'concentration': 7.5,
                                 'concentration_unit': '%v/v'
                             },
                             {
                                 'molecule_id': '123e4567-e89b-12d3-a456-426614174006',
                                 'concentration': 0.1,
                                 'concentration_unit': '%w/v'
                             }
                         ]
                     }
                 }
             }
         },
         responses={
             200: {
                 'description': 'Mixture successfully updated',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the mixture'},
                                 'name': {'type': 'string', 'description': 'Name of the mixture'},
                                 'description': {'type': 'string', 'description': 'Description of the mixture'},
                                 'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                 'updated_at': {'type': 'string', 'format': 'date-time', 'description': 'Last update timestamp'},
                                 'created_by': {'type': 'string', 'format': 'uuid', 'description': 'ID of the user who created the mixture'},
                                 'components': {
                                     'type': 'array',
                                     'description': 'List of components in the mixture',
                                     'items': {
                                         'type': 'object',
                                         'properties': {
                                             'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the component'},
                                             'mixture_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the mixture this component belongs to'},
                                             'molecule_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the molecule in this component'},
                                             'concentration': {'type': 'number', 'format': 'float', 'description': 'Concentration value'},
                                             'concentration_unit': {'type': 'string', 'description': 'Unit of concentration (e.g., mg/mL, %w/v)'},
                                             'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'},
                                             'molecule': {
                                                 'type': 'object',
                                                 'description': 'Details of the molecule in this component',
                                                 'properties': {
                                                     'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the molecule'},
                                                     'name': {'type': 'string', 'description': 'Name of the molecule'},
                                                     'molecular_formula': {'type': 'string', 'description': 'Chemical formula'},
                                                     'smiles': {'type': 'string', 'description': 'SMILES notation representing molecular structure'}
                                                 }
                                             }
                                         }
                                     }
                                 }
                             }
                         },
                         'example': {
                             'id': '123e4567-e89b-12d3-a456-426614174000',
                             'name': 'Improved Cryoprotective Solution A',
                             'description': 'Enhanced cryoprotective solution for improved cell preservation',
                             'created_at': '2025-01-15T10:30:00Z',
                             'updated_at': '2025-01-15T11:45:00Z',
                             'created_by': '123e4567-e89b-12d3-a456-426614174001',
                             'components': [
                                 {
                                     'id': '123e4567-e89b-12d3-a456-426614174007',
                                     'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                     'molecule_id': '123e4567-e89b-12d3-a456-426614174003',
                                     'concentration': 12.0,
                                     'concentration_unit': 'mg/mL',
                                     'created_at': '2025-01-15T11:45:00Z',
                                     'molecule': {
                                         'id': '123e4567-e89b-12d3-a456-426614174003',
                                         'name': 'Glycerol',
                                         'molecular_formula': 'C3H8O3',
                                         'smiles': 'C(C(CO)O)O'
                                     }
                                 },
                                 {
                                     'id': '123e4567-e89b-12d3-a456-426614174008',
                                     'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                     'molecule_id': '123e4567-e89b-12d3-a456-426614174005',
                                     'concentration': 7.5,
                                     'concentration_unit': '%v/v',
                                     'created_at': '2025-01-15T11:45:00Z',
                                     'molecule': {
                                         'id': '123e4567-e89b-12d3-a456-426614174005',
                                         'name': 'Dimethyl sulfoxide',
                                         'molecular_formula': 'C2H6OS',
                                         'smiles': 'CS(=O)C'
                                     }
                                 },
                                 {
                                     'id': '123e4567-e89b-12d3-a456-426614174009',
                                     'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                     'molecule_id': '123e4567-e89b-12d3-a456-426614174006',
                                     'concentration': 0.1,
                                     'concentration_unit': '%w/v',
                                     'created_at': '2025-01-15T11:45:00Z',
                                     'molecule': {
                                         'id': '123e4567-e89b-12d3-a456-426614174006',
                                         'name': 'Trehalose',
                                         'molecular_formula': 'C12H22O11',
                                         'smiles': 'C(C1C(C(C(C(O1)OC2C(OC(C2O)O)C(O)C(O)C(O)CO)O)O)O)O'
                                     }
                                 }
                             ]
                         }
                     }
                 }
             },
             400: {
                 'description': 'Invalid request parameters',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Validation error',
                             'message': 'Name is required. Components must have at least one item.',
                             'status': 400
                         }
                     }
                 }
             },
             401: {
                 'description': 'Authentication required',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Unauthorized',
                             'message': 'Authentication token is missing or invalid',
                             'status': 401
                         }
                     }
                 }
             },
             404: {
                 'description': 'Mixture not found',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Not found',
                             'message': 'Mixture with ID 123e4567-e89b-12d3-a456-426614174000 not found',
                             'status': 404
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error updating mixture',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(mixture_fields, code=200, description='Updated mixture')
    @auth_required('bearerAuth')
    @token_required
    @marshal_with(mixture_fields)
    def put(self, mixture_id: str, user=None) -> Tuple[Dict[str, Any], int]:
        """Update a mixture's information and components.
        
        Args:
            mixture_id (str): ID of the mixture to update
            user: User object from token_required decorator
            
        Returns:
            Tuple[Dict[str, Any], int]: Updated mixture object and HTTP status code
            
        Raises:
            400: If request data is invalid
            401: If authentication is missing or invalid
            404: If mixture not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            JSON object with mixture data following the MixtureSchema:
            - name (str): Name of the mixture (required)
            - description (str): Description of the mixture (optional)
            - components (list): List of components (required)
              - molecule_id (str): ID of the molecule (required)
              - concentration (float): Concentration value (required)
              - concentration_unit (str): Unit of concentration (required)
              
        Example Request:
            PUT /api/mixtures/{mixture_id}
            Authorization: Bearer <token>
            Content-Type: application/json
            
            {
                "name": "Improved Cryoprotective Solution A",
                "description": "Enhanced cryoprotective solution for improved cell preservation",
                "components": [
                    {
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174003",
                        "concentration": 12.0,
                        "concentration_unit": "mg/mL"
                    },
                    {
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174005",
                        "concentration": 7.5,
                        "concentration_unit": "%v/v"
                    },
                    {
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174006",
                        "concentration": 0.1,
                        "concentration_unit": "%w/v"
                    }
                ]
            }
            
        Example Response:
            Status: 200 OK
            Content-Type: application/json
            
            {
                "id": "123e4567-e89b-12d3-a456-426614174000",
                "name": "Improved Cryoprotective Solution A",
                "description": "Enhanced cryoprotective solution for improved cell preservation",
                "created_at": "2025-01-15T10:30:00Z",
                "updated_at": "2025-01-15T11:45:00Z",
                "created_by": "123e4567-e89b-12d3-a456-426614174001",
                "components": [
                    {
                        "id": "123e4567-e89b-12d3-a456-426614174007",
                        "mixture_id": "123e4567-e89b-12d3-a456-426614174000",
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174003",
                        "concentration": 12.0,
                        "concentration_unit": "mg/mL",
                        "created_at": "2025-01-15T11:45:00Z",
                        "molecule": {
                            "id": "123e4567-e89b-12d3-a456-426614174003",
                            "name": "Glycerol",
                            "molecular_formula": "C3H8O3",
                            "smiles": "C(C(CO)O)O"
                        }
                    },
                    {
                        "id": "123e4567-e89b-12d3-a456-426614174008",
                        "mixture_id": "123e4567-e89b-12d3-a456-426614174000",
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174005",
                        "concentration": 7.5,
                        "concentration_unit": "%v/v",
                        "created_at": "2025-01-15T11:45:00Z",
                        "molecule": {
                            "id": "123e4567-e89b-12d3-a456-426614174005",
                            "name": "Dimethyl sulfoxide",
                            "molecular_formula": "C2H6OS",
                            "smiles": "CS(=O)C"
                        }
                    },
                    {
                        "id": "123e4567-e89b-12d3-a456-426614174009",
                        "mixture_id": "123e4567-e89b-12d3-a456-426614174000",
                        "molecule_id": "123e4567-e89b-12d3-a456-426614174006",
                        "concentration": 0.1,
                        "concentration_unit": "%w/v",
                        "created_at": "2025-01-15T11:45:00Z",
                        "molecule": {
                            "id": "123e4567-e89b-12d3-a456-426614174006",
                            "name": "Trehalose",
                            "molecular_formula": "C12H22O11",
                            "smiles": "C(C1C(C(C(C(O1)OC2C(OC(C2O)O)C(O)C(O)C(O)CO)O)O)O)O"
                        }
                    }
                ]
            }
            
        Notes for Developers:
            - The endpoint completely replaces all components of the mixture
            - Existing components are deleted and new ones are created based on the request
            - All component IDs will change after an update, even if the molecule_id remains the same
            - The updated_at timestamp is automatically updated to the current time
            - The created_by field remains unchanged from the original creation
            - Molecule IDs must exist in the database before they can be used in components
        """
        try:
            # Validate request data
            schema = MixtureSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating mixture update data",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Check if mixture exists
            check_response = supabase.table("mixtures").select("id").eq("id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(check_response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Verifying mixture (ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            if not check_response.data:
                error_response, error_status = handle_error(
                    f"Mixture with ID {mixture_id} not found",
                    context="Updating mixture",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Update mixture basic info
            mixture_update = {
                "name": data["name"],
                "description": data.get("description", ""),
                "updated_at": datetime.now().isoformat()
            }
            
            response = supabase.table("mixtures").update(mixture_update).eq("id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Updating mixture info (ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            # Delete existing components
            response = supabase.table("mixture_components").delete().eq("mixture_id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Deleting mixture components (Mixture ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            # Add new components
            component_inserts = [{
                "mixture_id": mixture_id,
                "molecule_id": comp["molecule_id"],
                "concentration": comp["concentration"],
                "concentration_unit": comp["concentration_unit"],
                "created_by": user_id
            } for comp in data["components"]]
            
            response = supabase.table("mixture_components").insert(component_inserts).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Adding new mixture components (Mixture ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            # Get the updated mixture with components
            response = supabase.table("mixtures_with_components").select("*").eq("id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Retrieving updated mixture (ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            if not response.data:
                error_response, error_status = handle_error(
                    f"Mixture with ID {mixture_id} not found",
                    context="Retrieving updated mixture",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            current_app.logger.info(f"Mixture {mixture_id} updated successfully")
            return response.data[0], 200
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context=f"Updating mixture (ID: {mixture_id})",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)
    
    @doc(description='Delete a mixture and its components.',
         tags=['Mixtures'],
         params={'mixture_id': {'description': 'ID of the mixture to delete', 'type': 'string', 'required': True}},
         responses={
             200: {
                 'description': 'Mixture successfully deleted',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'message': {'type': 'string', 'description': 'Success message'}
                             }
                         },
                         'example': {
                             'message': 'Mixture with ID 123e4567-e89b-12d3-a456-426614174000 deleted successfully'
                         }
                     }
                 }
             },
             401: {
                 'description': 'Authentication required',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Unauthorized',
                             'message': 'Authentication token is missing or invalid',
                             'status': 401
                         }
                     }
                 }
             },
             404: {
                 'description': 'Mixture not found',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Not found',
                             'message': 'Mixture with ID 123e4567-e89b-12d3-a456-426614174000 not found',
                             'status': 404
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error deleting mixture',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(fields.Raw, code=200, description='Success message')
    @auth_required('bearerAuth')
    @token_required
    def delete(self, mixture_id: str, user=None) -> Tuple[Dict[str, str], int]:
        """Delete a mixture and its components.
        
        Args:
            mixture_id (str): ID of the mixture to delete
            user: User object from token_required decorator
            
        Returns:
            Tuple[Dict[str, str], int]: Success message and HTTP status code
            
        Raises:
            401: If authentication is missing or invalid
            404: If mixture not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Example Request:
            DELETE /api/mixtures/{mixture_id}
            Authorization: Bearer <token>
            
        Example Response:
            Status: 200 OK
            Content-Type: application/json
            
            {
                "message": "Mixture with ID 123e4567-e89b-12d3-a456-426614174000 deleted successfully"
            }
            
        Notes for Developers:
            - This operation permanently deletes the mixture and all its components
            - The deletion is performed using a cascading delete in the database
            - Any predictions or experiments associated with this mixture will also be deleted
            - This operation cannot be undone, so appropriate confirmation should be implemented in the client
        """
        try:
            supabase = get_supabase_client()
            
            # Check if mixture exists
            check_response = supabase.table("mixtures").select("id").eq("id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(check_response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Verifying mixture (ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            if not check_response.data:
                error_response, error_status = handle_error(
                    f"Mixture with ID {mixture_id} not found",
                    context="Deleting mixture",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Delete the mixture (components will be deleted via cascade)
            response = supabase.table("mixtures").delete().eq("id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context=f"Deleting mixture (ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
                abort(error_status, **error_response)
            
            current_app.logger.info(f"Mixture {mixture_id} deleted successfully")
            return {"message": f"Mixture with ID {mixture_id} deleted successfully"}, 200
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context=f"Deleting mixture (ID: {mixture_id})",
                log_level='error',
                return_response=True
            )
            abort(error_status, **error_response)


class PredictionResource(MethodResource):
    """Resource for creating, retrieving, updating, and deleting predictions.
    
    Provides endpoints for managing predictions associated with mixtures.
    Predictions represent calculated or estimated property values for mixtures
    using various calculation methods.
    """
    
    # Request parser for prediction operations
    prediction_parser = reqparse.RequestParser()
    prediction_parser.add_argument('property_name', type=str, required=True,
                                  help='Property name is required')
    prediction_parser.add_argument('calculation_method', type=str, required=True,
                                  help='Calculation method is required')
    prediction_parser.add_argument('value', required=True,
                                  help='Value is required')
    prediction_parser.add_argument('confidence', type=float, required=True,
                                  help='Confidence score is required')
    
    @doc(description='Get a paginated list of predictions for a mixture.',
         tags=['Predictions'],
         params={'mixture_id': {'description': 'ID of the mixture to get predictions for', 'type': 'string', 'required': True},
                'limit': {'description': 'Maximum number of records to return (max 500)', 'type': 'integer', 'default': 100},
                'offset': {'description': 'Number of records to skip', 'type': 'integer', 'default': 0}},
         security=[{'bearerAuth': []}],
         responses={
             200: {
                 'description': 'List of predictions for the mixture',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'array',
                             'items': {
                                 'type': 'object',
                                 'properties': {
                                     'id': {'type': 'string', 'format': 'uuid', 'description': 'Unique identifier for the prediction'},
                                     'mixture_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the mixture this prediction is for'},
                                     'property_type_id': {'type': 'string', 'format': 'uuid', 'description': 'ID of the property type'},
                                     'property_name': {'type': 'string', 'description': 'Name of the property'},
                                     'numeric_value': {'type': 'number', 'format': 'float', 'description': 'Numeric value (if property is numeric)'},
                                     'text_value': {'type': 'string', 'description': 'Text value (if property is text)'},
                                     'boolean_value': {'type': 'boolean', 'description': 'Boolean value (if property is boolean)'},
                                     'confidence': {'type': 'number', 'format': 'float', 'description': 'Confidence score for the prediction'},
                                     'calculation_method': {'type': 'string', 'description': 'Method used for calculation'},
                                     'created_at': {'type': 'string', 'format': 'date-time', 'description': 'Creation timestamp'}
                                 }
                             }
                         },
                         'example': [
                             {
                                 'id': '123e4567-e89b-12d3-a456-426614174020',
                                 'mixture_id': '123e4567-e89b-12d3-a456-426614174000',
                                 'property_type_id': '123e4567-e89b-12d3-a456-426614174011',
                                 'property_name': 'freezing_point',
                                 'numeric_value': -20.5,
                                 'text_value': None,
                                 'boolean_value': None,
                                 'confidence': 0.85,
                                 'calculation_method': 'molecular_dynamics',
                                 'created_at': '2025-01-15T11:00:00Z'
                             }
                         ]
                     }
                 }
             },
             400: {
                 'description': 'Invalid pagination parameters',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Invalid pagination parameters',
                             'message': 'Limit must be a positive integer not exceeding 500',
                             'status': 400
                         }
                     }
                 }
             },
             401: {
                 'description': 'Authentication required',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Unauthorized',
                             'message': 'Authentication token is missing or invalid',
                             'status': 401
                         }
                     }
                 }
             },
             403: {
                 'description': 'Insufficient permissions',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Forbidden',
                             'message': 'You do not have permission to access predictions for this mixture',
                             'status': 403
                         }
                     }
                 }
             },
             500: {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'example': {
                             'error': 'Database error',
                             'message': 'Error fetching predictions',
                             'status': 500
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(prediction_fields, code=200, description='List of predictions')
    @marshal_with(prediction_fields)
    @auth_required('bearerAuth')
    @token_required
    def get(self, mixture_id: str) -> Tuple[List[Dict[str, Any]], int]:
        """Get a paginated list of predictions for a mixture.
        
        Args:
            mixture_id (str): ID of the mixture to get predictions for
            
        Returns:
            Tuple[List[Dict[str, Any]], int]: List of prediction objects and HTTP status code
            
        Raises:
            400: If pagination parameters are invalid
            401: If authentication token is missing or invalid
            403: If user doesn't have permission to access predictions for this mixture
            500: If database error occurs

        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Authorization:
            Users can only access predictions for mixtures they have created or have been granted access to.
            Administrators can access all predictions.

        Performance Optimization:
        - Added pagination (limit, offset) to avoid loading all predictions at once.
        - This reduces memory usage and improves response time for mixtures with many predictions.
        """
        try:
            supabase = get_supabase_client()
            # Get pagination parameters from query args, with sensible defaults
            limit = min(int(request.args.get("limit", 100)), 500)
            offset = int(request.args.get("offset", 0))
            # Join with property_types and calculation_methods to get names
            response = supabase.from_("predictions").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value,
                confidence, created_at,
                property_types(name), calculation_methods(name)
            """).eq("mixture_id", mixture_id).range(offset, offset + limit - 1).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                return handle_error(
                    error_message,
                    context=f"Fetching predictions (Mixture ID: {mixture_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
            
            # Transform the response to match the prediction_fields structure
            predictions = []
            for pred in response.data:
                prediction = {
                    'id': pred['id'],
                    'mixture_id': pred['mixture_id'],
                    'property_type_id': pred['property_type_id'],
                    'property_name': pred['property_types']['name'],
                    'numeric_value': pred['numeric_value'],
                    'text_value': pred['text_value'],
                    'boolean_value': pred['boolean_value'],
                    'confidence': pred['confidence'],
                    'calculation_method': pred['calculation_methods']['name'],
                    'created_at': pred['created_at']
                }
                predictions.append(prediction)
            
            return predictions, 200
        except Exception as e:
            return handle_error(
                e,
                context=f"Fetching predictions (Mixture ID: {mixture_id})",
                log_level='error',
                return_response=True
            )
    
    @doc(description='Add a prediction for a mixture.',
         tags=['Predictions'],
         params={'mixture_id': {'description': 'ID of the mixture to add prediction for', 'type': 'string', 'required': True}},
         request_body=PredictionSchema)
    @apispec_marshal_with(prediction_fields, code=201, description='Created prediction')
    @auth_required('bearerAuth')
    @token_required
    @marshal_with(prediction_fields)
    def post(self, mixture_id: str, user=None) -> Tuple[Dict[str, Any], int]:
        """Add a prediction for a mixture.
        
        Args:
            mixture_id (str): ID of the mixture to add prediction for
            user: User object from token_required decorator
            
        Returns:
            Tuple[Dict[str, Any], int]: Created prediction object and HTTP status code
            
        Raises:
            400: If request data is invalid
            404: If property type or calculation method not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            JSON object with prediction data following the PredictionSchema:
            - property_name (str): Name of the property (required)
            - calculation_method (str): Method used for calculation (required)
            - value (mixed): Property value (required)
            - confidence (float): Confidence score (required)
        """
        try:
            # Validate request data
            schema = PredictionSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return handle_error(
                    err,
                    context="Validating prediction data",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Get property type ID
            response = supabase.table("property_types").select("id, data_type").eq("name", data["property_name"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                return handle_error(
                    error_message,
                    context=f"Fetching property type '{data['property_name']}'",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
            
            if not response.data:
                return handle_error(
                    f"Property type '{data['property_name']}' not found",
                    context="Adding prediction",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
            
            property_type = response.data[0]
            
            # Get calculation method ID
            response = supabase.table("calculation_methods").select("id").eq("name", data["calculation_method"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                return handle_error(
                    error_message,
                    context=f"Fetching calculation method '{data['calculation_method']}'",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
            
            if not response.data:
                return handle_error(
                    f"Calculation method '{data['calculation_method']}' not found",
                    context="Adding prediction",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
            
            calculation_method = response.data[0]
            
            # Prepare prediction insert
            prediction = {
                "mixture_id": mixture_id,
                "property_type_id": property_type["id"],
                "calculation_method_id": calculation_method["id"],
                "confidence": data["confidence"],
                "created_by": user_id
            }
            
            # Set the appropriate value field based on data type
            if property_type["data_type"] == "numeric":
                prediction["numeric_value"] = float(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "text":
                prediction["text_value"] = str(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "boolean":
                prediction["boolean_value"] = bool(data["value"]) if data["value"] is not None else None
            else:
                return handle_error(
                    f"Unknown data type '{property_type['data_type']}'",
                    context="Adding prediction",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
            
            # Insert prediction
            response = supabase.table("predictions").upsert(
                prediction,
                on_conflict="mixture_id,property_type_id,calculation_method_id"
            ).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                return handle_error(
                    error_message,
                    context="Inserting prediction",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
            
            # Get the created prediction with property and method names
            prediction_id = response.data[0]["id"]
            response = supabase.from_("predictions").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value,
                confidence, created_at,
                property_types(name), calculation_methods(name)
            """).eq("id", prediction_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                return handle_error(
                    error_message,
                    context=f"Retrieving created prediction (ID: {prediction_id})",
                    log_level='error',
                    return_response=True,
                    status_code=status_code
                )
            
            if not response.data:
                return handle_error(
                    f"Prediction with ID {prediction_id} not found",
                    context="Retrieving created prediction",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
            
            # Transform the response to match the prediction_fields structure
            pred = response.data[0]
            prediction = {
                'id': pred['id'],
                'mixture_id': pred['mixture_id'],
                'property_type_id': pred['property_type_id'],
                'property_name': pred['property_types']['name'],
                'numeric_value': pred['numeric_value'],
                'text_value': pred['text_value'],
                'boolean_value': pred['boolean_value'],
                'confidence': pred['confidence'],
                'calculation_method': pred['calculation_methods']['name'],
                'created_at': pred['created_at']
            }
            
            return prediction, 201
        except Exception as e:
            return handle_error(
                e,
                context="Adding prediction",
                log_level='error',
                return_response=True
            )
    
    @doc(description='Update an existing prediction.',
         tags=['Predictions'],
         params={'mixture_id': {'description': 'ID of the mixture the prediction belongs to', 'type': 'string', 'required': True},
                'prediction_id': {'description': 'ID of the prediction to update (optional)', 'type': 'string'}},
         request_body=PredictionSchema)
    @apispec_marshal_with(prediction_fields, code=200, description='Updated prediction')
    @auth_required('bearerAuth')
    @token_required
    @marshal_with(prediction_fields)
    def put(self, mixture_id: str, prediction_id: Optional[str]=None, user=None) -> Tuple[Dict[str, Any], int]:
        """Update an existing prediction.
        
        Args:
            mixture_id (str): ID of the mixture the prediction belongs to
            prediction_id (Optional[str]): Optional ID of the prediction to update
            user: User object from token_required decorator
            
        Returns:
            Tuple[Dict[str, Any], int]: Updated prediction object and HTTP status code
            
        Raises:
            400: If request data is invalid
            404: If prediction, property type, or calculation method not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            JSON object with prediction data following the PredictionSchema:
            - property_name (str): Name of the property (required)
            - calculation_method (str): Method used for calculation (required)
            - value (mixed): Property value (required)
            - confidence (float): Confidence score (required)
            
        Notes:
            If prediction_id is not provided, the system will attempt to find the prediction
            using the property_name and calculation_method from the request body.
        """
        try:
            # Validate request data
            schema = PredictionSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return handle_error(
                    err,
                    context="Validating prediction update data",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # If prediction_id is not provided, try to find by property name and calculation method
            if prediction_id is None:
                # Get property type ID
                property_response = supabase.table("property_types").select("id").eq("name", data["property_name"]).execute()
                
                error_message, status_code = handle_supabase_error(property_response)
                if error_message:
                    error_response, error_status = handle_error(
                        error_message,
                        context="Fetching property type ID for prediction update",
                        log_level='error',
                        return_response=True
                    )
                    abort(error_status, **error_response)
                
                if not property_response.data:
                    error_response, error_status = handle_error(
                        f"Property type '{data['property_name']}' not found",
                        context="Fetching property type ID for prediction update",
                        log_level='error',
                        return_response=True,
                        status_code=404
                    )
                    abort(error_status, **error_response)
                
                property_type_id = property_response.data[0]["id"]
                
                # Get calculation method ID
                method_response = supabase.table("calculation_methods").select("id").eq("name", data["calculation_method"]).execute()
                
                error_message, status_code = handle_supabase_error(method_response)
                if error_message:
                    error_response, error_status = handle_error(
                        error_message,
                        context="Fetching calculation method ID for prediction update",
                        log_level='error',
                        return_response=True
                    )
                    abort(error_status, **error_response)
                
                if not method_response.data:
                    error_response, error_status = handle_error(
                        f"Calculation method '{data['calculation_method']}' not found",
                        context="Fetching calculation method ID for prediction update",
                        log_level='error',
                        return_response=True,
                        status_code=404
                    )
                    abort(error_status, **error_response)
                
                calculation_method_id = method_response.data[0]["id"]
                
                # Find the prediction
                pred_response = supabase.table("predictions").select("id").eq("mixture_id", mixture_id).eq("property_type_id", property_type_id).eq("calculation_method_id", calculation_method_id).execute()
                
                error_message, status_code = handle_supabase_error(pred_response)
                if error_message:
                    error_response, error_status = handle_error(
                        error_message,
                        context="Finding prediction for update",
                        log_level='error',
                        return_response=True
                    )
                    abort(error_status, **error_response)
                
                if not pred_response.data:
                    error_response, error_status = handle_error(
                        f"Prediction for mixture {mixture_id} with property '{data['property_name']}' and calculation method '{data['calculation_method']}' not found",
                        context="Finding prediction for update",
                        log_level='error',
                        return_response=True,
                        status_code=404
                    )
                    abort(error_status, **error_response)
                
                prediction_id = pred_response.data[0]["id"]
            
            # Get property type ID and data type
            response = supabase.table("property_types").select("id, data_type").eq("name", data["property_name"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context="Fetching property type for prediction update",
                    log_level='error',
                    return_response=True
                )
                abort(error_status, **error_response)
            
            if not response.data:
                error_response, error_status = handle_error(
                    f"Property type '{data['property_name']}' not found",
                    context="Fetching property type for prediction update",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            property_type = response.data[0]
            
            # Get calculation method ID
            response = supabase.table("calculation_methods").select("id").eq("name", data["calculation_method"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context="Fetching calculation method for prediction update",
                    log_level='error',
                    return_response=True
                )
                abort(error_status, **error_response)
            
            if not response.data:
                error_response, error_status = handle_error(
                    f"Calculation method '{data['calculation_method']}' not found",
                    context="Fetching calculation method for prediction update",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            calculation_method = response.data[0]
            
            # Prepare prediction update
            prediction = {
                "property_type_id": property_type["id"],
                "calculation_method_id": calculation_method["id"],
                "confidence": data["confidence"],
                "updated_at": datetime.now().isoformat()
            }
            
            # Set the appropriate value field based on data type
            if property_type["data_type"] == "numeric":
                prediction["numeric_value"] = float(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "text":
                prediction["text_value"] = str(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "boolean":
                prediction["boolean_value"] = bool(data["value"]) if data["value"] is not None else None
            else:
                error_response, error_status = handle_error(
                    f"Unknown data type '{property_type['data_type']}'",
                    context="Updating prediction value",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Update prediction
            response = supabase.table("predictions").update(prediction).eq("id", prediction_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context="Updating prediction in database",
                    log_level='error',
                    return_response=True
                )
                abort(error_status, **error_response)
            
            # Get the updated prediction with property and method names
            response = supabase.from_("predictions").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value,
                confidence, created_at,
                property_types(name), calculation_methods(name)
            """).eq("id", prediction_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context="Retrieving updated prediction",
                    log_level='error',
                    return_response=True
                )
                abort(error_status, **error_response)
            
            if not response.data:
                error_response, error_status = handle_error(
                    f"Prediction with ID {prediction_id} not found",
                    context="Retrieving updated prediction",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Transform the response to match the prediction_fields structure
            pred = response.data[0]
            prediction = {
                'id': pred['id'],
                'mixture_id': pred['mixture_id'],
                'property_type_id': pred['property_type_id'],
                'property_name': pred['property_types']['name'],
                'numeric_value': pred['numeric_value'],
                'text_value': pred['text_value'],
                'boolean_value': pred['boolean_value'],
                'confidence': pred['confidence'],
                'calculation_method': pred['calculation_methods']['name'],
                'created_at': pred['created_at']
            }
            
            current_app.logger.info(f"Prediction {prediction_id} updated successfully")
            return _handle_json_serialization(prediction), 200
        except Exception as e:
            current_app.logger.error(f"Error updating prediction: {str(e)}", exc_info=True)
            error_response, error_status = handle_error(
                f"Error updating prediction: {str(e)}",
                context="Updating prediction",
                log_level='error',
                return_response=True,
                status_code=500
            )
            abort(error_status, **error_response)
    
    @doc(description='Delete a prediction.',
         tags=['Predictions'],
         params={'mixture_id': {'description': 'ID of the mixture the prediction belongs to', 'type': 'string', 'required': True},
                'prediction_id': {'description': 'ID of the prediction to delete (optional)', 'type': 'string'}})
    @apispec_marshal_with(fields.Raw, code=200, description='Success message')
    @auth_required('bearerAuth')
    @token_required
    def delete(self, mixture_id: str, prediction_id: Optional[str]=None) -> Tuple[Dict[str, str], int]:
        """Delete a prediction.
        
        Args:
            mixture_id (str): ID of the mixture the prediction belongs to
            prediction_id (Optional[str]): Optional ID of the prediction to delete
            
        Returns:
            Tuple[Dict[str, str], int]: Success message and HTTP status code
            
        Raises:
            400: If neither prediction_id nor property_name and calculation_method are provided
            404: If prediction not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Notes:
            If prediction_id is not provided, the request body must include property_name
            and calculation_method to identify the prediction to delete.
        """
        try:
            supabase = get_supabase_client()
            
            # If prediction_id is not provided, try to find by property name and calculation method
            if prediction_id is None and request.json:
                data = request.json
                property_name = data.get("property_name")
                calculation_method = data.get("calculation_method")
                
                if not property_name or not calculation_method:
                    error_response, error_status = handle_error(
                        "Either prediction_id or both property_name and calculation_method must be provided",
                        context="Deleting prediction",
                        log_level='error',
                        return_response=True,
                        status_code=400
                    )
                    abort(error_status, **error_response)
                
                # Get property type ID
                property_response = supabase.table("property_types").select("id").eq("name", property_name).execute()
                
                error_message, status_code = handle_supabase_error(property_response)
                if error_message:
                    error_response, error_status = handle_error(
                        error_message,
                        context="Fetching property type ID for prediction delete",
                        log_level='error',
                        return_response=True
                    )
                    abort(error_status, **error_response)
                
                if not property_response.data:
                    error_response, error_status = handle_error(
                        f"Property type '{property_name}' not found",
                        context="Fetching property type ID for prediction delete",
                        log_level='error',
                        return_response=True,
                        status_code=404
                    )
                    abort(error_status, **error_response)
                
                property_type_id = property_response.data[0]["id"]
                
                # Get calculation method ID
                method_response = supabase.table("calculation_methods").select("id").eq("name", calculation_method).execute()
                
                error_message, status_code = handle_supabase_error(method_response)
                if error_message:
                    error_response, error_status = handle_error(
                        error_message,
                        context="Fetching calculation method ID for prediction delete",
                        log_level='error',
                        return_response=True
                    )
                    abort(error_status, **error_response)
                
                if not method_response.data:
                    error_response, error_status = handle_error(
                        f"Calculation method '{calculation_method}' not found",
                        context="Fetching calculation method ID for prediction delete",
                        log_level='error',
                        return_response=True,
                        status_code=404
                    )
                    abort(error_status, **error_response)
                
                calculation_method_id = method_response.data[0]["id"]
                
                # Find the prediction
                pred_response = supabase.table("predictions").select("id").eq("mixture_id", mixture_id).eq("property_type_id", property_type_id).eq("calculation_method_id", calculation_method_id).execute()
                
                error_message, status_code = handle_supabase_error(pred_response)
                if error_message:
                    error_response, error_status = handle_error(
                        error_message,
                        context="Finding prediction for delete",
                        log_level='error',
                        return_response=True
                    )
                    abort(error_status, **error_response)
                
                if not pred_response.data:
                    error_response, error_status = handle_error(
                        f"Prediction for mixture {mixture_id} with property '{property_name}' and calculation method '{calculation_method}' not found",
                        context="Finding prediction for delete",
                        log_level='error',
                        return_response=True,
                        status_code=404
                    )
                    abort(error_status, **error_response)
                
                prediction_id = pred_response.data[0]["id"]
            elif prediction_id is None:
                error_response, error_status = handle_error(
                    "Either prediction_id or both property_name and calculation_method must be provided",
                    context="Deleting prediction",
                    log_level='error',
                    return_response=True,
                    status_code=400
                )
                abort(error_status, **error_response)
            
            # Check if prediction exists
            check_response = supabase.table("predictions").select("id").eq("id", prediction_id).execute()
            
            error_message, status_code = handle_supabase_error(check_response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context="Checking prediction existence for delete",
                    log_level='error',
                    return_response=True
                )
                abort(error_status, **error_response)
            
            if not check_response.data:
                error_response, error_status = handle_error(
                    f"Prediction with ID {prediction_id} not found",
                    context="Checking prediction existence for delete",
                    log_level='error',
                    return_response=True,
                    status_code=404
                )
                abort(error_status, **error_response)
            
            # Delete the prediction
            response = supabase.table("predictions").delete().eq("id", prediction_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                error_response, error_status = handle_error(
                    error_message,
                    context="Deleting prediction in database",
                    log_level='error',
                    return_response=True
                )
                abort(error_status, **error_response)
            
            current_app.logger.info(f"Prediction {prediction_id} deleted successfully")
            return {"message": f"Prediction with ID {prediction_id} deleted successfully"}, 200
        except Exception as e:
            current_app.logger.error(f"Error deleting prediction: {str(e)}", exc_info=True)
            error_response, error_status = handle_error(
                f"Error deleting prediction: {str(e)}",
                context="Deleting prediction",
                log_level='error',
                return_response=True,
                status_code=500
            )
            abort(error_status, **error_response)


class ExperimentResource(MethodResource):
    """Resource for creating, retrieving, updating, and deleting experiments.
    
    Provides endpoints for managing experimental data associated with mixtures.
    Experiments represent measured property values for mixtures under specific
    experimental conditions.
    """
    
    @doc(description='Get experiments for a mixture.',
         tags=['Experiments'],
         params={'mixture_id': {'description': 'ID of the mixture to get experiments for', 'type': 'string', 'required': True}},
         security=[{'bearerAuth': []}])
    @apispec_marshal_with(experiment_fields, code=200, description='List of experiments')
    @marshal_with(experiment_fields)
    @auth_required('bearerAuth')
    @token_required
    def get(self, mixture_id: str) -> Tuple[List[Dict[str, Any]], int]:
        """Get experiments for a mixture.
        
        Args:
            mixture_id (str): ID of the mixture to get experiments for
            
        Returns:
            Tuple[List[Dict[str, Any]], int]: List of experiment objects and HTTP status code
            
        Raises:
            401: If authentication token is missing or invalid
            403: If user doesn't have permission to access experiments for this mixture
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Authorization:
            Users can only access experiments for mixtures they have created or have been granted access to.
            Administrators can access all experiments.
        """
        try:
            supabase = get_supabase_client()
            
            # Join with property_types to get names
            response = supabase.from_("experiments").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value, 
                experimental_conditions, date_performed, created_at, 
                property_types(name)
            """).eq("mixture_id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Transform the response to match the experiment_fields structure
            experiments = []
            for exp in response.data:
                experiment = {
                    'id': exp['id'],
                    'mixture_id': exp['mixture_id'],
                    'property_type_id': exp['property_type_id'],
                    'property_name': exp['property_types']['name'],
                    'numeric_value': exp['numeric_value'],
                    'text_value': exp['text_value'],
                    'boolean_value': exp['boolean_value'],
                    'experimental_conditions': exp['experimental_conditions'],
                    'date_performed': exp['date_performed'],
                    'created_at': exp['created_at']
                }
                experiments.append(experiment)
            
            return _handle_json_serialization(experiments), 200
        except Exception as e:
            current_app.logger.error(f"Error fetching experiments: {str(e)}", exc_info=True)
            abort(500, message=f"Error fetching experiments: {str(e)}")
    
    @doc(description='Record an experiment for a mixture.',
         tags=['Experiments'],
         params={'mixture_id': {'description': 'ID of the mixture to record experiment for', 'type': 'string', 'required': True}},
         request_body=ExperimentSchema)
    @apispec_marshal_with(experiment_fields, code=201, description='Created experiment')
    @auth_required('bearerAuth')
    @token_required
    @marshal_with(experiment_fields)
    def post(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """Record an experiment for a mixture.
        
        Args:
            mixture_id (str): ID of the mixture to record experiment for
            
        Returns:
            Tuple[Dict[str, Any], int]: Created experiment object and HTTP status code
            
        Raises:
            400: If request data is invalid
            404: If property type not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            JSON object with experiment data following the ExperimentSchema:
            - property_name (str): Name of the property (required)
            - value (mixed): Property value (required)
            - experimental_conditions (str): Description of experimental conditions (optional)
            - date_performed (date): Date when experiment was performed (required)
        """
        try:
            # Validate request data
            schema = ExperimentSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Get property type ID
            response = supabase.table("property_types").select("id, data_type").eq("name", data["property_name"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Property type '{data['property_name']}' not found")
            
            property_type = response.data[0]
            
            # Prepare experiment insert
            experiment = {
                "mixture_id": mixture_id,
                "property_type_id": property_type["id"],
                "experimental_conditions": data.get("experimental_conditions", ""),
                "date_performed": data["date_performed"].isoformat(),
                "created_by": user_id
            }
            
            # Set the appropriate value field based on data type
            if property_type["data_type"] == "numeric":
                experiment["numeric_value"] = float(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "text":
                experiment["text_value"] = str(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "boolean":
                experiment["boolean_value"] = bool(data["value"]) if data["value"] is not None else None
            else:
                abort(400, message=f"Unknown data type '{property_type['data_type']}'")
            
            # Insert experiment
            response = supabase.table("experiments").insert(experiment).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Get the created experiment with property name
            experiment_id = response.data[0]["id"]
            response = supabase.from_("experiments").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value, 
                experimental_conditions, date_performed, created_at, 
                property_types(name)
            """).eq("id", experiment_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Experiment with ID {experiment_id} not found")
            
            # Transform the response to match the experiment_fields structure
            exp = response.data[0]
            experiment = {
                'id': exp['id'],
                'mixture_id': exp['mixture_id'],
                'property_type_id': exp['property_type_id'],
                'property_name': exp['property_types']['name'],
                'numeric_value': exp['numeric_value'],
                'text_value': exp['text_value'],
                'boolean_value': exp['boolean_value'],
                'experimental_conditions': exp['experimental_conditions'],
                'date_performed': exp['date_performed'],
                'created_at': exp['created_at']
            }
            
            return _handle_json_serialization(experiment), 201
        except Exception as e:
            current_app.logger.error(f"Error recording experiment: {str(e)}", exc_info=True)
            abort(500, message=f"Error recording experiment: {str(e)}")
    
    @doc(description='Update an existing experiment.',
         tags=['Experiments'],
         params={'mixture_id': {'description': 'ID of the mixture the experiment belongs to', 'type': 'string', 'required': True},
                'experiment_id': {'description': 'ID of the experiment to update (optional)', 'type': 'string'}},
         request_body=ExperimentSchema)
    @apispec_marshal_with(experiment_fields, code=200, description='Updated experiment')
    @auth_required('bearerAuth')
    @token_required
    @marshal_with(experiment_fields)
    def put(self, mixture_id: str, experiment_id: Optional[str]=None) -> Tuple[Dict[str, Any], int]:
        """Update an existing experiment.
        
        Args:
            mixture_id (str): ID of the mixture the experiment belongs to
            experiment_id (Optional[str]): Optional ID of the experiment to update
            
        Returns:
            Tuple[Dict[str, Any], int]: Updated experiment object and HTTP status code
            
        Raises:
            400: If request data is invalid
            404: If experiment or property type not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Request Body:
            JSON object with experiment data following the ExperimentSchema:
            - property_name (str): Name of the property (required)
            - value (mixed): Property value (required)
            - experimental_conditions (str): Description of experimental conditions (optional)
            - date_performed (date): Date when experiment was performed (required)
            
        Notes:
            If experiment_id is not provided, the system will attempt to find the experiment
            using the property_name from the request body.
        """
        try:
            # Validate request data
            schema = ExperimentSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # If experiment_id is not provided, try to find by property name
            if experiment_id is None:
                # Get property type ID
                property_response = supabase.table("property_types").select("id").eq("name", data["property_name"]).execute()
                
                error_message, status_code = handle_supabase_error(property_response)
                if error_message:
                    abort(status_code, message=error_message)
                
                if not property_response.data:
                    abort(404, message=f"Property type '{data['property_name']}' not found")
                
                property_type_id = property_response.data[0]["id"]
                
                # Find the experiment
                exp_response = supabase.table("experiments").select("id").eq("mixture_id", mixture_id).eq("property_type_id", property_type_id).execute()
                
                error_message, status_code = handle_supabase_error(exp_response)
                if error_message:
                    abort(status_code, message=error_message)
                
                if not exp_response.data:
                    abort(404, message=f"Experiment for mixture {mixture_id} with property '{data['property_name']}' not found")
                
                experiment_id = exp_response.data[0]["id"]
            
            # Get property type ID and data type
            response = supabase.table("property_types").select("id, data_type").eq("name", data["property_name"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Property type '{data['property_name']}' not found")
            
            property_type = response.data[0]
            
            # Prepare experiment update
            experiment = {
                "property_type_id": property_type["id"],
                "experimental_conditions": data.get("experimental_conditions", ""),
                "date_performed": data["date_performed"].isoformat(),
                "updated_at": datetime.now().isoformat()
            }
            
            # Set the appropriate value field based on data type
            if property_type["data_type"] == "numeric":
                experiment["numeric_value"] = float(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "text":
                experiment["text_value"] = str(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "boolean":
                experiment["boolean_value"] = bool(data["value"]) if data["value"] is not None else None
            else:
                abort(400, message=f"Unknown data type '{property_type['data_type']}'")
            
            # Update experiment
            response = supabase.table("experiments").update(experiment).eq("id", experiment_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Get the updated experiment with property name
            response = supabase.from_("experiments").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value,
                experimental_conditions, date_performed, created_at,
                property_types(name)
            """).eq("id", experiment_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Experiment with ID {experiment_id} not found")
            
            # Transform the response to match the experiment_fields structure
            exp = response.data[0]
            experiment = {
                'id': exp['id'],
                'mixture_id': exp['mixture_id'],
                'property_type_id': exp['property_type_id'],
                'property_name': exp['property_types']['name'],
                'numeric_value': exp['numeric_value'],
                'text_value': exp['text_value'],
                'boolean_value': exp['boolean_value'],
                'experimental_conditions': exp['experimental_conditions'],
                'date_performed': exp['date_performed'],
                'created_at': exp['created_at']
            }
            
            current_app.logger.info(f"Experiment {experiment_id} updated successfully")
            return _handle_json_serialization(experiment), 200
        except Exception as e:
            current_app.logger.error(f"Error updating experiment: {str(e)}", exc_info=True)
            abort(500, message=f"Error updating experiment: {str(e)}")
    
    @doc(description='Delete an experiment.',
         tags=['Experiments'],
         params={'mixture_id': {'description': 'ID of the mixture the experiment belongs to', 'type': 'string', 'required': True},
                'experiment_id': {'description': 'ID of the experiment to delete (optional)', 'type': 'string'}})
    @apispec_marshal_with(fields.Raw, code=200, description='Success message')
    @auth_required('bearerAuth')
    @token_required
    def delete(self, mixture_id: str, experiment_id: Optional[str]=None) -> Tuple[Dict[str, str], int]:
        """Delete an experiment.
        
        Args:
            mixture_id (str): ID of the mixture the experiment belongs to
            experiment_id (Optional[str]): Optional ID of the experiment to delete
            
        Returns:
            Tuple[Dict[str, str], int]: Success message and HTTP status code
            
        Raises:
            400: If neither experiment_id nor property_name is provided
            404: If experiment not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Notes:
            If experiment_id is not provided, the request body must include property_name
            to identify the experiment to delete.
        """
        try:
            supabase = get_supabase_client()
            
            # If experiment_id is not provided, try to find by property name
            if experiment_id is None and request.json:
                data = request.json
                property_name = data.get("property_name")
                
                if not property_name:
                    abort(400, message="Either experiment_id or property_name must be provided")
                
                # Get property type ID
                property_response = supabase.table("property_types").select("id").eq("name", property_name).execute()
                
                error_message, status_code = handle_supabase_error(property_response)
                if error_message:
                    abort(status_code, message=error_message)
                
                if not property_response.data:
                    abort(404, message=f"Property type '{property_name}' not found")
                
                property_type_id = property_response.data[0]["id"]
                
                # Find the experiment
                exp_response = supabase.table("experiments").select("id").eq("mixture_id", mixture_id).eq("property_type_id", property_type_id).execute()
                
                error_message, status_code = handle_supabase_error(exp_response)
                if error_message:
                    abort(status_code, message=error_message)
                
                if not exp_response.data:
                    abort(404, message=f"Experiment for mixture {mixture_id} with property '{property_name}' not found")
                
                experiment_id = exp_response.data[0]["id"]
            elif experiment_id is None:
                abort(400, message="Either experiment_id or property_name must be provided")
            
            # Check if experiment exists
            check_response = supabase.table("experiments").select("id").eq("id", experiment_id).execute()
            
            error_message, status_code = handle_supabase_error(check_response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not check_response.data:
                abort(404, message=f"Experiment with ID {experiment_id} not found")
            
            # Delete the experiment
            response = supabase.table("experiments").delete().eq("id", experiment_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            current_app.logger.info(f"Experiment {experiment_id} deleted successfully")
            return {"message": f"Experiment with ID {experiment_id} deleted successfully"}, 200
        except Exception as e:
            current_app.logger.error(f"Error deleting experiment: {str(e)}", exc_info=True)
            abort(500, message=f"Error deleting experiment: {str(e)}")


class ComparisonResource(MethodResource):
    """Resource for comparing predictions with experiments.
    
    Provides an endpoint for comparing predicted property values with
    experimental measurements for a mixture, calculating differences
    and percent errors.
    """
    
    @doc(description='Compare prediction with experiment for a mixture.',
         tags=['Comparisons'],
         params={'mixture_id': {'description': 'ID of the mixture to compare data for', 'type': 'string', 'required': True},
                'property_name': {'description': 'Name of the property to compare', 'type': 'string', 'required': True}},
         security=[{'bearerAuth': []}])
    @apispec_marshal_with(comparison_fields, code=200, description='Comparison results')
    @marshal_with(comparison_fields)
    @auth_required('bearerAuth')
    @token_required
    def get(self, mixture_id: str) -> Tuple[Dict[str, Any], int]:
        """Compare prediction with experiment for a mixture.
        
        Args:
            mixture_id (str): ID of the mixture to compare data for
            
        Returns:
            Tuple[Dict[str, Any], int]: Comparison results and HTTP status code
            
        Raises:
            400: If property_name parameter is missing or invalid
            401: If authentication token is missing or invalid
            403: If user doesn't have permission to access comparison data for this mixture
            404: If property type not found
            500: If database error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Authorization:
            Users can only access comparison data for mixtures they have created or have been granted access to.
            Administrators can access all comparison data.
            
        Query Parameters:
            property_name (str): Name of the property to compare
            
        Notes:
            This endpoint compares the predicted value with the experimental value
            for the specified property, calculating the absolute difference and
            percent error.
        """
        try:
            # Validate query parameters
            schema = ComparisonQuerySchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            supabase = get_supabase_client()
            
            # Get property type ID
            response = supabase.table("property_types").select("id").eq("name", args["property_name"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Property type '{args['property_name']}' not found")
            
            property_type = response.data[0]
            
            # Call the database function to compare
            response = supabase.rpc(
                "compare_prediction_with_experiment",
                {
                    "p_mixture_id": mixture_id,
                    "p_property_type_id": property_type["id"]
                }
            ).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Handle potential JSON serialization issues
            try:
                # If response.data is a string, try to parse it as JSON
                if isinstance(response.data, str):
                    import json
                    result = json.loads(response.data)
                else:
                    result = response.data
                
                # Ensure all fields are properly formatted for serialization
                if isinstance(result, dict):
                    # Convert any None values to appropriate types based on field definitions
                    if 'difference' in result and result['difference'] is None:
                        result['difference'] = 0.0
                    if 'percent_error' in result and result['percent_error'] is None:
                        result['percent_error'] = 0.0
                
                return _handle_json_serialization(result), 200
            except Exception as e:
                current_app.logger.error(f"JSON serialization error: {str(e)}", exc_info=True)
                abort(500, message=f"Error serializing comparison result: {str(e)}")
        except Exception as e:
            current_app.logger.error(f"Error comparing prediction with experiment: {str(e)}", exc_info=True)
            abort(500, message=f"Error comparing prediction with experiment: {str(e)}")
        
        
class PropertyComparisonResource(MethodResource):
    """Resource for comparing properties of molecules and mixtures side-by-side.
    
    Provides an endpoint for comparing properties of multiple entities (molecules
    and mixtures) in a side-by-side view, highlighting similarities and differences.
    
    Endpoint:
        POST /api/compare-properties
        
    Request JSON:
        { "ids": ["mol_id1", "mix_id2", ...] }
        
    Response:
        {
            "comparison": [ ... ],
            "differences": { ... },
            "properties_compared": [ ... ]
        }
    """
    # Define response fields
    property_comparison_fields = {
        'comparison': fields.Raw,
        'differences': fields.Raw,
        'properties_compared': fields.List(fields.String)
    }
    
    # Request parser
    comparison_parser = reqparse.RequestParser()
    comparison_parser.add_argument('ids', type=list, required=True,
                                 location='json',
                                 help='List of entity IDs to compare is required')
    
    @doc(description='Compare properties of multiple entities side-by-side.',
         tags=['Comparisons'],
         request_body={'description': 'List of entity IDs to compare',
                      'content': {'application/json': {'schema': {'type': 'object', 'properties': {'ids': {'type': 'array', 'items': {'type': 'string'}}}}}}})
    @apispec_marshal_with(property_comparison_fields, code=200, description='Comparison results')
    @auth_required('bearerAuth')
    @marshal_with(property_comparison_fields)
    def post(self) -> Tuple[Dict[str, Any], int]:
        """Compare properties of multiple entities side-by-side.
        
        Args:
            Request JSON with 'ids' list of entity IDs to compare
            
        Returns:
            Tuple[Dict[str, Any], int]: Comparison data with differences and properties compared, and HTTP status code
            
        Raises:
            400: If request data is invalid
            500: If comparison error occurs
            
        Request Body:
            JSON object with:
            - ids (list): List of entity IDs (molecules or mixtures) to compare
            
        Notes:
            This endpoint compares properties across multiple entities, organizing them
            in a side-by-side view and highlighting differences. It can compare both
            molecules and mixtures in the same request.
        """
        try:
            data = request.get_json(force=True)
            schema = PropertyComparisonRequestSchema()
            try:
                validated = schema.load(data)
            except Exception as err:
                return handle_error(
                    {"error": "Validation failed", "messages": err.messages},
                    context="Property comparison",
                    log_level='warning',
                    return_response=True,
                    status_code=400
                )

            ids = validated["ids"]
            result = compare_entities(ids)
            return result, 200
        except Exception as e:
            return handle_error(
                e,
                context="Property comparison",
                log_level='error',
                return_response=True
            )