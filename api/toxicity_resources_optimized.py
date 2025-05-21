"""
CryoProtect Analyzer API - Optimized Toxicity Resources

This module provides RESTful API endpoints for accessing and analyzing toxicity data
using the optimized toxicity database schema. This implementation leverages materialized
views, efficient caching, and specialized endpoints for improved performance.

Endpoints include:
- GET /api/toxicity/molecule/{id}: Get toxicity data for a specific molecule
- GET /api/toxicity/endpoints: Get available toxicity endpoints
- GET /api/toxicity/assays: Get available toxicity assays
- GET /api/toxicity/ld50/molecule/{id}: Get LD50 data for a molecule
- GET /api/toxicity/tox21/molecule/{id}: Get Tox21 data for a molecule
- GET /api/toxicity/classifications/molecule/{id}: Get hazard classifications for a molecule
- GET /api/toxicity/summary/molecule/{id}: Get toxicity summary for a molecule
- GET /api/toxicity/scores/molecule/{id}: Get toxicity scores for a molecule
- GET /api/toxicity/scores/mixture/{id}: Get toxicity scores for a mixture
- GET /api/toxicity/unified/molecule/{id}: Get unified score for a molecule
- GET /api/toxicity/unified/mixture/{id}: Get unified score for a mixture
- GET /api/toxicity/unified/contexts: Get available application contexts for unified scoring
- GET /api/toxicity/similar/{id}: Get molecules with similar toxicity profiles
- POST /api/toxicity/bulk/molecules: Get toxicity data for multiple molecules
- POST /api/toxicity/unified/batch: Calculate unified scores for multiple molecules/mixtures

Implementation Notes:
- Uses materialized views for common queries
- Implements efficient caching with ETags
- Provides bulk endpoints for retrieving multiple records
- Leverages database functions for complex calculations
"""

import logging
import hashlib
import time
import functools
from typing import Dict, List, Optional, Union, Any
from flask import request, jsonify, make_response, Response
from flask_restful import Resource, reqparse, marshal_with, fields
from marshmallow import Schema, fields as ma_fields, validate, ValidationError

from api.models import (
    toxicity_data_fields, toxicity_assay_fields, toxicity_endpoint_fields,
    toxicity_score_fields, unified_score_fields, ToxicityData, ToxicityAssay,
    ToxicityEndpoint, ToxicityScore, ToxicityQuerySchema, UnifiedScoreQuerySchema
)
from api.api_decorators import rate_limit
from api.utils import token_required, handle_error, get_user_id
from api.unified_scoring import (
    calculate_unified_molecule_score, calculate_unified_mixture_score,
    batch_calculate_unified_scores, get_available_application_contexts
)
from chemical_data.toxicity.toxicity_scorer import ToxicityScorer

# Set up logging
logger = logging.getLogger(__name__)

# Initialize toxicity scorer
toxicity_scorer = ToxicityScorer()

# Cache control decorator
def cache_control(max_age=3600):
    """
    Decorator for adding cache control headers to responses.
    
    Args:
        max_age: Cache lifetime in seconds
    """
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            response = make_response(f(*args, **kwargs))
            response.headers['Cache-Control'] = f'public, max-age={max_age}'
            return response
        return wrapper
    return decorator

# ETag generator
def generate_etag(data):
    """
    Generate an ETag for the given data.
    
    Args:
        data: Data to generate ETag from
        
    Returns:
        str: ETag string
    """
    return hashlib.md5(str(data).encode('utf-8')).hexdigest()

# Response wrapper with ETag support
def etag_response(data, etag=None):
    """
    Create a response with ETag support.
    
    Args:
        data: Data to return
        etag: Optional ETag string
        
    Returns:
        Response: Flask response object
    """
    response = make_response(jsonify(data))
    if etag:
        response.headers['ETag'] = etag
    return response

class ToxicitySummaryResource(Resource):
    """Resource for accessing toxicity summary for a molecule."""
    
    @token_required
    @rate_limit(limit=20, period=60)
    @cache_control(max_age=3600)  # 1 hour cache
    def get(self, molecule_id: str):
        """
        Get toxicity summary for a specific molecule.
        
        Uses the optimized toxicity_summary materialized view for efficient retrieval.
        Includes ETag support for client-side caching.
        
        Args:
            molecule_id: UUID of the molecule
            
        Returns:
            JSON response with toxicity summary
        """
        try:
            # Get last update time for ETag calculation
            last_update_query = """
                SELECT MAX(updated_at) FROM toxicity_data WHERE molecule_id = %s
                UNION ALL
                SELECT MAX(updated_at) FROM toxicity_classification WHERE molecule_id = %s
            """
            
            db = get_supabase_client()
            last_update_response = db.rpc('exec_sql', {
                'query': last_update_query,
                'params': [molecule_id, molecule_id]
            }).execute()
            
            if last_update_response.error:
                logger.error(f"Error getting last update time: {last_update_response.error}")
                return handle_error(last_update_response.error)
                
            # Generate ETag
            last_update = last_update_response.data[0]['max'] if last_update_response.data else None
            etag = generate_etag(f"{molecule_id}:{last_update}")
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Query the materialized view
            query = """
                SELECT * FROM toxicity_summary WHERE molecule_id = %s
            """
            
            response = db.rpc('exec_sql', {
                'query': query,
                'params': [molecule_id]
            }).execute()
            
            if response.error:
                logger.error(f"Error getting toxicity summary: {response.error}")
                return handle_error(response.error)
                
            if not response.data:
                return {"message": "No toxicity data found for this molecule"}, 404
                
            # Return the toxicity summary
            result = {
                "molecule_id": molecule_id,
                "toxicity_summary": response.data[0],
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting toxicity summary: {str(e)}")
            return handle_error(e)

class ToxicityDataResource(Resource):
    """Resource for accessing toxicity data for a specific molecule."""
    
    @token_required
    @rate_limit(limit=20, period=60)
    @cache_control(max_age=1800)  # 30 minutes cache
    def get(self, molecule_id: str):
        """
        Get toxicity data for a specific molecule.
        
        Args:
            molecule_id: UUID of the molecule
            
        Returns:
            JSON response with toxicity data
        """
        try:
            # Get query parameters
            parser = reqparse.RequestParser()
            parser.add_argument('source', type=str, help='Filter by data source')
            parser.add_argument('toxicity_type', type=str, help='Filter by toxicity type')
            parser.add_argument('species', type=str, help='Filter by species')
            parser.add_argument('is_predicted', type=bool, help='Filter by prediction status')
            args = parser.parse_args()
            
            # Generate ETag based on query params
            etag_data = f"{molecule_id}:{args}"
            etag = generate_etag(etag_data)
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Build query
            query = "SELECT * FROM toxicity_data WHERE molecule_id = %s"
            params = [molecule_id]
            
            # Add filters
            if args.get('source'):
                query += " AND source = %s"
                params.append(args.get('source'))
                
            if args.get('toxicity_type'):
                query += " AND toxicity_type = %s"
                params.append(args.get('toxicity_type'))
                
            if args.get('species'):
                query += " AND species = %s"
                params.append(args.get('species'))
                
            if args.get('is_predicted') is not None:
                query += " AND is_predicted = %s"
                params.append(args.get('is_predicted'))
            
            # Execute query
            db = get_supabase_client()
            response = db.rpc('exec_sql', {
                'query': query,
                'params': params
            }).execute()
            
            if response.error:
                logger.error(f"Error getting toxicity data: {response.error}")
                return handle_error(response.error)
                
            # Check if data exists
            if not response.data:
                return {"message": "No toxicity data found for this molecule"}, 404
                
            # Return the toxicity data
            result = {
                "molecule_id": molecule_id,
                "toxicity_data": response.data,
                "count": len(response.data)
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting toxicity data: {str(e)}")
            return handle_error(e)

class ToxicityEndpointsResource(Resource):
    """Resource for accessing available toxicity endpoints."""
    
    @token_required
    @rate_limit(limit=10, period=60)
    @cache_control(max_age=86400)  # 24 hours cache (endpoints rarely change)
    def get(self):
        """
        Get available toxicity endpoints.
        
        Returns:
            JSON response with toxicity endpoints
        """
        try:
            endpoints = ToxicityEndpoint.get_all()
            
            # Generate ETag
            etag = generate_etag(str(endpoints))
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Return the endpoints
            result = {
                "endpoints": endpoints,
                "count": len(endpoints)
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting toxicity endpoints: {str(e)}")
            return handle_error(e)

class ToxicityAssaysResource(Resource):
    """Resource for accessing available toxicity assays."""
    
    @token_required
    @rate_limit(limit=10, period=60)
    @cache_control(max_age=86400)  # 24 hours cache (assays rarely change)
    def get(self):
        """
        Get available toxicity assays.
        
        Returns:
            JSON response with toxicity assays
        """
        try:
            # Get query parameters
            parser = reqparse.RequestParser()
            parser.add_argument('endpoint', type=str, help='Filter by toxicological endpoint')
            args = parser.parse_args()
            
            # Generate ETag based on query params
            etag_data = f"assays:{args.get('endpoint')}"
            etag = generate_etag(etag_data)
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Get assays
            assays = ToxicityAssay.get_all(endpoint=args.get('endpoint'))
            
            # Return the assays
            result = {
                "assays": assays,
                "count": len(assays)
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting toxicity assays: {str(e)}")
            return handle_error(e)

class ToxicityLD50Resource(Resource):
    """Resource for accessing LD50 toxicity data for a molecule."""
    
    @token_required
    @rate_limit(limit=20, period=60)
    @cache_control(max_age=3600)  # 1 hour cache
    def get(self, molecule_id: str):
        """
        Get LD50 toxicity data for a specific molecule.
        
        Uses the optimized ld50_summary materialized view for efficient retrieval.
        
        Args:
            molecule_id: UUID of the molecule
            
        Returns:
            JSON response with LD50 toxicity data
        """
        try:
            # Get query parameters
            parser = reqparse.RequestParser()
            parser.add_argument('species', type=str, help='Filter by species')
            parser.add_argument('route', type=str, help='Filter by route of administration')
            parser.add_argument('is_predicted', type=bool, help='Filter by prediction status')
            args = parser.parse_args()
            
            # Generate ETag based on query params
            etag_data = f"{molecule_id}:ld50:{args}"
            etag = generate_etag(etag_data)
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Build query
            query = "SELECT * FROM ld50_summary WHERE molecule_id = %s"
            params = [molecule_id]
            
            # Add filters
            if args.get('species'):
                query += " AND species = %s"
                params.append(args.get('species'))
                
            if args.get('route'):
                query += " AND route_of_administration = %s"
                params.append(args.get('route'))
                
            if args.get('is_predicted') is not None:
                query += " AND is_predicted = %s"
                params.append(args.get('is_predicted'))
                
            # Add ordering
            query += " ORDER BY ld50_value ASC"
            
            # Execute query
            db = get_supabase_client()
            response = db.rpc('exec_sql', {
                'query': query,
                'params': params
            }).execute()
            
            if response.error:
                logger.error(f"Error getting LD50 data: {response.error}")
                return handle_error(response.error)
                
            # Check if data exists
            if not response.data:
                return {"message": "No LD50 data found for this molecule"}, 404
                
            # Return the LD50 data
            result = {
                "molecule_id": molecule_id,
                "ld50_data": response.data,
                "count": len(response.data)
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting LD50 data: {str(e)}")
            return handle_error(e)

class ToxicityTox21Resource(Resource):
    """Resource for accessing Tox21 assay data for a molecule."""
    
    @token_required
    @rate_limit(limit=20, period=60)
    @cache_control(max_age=3600)  # 1 hour cache
    def get(self, molecule_id: str):
        """
        Get Tox21 assay data for a specific molecule.
        
        Uses the optimized tox21_activity_summary materialized view for efficient retrieval.
        
        Args:
            molecule_id: UUID of the molecule
            
        Returns:
            JSON response with Tox21 assay data
        """
        try:
            # Get query parameters
            parser = reqparse.RequestParser()
            parser.add_argument('assay_id', type=str, help='Filter by assay ID')
            parser.add_argument('target_family', type=str, help='Filter by target family')
            parser.add_argument('activity', type=str, help='Filter by activity outcome (Active/Inactive)')
            args = parser.parse_args()
            
            # Generate ETag based on query params
            etag_data = f"{molecule_id}:tox21:{args}"
            etag = generate_etag(etag_data)
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Build query
            query = "SELECT * FROM tox21_activity_summary WHERE molecule_id = %s"
            params = [molecule_id]
            
            # Add filters
            if args.get('assay_id'):
                query += " AND assay_id = %s"
                params.append(args.get('assay_id'))
                
            if args.get('target_family'):
                query += " AND intended_target_family = %s"
                params.append(args.get('target_family'))
                
            if args.get('activity'):
                query += " AND activity_outcome = %s"
                params.append(args.get('activity'))
                
            # Execute query
            db = get_supabase_client()
            response = db.rpc('exec_sql', {
                'query': query,
                'params': params
            }).execute()
            
            if response.error:
                logger.error(f"Error getting Tox21 data: {response.error}")
                return handle_error(response.error)
                
            # Check if data exists
            if not response.data:
                return {"message": "No Tox21 data found for this molecule"}, 404
                
            # Return the Tox21 data
            result = {
                "molecule_id": molecule_id,
                "tox21_data": response.data,
                "count": len(response.data),
                "active_count": sum(1 for item in response.data if item.get('activity_outcome') == 'Active')
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting Tox21 data: {str(e)}")
            return handle_error(e)

class ToxicityClassificationResource(Resource):
    """Resource for accessing hazard classifications for a molecule."""
    
    @token_required
    @rate_limit(limit=20, period=60)
    @cache_control(max_age=3600)  # 1 hour cache
    def get(self, molecule_id: str):
        """
        Get hazard classifications for a specific molecule.
        
        Uses the optimized hazard_classification_summary materialized view for efficient retrieval.
        
        Args:
            molecule_id: UUID of the molecule
            
        Returns:
            JSON response with hazard classifications
        """
        try:
            # Get query parameters
            parser = reqparse.RequestParser()
            parser.add_argument('system', type=str, help='Filter by classification system')
            parser.add_argument('hazard_class', type=str, help='Filter by hazard class')
            parser.add_argument('is_predicted', type=bool, help='Filter by prediction status')
            args = parser.parse_args()
            
            # Generate ETag based on query params
            etag_data = f"{molecule_id}:classification:{args}"
            etag = generate_etag(etag_data)
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Build query
            query = "SELECT * FROM hazard_classification_summary WHERE molecule_id = %s"
            params = [molecule_id]
            
            # Add filters
            if args.get('system'):
                query += " AND classification_system = %s"
                params.append(args.get('system'))
                
            if args.get('hazard_class'):
                query += " AND hazard_class = %s"
                params.append(args.get('hazard_class'))
                
            if args.get('is_predicted') is not None:
                query += " AND is_predicted = %s"
                params.append(args.get('is_predicted'))
                
            # Execute query
            db = get_supabase_client()
            response = db.rpc('exec_sql', {
                'query': query,
                'params': params
            }).execute()
            
            if response.error:
                logger.error(f"Error getting hazard classifications: {response.error}")
                return handle_error(response.error)
                
            # Check if data exists
            if not response.data:
                return {"message": "No hazard classifications found for this molecule"}, 404
                
            # Return the hazard classifications
            result = {
                "molecule_id": molecule_id,
                "classifications": response.data,
                "count": len(response.data)
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting hazard classifications: {str(e)}")
            return handle_error(e)

class ToxicityScoreResource(Resource):
    """Resource for accessing toxicity scores for a molecule."""
    
    @token_required
    @rate_limit(limit=20, period=60)
    @cache_control(max_age=7200)  # 2 hour cache
    def get(self, molecule_id: str):
        """
        Get toxicity scores for a specific molecule.
        
        Args:
            molecule_id: UUID of the molecule
            
        Returns:
            JSON response with toxicity scores
        """
        try:
            # Generate ETag
            etag = generate_etag(f"{molecule_id}:score")
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Use database function for efficient calculation
            db = get_supabase_client()
            response = db.rpc('calculate_toxicity_score', {
                'molecule_uuid': molecule_id
            }).execute()
            
            if response.error:
                logger.error(f"Error calculating toxicity score: {response.error}")
                return handle_error(response.error)
            
            # Get detailed toxicity data
            scores = ToxicityScore.get_for_molecule(molecule_id)
            
            if not scores and not response.data:
                # Try to calculate scores if they don't exist
                calculation_result = toxicity_scorer._calculate_molecule_toxicity_score(molecule_id)
                if calculation_result.get("success"):
                    scores = ToxicityScore.get_for_molecule(molecule_id)
            
            # Return the toxicity score
            result = {
                "molecule_id": molecule_id,
                "toxicity_score": response.data,
                "scores": scores,
                "score_count": len(scores) if scores else 0
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting toxicity score: {str(e)}")
            return handle_error(e)

class MixtureToxicityScoreResource(Resource):
    """Resource for accessing toxicity scores for a mixture."""
    
    @token_required
    @rate_limit(limit=20, period=60)
    @cache_control(max_age=7200)  # 2 hour cache
    def get(self, mixture_id: str):
        """
        Get toxicity scores for a specific mixture.
        
        Args:
            mixture_id: UUID of the mixture
            
        Returns:
            JSON response with toxicity scores
        """
        try:
            # Generate ETag
            etag = generate_etag(f"{mixture_id}:mixture_score")
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Get toxicity scores for mixture components
            component_scores = ToxicityScore.get_for_mixture(mixture_id)
            
            if not component_scores:
                return {"message": "No toxicity scores found for this mixture's components"}, 404
                
            # Calculate aggregate score for the mixture
            aggregate_score = ToxicityScore.calculate_mixture_toxicity_score(mixture_id)
                
            # Return the mixture toxicity score
            result = {
                "mixture_id": mixture_id,
                "aggregate_score": aggregate_score,
                "component_scores": component_scores,
                "component_count": len(component_scores)
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting mixture toxicity scores: {str(e)}")
            return handle_error(e)

class SimilarToxicityResource(Resource):
    """Resource for finding molecules with similar toxicity profiles."""
    
    @token_required
    @rate_limit(limit=10, period=60)
    def get(self, molecule_id: str):
        """
        Get molecules with similar toxicity profiles.
        
        Uses the database function find_similar_toxicity_profiles for efficient retrieval.
        
        Args:
            molecule_id: UUID of the reference molecule
            
        Returns:
            JSON response with similar molecules
        """
        try:
            # Get query parameters
            parser = reqparse.RequestParser()
            parser.add_argument('limit', type=int, default=10, help='Maximum number of results to return')
            args = parser.parse_args()
            
            # Use database function for finding similar molecules
            db = get_supabase_client()
            response = db.rpc('find_similar_toxicity_profiles', {
                'molecule_uuid': molecule_id,
                'limit_count': args.get('limit', 10)
            }).execute()
            
            if response.error:
                logger.error(f"Error finding similar molecules: {response.error}")
                return handle_error(response.error)
                
            # Check if data exists
            if not response.data:
                return {"message": "No similar molecules found"}, 404
                
            # Return the similar molecules
            result = {
                "reference_molecule_id": molecule_id,
                "similar_molecules": response.data,
                "count": len(response.data)
            }
            
            return jsonify(result), 200
            
        except Exception as e:
            logger.error(f"Error finding similar molecules: {str(e)}")
            return handle_error(e)

class BulkToxicityResource(Resource):
    """Resource for retrieving toxicity data for multiple molecules."""
    
    @token_required
    @rate_limit(limit=5, period=60)  # Lower rate limit due to higher resource usage
    def post(self):
        """
        Get toxicity data for multiple molecules.
        
        Request Body:
            molecule_ids: List of molecule IDs
            data_type: Type of toxicity data to retrieve
            
        Returns:
            JSON response with toxicity data for each molecule
        """
        try:
            # Parse request body
            parser = reqparse.RequestParser()
            parser.add_argument('molecule_ids', type=list, required=True, location='json',
                              help='List of molecule IDs is required')
            parser.add_argument('data_type', type=str, default='summary', location='json',
                              choices=('summary', 'ld50', 'tox21', 'classification', 'score'),
                              help='Type of toxicity data to retrieve')
            args = parser.parse_args()
            
            molecule_ids = args['molecule_ids']
            data_type = args['data_type']
            
            if not molecule_ids:
                return {"error": "No molecule IDs provided"}, 400
                
            if len(molecule_ids) > 100:
                return {"error": "Maximum of 100 molecule IDs allowed"}, 400
            
            # Determine which materialized view to query
            view_map = {
                'summary': 'toxicity_summary',
                'ld50': 'ld50_summary',
                'tox21': 'tox21_activity_summary',
                'classification': 'hazard_classification_summary',
                'score': None  # Special case, uses function
            }
            
            view_name = view_map.get(data_type)
            
            # Special case for score
            if data_type == 'score':
                # Calculate scores for all molecules
                results = {}
                db = get_supabase_client()
                
                for mol_id in molecule_ids:
                    try:
                        score_response = db.rpc('calculate_toxicity_score', {
                            'molecule_uuid': mol_id
                        }).execute()
                        
                        if not score_response.error:
                            results[mol_id] = score_response.data
                    except Exception as e:
                        logger.warning(f"Error calculating score for molecule {mol_id}: {str(e)}")
                        # Skip this molecule and continue with others
                
                return {
                    "molecule_count": len(molecule_ids),
                    "results_count": len(results),
                    "results": results
                }, 200
            
            # For other data types
            if view_name:
                # Build query to get data for all molecules at once
                placeholders = ', '.join(['%s'] * len(molecule_ids))
                query = f"SELECT * FROM {view_name} WHERE molecule_id IN ({placeholders})"
                
                # Execute query
                db = get_supabase_client()
                response = db.rpc('exec_sql', {
                    'query': query,
                    'params': molecule_ids
                }).execute()
                
                if response.error:
                    logger.error(f"Error getting bulk toxicity data: {response.error}")
                    return handle_error(response.error)
                
                # Group results by molecule_id
                results = {}
                for item in response.data:
                    mol_id = item['molecule_id']
                    if mol_id not in results:
                        results[mol_id] = []
                    results[mol_id].append(item)
                
                return {
                    "molecule_count": len(molecule_ids),
                    "results_count": len(results),
                    "results": results
                }, 200
            
            return {"error": "Invalid data type specified"}, 400
            
        except Exception as e:
            logger.error(f"Error getting bulk toxicity data: {str(e)}")
            return handle_error(e)

class UnifiedScoreMoleculeResource(Resource):
    """Resource for accessing unified scores for a molecule."""
    
    @token_required
    @rate_limit(limit=10, period=60)
    @cache_control(max_age=3600)  # 1 hour cache
    def get(self, molecule_id: str):
        """
        Get unified score for a specific molecule.
        
        Args:
            molecule_id: UUID of the molecule
            
        Returns:
            JSON response with unified score
        """
        try:
            # Get query parameters
            args = request.args
            application_context = args.get('context', 'general')
            algorithm = args.get('algorithm', 'random_forest')
            recalculate = args.get('recalculate', 'false').lower() == 'true'
            
            # Generate ETag based on query params
            etag_data = f"{molecule_id}:unified:{application_context}:{algorithm}:{recalculate}"
            etag = generate_etag(etag_data)
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag and not recalculate:
                return '', 304  # Not Modified
            
            # Calculate unified score
            result = calculate_unified_molecule_score(
                molecule_id=molecule_id,
                application_context=application_context,
                algorithm=algorithm,
                recalculate=recalculate
            )
            
            if "error" in result:
                return {"message": result["error"]}, 404
            
            # Add ETag header
            response = make_response(jsonify(result))
            response.headers['ETag'] = etag
            return response
            
        except Exception as e:
            logger.error(f"Error calculating unified score: {str(e)}")
            return handle_error(e)

class UnifiedScoreMixtureResource(Resource):
    """Resource for accessing unified scores for a mixture."""
    
    @token_required
    @rate_limit(limit=10, period=60)
    @cache_control(max_age=3600)  # 1 hour cache
    def get(self, mixture_id: str):
        """
        Get unified score for a specific mixture.
        
        Args:
            mixture_id: UUID of the mixture
            
        Returns:
            JSON response with unified score
        """
        try:
            # Get query parameters
            args = request.args
            application_context = args.get('context', 'general')
            algorithm = args.get('algorithm', 'random_forest')
            recalculate = args.get('recalculate', 'false').lower() == 'true'
            
            # Generate ETag based on query params
            etag_data = f"{mixture_id}:unified:{application_context}:{algorithm}:{recalculate}"
            etag = generate_etag(etag_data)
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag and not recalculate:
                return '', 304  # Not Modified
            
            # Calculate unified score
            result = calculate_unified_mixture_score(
                mixture_id=mixture_id,
                application_context=application_context,
                algorithm=algorithm,
                recalculate=recalculate
            )
            
            if "error" in result:
                return {"message": result["error"]}, 404
            
            # Add ETag header
            response = make_response(jsonify(result))
            response.headers['ETag'] = etag
            return response
            
        except Exception as e:
            logger.error(f"Error calculating unified score: {str(e)}")
            return handle_error(e)

class UnifiedScoreContextsResource(Resource):
    """Resource for accessing available application contexts for unified scoring."""
    
    @token_required
    @rate_limit(limit=5, period=60)
    @cache_control(max_age=86400)  # 24 hours cache (contexts rarely change)
    def get(self):
        """
        Get available application contexts for unified scoring.
        
        Returns:
            JSON response with application contexts
        """
        try:
            contexts = get_available_application_contexts()
            
            # Generate ETag
            etag = generate_etag(str(contexts))
            
            # Check If-None-Match header
            if request.headers.get('If-None-Match') == etag:
                return '', 304  # Not Modified
            
            # Return the contexts
            result = {
                "contexts": contexts
            }
            
            # Add ETag header
            return etag_response(result, etag)
            
        except Exception as e:
            logger.error(f"Error getting application contexts: {str(e)}")
            return handle_error(e)

class BatchUnifiedScoreResource(Resource):
    """Resource for calculating unified scores for multiple molecules/mixtures."""
    
    @token_required
    @rate_limit(limit=5, period=60)
    def post(self):
        """
        Calculate unified scores for multiple molecules/mixtures.
        
        Request Body:
            entity_ids: List of molecule or mixture IDs
            entity_type: Type of entities ("molecule" or "mixture")
            application_context: Application context for weighting
            algorithm: Algorithm to use for predictive models
            recalculate: Whether to recalculate scores even if they already exist
            
        Returns:
            JSON response with results for each entity
        """
        try:
            # Parse request body
            parser = reqparse.RequestParser()
            parser.add_argument('entity_ids', type=list, required=True, location='json',
                               help='List of molecule or mixture IDs is required')
            parser.add_argument('entity_type', type=str, required=True, location='json',
                               choices=('molecule', 'mixture'),
                               help='Entity type must be "molecule" or "mixture"')
            parser.add_argument('application_context', type=str, default='general', location='json',
                               help='Application context for weighting')
            parser.add_argument('algorithm', type=str, default='random_forest', location='json',
                               help='Algorithm to use for predictive models')
            parser.add_argument('recalculate', type=bool, default=False, location='json',
                               help='Whether to recalculate scores even if they already exist')
            args = parser.parse_args()
            
            # Check entity_ids count
            if len(args['entity_ids']) > 50:
                return {"error": "Maximum of 50 entity IDs allowed"}, 400
            
            # Calculate batch scores
            result = batch_calculate_unified_scores(
                entity_ids=args['entity_ids'],
                entity_type=args['entity_type'],
                application_context=args['application_context'],
                algorithm=args['algorithm'],
                recalculate=args['recalculate']
            )
            
            return result, 200
            
        except Exception as e:
            logger.error(f"Error calculating batch unified scores: {str(e)}")
            return handle_error(e)

def register_toxicity_resources(api):
    """Register toxicity resources with the API."""
    # Summary endpoint
    api.add_resource(ToxicitySummaryResource, '/api/toxicity/summary/molecule/<string:molecule_id>')
    
    # Basic data endpoints
    api.add_resource(ToxicityDataResource, '/api/toxicity/molecule/<string:molecule_id>')
    api.add_resource(ToxicityEndpointsResource, '/api/toxicity/endpoints')
    api.add_resource(ToxicityAssaysResource, '/api/toxicity/assays')
    
    # Specialized data endpoints
    api.add_resource(ToxicityLD50Resource, '/api/toxicity/ld50/molecule/<string:molecule_id>')
    api.add_resource(ToxicityTox21Resource, '/api/toxicity/tox21/molecule/<string:molecule_id>')
    api.add_resource(ToxicityClassificationResource, '/api/toxicity/classifications/molecule/<string:molecule_id>')
    
    # Score endpoints
    api.add_resource(ToxicityScoreResource, '/api/toxicity/scores/molecule/<string:molecule_id>')
    api.add_resource(MixtureToxicityScoreResource, '/api/toxicity/scores/mixture/<string:mixture_id>')
    
    # Unified score endpoints
    api.add_resource(UnifiedScoreMoleculeResource, '/api/toxicity/unified/molecule/<string:molecule_id>')
    api.add_resource(UnifiedScoreMixtureResource, '/api/toxicity/unified/mixture/<string:mixture_id>')
    api.add_resource(UnifiedScoreContextsResource, '/api/toxicity/unified/contexts')
    
    # Similar compounds endpoint
    api.add_resource(SimilarToxicityResource, '/api/toxicity/similar/<string:molecule_id>')
    
    # Bulk endpoints
    api.add_resource(BulkToxicityResource, '/api/toxicity/bulk/molecules')
    api.add_resource(BatchUnifiedScoreResource, '/api/toxicity/unified/batch')