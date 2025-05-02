"""
CryoProtect Analyzer API - Toxicity Resources

This module provides RESTful API endpoints for accessing and analyzing toxicity data
from the Tox21 database, as well as unified scoring that combines toxicity with
efficacy and glass transition temperature (Tg) data.

Endpoints include:
- GET /api/toxicity/molecule/{id}: Get toxicity data for a specific molecule
- GET /api/toxicity/endpoints: Get available toxicity endpoints
- GET /api/toxicity/assays: Get available toxicity assays
- GET /api/toxicity/scores/molecule/{id}: Get toxicity scores for a molecule
- GET /api/toxicity/scores/mixture/{id}: Get toxicity scores for a mixture
- GET /api/toxicity/unified/molecule/{id}: Get unified score for a molecule
- GET /api/toxicity/unified/mixture/{id}: Get unified score for a mixture
- GET /api/toxicity/unified/contexts: Get available application contexts for unified scoring
- POST /api/toxicity/unified/batch: Calculate unified scores for multiple molecules/mixtures

Scientific Context:
    Toxicity assessment is a critical component of cryoprotectant evaluation.
    This module leverages Tox21 data, which includes high-throughput screening
    results across multiple toxicological endpoints. The unified scoring system
    balances efficacy, toxicity, and physical properties to provide a comprehensive
    assessment tailored to different application contexts.

References:
    - Huang, R., et al. (2016). The Tox21 10K Compound Library: Collaborative Chemistry
      Advancing Toxicology. Chemical Research in Toxicology, 29(8), 1225-1233.
    - Best, B. P. (2015). Cryoprotectant toxicity: facts, issues, and questions.
      Rejuvenation research, 18(5), 422-436.
"""

import logging
from typing import Dict, List, Optional, Union, Any
from flask import request, jsonify
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

class ToxicityDataResource(Resource):
    """Resource for accessing toxicity data for a specific molecule."""
    
    @token_required
    @rate_limit(limit=20, period=60)
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
            parser.add_argument('assay_id', type=str, help='Filter by assay ID')
            parser.add_argument('endpoint', type=str, help='Filter by toxicological endpoint')
            parser.add_argument('active_only', type=bool, default=False, help='Return only active results')
            args = parser.parse_args()
            
            # Get toxicity data
            toxicity_data = ToxicityData.get_for_molecule(
                molecule_id=molecule_id,
                assay_id=args.get('assay_id'),
                endpoint=args.get('endpoint'),
                active_only=args.get('active_only', False)
            )
            
            if not toxicity_data:
                return {"message": "No toxicity data found for this molecule"}, 404
                
            return {
                "molecule_id": molecule_id,
                "toxicity_data": toxicity_data,
                "count": len(toxicity_data)
            }, 200
            
        except Exception as e:
            logger.error(f"Error getting toxicity data: {str(e)}")
            return handle_error(e)

class ToxicityEndpointsResource(Resource):
    """Resource for accessing available toxicity endpoints."""
    
    @token_required
    @rate_limit(limit=10, period=60)
    def get(self):
        """
        Get available toxicity endpoints.
        
        Returns:
            JSON response with toxicity endpoints
        """
        try:
            endpoints = ToxicityEndpoint.get_all()
            
            return {
                "endpoints": endpoints,
                "count": len(endpoints)
            }, 200
            
        except Exception as e:
            logger.error(f"Error getting toxicity endpoints: {str(e)}")
            return handle_error(e)

class ToxicityAssaysResource(Resource):
    """Resource for accessing available toxicity assays."""
    
    @token_required
    @rate_limit(limit=10, period=60)
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
            
            # Get assays
            assays = ToxicityAssay.get_all(endpoint=args.get('endpoint'))
            
            return {
                "assays": assays,
                "count": len(assays)
            }, 200
            
        except Exception as e:
            logger.error(f"Error getting toxicity assays: {str(e)}")
            return handle_error(e)

class ToxicityScoreResource(Resource):
    """Resource for accessing toxicity scores for a molecule."""
    
    @token_required
    @rate_limit(limit=20, period=60)
    def get(self, molecule_id: str):
        """
        Get toxicity scores for a specific molecule.
        
        Args:
            molecule_id: UUID of the molecule
            
        Returns:
            JSON response with toxicity scores
        """
        try:
            # Get toxicity scores
            scores = ToxicityScore.get_for_molecule(molecule_id)
            
            if not scores:
                # Try to calculate scores if they don't exist
                calculation_result = toxicity_scorer._calculate_molecule_toxicity_score(molecule_id)
                if calculation_result.get("success"):
                    scores = ToxicityScore.get_for_molecule(molecule_id)
                
            if not scores:
                return {"message": "No toxicity scores found for this molecule"}, 404
                
            # Get endpoint details for context
            endpoint_ids = [score.get("endpoint_id") for score in scores if score.get("endpoint_id")]
            endpoints = {}
            if endpoint_ids:
                endpoint_data = ToxicityEndpoint.get_by_ids(endpoint_ids)
                endpoints = {endpoint["id"]: endpoint for endpoint in endpoint_data}
            
            # Enhance scores with endpoint details
            for score in scores:
                if score.get("endpoint_id") and score["endpoint_id"] in endpoints:
                    score["endpoint"] = endpoints[score["endpoint_id"]]
                
            return {
                "molecule_id": molecule_id,
                "scores": scores,
                "count": len(scores)
            }, 200
            
        except Exception as e:
            logger.error(f"Error getting toxicity scores: {str(e)}")
            return handle_error(e)

class MixtureToxicityScoreResource(Resource):
    """Resource for accessing toxicity scores for a mixture."""
    
    @token_required
    @rate_limit(limit=20, period=60)
    def get(self, mixture_id: str):
        """
        Get toxicity scores for a specific mixture.
        
        Args:
            mixture_id: UUID of the mixture
            
        Returns:
            JSON response with toxicity scores
        """
        try:
            # Get toxicity scores for mixture components
            component_scores = ToxicityScore.get_for_mixture(mixture_id)
            
            if not component_scores:
                return {"message": "No toxicity scores found for this mixture's components"}, 404
                
            # Calculate aggregate score for the mixture
            aggregate_score = ToxicityScore.calculate_mixture_toxicity_score(mixture_id)
                
            return {
                "mixture_id": mixture_id,
                "aggregate_score": aggregate_score,
                "component_scores": component_scores,
                "component_count": len(component_scores)
            }, 200
            
        except Exception as e:
            logger.error(f"Error getting mixture toxicity scores: {str(e)}")
            return handle_error(e)

class UnifiedScoreMoleculeResource(Resource):
    """Resource for accessing unified scores for a molecule."""
    
    @token_required
    @rate_limit(limit=10, period=60)
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
            
            # Calculate unified score
            result = calculate_unified_molecule_score(
                molecule_id=molecule_id,
                application_context=application_context,
                algorithm=algorithm,
                recalculate=recalculate
            )
            
            if "error" in result:
                return {"message": result["error"]}, 404
                
            return result, 200
            
        except Exception as e:
            logger.error(f"Error calculating unified score: {str(e)}")
            return handle_error(e)

class UnifiedScoreMixtureResource(Resource):
    """Resource for accessing unified scores for a mixture."""
    
    @token_required
    @rate_limit(limit=10, period=60)
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
            
            # Calculate unified score
            result = calculate_unified_mixture_score(
                mixture_id=mixture_id,
                application_context=application_context,
                algorithm=algorithm,
                recalculate=recalculate
            )
            
            if "error" in result:
                return {"message": result["error"]}, 404
                
            return result, 200
            
        except Exception as e:
            logger.error(f"Error calculating unified score: {str(e)}")
            return handle_error(e)

class UnifiedScoreContextsResource(Resource):
    """Resource for accessing available application contexts for unified scoring."""
    
    @token_required
    @rate_limit(limit=5, period=60)
    def get(self):
        """
        Get available application contexts for unified scoring.
        
        Returns:
            JSON response with application contexts
        """
        try:
            contexts = get_available_application_contexts()
            
            return {
                "contexts": contexts
            }, 200
            
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
    api.add_resource(ToxicityDataResource, '/api/toxicity/molecule/<string:molecule_id>')
    api.add_resource(ToxicityEndpointsResource, '/api/toxicity/endpoints')
    api.add_resource(ToxicityAssaysResource, '/api/toxicity/assays')
    api.add_resource(ToxicityScoreResource, '/api/toxicity/scores/molecule/<string:molecule_id>')
    api.add_resource(MixtureToxicityScoreResource, '/api/toxicity/scores/mixture/<string:mixture_id>')
    api.add_resource(UnifiedScoreMoleculeResource, '/api/toxicity/unified/molecule/<string:molecule_id>')
    api.add_resource(UnifiedScoreMixtureResource, '/api/toxicity/unified/mixture/<string:mixture_id>')
    api.add_resource(UnifiedScoreContextsResource, '/api/toxicity/unified/contexts')
    api.add_resource(BatchUnifiedScoreResource, '/api/toxicity/unified/batch')