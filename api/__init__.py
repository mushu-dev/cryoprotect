"""
CryoProtect Analyzer API Package

This package contains the API resources, models, and utilities for the CryoProtect Analyzer.

The API follows standardized patterns for:
- Response formats: All endpoints return responses in a consistent format
- Error handling: Errors are reported in a standardized way
- HTTP status codes: Proper status codes are used for all responses
- Documentation: All endpoints are documented in OpenAPI format
"""

import logging
from flask import Blueprint, jsonify, request, current_app
from flask_restful import Api

# Set up logging
logger = logging.getLogger(__name__)

# Create a blueprint for the API
api_bp = Blueprint('api', __name__, url_prefix='/api/v1')
api = Api(api_bp)

# Import API standardization utilities
from api.api_standards import (
    create_standard_response,
    create_error_response,
    create_success_response,
    jsonify_standard_response,
    HTTP_STATUS_CODES
)

# Import API decorators
from api.api_decorators import (
    standardize_response,
    validate_request_schema,
    rate_limit,
    csrf_protect,
    csrf_exempt
)

# Import API documentation utilities
from api.api_docs import (
    init_docs,
    register_resource,
    document_endpoint,
    auth_required,
    paginated_response,
    generate_openapi_spec
)

# Import resources to register them with the API
from api.resources import (
    MoleculeResource, MoleculeListResource,
    MixtureResource, MixtureListResource,
    PredictionResource, ExperimentResource,
    ComparisonResource,
    PropertyComparisonResource
)
from api.lab_verification_resources import LabVerificationResource, VerificationStatsResource
# Import Experiment resources
from api.experiment_resources import (
    ExperimentList, ExperimentDetail, ExperimentResultList,
    ExperimentResultDetail, ExperimentTimeSeriesList,
    ExperimentAnalysis, ExperimentExport, ExperimentImport,
    ExperimentSearch, experiment_bp
)
# Import Batch Operations resource
from api.batch_resources import BatchOperationResource

# Import user profile resource
from api.user_profile_resources import UserProfileResource

# Import RDKit resources
from api.rdkit_resources import (
    MoleculePropertyResource, MoleculeVisualizationResource,
    SubstructureSearchResource, SimilaritySearchResource,
    MoleculePropertyCalculationResource
)

# Import Enhanced RDKit resources
from api.rdkit_enhanced_resources import register_resources as register_rdkit_enhanced_resources

# Import Scoring resources
from api.scoring_resources import (
    MoleculeScoreResource, MoleculeIdScoreResource,
    MixtureScoreResource
)

# Import Mixture Analysis resources
from api.mixture_analysis_resources import register_resources

# Import Protocol Designer resources
from api.protocol_designer_resources import register_resources as register_protocol_resources

# Import Predictive Models resources
from api.predictive_models_resources import register_resources as register_predictive_models_resources

# Import Dashboard resources
from api.dashboard_resources import register_resources as register_dashboard_resources

# Import Team resources
from api.team_resources import register_resources as register_team_resources

# Import Export and Sharing resources
from api.export_api_resources import register_resources as register_export_resources

# Import System resources
from api.system_resources import register_resources as register_system_resources

# Import Toxicity resources
from api.toxicity_resources import register_toxicity_resources

# Register resources with the API
api.add_resource(MoleculeListResource, '/molecules')
api.add_resource(MoleculeResource, '/molecules/<string:molecule_id>')
api.add_resource(MixtureListResource, '/mixtures')
api.add_resource(BatchOperationResource, '/batch')
api.add_resource(MixtureResource, '/mixtures/<string:mixture_id>')
api.add_resource(PredictionResource,
                '/mixtures/<string:mixture_id>/predictions',
                '/mixtures/<string:mixture_id>/predictions/<string:prediction_id>')
api.add_resource(ExperimentResource,
                '/mixtures/<string:mixture_id>/experiments',
                '/mixtures/<string:mixture_id>/experiments/<string:experiment_id>')
# Register ComparisonResource with a single endpoint to avoid duplication
api.add_resource(
    ComparisonResource,
    '/mixtures/<string:mixture_id>/compare',
    endpoint='comparison_resource'
)
api.add_resource(PropertyComparisonResource, '/compare-properties')
 
# Register user profile resource
api.add_resource(UserProfileResource, '/user_profile')

# Register RDKit resources
api.add_resource(MoleculePropertyResource, '/rdkit/properties')
api.add_resource(MoleculeVisualizationResource, '/rdkit/visualization')
api.add_resource(SubstructureSearchResource, '/rdkit/substructure')
api.add_resource(SimilaritySearchResource, '/rdkit/similarity')
api.add_resource(MoleculePropertyCalculationResource, '/molecules/<string:molecule_id>/calculate-properties')

# Register Scoring resources
api.add_resource(MoleculeScoreResource, '/scoring/molecules')
api.add_resource(MoleculeIdScoreResource, '/molecules/<string:molecule_id>/score')
api.add_resource(MixtureScoreResource, '/mixtures/<string:mixture_id>/score')

# Register Mixture Analysis resources
register_resources(api)

# Register Protocol Designer resources
register_protocol_resources(api)

# Register Predictive Models resources
register_predictive_models_resources(api)

# Register Dashboard resources
register_dashboard_resources(api)

# Register Team resources
register_team_resources(api)

# Register Export and Sharing resources
register_export_resources(api)

# Register Enhanced RDKit resources
register_rdkit_enhanced_resources(api)

# Register System resources
register_system_resources(api)

# Register Toxicity resources
register_toxicity_resources(api)

# Register Lab Verification resources
api.add_resource(LabVerificationResource,
                '/experiments/<string:experiment_id>/verification',
                endpoint='experiment_verification')
api.add_resource(LabVerificationResource,
                '/verifications/<string:verification_id>',
                endpoint='verification_update')
api.add_resource(VerificationStatsResource,
                '/verification/stats',
                endpoint='verification_stats')

# Import RBAC routes
from api.rbac_routes import rbac_bp

# Add error handlers to the API blueprint
@api_bp.errorhandler(400)
def handle_bad_request(error):
    """Handle bad request errors."""
    return jsonify_standard_response(
        *create_error_response(
            error=str(error),
            status_code=400,
            context="Bad request"
        )
    )

@api_bp.errorhandler(401)
def handle_unauthorized(error):
    """Handle unauthorized errors."""
    return jsonify_standard_response(
        *create_error_response(
            error=str(error),
            status_code=401,
            context="Unauthorized"
        )
    )

@api_bp.errorhandler(403)
def handle_forbidden(error):
    """Handle forbidden errors."""
    return jsonify_standard_response(
        *create_error_response(
            error=str(error),
            status_code=403,
            context="Forbidden"
        )
    )

@api_bp.errorhandler(404)
def handle_not_found(error):
    """Handle not found errors."""
    return jsonify_standard_response(
        *create_error_response(
            error=str(error),
            status_code=404,
            context="Not found"
        )
    )

@api_bp.errorhandler(429)
def handle_too_many_requests(error):
    """Handle too many requests errors."""
    return jsonify_standard_response(
        *create_error_response(
            error=str(error),
            status_code=429,
            context="Too many requests"
        )
    )

@api_bp.errorhandler(500)
def handle_server_error(error):
    """Handle server errors."""
    logger.error(f"Server error: {str(error)}", exc_info=True)
    return jsonify_standard_response(
        *create_error_response(
            error=str(error),
            status_code=500,
            context="Server error"
        )
    )

def init_app(app):
    """Initialize the API with the Flask app."""
    # Register blueprints
    app.register_blueprint(api_bp)
    app.register_blueprint(rbac_bp, url_prefix='/api/v1/rbac')
    app.register_blueprint(experiment_bp, url_prefix='/api/v1/experiments')
    
    # Initialize API documentation
    docs = init_docs(app)

    # Register all main API resources with documentation
    register_resource(docs, MoleculeListResource, 'moleculelistresource')
    register_resource(docs, MoleculeResource, 'moleculeresource')
    register_resource(docs, MixtureListResource, 'mixturelistresource')
    register_resource(docs, MixtureResource, 'mixtureresource')
    register_resource(docs, BatchOperationResource, 'batchoperationresource')
    register_resource(docs, PredictionResource, 'predictionresource')
    register_resource(docs, ExperimentResource, 'experimentresource')
    register_resource(docs, ComparisonResource, 'comparison_resource')
    register_resource(docs, PropertyComparisonResource, 'propertycomparisonresource')
    register_resource(docs, UserProfileResource, 'userprofileresource')
    register_resource(docs, MoleculePropertyResource, 'moleculepropertyresource')
    register_resource(docs, MoleculeVisualizationResource, 'moleculevisualizationresource')
    register_resource(docs, SubstructureSearchResource, 'substructuresearchresource')
    register_resource(docs, SimilaritySearchResource, 'similaritysearchresource')
    register_resource(docs, MoleculePropertyCalculationResource, 'moleculepropertycalculationresource')
    register_resource(docs, MoleculeScoreResource, 'moleculescoreresource')
    register_resource(docs, MoleculeIdScoreResource, 'moleculeidscoreresource')
    register_resource(docs, MixtureScoreResource, 'mixturescoreresource')
    register_resource(docs, LabVerificationResource, 'experiment_verification')
    register_resource(docs, LabVerificationResource, 'verification_update')
    register_resource(docs, VerificationStatsResource, 'verification_stats')
    
    # Register experiment resources
    register_resource(docs, ExperimentList, 'experimentlist')
    register_resource(docs, ExperimentDetail, 'experimentdetail')
    register_resource(docs, ExperimentResultList, 'experimentresultlist')
    register_resource(docs, ExperimentResultDetail, 'experimentresultdetail')
    register_resource(docs, ExperimentTimeSeriesList, 'experimenttimeserieslist')
    register_resource(docs, ExperimentAnalysis, 'experimentanalysis')
    register_resource(docs, ExperimentExport, 'experimentexport')
    register_resource(docs, ExperimentImport, 'experimentimport')
    register_resource(docs, ExperimentSearch, 'experimentsearch')
    
    # Register additional resources from modules with their own register_resources functions if needed

    # Register global error handlers
    @app.errorhandler(Exception)
    def handle_exception(error):
        """Handle uncaught exceptions."""
        logger.error(f"Uncaught exception: {str(error)}", exc_info=True)
        return jsonify_standard_response(
            *create_error_response(
                error=error,
                context="Uncaught exception"
            )
        )

    # Log API initialization
    logger.info("API initialized with standardized response formats and error handling")