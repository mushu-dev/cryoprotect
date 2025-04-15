"""
CryoProtect Analyzer API Package

This package contains the API resources, models, and utilities for the CryoProtect Analyzer.
"""

from flask import Blueprint
from flask_restful import Api

# Create a blueprint for the API
api_bp = Blueprint('api', __name__, url_prefix='/api/v1')
api = Api(api_bp)

# Import resources to register them with the API
from api.resources import (
    MoleculeResource, MoleculeListResource,
    MixtureResource, MixtureListResource,
    PredictionResource, ExperimentResource,
    ComparisonResource
)

# Import RDKit resources
from api.rdkit_resources import (
    MoleculePropertyResource, MoleculeVisualizationResource,
    SubstructureSearchResource, SimilaritySearchResource,
    MoleculePropertyCalculationResource
)

# Register resources with the API
api.add_resource(MoleculeListResource, '/molecules')
api.add_resource(MoleculeResource, '/molecules/<string:molecule_id>')
api.add_resource(MixtureListResource, '/mixtures')
api.add_resource(MixtureResource, '/mixtures/<string:mixture_id>')
api.add_resource(PredictionResource, '/mixtures/<string:mixture_id>/predictions')
api.add_resource(ExperimentResource, '/mixtures/<string:mixture_id>/experiments')
api.add_resource(ComparisonResource, '/mixtures/<string:mixture_id>/comparisons')

# Register RDKit resources
api.add_resource(MoleculePropertyResource, '/rdkit/properties')
api.add_resource(MoleculeVisualizationResource, '/rdkit/visualization')
api.add_resource(SubstructureSearchResource, '/rdkit/substructure')
api.add_resource(SimilaritySearchResource, '/rdkit/similarity')
api.add_resource(MoleculePropertyCalculationResource, '/molecules/<string:molecule_id>/calculate-properties')

def init_app(app):
    """Initialize the API with the Flask app."""
    app.register_blueprint(api_bp)