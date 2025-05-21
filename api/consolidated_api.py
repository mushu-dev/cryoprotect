"""
CryoProtect Analyzer API - Consolidated API Resources Registration

This module registers the API resources for handling consolidated molecules.
It provides a function to register all consolidated molecule resources with
the Flask-RESTful API.
"""

from flask_restful import Api
from api.api_docs import register_resource

from api.consolidated_molecule_resource import (
    ConsolidatedMoleculeResource as OriginalConsolidatedMoleculeResource,
    ConsolidatedMoleculeBatchResource,
    PrimaryMoleculeResource,
    ConsolidatedMoleculesListResource,
    MoleculeConsolidationResource,
    MoleculePropertyMigrationResource
)

from api.differentiation_resources import (
    DifferentiationGroupListResource,
    DifferentiationGroupResource,
    MoleculeDifferentiationResource
)

def register_consolidated_resources(api: Api):
    """
    Register consolidated molecule API resources.
    
    Args:
        api: Flask-RESTful API instance
    """
    # Register consolidated molecule resources
    api.add_resource(
        OriginalConsolidatedMoleculeResource,
        '/consolidated/molecules/<string:molecule_id>',
        endpoint='consolidated_molecule'
    )
    
    api.add_resource(
        ConsolidatedMoleculeBatchResource,
        '/consolidated/batch',
        endpoint='consolidated_batch'
    )
    
    api.add_resource(
        PrimaryMoleculeResource,
        '/molecules/<string:molecule_id>/primary',
        endpoint='primary_molecule'
    )
    
    api.add_resource(
        ConsolidatedMoleculesListResource,
        '/consolidated',
        endpoint='consolidated_list'
    )
    
    # Register differentiation resources
    api.add_resource(
        DifferentiationGroupListResource,
        '/differentiation/groups',
        endpoint='differentiation_groups'
    )
    
    api.add_resource(
        DifferentiationGroupResource,
        '/differentiation/groups/<string:group_id>',
        endpoint='differentiation_group'
    )
    
    api.add_resource(
        MoleculeDifferentiationResource,
        '/molecules/<string:molecule_id>/differentiation',
        endpoint='molecule_differentiation'
    )
    
    # Register new consolidated molecule endpoints for CHEMBL verification
    api.add_resource(
        MoleculeConsolidationResource,
        '/molecule-consolidation',
        endpoint='molecule_consolidation'
    )
    
    api.add_resource(
        MoleculePropertyMigrationResource,
        '/molecule-property-migration',
        endpoint='molecule_property_migration'
    )

def register_consolidated_docs(docs):
    """
    Register consolidated molecule API resources for documentation.
    
    Args:
        docs: FlaskApiSpec instance
    """
    # Register consolidated molecule resources
    register_resource(docs, OriginalConsolidatedMoleculeResource, 'consolidated_molecule')
    register_resource(docs, ConsolidatedMoleculeBatchResource, 'consolidated_batch')
    register_resource(docs, PrimaryMoleculeResource, 'primary_molecule')
    register_resource(docs, ConsolidatedMoleculesListResource, 'consolidated_list')
    
    # Register differentiation resources
    register_resource(docs, DifferentiationGroupListResource, 'differentiation_groups')
    register_resource(docs, DifferentiationGroupResource, 'differentiation_group')
    register_resource(docs, MoleculeDifferentiationResource, 'molecule_differentiation')
    
    # Register new consolidated molecule endpoints for CHEMBL verification
    register_resource(docs, MoleculeConsolidationResource, 'molecule_consolidation')
    register_resource(docs, MoleculePropertyMigrationResource, 'molecule_property_migration')