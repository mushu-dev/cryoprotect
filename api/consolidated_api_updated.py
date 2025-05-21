"""
CryoProtect Analyzer API - Enhanced Consolidated API Resources Registration

This module registers the enhanced API resources for handling consolidated molecules.
It provides functions to register all consolidated molecule resources with
the Flask-RESTful API and the API documentation.
"""

from flask_restful import Api
from api.api_docs import register_resource

from api.consolidated_molecule_resource import (
    ConsolidatedMoleculeResource,
    ConsolidatedMoleculeBatchResource,
    PrimaryMoleculeResource,
    ConsolidatedMoleculesListResource,
    ConsolidatedMoleculeAuditResource,
    ConsolidatedMoleculePropertyMigrationResource,
    ConsolidatedMoleculeSearchResource
)

from api.differentiation_resources import (
    DifferentiationGroupListResource,
    DifferentiationGroupResource,
    MoleculeDifferentiationResource
)

def register_enhanced_consolidated_resources(api: Api):
    """
    Register enhanced consolidated molecule API resources.
    
    Args:
        api: Flask-RESTful API instance
    """
    # Register consolidated molecule resources
    api.add_resource(
        ConsolidatedMoleculeResource,
        '/consolidated/molecules/<string:molecule_id>',
        endpoint='consolidated_molecule_enhanced'
    )
    
    api.add_resource(
        ConsolidatedMoleculeBatchResource,
        '/consolidated/batch',
        endpoint='consolidated_batch_enhanced'
    )
    
    api.add_resource(
        PrimaryMoleculeResource,
        '/molecules/<string:molecule_id>/primary',
        endpoint='primary_molecule_enhanced'
    )
    
    api.add_resource(
        ConsolidatedMoleculesListResource,
        '/consolidated',
        endpoint='consolidated_list_enhanced'
    )
    
    # Register new enhanced resources
    api.add_resource(
        ConsolidatedMoleculeAuditResource,
        '/consolidated/molecules/<string:molecule_id>/audit',
        endpoint='consolidated_molecule_audit'
    )
    
    api.add_resource(
        ConsolidatedMoleculePropertyMigrationResource,
        '/consolidated/molecules/<string:molecule_id>/migrate-properties',
        endpoint='consolidated_molecule_property_migration'
    )
    
    api.add_resource(
        ConsolidatedMoleculeSearchResource,
        '/consolidated/search',
        endpoint='consolidated_molecule_search'
    )
    
    # Register differentiation resources
    api.add_resource(
        DifferentiationGroupListResource,
        '/differentiation/groups',
        endpoint='differentiation_groups_enhanced'
    )
    
    api.add_resource(
        DifferentiationGroupResource,
        '/differentiation/groups/<string:group_id>',
        endpoint='differentiation_group_enhanced'
    )
    
    api.add_resource(
        MoleculeDifferentiationResource,
        '/molecules/<string:molecule_id>/differentiation',
        endpoint='molecule_differentiation_enhanced'
    )

def register_enhanced_consolidated_docs(docs):
    """
    Register enhanced consolidated molecule API resources for documentation.
    
    Args:
        docs: FlaskApiSpec instance
    """
    # Register consolidated molecule resources
    register_resource(docs, ConsolidatedMoleculeResource, 'consolidated_molecule_enhanced')
    register_resource(docs, ConsolidatedMoleculeBatchResource, 'consolidated_batch_enhanced')
    register_resource(docs, PrimaryMoleculeResource, 'primary_molecule_enhanced')
    register_resource(docs, ConsolidatedMoleculesListResource, 'consolidated_list_enhanced')
    
    # Register new enhanced resources
    register_resource(docs, ConsolidatedMoleculeAuditResource, 'consolidated_molecule_audit')
    register_resource(docs, ConsolidatedMoleculePropertyMigrationResource, 'consolidated_molecule_property_migration')
    register_resource(docs, ConsolidatedMoleculeSearchResource, 'consolidated_molecule_search')
    
    # Register differentiation resources
    register_resource(docs, DifferentiationGroupListResource, 'differentiation_groups_enhanced')
    register_resource(docs, DifferentiationGroupResource, 'differentiation_group_enhanced')
    register_resource(docs, MoleculeDifferentiationResource, 'molecule_differentiation_enhanced')