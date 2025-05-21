"""
CryoProtect Analyzer API - Consolidated Resources Registration

This module registers the updated API resources for consolidated molecule handling.
"""

from flask_restful import Api
from api.updated_resources import (
    UpdatedMoleculeResource,
    ConsolidationStatusResource,
    SecondaryMoleculesResource,
    DifferentiationGroupResource,
    DifferentiationGroupListResource,
    UpdatedBatchOperationResource
)
from api.api_docs import register_resource

def register_updated_resources(api: Api):
    """
    Register updated API resources with consolidated molecule handling.
    
    Args:
        api: Flask-RESTful API instance
    """
    # Register updated molecule resource
    api.add_resource(
        UpdatedMoleculeResource,
        '/api/v1/molecules-v2/<string:molecule_id>',
        endpoint='molecule_resource_v2'
    )
    
    # Register consolidation status resource
    api.add_resource(
        ConsolidationStatusResource,
        '/api/v1/molecules/<string:molecule_id>/consolidation-status',
        endpoint='consolidation_status'
    )
    
    # Register secondary molecules resource
    api.add_resource(
        SecondaryMoleculesResource,
        '/api/v1/molecules/<string:molecule_id>/secondary-molecules',
        endpoint='secondary_molecules'
    )
    
    # Register differentiation group resources
    api.add_resource(
        DifferentiationGroupResource,
        '/api/v1/differentiation-groups/<string:group_id>',
        endpoint='differentiation_group'
    )
    
    api.add_resource(
        DifferentiationGroupListResource,
        '/api/v1/differentiation-groups',
        endpoint='differentiation_groups'
    )
    
    # Register updated batch operation resource
    api.add_resource(
        UpdatedBatchOperationResource,
        '/api/v1/batch-v2',
        endpoint='batch_operation_v2'
    )

def register_updated_docs(docs):
    """
    Register updated API resources for documentation.
    
    Args:
        docs: FlaskApiSpec instance
    """
    # Register updated molecule resource
    register_resource(docs, UpdatedMoleculeResource, 'molecule_resource_v2')
    
    # Register consolidation status resource
    register_resource(docs, ConsolidationStatusResource, 'consolidation_status')
    
    # Register secondary molecules resource
    register_resource(docs, SecondaryMoleculesResource, 'secondary_molecules')
    
    # Register differentiation group resources
    register_resource(docs, DifferentiationGroupResource, 'differentiation_group')
    register_resource(docs, DifferentiationGroupListResource, 'differentiation_groups')
    
    # Register updated batch operation resource
    register_resource(docs, UpdatedBatchOperationResource, 'batch_operation_v2')