"""
CryoProtect Analyzer API - Differentiation Resources

This module provides API resources for handling differentiated molecules.
These endpoints allow for querying molecules within differentiation groups
and retrieving information about molecule differentiation.
"""

from flask import request, g, current_app
from flask_restful import Resource, reqparse, abort, fields, marshal_with
from marshmallow import ValidationError
from typing import Dict, List, Optional, Union, Any, Tuple

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

from api.consolidated_utils import (
    get_differentiation_group,
    get_differentiation_group_members,
    get_all_differentiation_groups,
    get_differentiation_description,
    enrich_molecule_data
)
from api.api_decorators import standardize_response
from api.models import molecule_fields

# Define response fields for differentiation groups
differentiation_group_fields = {
    'id': fields.String,
    'name': fields.String,
    'member_count': fields.Integer,
    'members': fields.List(fields.String)
}

class DifferentiationGroupListResource(MethodResource):
    """
    Resource for listing all differentiation groups.
    
    This endpoint returns information about all differentiation groups
    in the system, including their members and descriptions.
    """
    
    @doc(description='Get a list of all differentiation groups',
         tags=['Molecules', 'Differentiation'])
    @standardize_response
    def get(self):
        """Get a list of all differentiation groups."""
        try:
            # Get all differentiation groups
            groups = get_all_differentiation_groups()
            
            return {
                'differentiation_groups': groups,
                'count': len(groups)
            }, 200
            
        except Exception as e:
            current_app.logger.error(f"Error getting differentiation groups: {str(e)}")
            abort(500, message=f"Error getting differentiation groups: {str(e)}")

class DifferentiationGroupResource(MethodResource):
    """
    Resource for getting information about a specific differentiation group.
    
    This endpoint returns detailed information about a differentiation group,
    including all molecules that belong to the group.
    """
    
    @doc(description='Get information about a specific differentiation group',
         tags=['Molecules', 'Differentiation'],
         params={
             'group_id': {
                 'description': 'ID of the differentiation group',
                 'type': 'string',
                 'required': True
             }
         })
    @standardize_response
    def get(self, group_id):
        """Get information about a specific differentiation group."""
        try:
            # Get all members of the differentiation group
            members = get_differentiation_group_members(group_id)
            
            if not members:
                abort(404, message=f"Differentiation group with ID {group_id} not found")
            
            # Get molecule information for each member
            molecules = []
            
            # Import here to avoid circular imports
            from api.utils import get_supabase_client
            
            supabase = get_supabase_client()
            
            # Batch query for all molecules in the group
            molecule_data = (
                supabase.table('molecules')
                .select('*')
                .in_('id', members)
                .execute()
            )
            
            if hasattr(molecule_data, 'data'):
                for molecule in molecule_data.data:
                    # Enrich with differentiation information
                    enriched = enrich_molecule_data(molecule)
                    molecules.append(enriched)
            
            # Create group information
            group_info = {
                'id': group_id,
                'name': f"Differentiation Group {group_id}",
                'description': "Group of molecules with similar names but different structures",
                'member_count': len(members),
                'members': members,
                'molecules': molecules
            }
            
            return group_info, 200
            
        except Exception as e:
            current_app.logger.error(f"Error getting differentiation group: {str(e)}")
            abort(500, message=f"Error getting differentiation group: {str(e)}")

class MoleculeDifferentiationResource(MethodResource):
    """
    Resource for getting differentiation information for a specific molecule.
    
    This endpoint returns information about a molecule's differentiation,
    including the differentiation group it belongs to and other molecules
    in the same group.
    """
    
    @doc(description='Get differentiation information for a molecule',
         tags=['Molecules', 'Differentiation'],
         params={
             'molecule_id': {
                 'description': 'ID of the molecule',
                 'type': 'string',
                 'required': True
             }
         })
    @standardize_response
    def get(self, molecule_id):
        """Get differentiation information for a molecule."""
        try:
            # Get differentiation group for the molecule
            diff_group = get_differentiation_group(molecule_id)
            
            if not diff_group:
                abort(404, message=f"Molecule with ID {molecule_id} is not part of a differentiation group")
            
            # Get all members of the differentiation group
            members = get_differentiation_group_members(diff_group)
            
            # Get differentiation description
            description = get_differentiation_description(molecule_id)
            
            # Create differentiation information
            diff_info = {
                'molecule_id': molecule_id,
                'differentiation_group': diff_group,
                'differentiation_description': description,
                'group_members': members,
                'similar_molecules_count': len(members) - 1  # Exclude the molecule itself
            }
            
            return diff_info, 200
            
        except Exception as e:
            current_app.logger.error(f"Error getting molecule differentiation: {str(e)}")
            abort(500, message=f"Error getting molecule differentiation: {str(e)}")