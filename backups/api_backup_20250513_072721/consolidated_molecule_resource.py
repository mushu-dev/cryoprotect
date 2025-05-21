"""
CryoProtect Analyzer API - Consolidated Molecule Resources

This module provides API resources for handling consolidated molecules.
These endpoints demonstrate how to work with the consolidated molecule pattern
established in Phase 2 of the data quality enhancement.
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

from api.consolidated_decorators import (
    handle_consolidated_molecules,
    handle_batch_consolidated_molecules
)
from api.consolidated_utils import (
    is_consolidated,
    get_primary_molecule,
    get_consolidated_molecules,
    enrich_molecule_data
)
from api.api_decorators import standardize_response
from api.models import molecule_fields

class ConsolidatedMoleculeResource(MethodResource):
    """
    Resource for handling consolidated molecules.
    
    This endpoint demonstrates how to handle consolidated molecules
    in API responses, redirecting to the primary molecule when needed
    and enriching responses with consolidated information.
    """
    
    @doc(description='Get a molecule with consolidated molecule handling',
         tags=['Molecules', 'Consolidated'])
    @handle_consolidated_molecules
    @standardize_response
    def get(self, molecule_id):
        """
        Get a molecule with consolidated molecule handling.
        
        If the requested molecule is consolidated, the endpoint will
        automatically handle the redirection to the primary molecule.
        """
        try:
            # Import here to avoid circular imports
            from api.utils import get_supabase_client
            
            supabase = get_supabase_client()
            
            # Query for the molecule
            result = (
                supabase.table('molecules')
                .select('*')
                .eq('id', molecule_id)
                .single()
                .execute()
            )
            
            if not hasattr(result, 'data') or not result.data:
                abort(404, message=f"Molecule with ID {molecule_id} not found")
            
            # Check if original_molecule_id is different from molecule_id
            # This indicates that a consolidated molecule was requested
            if hasattr(g, 'original_molecule_id') and g.original_molecule_id != molecule_id:
                # Add a note about the redirection
                result.data['redirection_note'] = (
                    f"This response shows the primary molecule. The requested molecule "
                    f"({g.original_molecule_id}) has been consolidated into this molecule."
                )
            
            return result.data, 200
            
        except Exception as e:
            current_app.logger.error(f"Error getting molecule: {str(e)}")
            abort(500, message=f"Error getting molecule: {str(e)}")

class ConsolidatedMoleculeBatchResource(MethodResource):
    """
    Resource for batch operations on consolidated molecules.
    
    This endpoint demonstrates how to handle batch operations with
    consolidated molecules, ensuring proper redirection and enrichment.
    """
    
    # Parser for batch molecule requests
    batch_parser = reqparse.RequestParser()
    batch_parser.add_argument('molecule_ids', type=list, required=True,
                           help='List of molecule IDs is required',
                           location='json')
    
    @doc(description='Get multiple molecules with consolidated molecule handling',
         tags=['Molecules', 'Consolidated', 'Batch'])
    @handle_batch_consolidated_molecules
    @standardize_response
    def post(self):
        """
        Get multiple molecules with consolidated molecule handling.
        
        This endpoint accepts a list of molecule IDs and returns their data,
        handling any consolidated molecules by redirecting to their primary molecules.
        """
        try:
            # Import here to avoid circular imports
            from api.utils import get_supabase_client
            
            # Parse the request
            args = self.batch_parser.parse_args()
            molecule_ids = args['molecule_ids']
            
            if not molecule_ids:
                return {'molecules': [], 'count': 0}, 200
            
            # Get supabase client
            supabase = get_supabase_client()
            
            # Query for all requested molecules
            result = (
                supabase.table('molecules')
                .select('*')
                .in_('id', molecule_ids)
                .execute()
            )
            
            if not hasattr(result, 'data'):
                return {'molecules': [], 'count': 0}, 200
            
            # Add metadata about any redirections that occurred
            meta = {}
            if hasattr(g, 'consolidated_map') and g.consolidated_map:
                meta['consolidated_redirections'] = g.consolidated_map
            
            return {
                'molecules': result.data,
                'count': len(result.data),
                'meta': meta
            }, 200
            
        except Exception as e:
            current_app.logger.error(f"Error getting batch molecules: {str(e)}")
            abort(500, message=f"Error getting batch molecules: {str(e)}")

class PrimaryMoleculeResource(MethodResource):
    """
    Resource for getting the primary molecule for a given molecule.
    
    This endpoint returns the primary molecule for a given molecule ID,
    which is useful for clients that need to explicitly handle consolidation.
    """
    
    @doc(description='Get the primary molecule for a given molecule ID',
         tags=['Molecules', 'Consolidated'])
    @standardize_response
    def get(self, molecule_id):
        """
        Get the primary molecule for a given molecule ID.
        
        If the molecule is already a primary molecule, it returns the molecule itself.
        If the molecule is consolidated, it returns the primary molecule.
        """
        try:
            # Check if the molecule is consolidated
            if is_consolidated(molecule_id):
                # Get the primary molecule ID
                primary_id = get_primary_molecule(molecule_id)
                
                # Import here to avoid circular imports
                from api.utils import get_supabase_client
                
                supabase = get_supabase_client()
                
                # Query for the primary molecule
                result = (
                    supabase.table('molecules')
                    .select('*')
                    .eq('id', primary_id)
                    .single()
                    .execute()
                )
                
                if not hasattr(result, 'data') or not result.data:
                    abort(404, message=f"Primary molecule with ID {primary_id} not found")
                
                # Add consolidation information
                result.data['is_consolidated'] = False  # Primary molecules are not consolidated
                result.data['consolidated_molecules'] = get_consolidated_molecules(primary_id)
                result.data['consolidation_info'] = {
                    'requested_molecule': molecule_id,
                    'is_primary': False,
                    'primary_molecule': primary_id
                }
                
                return result.data, 200
            else:
                # The molecule is already a primary, get its details
                from api.utils import get_supabase_client
                
                supabase = get_supabase_client()
                
                # Query for the molecule
                result = (
                    supabase.table('molecules')
                    .select('*')
                    .eq('id', molecule_id)
                    .single()
                    .execute()
                )
                
                if not hasattr(result, 'data') or not result.data:
                    abort(404, message=f"Molecule with ID {molecule_id} not found")
                
                # Add consolidation information
                result.data['is_consolidated'] = False  # This is a primary molecule
                result.data['consolidated_molecules'] = get_consolidated_molecules(molecule_id)
                result.data['consolidation_info'] = {
                    'requested_molecule': molecule_id,
                    'is_primary': True,
                    'primary_molecule': molecule_id
                }
                
                return result.data, 200
                
        except Exception as e:
            current_app.logger.error(f"Error getting primary molecule: {str(e)}")
            abort(500, message=f"Error getting primary molecule: {str(e)}")

class ConsolidatedMoleculesListResource(MethodResource):
    """
    Resource for listing all consolidated molecules.
    
    This endpoint returns all consolidated molecule relationships,
    which is useful for clients that need to understand the consolidation structure.
    """
    
    @doc(description='Get a list of all consolidated molecule relationships',
         tags=['Molecules', 'Consolidated'])
    @standardize_response
    def get(self):
        """
        Get a list of all consolidated molecule relationships.
        
        Returns all primary molecules and their secondary (consolidated) molecules.
        """
        try:
            # Import here to avoid circular imports
            from api.utils import get_supabase_client
            
            supabase = get_supabase_client()
            
            # Query for all molecules that are consolidated (have a consolidated_to value)
            result = (
                supabase.table('molecules')
                .select('id, consolidated_to')
                .not_.is_('consolidated_to', 'null')
                .execute()
            )
            
            if not hasattr(result, 'data'):
                return {'consolidated_relationships': [], 'count': 0}, 200
            
            # Organize by primary molecule
            relationships = {}
            for molecule in result.data:
                primary_id = molecule['consolidated_to']
                secondary_id = molecule['id']
                
                if primary_id not in relationships:
                    relationships[primary_id] = []
                
                relationships[primary_id].append(secondary_id)
            
            # Format as a list of primary -> secondary relationships
            formatted_relationships = []
            for primary_id, secondary_ids in relationships.items():
                formatted_relationships.append({
                    'primary_molecule_id': primary_id,
                    'consolidated_molecule_ids': secondary_ids,
                    'count': len(secondary_ids)
                })
            
            return {
                'consolidated_relationships': formatted_relationships,
                'count': len(formatted_relationships),
                'total_consolidated_molecules': len(result.data)
            }, 200
            
        except Exception as e:
            current_app.logger.error(f"Error getting consolidated relationships: {str(e)}")
            abort(500, message=f"Error getting consolidated relationships: {str(e)}")