"""
CryoProtect Analyzer API - Updated Resources

This module contains updated API resource classes that use the consolidated molecule
middleware for proper handling of consolidated molecules.
"""

from flask import request, g, current_app
from flask_restful import Resource, abort
from typing import Dict, List, Optional, Union, Any, Tuple
import uuid

from api.middleware import (
    handle_consolidated_molecules,
    molecule_batch_middleware,
    get_molecule_with_consolidated_info,
    get_differentiation_group_members
)
from api.api_standards import (
    create_standard_response,
    create_error_response, 
    create_success_response
)
from api.models import molecule_fields
from api.utils import get_supabase_client, token_required

class UpdatedMoleculeResource(Resource):
    """
    Resource for retrieving and manipulating molecules with consolidated molecule handling.
    
    This resource class enhances the original MoleculeResource with consolidated
    molecule handling, ensuring that consolidated molecules are properly redirected
    to their primary molecules and that responses include appropriate metadata.
    """
    
    @handle_consolidated_molecules
    def get(self, molecule_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Get a molecule by ID with consolidated molecule handling.
        
        If the requested molecule is consolidated, this endpoint automatically
        redirects to the primary molecule and includes information about the
        consolidation in the response.
        
        Args:
            molecule_id: The molecule ID to retrieve
            
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            # Get molecule with consolidated info
            molecule_data = get_molecule_with_consolidated_info(molecule_id)
            
            if not molecule_data:
                return create_error_response(
                    error=f"Molecule with ID {molecule_id} not found",
                    status_code=404,
                    context="Molecule retrieval"
                )
            
            # Create success response
            return create_success_response(
                data=molecule_data,
                message="Molecule retrieved successfully",
                status_code=200
            )
            
        except Exception as e:
            current_app.logger.error(f"Error retrieving molecule {molecule_id}: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Molecule retrieval"
            )

class ConsolidationStatusResource(Resource):
    """
    Resource for retrieving molecule consolidation status.
    
    This resource provides detailed information about a molecule's consolidation
    status, including references to primary and secondary molecules.
    """
    
    def get(self, molecule_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Get consolidation status for a molecule.
        
        Args:
            molecule_id: The molecule ID to check
            
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            # Validate UUID format
            try:
                uuid.UUID(molecule_id)
            except ValueError:
                return create_error_response(
                    error=f"Invalid molecule ID format: {molecule_id}",
                    status_code=400,
                    context="Consolidation status"
                )
            
            # Get connection
            supabase = get_supabase_client()
            
            # Check if molecule exists
            molecule_result = (
                supabase.table('molecules')
                .select('id, name, consolidated_to')
                .eq('id', molecule_id)
                .single()
                .execute()
            )
            
            if not hasattr(molecule_result, 'data') or not molecule_result.data:
                return create_error_response(
                    error=f"Molecule with ID {molecule_id} not found",
                    status_code=404,
                    context="Consolidation status"
                )
            
            molecule = molecule_result.data
            
            # Check if this is a consolidated (secondary) molecule
            if molecule['consolidated_to']:
                # Get primary molecule
                primary_result = (
                    supabase.table('molecules')
                    .select('id, name')
                    .eq('id', molecule['consolidated_to'])
                    .single()
                    .execute()
                )
                
                primary_molecule = primary_result.data if hasattr(primary_result, 'data') else None
                
                # Prepare response
                consolidation_info = {
                    'is_consolidated': True,
                    'status': 'CONSOLIDATED',
                    'molecule_id': molecule_id,
                    'molecule_name': molecule['name'],
                    'primary_molecule_id': molecule['consolidated_to'],
                    'primary_molecule_name': primary_molecule['name'] if primary_molecule else None
                }
                
                return create_success_response(
                    data=consolidation_info,
                    message="Molecule is consolidated",
                    status_code=200
                )
            
            # Check if this is a primary molecule with secondaries
            secondary_result = (
                supabase.table('molecules')
                .select('id, name')
                .eq('consolidated_to', molecule_id)
                .execute()
            )
            
            secondary_molecules = []
            if hasattr(secondary_result, 'data') and secondary_result.data:
                secondary_molecules = [
                    {
                        'id': m['id'],
                        'name': m['name']
                    }
                    for m in secondary_result.data
                ]
            
            # Prepare response
            consolidation_info = {
                'is_consolidated': False,
                'status': 'PRIMARY' if secondary_molecules else 'INDEPENDENT',
                'molecule_id': molecule_id,
                'molecule_name': molecule['name'],
                'secondary_molecules': secondary_molecules,
                'secondary_count': len(secondary_molecules)
            }
            
            return create_success_response(
                data=consolidation_info,
                message="Molecule consolidation status retrieved",
                status_code=200
            )
            
        except Exception as e:
            current_app.logger.error(f"Error retrieving consolidation status for {molecule_id}: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Consolidation status"
            )

class SecondaryMoleculesResource(Resource):
    """
    Resource for retrieving secondary molecules for a primary molecule.
    
    This resource provides a list of all secondary (consolidated) molecules
    that point to a given primary molecule.
    """
    
    def get(self, molecule_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Get secondary molecules for a primary molecule.
        
        Args:
            molecule_id: The primary molecule ID
            
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            # Validate UUID format
            try:
                uuid.UUID(molecule_id)
            except ValueError:
                return create_error_response(
                    error=f"Invalid molecule ID format: {molecule_id}",
                    status_code=400,
                    context="Secondary molecules"
                )
            
            # Get connection
            supabase = get_supabase_client()
            
            # Check if molecule exists
            molecule_result = (
                supabase.table('molecules')
                .select('id, name')
                .eq('id', molecule_id)
                .single()
                .execute()
            )
            
            if not hasattr(molecule_result, 'data') or not molecule_result.data:
                return create_error_response(
                    error=f"Molecule with ID {molecule_id} not found",
                    status_code=404,
                    context="Secondary molecules"
                )
            
            # Check if this is a consolidated molecule
            consolidated_result = (
                supabase.table('molecules')
                .select('consolidated_to')
                .eq('id', molecule_id)
                .single()
                .execute()
            )
            
            if hasattr(consolidated_result, 'data') and consolidated_result.data and consolidated_result.data['consolidated_to']:
                return create_error_response(
                    error=f"Molecule with ID {molecule_id} is a secondary molecule, not a primary",
                    status_code=400,
                    context="Secondary molecules"
                )
            
            # Get secondary molecules
            secondary_result = (
                supabase.table('molecules')
                .select('id, name, molecular_formula, smiles, inchi_key, created_at, updated_at')
                .eq('consolidated_to', molecule_id)
                .execute()
            )
            
            secondary_molecules = []
            if hasattr(secondary_result, 'data'):
                secondary_molecules = secondary_result.data
            
            # Prepare response
            return create_success_response(
                data={
                    'primary_molecule_id': molecule_id,
                    'primary_molecule_name': molecule_result.data['name'],
                    'secondary_molecules': secondary_molecules,
                    'count': len(secondary_molecules)
                },
                message="Secondary molecules retrieved successfully",
                status_code=200
            )
            
        except Exception as e:
            current_app.logger.error(f"Error retrieving secondary molecules for {molecule_id}: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Secondary molecules"
            )

class DifferentiationGroupResource(Resource):
    """
    Resource for retrieving information about a differentiation group.
    
    This resource provides detailed information about a differentiation group,
    including all molecules in the group and their differentiation descriptions.
    """
    
    def get(self, group_id: str) -> Tuple[Dict[str, Any], int]:
        """
        Get information about a differentiation group.
        
        Args:
            group_id: The differentiation group ID
            
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            # Validate input
            if not group_id or not isinstance(group_id, str):
                return create_error_response(
                    error=f"Invalid differentiation group ID: {group_id}",
                    status_code=400,
                    context="Differentiation group"
                )
            
            # Get molecules in the differentiation group
            members = get_differentiation_group_members(group_id)
            
            if not members:
                return create_error_response(
                    error=f"Differentiation group with ID {group_id} not found or has no members",
                    status_code=404,
                    context="Differentiation group"
                )
            
            # Prepare response
            return create_success_response(
                data={
                    'group_id': group_id,
                    'members': members,
                    'count': len(members)
                },
                message="Differentiation group retrieved successfully",
                status_code=200
            )
            
        except Exception as e:
            current_app.logger.error(f"Error retrieving differentiation group {group_id}: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Differentiation group"
            )

class DifferentiationGroupListResource(Resource):
    """
    Resource for listing all differentiation groups.
    
    This resource provides a list of all differentiation groups in the system,
    along with basic information about each group.
    """
    
    def get(self) -> Tuple[Dict[str, Any], int]:
        """
        List all differentiation groups.
        
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            # Get connection
            supabase = get_supabase_client()
            
            # Get property type ID for differentiationGroup
            property_type_result = (
                supabase.table('property_types')
                .select('id')
                .eq('name', 'differentiationGroup')
                .single()
                .execute()
            )
            
            if not hasattr(property_type_result, 'data') or not property_type_result.data:
                return create_success_response(
                    data={
                        'groups': [],
                        'count': 0
                    },
                    message="No differentiation groups found",
                    status_code=200
                )
            
            property_type_id = property_type_result.data['id']
            
            # Get distinct differentiation groups
            groups_result = (
                supabase.table('molecular_properties')
                .select('property_value, count(*)')
                .eq('property_type_id', property_type_id)
                .group_by('property_value')
                .execute()
            )
            
            groups = []
            if hasattr(groups_result, 'data'):
                for group in groups_result.data:
                    # Get a sample molecule from the group
                    sample_result = (
                        supabase.table('molecular_properties')
                        .select('molecule_id')
                        .eq('property_type_id', property_type_id)
                        .eq('property_value', group['property_value'])
                        .limit(1)
                        .execute()
                    )
                    
                    sample_molecule_id = None
                    if hasattr(sample_result, 'data') and sample_result.data:
                        sample_molecule_id = sample_result.data[0]['molecule_id']
                    
                    groups.append({
                        'id': group['property_value'],
                        'molecule_count': group['count'],
                        'sample_molecule_id': sample_molecule_id
                    })
            
            # Prepare response
            return create_success_response(
                data={
                    'groups': groups,
                    'count': len(groups)
                },
                message="Differentiation groups retrieved successfully",
                status_code=200
            )
            
        except Exception as e:
            current_app.logger.error(f"Error retrieving differentiation groups: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Differentiation groups"
            )

class UpdatedBatchOperationResource(Resource):
    """
    Resource for performing batch operations with consolidated molecule handling.
    
    This resource enhances the original BatchOperationResource with consolidated
    molecule handling, ensuring that consolidated molecules are properly redirected
    to their primary molecules in batch operations.
    """
    
    @molecule_batch_middleware
    def post(self) -> Tuple[Dict[str, Any], int]:
        """
        Perform batch operations with consolidated molecule handling.
        
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            if not request.is_json:
                return create_error_response(
                    error="Request must be JSON",
                    status_code=400,
                    context="Batch operation"
                )
            
            # Get request data
            request_data = request.get_json()
            
            # Check if this is a batch molecule request
            if 'molecule_ids' in request_data:
                return self._handle_batch_molecules(request_data['molecule_ids'])
            
            # Default response for unsupported batch operation
            return create_error_response(
                error="Unsupported batch operation",
                status_code=400,
                context="Batch operation"
            )
            
        except Exception as e:
            current_app.logger.error(f"Error in batch operation: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Batch operation"
            )
    
    def _handle_batch_molecules(self, molecule_ids: List[str]) -> Tuple[Dict[str, Any], int]:
        """
        Handle batch molecule retrieval.
        
        Args:
            molecule_ids: List of molecule IDs to retrieve
            
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            if not molecule_ids or not isinstance(molecule_ids, list):
                return create_error_response(
                    error="molecule_ids must be a non-empty list",
                    status_code=400,
                    context="Batch molecule retrieval"
                )
            
            # Get connection
            supabase = get_supabase_client()
            
            # Query for all molecules
            result = (
                supabase.table('molecules')
                .select('*')
                .in_('id', molecule_ids)
                .execute()
            )
            
            molecules = []
            if hasattr(result, 'data'):
                molecules = result.data
            
            # Prepare response
            return create_success_response(
                data={
                    'molecules': molecules,
                    'count': len(molecules)
                },
                message="Batch molecule retrieval successful",
                status_code=200
            )
            
        except Exception as e:
            current_app.logger.error(f"Error in batch molecule retrieval: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Batch molecule retrieval"
            )