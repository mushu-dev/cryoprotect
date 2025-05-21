"""
CryoProtect Analyzer API - Consolidated Molecule Resources

This module provides API resources for handling consolidated molecules.
These endpoints allow for working with the consolidated molecule system,
including retrieving primary/duplicate relationships, updating consolidated
molecule relationships, and migrating properties.
"""

from flask import request, g, current_app, jsonify
from flask_restful import Resource, reqparse, abort, fields, marshal_with
from marshmallow import ValidationError
from typing import Dict, List, Optional, Union, Any, Tuple
from datetime import datetime
import uuid
import json

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

# Import helper functions - note we're using the updated version
from api.consolidated_utils_updated import (
    is_consolidated,
    get_primary_molecule,
    get_consolidated_molecules,
    enrich_molecule_data,
    get_molecule_audit_history
)
from api.api_decorators import standardize_response

class ConsolidatedMoleculeResource(MethodResource):
    """
    Resource for handling individual consolidated molecules.
    
    This endpoint provides access to a specific consolidated molecule,
    handling redirection to primary molecules and enrichment of responses.
    """
    
    @doc(description='Get a molecule with consolidated molecule handling',
         tags=['Molecules', 'Consolidated'])
    @standardize_response
    def get(self, molecule_id):
        """
        Get a molecule with consolidated molecule handling.
        
        If the requested molecule is consolidated, the endpoint will
        provide information about the primary molecule and consolidation status.
        """
        try:
            from api.utils import get_supabase_client
            
            supabase = get_supabase_client()
            
            # Get molecule data - use consolidated_molecules table
            result = (
                supabase.table('consolidated_molecules')
                .select('*')
                .eq('id', molecule_id)
                .single()
                .execute()
            )
            
            if not hasattr(result, 'data') or not result.data:
                abort(404, message=f"Molecule with ID {molecule_id} not found")
            
            molecule_data = result.data
            
            # Check if this is a duplicate molecule
            if molecule_data.get('molecule_status') == 'duplicate' and molecule_data.get('primary_molecule_id'):
                # Add a note about primary molecule
                molecule_data['primary_note'] = (
                    f"This molecule is a duplicate. The primary molecule is "
                    f"{molecule_data.get('primary_molecule_id')} "
                    f"({molecule_data.get('primary_molecule_name')})"
                )
                
                # Include audit history for this consolidation
                molecule_data['audit_history'] = get_molecule_audit_history(molecule_id)
            
            # If this is a primary molecule, get its duplicates
            elif molecule_data.get('molecule_status') == 'primary':
                duplicates = get_consolidated_molecules(molecule_id)
                
                if duplicates:
                    # Get duplicate molecule details
                    dup_result = (
                        supabase.table('consolidated_molecules')
                        .select('id, name')
                        .in_('id', duplicates)
                        .execute()
                    )
                    
                    if hasattr(dup_result, 'data'):
                        molecule_data['duplicate_molecules'] = dup_result.data
                        molecule_data['duplicate_count'] = len(dup_result.data)
                else:
                    molecule_data['duplicate_molecules'] = []
                    molecule_data['duplicate_count'] = 0
            
            return molecule_data, 200
            
        except Exception as e:
            current_app.logger.error(f"Error getting consolidated molecule: {str(e)}")
            abort(500, message=f"Error getting consolidated molecule: {str(e)}")
    
    @doc(description='Update a consolidated molecule relationship',
         tags=['Molecules', 'Consolidated'])
    @standardize_response
    def put(self, molecule_id):
        """
        Update a consolidated molecule relationship.
        
        This endpoint allows setting a molecule as a duplicate of another molecule,
        or changing its primary molecule.
        """
        try:
            from api.utils import get_supabase_client, token_required
            
            # This operation requires authentication
            current_user = token_required()
            if not current_user:
                abort(401, message="Unauthorized")
            
            # Parse request body
            data = request.get_json()
            if not data:
                abort(400, message="No data provided")
            
            primary_id = data.get('primary_molecule_id')
            if not primary_id:
                abort(400, message="primary_molecule_id is required")
            
            # Check if molecule and primary molecule exist
            supabase = get_supabase_client()
            
            # Check molecule exists
            molecule_result = (
                supabase.table('consolidated_molecules')
                .select('id, name, molecule_status, primary_molecule_id')
                .eq('id', molecule_id)
                .single()
                .execute()
            )
            
            if not hasattr(molecule_result, 'data') or not molecule_result.data:
                abort(404, message=f"Molecule with ID {molecule_id} not found")
            
            # Check primary molecule exists
            primary_result = (
                supabase.table('consolidated_molecules')
                .select('id, name, molecule_status')
                .eq('id', primary_id)
                .single()
                .execute()
            )
            
            if not hasattr(primary_result, 'data') or not primary_result.data:
                abort(404, message=f"Primary molecule with ID {primary_id} not found")
            
            # Check that primary molecule is not already a duplicate
            primary_molecule = primary_result.data
            if primary_molecule.get('molecule_status') == 'duplicate':
                abort(400, message="Cannot set a duplicate molecule as primary. Use its primary instead.")
            
            # Check that we're not creating a circular reference
            if molecule_id == primary_id:
                abort(400, message="Cannot set a molecule as its own primary")
            
            # Update the molecule to be a duplicate
            update_result = (
                supabase.table('consolidated_molecules')
                .update({
                    'molecule_status': 'duplicate',
                    'is_consolidated': True,
                    'primary_molecule_id': primary_id,
                    'primary_molecule_name': primary_molecule.get('name'),
                    'updated_at': datetime.now().isoformat()
                })
                .eq('id', molecule_id)
                .execute()
            )
            
            if not hasattr(update_result, 'data'):
                abort(500, message="Failed to update molecule")
            
            # Update the primary molecule status if needed
            if primary_molecule.get('molecule_status') == 'original':
                primary_update = (
                    supabase.table('consolidated_molecules')
                    .update({
                        'molecule_status': 'primary',
                        'is_consolidated': True,
                        'updated_at': datetime.now().isoformat()
                    })
                    .eq('id', primary_id)
                    .execute()
                )
            
            # Record this change in the audit log
            audit_id = str(uuid.uuid4())
            audit_data = {
                'id': audit_id,
                'table_name': 'consolidated_molecules',
                'record_id': molecule_id,
                'operation': 'consolidation',
                'old_value': json.dumps(molecule_result.data, default=str),
                'new_value': json.dumps({
                    'molecule_status': 'duplicate',
                    'primary_molecule_id': primary_id,
                    'primary_molecule_name': primary_molecule.get('name')
                }, default=str),
                'user_id': current_user.get('id'),
                'timestamp': datetime.now().isoformat()
            }
            
            audit_result = (
                supabase.table('scientific_data_audit')
                .insert(audit_data)
                .execute()
            )
            
            # Return updated molecule data
            updated_result = (
                supabase.table('consolidated_molecules')
                .select('*')
                .eq('id', molecule_id)
                .single()
                .execute()
            )
            
            if hasattr(updated_result, 'data'):
                return {
                    'message': 'Molecule successfully consolidated',
                    'molecule': updated_result.data,
                    'audit_id': audit_id
                }, 200
            else:
                return {
                    'message': 'Molecule successfully consolidated, but could not retrieve updated data',
                    'audit_id': audit_id
                }, 200
            
        except Exception as e:
            current_app.logger.error(f"Error updating consolidated molecule: {str(e)}")
            abort(500, message=f"Error updating consolidated molecule: {str(e)}")

class MoleculeConsolidationResource(MethodResource):
    """
    Resource for searching potential molecule duplicates and managing consolidation.
    """
    
    @doc(description='Find potential duplicate molecules by InChIKey',
         tags=['Molecules', 'Consolidated'])
    @standardize_response
    def get(self):
        """
        Find potential duplicate molecules by InChIKey.
        
        This endpoint finds all molecules that share the same InChIKey,
        which are potential candidates for consolidation.
        """
        try:
            from api.utils import get_supabase_client
            
            # Parse query parameters
            inchikey = request.args.get('inchikey')
            limit = request.args.get('limit', 100)
            try:
                limit = int(limit)
            except:
                limit = 100
            
            if not inchikey:
                abort(400, message="inchikey parameter is required")
            
            supabase = get_supabase_client()
            
            # Find molecules with the same InChIKey
            result = (
                supabase.table('consolidated_molecules')
                .select('*')
                .eq('inchikey', inchikey)
                .limit(limit)
                .execute()
            )
            
            if not hasattr(result, 'data'):
                return {'molecules': [], 'count': 0}, 200
            
            molecules = result.data
            
            # Group by molecule_status
            grouped = {
                'primary': [],
                'duplicate': [],
                'original': []
            }
            
            for mol in molecules:
                status = mol.get('molecule_status', 'original')
                if status in grouped:
                    grouped[status].append(mol)
                else:
                    grouped['original'].append(mol)
            
            return {
                'molecules': molecules,
                'count': len(molecules),
                'grouped_by_status': grouped,
                'inchikey': inchikey
            }, 200
            
        except Exception as e:
            current_app.logger.error(f"Error finding potential duplicates: {str(e)}")
            abort(500, message=f"Error finding potential duplicates: {str(e)}")
    
    @doc(description='Consolidate a batch of molecules',
         tags=['Molecules', 'Consolidated'])
    @standardize_response
    def post(self):
        """
        Consolidate a batch of molecules.
        
        This endpoint allows consolidating multiple molecules under a single
        primary molecule in one operation.
        """
        try:
            from api.utils import get_supabase_client, token_required
            
            # This operation requires authentication
            current_user = token_required()
            if not current_user:
                abort(401, message="Unauthorized")
            
            # Parse request body
            data = request.get_json()
            if not data:
                abort(400, message="No data provided")
            
            primary_id = data.get('primary_molecule_id')
            if not primary_id:
                abort(400, message="primary_molecule_id is required")
            
            duplicate_ids = data.get('duplicate_molecule_ids', [])
            if not duplicate_ids or not isinstance(duplicate_ids, list):
                abort(400, message="duplicate_molecule_ids must be a non-empty list")
            
            # Check if primary molecule exists
            supabase = get_supabase_client()
            
            primary_result = (
                supabase.table('consolidated_molecules')
                .select('id, name, molecule_status')
                .eq('id', primary_id)
                .single()
                .execute()
            )
            
            if not hasattr(primary_result, 'data') or not primary_result.data:
                abort(404, message=f"Primary molecule with ID {primary_id} not found")
            
            primary_molecule = primary_result.data
            
            # Ensure primary molecule is not a duplicate
            if primary_molecule.get('molecule_status') == 'duplicate':
                abort(400, message="Cannot set a duplicate molecule as primary")
            
            # Remove primary ID from duplicates list if present
            if primary_id in duplicate_ids:
                duplicate_ids.remove(primary_id)
                
            if not duplicate_ids:
                abort(400, message="No valid duplicate molecules provided")
            
            # Update each duplicate molecule
            updated_ids = []
            failed_ids = []
            
            for dup_id in duplicate_ids:
                try:
                    # Update the molecule
                    update_result = (
                        supabase.table('consolidated_molecules')
                        .update({
                            'molecule_status': 'duplicate',
                            'is_consolidated': True,
                            'primary_molecule_id': primary_id,
                            'primary_molecule_name': primary_molecule.get('name'),
                            'updated_at': datetime.now().isoformat()
                        })
                        .eq('id', dup_id)
                        .execute()
                    )
                    
                    if hasattr(update_result, 'data'):
                        updated_ids.append(dup_id)
                        
                        # Add audit entry
                        audit_id = str(uuid.uuid4())
                        audit_data = {
                            'id': audit_id,
                            'table_name': 'consolidated_molecules',
                            'record_id': dup_id,
                            'operation': 'batch_consolidation',
                            'new_value': json.dumps({
                                'primary_molecule_id': primary_id,
                                'primary_molecule_name': primary_molecule.get('name')
                            }, default=str),
                            'user_id': current_user.get('id'),
                            'timestamp': datetime.now().isoformat()
                        }
                        
                        supabase.table('scientific_data_audit').insert(audit_data).execute()
                    else:
                        failed_ids.append(dup_id)
                except Exception as e:
                    current_app.logger.error(f"Error consolidating molecule {dup_id}: {str(e)}")
                    failed_ids.append(dup_id)
            
            # Update primary molecule status if needed
            if primary_molecule.get('molecule_status') == 'original' and updated_ids:
                primary_update = (
                    supabase.table('consolidated_molecules')
                    .update({
                        'molecule_status': 'primary',
                        'is_consolidated': True,
                        'updated_at': datetime.now().isoformat()
                    })
                    .eq('id', primary_id)
                    .execute()
                )
            
            return {
                'message': f"Successfully consolidated {len(updated_ids)} molecules",
                'primary_molecule_id': primary_id,
                'primary_molecule_name': primary_molecule.get('name'),
                'updated_molecule_ids': updated_ids,
                'failed_molecule_ids': failed_ids
            }, 200
            
        except Exception as e:
            current_app.logger.error(f"Error batch consolidating molecules: {str(e)}")
            abort(500, message=f"Error batch consolidating molecules: {str(e)}")

class MoleculePropertyMigrationResource(MethodResource):
    """
    Resource for migrating properties between consolidated molecules.
    """
    
    @doc(description='Migrate properties from a duplicate to its primary molecule',
         tags=['Molecules', 'Consolidated'])
    @standardize_response
    def post(self):
        """
        Migrate properties from a duplicate to its primary molecule.
        
        This endpoint allows migrating specific properties from a duplicate
        molecule to its primary molecule, or migrating all properties that
        don't already exist in the primary.
        """
        try:
            from api.utils import get_supabase_client, token_required
            
            # This operation requires authentication
            current_user = token_required()
            if not current_user:
                abort(401, message="Unauthorized")
            
            # Parse request body
            data = request.get_json()
            if not data:
                abort(400, message="No data provided")
            
            source_id = data.get('source_molecule_id')
            if not source_id:
                abort(400, message="source_molecule_id is required")
            
            target_id = data.get('target_molecule_id')
            property_ids = data.get('property_ids', [])  # Optional
            
            # If target_id is not provided, get primary molecule ID
            if not target_id:
                target_id = get_primary_molecule(source_id)
                
                if target_id == source_id:
                    abort(400, message="Source molecule is not a duplicate. Please specify target_molecule_id")
            
            # Check that molecules exist
            supabase = get_supabase_client()
            
            source_result = (
                supabase.table('consolidated_molecules')
                .select('id, name, molecule_status, primary_molecule_id')
                .eq('id', source_id)
                .single()
                .execute()
            )
            
            if not hasattr(source_result, 'data') or not source_result.data:
                abort(404, message=f"Source molecule with ID {source_id} not found")
            
            target_result = (
                supabase.table('consolidated_molecules')
                .select('id, name')
                .eq('id', target_id)
                .single()
                .execute()
            )
            
            if not hasattr(target_result, 'data') or not target_result.data:
                abort(404, message=f"Target molecule with ID {target_id} not found")
            
            # Get properties from source molecule
            if property_ids and len(property_ids) > 0:
                # Get specific properties
                props_result = (
                    supabase.table('molecular_properties')
                    .select('*')
                    .eq('molecule_id', source_id)
                    .in_('id', property_ids)
                    .execute()
                )
            else:
                # Get all properties
                props_result = (
                    supabase.table('molecular_properties')
                    .select('*')
                    .eq('molecule_id', source_id)
                    .execute()
                )
            
            if not hasattr(props_result, 'data') or not props_result.data:
                return {
                    'message': "No properties found to migrate",
                    'source_molecule_id': source_id,
                    'target_molecule_id': target_id,
                    'migrated_count': 0
                }, 200
            
            source_properties = props_result.data
            
            # Get existing property types for target molecule
            target_props_result = (
                supabase.table('molecular_properties')
                .select('property_type_id')
                .eq('molecule_id', target_id)
                .execute()
            )
            
            existing_property_types = []
            if hasattr(target_props_result, 'data'):
                existing_property_types = [prop.get('property_type_id') for prop in target_props_result.data]
            
            # Migrate properties
            migrated = []
            skipped = []
            
            for prop in source_properties:
                property_type_id = prop.get('property_type_id')
                
                # Skip if target already has this property type
                if property_type_id in existing_property_types:
                    skipped.append(prop.get('id'))
                    continue
                
                # Create new property for target
                new_prop_id = str(uuid.uuid4())
                new_prop = {
                    'id': new_prop_id,
                    'molecule_id': target_id,
                    'property_type_id': property_type_id,
                    'numeric_value': prop.get('numeric_value'),
                    'text_value': prop.get('text_value'),
                    'property_name': prop.get('property_name'),
                    'property_type': prop.get('property_type'),
                    'source': f"{prop.get('source', 'unknown')}_migrated",
                    'created_at': datetime.now().isoformat(),
                    'updated_at': datetime.now().isoformat()
                }
                
                insert_result = (
                    supabase.table('molecular_properties')
                    .insert(new_prop)
                    .execute()
                )
                
                if hasattr(insert_result, 'data'):
                    migrated.append({
                        'source_property_id': prop.get('id'),
                        'target_property_id': new_prop_id,
                        'property_name': prop.get('property_name')
                    })
                    
                    # Add audit entry
                    audit_id = str(uuid.uuid4())
                    audit_data = {
                        'id': audit_id,
                        'table_name': 'molecular_properties',
                        'record_id': new_prop_id,
                        'operation': 'property_migration',
                        'old_value': json.dumps({
                            'source_molecule_id': source_id,
                            'source_property_id': prop.get('id')
                        }, default=str),
                        'new_value': json.dumps({
                            'target_molecule_id': target_id,
                            'property_type_id': property_type_id
                        }, default=str),
                        'user_id': current_user.get('id'),
                        'timestamp': datetime.now().isoformat()
                    }
                    
                    supabase.table('scientific_data_audit').insert(audit_data).execute()
            
            return {
                'message': f"Successfully migrated {len(migrated)} properties",
                'source_molecule_id': source_id,
                'target_molecule_id': target_id,
                'migrated_properties': migrated,
                'skipped_properties': skipped,
                'migrated_count': len(migrated),
                'skipped_count': len(skipped)
            }, 200
            
        except Exception as e:
            current_app.logger.error(f"Error migrating properties: {str(e)}")
            abort(500, message=f"Error migrating properties: {str(e)}")

def register_resources(api):
    """Register consolidated molecule resources with the API."""
    api.add_resource(ConsolidatedMoleculeResource, '/consolidated-molecules/<uuid:molecule_id>')
    api.add_resource(MoleculeConsolidationResource, '/molecule-consolidation')
    api.add_resource(MoleculePropertyMigrationResource, '/molecule-property-migration')