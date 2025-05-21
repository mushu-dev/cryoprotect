"""
CryoProtect Analyzer API - Audit Resources

This module provides API resources for retrieving audit information,
including the consolidated molecule audit trail.
"""

from flask import request, g, current_app
from flask_restful import Resource, abort
from typing import Dict, List, Optional, Union, Any, Tuple

from api.api_standards import (
    create_standard_response,
    create_error_response, 
    create_success_response
)
from api.utils import get_supabase_client, token_required

class ConsolidationAuditResource(Resource):
    """
    Resource for retrieving consolidated molecule audit trail.
    
    This resource provides access to the audit trail of molecule consolidation
    operations, including when molecules are consolidated, deconsolidated, or
    have their primary molecule changed.
    """
    
    @token_required
    def get(self, user=None) -> Tuple[Dict[str, Any], int]:
        """
        Get consolidated molecule audit trail.
        
        Args:
            user: User object from token_required decorator
            
        Returns:
            Tuple of (response_data, status_code)
        """
        try:
            # Check admin permission
            if not user or 'admin' not in user.get('app_metadata', {}).get('roles', []):
                return create_error_response(
                    error="Administrative privileges required",
                    status_code=403,
                    context="Consolidation audit"
                )
            
            # Parse pagination parameters
            page = request.args.get('page', 1, type=int)
            per_page = request.args.get('per_page', 20, type=int)
            
            # Validate and constrain pagination values
            page = max(1, page)
            per_page = max(1, min(100, per_page))
            
            # Parse filter parameters
            molecule_id = request.args.get('molecule_id')
            operation_type = request.args.get('operation_type')
            performed_by = request.args.get('performed_by')
            start_date = request.args.get('start_date')
            end_date = request.args.get('end_date')
            
            # Get connection
            supabase = get_supabase_client()
            
            # Build query
            query = (
                supabase.table('molecule_consolidation_audit')
                .select('*')
            )
            
            # Apply filters
            if molecule_id:
                query = query.or_(f'primary_molecule_id.eq.{molecule_id},secondary_molecule_id.eq.{molecule_id}')
            
            if operation_type:
                query = query.eq('operation_type', operation_type)
            
            if performed_by:
                query = query.eq('performed_by', performed_by)
            
            if start_date:
                query = query.gte('performed_at', start_date)
            
            if end_date:
                query = query.lte('performed_at', end_date)
            
            # Apply pagination
            query = query.order('performed_at', 'desc')
            query = query.range((page - 1) * per_page, page * per_page - 1)
            
            # Execute query
            result = query.execute()
            
            # Get total count for pagination
            count_query = (
                supabase.table('molecule_consolidation_audit')
                .select('id', count='exact')
            )
            
            # Apply the same filters to count query
            if molecule_id:
                count_query = count_query.or_(f'primary_molecule_id.eq.{molecule_id},secondary_molecule_id.eq.{molecule_id}')
            
            if operation_type:
                count_query = count_query.eq('operation_type', operation_type)
            
            if performed_by:
                count_query = count_query.eq('performed_by', performed_by)
            
            if start_date:
                count_query = count_query.gte('performed_at', start_date)
            
            if end_date:
                count_query = count_query.lte('performed_at', end_date)
            
            # Execute count query
            count_result = count_query.execute()
            total_count = count_result.count if hasattr(count_result, 'count') else 0
            
            # Build pagination info
            total_pages = (total_count + per_page - 1) // per_page if per_page > 0 else 0
            
            pagination = {
                'page': page,
                'per_page': per_page,
                'total_items': total_count,
                'total_pages': total_pages,
                'has_next': page < total_pages,
                'has_prev': page > 1
            }
            
            # Add molecule names to audit records
            audit_records = []
            if hasattr(result, 'data'):
                for record in result.data:
                    # Get primary molecule name
                    primary_result = (
                        supabase.table('molecules')
                        .select('name')
                        .eq('id', record['primary_molecule_id'])
                        .single()
                        .execute()
                    )
                    
                    primary_name = primary_result.data['name'] if hasattr(primary_result, 'data') and primary_result.data else 'Unknown'
                    
                    # Get secondary molecule name
                    secondary_result = (
                        supabase.table('molecules')
                        .select('name')
                        .eq('id', record['secondary_molecule_id'])
                        .single()
                        .execute()
                    )
                    
                    secondary_name = secondary_result.data['name'] if hasattr(secondary_result, 'data') and secondary_result.data else 'Unknown'
                    
                    # Add names to record
                    record['primary_molecule_name'] = primary_name
                    record['secondary_molecule_name'] = secondary_name
                    
                    audit_records.append(record)
            
            # Prepare response
            return create_success_response(
                data={
                    'audit_records': audit_records,
                    'count': len(audit_records)
                },
                message="Consolidation audit records retrieved successfully",
                status_code=200,
                pagination=pagination
            )
            
        except Exception as e:
            current_app.logger.error(f"Error retrieving consolidation audit records: {str(e)}", exc_info=True)
            return create_error_response(
                error=str(e),
                status_code=500,
                context="Consolidation audit"
            )

def register_audit_resources(api):
    """
    Register audit resources with the API.
    
    Args:
        api: Flask-RESTful API instance
    """
    api.add_resource(
        ConsolidationAuditResource,
        '/api/v1/admin/consolidation-audit',
        endpoint='consolidation_audit'
    )

def register_audit_docs(docs):
    """
    Register audit resources for documentation.
    
    Args:
        docs: FlaskApiSpec instance
    """
    from api.api_docs import register_resource
    register_resource(docs, ConsolidationAuditResource, 'consolidation_audit')