"""
CryoProtect Analyzer - Session Management Routes

This module provides API routes for session management,
allowing users to view and manage their active sessions.
"""

import logging
from flask import Blueprint, jsonify, request, g
from typing import Dict, Any, List

from .jwt_auth import jwt_required, get_current_user
from .session_management import get_session_manager
from .session_utils import get_user_sessions, revoke_session_by_id, get_session_audit_logs

# Setup logging
logger = logging.getLogger(__name__)

# Create Blueprint
session_bp = Blueprint('sessions', __name__)

@session_bp.route('/sessions', methods=['GET'])
@jwt_required
def get_sessions():
    """
    Get all active sessions for the current user.
    
    Returns:
        JSON response with sessions data
    """
    try:
        user = get_current_user()
        if not user:
            return jsonify({'message': 'User not authenticated'}), 401
        
        user_id = user.get('id')
        sessions = get_user_sessions(user_id)
        
        # Remove sensitive information
        for session in sessions:
            if 'refresh_token' in session:
                del session['refresh_token']
            if 'refresh_token_hash' in session:
                del session['refresh_token_hash']
            if 'previous_refresh_token_hash' in session:
                del session['previous_refresh_token_hash']
        
        return jsonify({
            'message': 'Sessions retrieved successfully',
            'sessions': sessions
        }), 200
    except Exception as e:
        logger.error(f"Error retrieving sessions: {str(e)}")
        return jsonify({'message': f'Error retrieving sessions: {str(e)}'}), 500

@session_bp.route('/sessions/<session_id>', methods=['DELETE'])
@jwt_required
def revoke_session(session_id):
    """
    Revoke a specific session.
    
    Args:
        session_id: The ID of the session to revoke
        
    Returns:
        JSON response with success/failure message
    """
    try:
        user = get_current_user()
        if not user:
            return jsonify({'message': 'User not authenticated'}), 401
        
        user_id = user.get('id')
        
        # Check if trying to revoke current session
        current_refresh_token = request.cookies.get('refresh_token')
        if current_refresh_token:
            session_manager = get_session_manager()
            is_valid, current_session = session_manager.validate_session(current_refresh_token)
            
            if is_valid and current_session and current_session['id'] == session_id:
                return jsonify({
                    'message': 'Cannot revoke current session. Use logout instead.'
                }), 400
        
        # Revoke the session
        success = revoke_session_by_id(user_id, session_id)
        
        if success:
            return jsonify({
                'message': 'Session revoked successfully'
            }), 200
        else:
            return jsonify({
                'message': 'Failed to revoke session or session not found'
            }), 404
    except Exception as e:
        logger.error(f"Error revoking session: {str(e)}")
        return jsonify({'message': f'Error revoking session: {str(e)}'}), 500

@session_bp.route('/sessions/all', methods=['DELETE'])
@jwt_required
def revoke_all_sessions():
    """
    Revoke all sessions for the current user except the current one.
    
    Returns:
        JSON response with success/failure message
    """
    try:
        user = get_current_user()
        if not user:
            return jsonify({'message': 'User not authenticated'}), 401
        
        user_id = user.get('id')
        
        # Get current session ID to exclude it
        current_session_id = None
        current_refresh_token = request.cookies.get('refresh_token')
        
        if current_refresh_token:
            session_manager = get_session_manager()
            is_valid, current_session = session_manager.validate_session(current_refresh_token)
            
            if is_valid and current_session:
                current_session_id = current_session['id']
        
        # Get all active sessions
        sessions = get_user_sessions(user_id)
        
        # Revoke all sessions except the current one
        revoked_count = 0
        for session in sessions:
            if current_session_id and session['id'] == current_session_id:
                continue
            
            if revoke_session_by_id(user_id, session['id']):
                revoked_count += 1
        
        return jsonify({
            'message': f'Successfully revoked {revoked_count} sessions',
            'revoked_count': revoked_count
        }), 200
    except Exception as e:
        logger.error(f"Error revoking all sessions: {str(e)}")
        return jsonify({'message': f'Error revoking all sessions: {str(e)}'}), 500

@session_bp.route('/sessions/audit', methods=['GET'])
@jwt_required
def get_session_audit():
    """
    Get audit logs for the current user's sessions.
    
    Query parameters:
        limit: Maximum number of logs to return (default: 50)
        offset: Offset for pagination (default: 0)
        
    Returns:
        JSON response with audit logs
    """
    try:
        user = get_current_user()
        if not user:
            return jsonify({'message': 'User not authenticated'}), 401
        
        user_id = user.get('id')
        
        # Get pagination parameters
        limit = request.args.get('limit', 50, type=int)
        offset = request.args.get('offset', 0, type=int)
        
        # Get audit logs
        logs = get_session_audit_logs(user_id, limit, offset)
        
        return jsonify({
            'message': 'Audit logs retrieved successfully',
            'logs': logs,
            'count': len(logs),
            'limit': limit,
            'offset': offset
        }), 200
    except Exception as e:
        logger.error(f"Error retrieving audit logs: {str(e)}")
        return jsonify({'message': f'Error retrieving audit logs: {str(e)}'}), 500

def register_session_routes(app):
    """
    Register session routes with the Flask app.
    
    Args:
        app: Flask application instance
    """
    app.register_blueprint(session_bp, url_prefix='/api')