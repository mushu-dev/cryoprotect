"""
CryoProtect Analyzer - Session Utilities

This module provides utility functions for session management,
including scheduled cleanup tasks and session verification.
"""

import logging
import threading
import time
from datetime import datetime, timedelta
from typing import Dict, Any, Optional, List, Tuple, Union
from flask import request, g

# Setup logging
logger = logging.getLogger(__name__)

def verify_session_status():
    """
    Verify the status of the current session.
    This function is called on protected routes to ensure the session is valid.
    
    Returns:
        Tuple of (is_valid, error_message)
    """
    from .session_management import get_session_manager
    
    # Get refresh token from cookies or request
    refresh_token = request.cookies.get('refresh_token')
    
    if not refresh_token:
        # No refresh token, can't verify session
        return True, None
    
    # Validate the session
    session_manager = get_session_manager()
    is_valid, session = session_manager.validate_session(refresh_token)
    
    if not is_valid:
        return False, "Session is invalid or expired"
    
    return True, None

def start_session_cleanup_thread(interval_seconds: int = 3600):
    """
    Start a background thread to clean up expired sessions.
    
    Args:
        interval_seconds: Interval between cleanup runs in seconds (default: 1 hour)
    """
    def cleanup_task():
        from .session_management import get_session_manager
        
        while True:
            try:
                logger.info("Running session cleanup task")
                session_manager = get_session_manager()
                cleaned = session_manager.clean_expired_sessions()
                logger.info(f"Cleaned up {cleaned} expired sessions")
            except Exception as e:
                logger.error(f"Error in session cleanup task: {str(e)}")
            
            # Sleep until next run
            time.sleep(interval_seconds)
    
    # Create and start the thread
    cleanup_thread = threading.Thread(target=cleanup_task, daemon=True)
    cleanup_thread.start()
    logger.info(f"Session cleanup thread started with interval {interval_seconds} seconds")

def get_user_sessions(user_id: str) -> List[Dict]:
    """
    Get all active sessions for a user.
    
    Args:
        user_id: The user ID
        
    Returns:
        List of active session data
    """
    from .session_management import get_session_manager
    session_manager = get_session_manager()
    return session_manager.get_active_sessions(user_id)

def revoke_session_by_id(user_id: str, session_id: str) -> bool:
    """
    Revoke a specific session for a user.
    
    Args:
        user_id: The user ID
        session_id: The session ID to revoke
        
    Returns:
        True if successful, False otherwise
    """
    from .session_management import get_session_manager
    session_manager = get_session_manager()
    
    # Get the session to verify ownership
    response = session_manager.supabase.table('sessions').select('*').eq('id', session_id).execute()
    
    if response.error or not response.data:
        logger.warning(f"Session not found or error: {session_id}")
        return False
    
    session = response.data[0]
    
    # Verify the session belongs to the user
    if session['user_id'] != user_id:
        logger.warning(f"Session {session_id} does not belong to user {user_id}")
        return False
    
    # Revoke the session
    return session_manager.revoke_session(session_id, reason='user_revoked')

def get_session_audit_logs(user_id: str, limit: int = 50, offset: int = 0) -> List[Dict]:
    """
    Get audit logs for a user's sessions.
    
    Args:
        user_id: The user ID
        limit: Maximum number of logs to return
        offset: Offset for pagination
        
    Returns:
        List of audit log entries
    """
    from .session_management import get_session_manager
    session_manager = get_session_manager()
    
    try:
        response = session_manager.supabase.table('session_audit_logs') \
            .select('*') \
            .eq('user_id', user_id) \
            .order('created_at', desc=True) \
            .limit(limit) \
            .offset(offset) \
            .execute()
        
        if response.error:
            logger.error(f"Error retrieving audit logs: {response.error.message}")
            return []
        
        return response.data
    except Exception as e:
        logger.error(f"Error retrieving audit logs: {str(e)}")
        return []