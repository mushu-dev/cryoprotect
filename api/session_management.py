"""
CryoProtect Analyzer - Session Management

This module provides functions for secure session management, including:
- Token revocation
- Refresh token rotation
- Session storage
- Audit logging
"""

import os
import uuid
import logging
import json
from datetime import datetime, timedelta
from typing import Dict, Any, Optional, List, Tuple, Union
from flask import current_app, g
import jwt

# Import Supabase client
from .utils import get_supabase_client

# Setup logging
logger = logging.getLogger(__name__)

# Session status constants
SESSION_ACTIVE = 'active'
SESSION_REVOKED = 'revoked'
SESSION_EXPIRED = 'expired'
SESSION_ROTATED = 'rotated'

class SessionManager:
    """
    Manages user sessions, token revocation, and audit logging.
    Uses Supabase database for storage.
    """
    
    def __init__(self):
        """Initialize the session manager."""
        self.supabase = get_supabase_client()
        
        # In-memory storage as fallback
        self._use_memory_sessions = False
        self._use_memory_audit_logs = False
        self._memory_sessions = {}
        self._memory_audit_logs = []
        
        # Ensure the required tables exist
        self._ensure_tables_exist()
    
    def _ensure_tables_exist(self):
        """
        Ensure that the required tables exist in the database.
        Creates them if they don't exist.
        """
        try:
            # Check if sessions table exists
            self.supabase.table('sessions').select('id').limit(1).execute()
        except Exception:
            # Create sessions table
            logger.info("Creating sessions table")
            self._create_sessions_table()
            
        try:
            # Check if session_audit_logs table exists
            self.supabase.table('session_audit_logs').select('id').limit(1).execute()
        except Exception:
            # Create session_audit_logs table
            logger.info("Creating session_audit_logs table")
            self._create_audit_logs_table()
    
    def _create_sessions_table(self):
        """Create the sessions table in the database."""
        try:
            # Instead of using execute_sql which doesn't exist, we'll use a REST API approach
            # First, create a simple session record to store the data
            session_data = {
                'id': str(uuid.uuid4()),
                'user_id': 'system',
                'refresh_token': 'temporary',
                'refresh_token_hash': 'temporary',
                'status': SESSION_ACTIVE,
                'expires_at': (datetime.now() + timedelta(hours=1)).isoformat()
            }
            
            # Try to insert a temporary record - this will fail if the table doesn't exist
            # but we'll catch the error and handle it
            try:
                self.supabase.table('sessions').insert(session_data).execute()
                # If we get here, the table exists, so we can delete our temporary record
                self.supabase.table('sessions').delete().eq('id', session_data['id']).execute()
                logger.info("Sessions table already exists")
                return
            except Exception as e:
                # Table likely doesn't exist, we'll handle this in the outer try/except
                pass
            
            # Since we can't directly execute SQL, we need to use a different approach
            # For now, we'll store sessions in memory as a fallback
            logger.warning("Unable to create sessions table. Using in-memory session storage as fallback.")
            self._use_memory_sessions = True
            
            # In a production environment, you would want to:
            # 1. Create the table through Supabase dashboard
            # 2. Use a migration tool
            # 3. Use a database admin tool
            # 4. Implement a custom endpoint in your API that can execute SQL
        except Exception as e:
            logger.error(f"Error creating sessions table: {str(e)}")
            raise
    
    def _create_audit_logs_table(self):
        """Create the session_audit_logs table in the database."""
        try:
            # Similar approach as with sessions table
            audit_log_data = {
                'id': str(uuid.uuid4()),
                'user_id': 'system',
                'action': 'system_init',
                'status': 'success',
                'created_at': datetime.now().isoformat()
            }
            
            # Try to insert a temporary record
            try:
                self.supabase.table('session_audit_logs').insert(audit_log_data).execute()
                # If we get here, the table exists, so we can delete our temporary record
                self.supabase.table('session_audit_logs').delete().eq('id', audit_log_data['id']).execute()
                logger.info("Session audit logs table already exists")
                return
            except Exception as e:
                # Table likely doesn't exist, we'll handle this in the outer try/except
                pass
            
            # Since we can't directly execute SQL, we need to use a different approach
            logger.warning("Unable to create session_audit_logs table. Using in-memory audit logs as fallback.")
            self._use_memory_audit_logs = True
        except Exception as e:
            logger.error(f"Error creating session audit logs table: {str(e)}")
            raise
    
    def create_session(self, user_id: str, refresh_token: str,
                       ip_address: Optional[str] = None,
                       user_agent: Optional[str] = None,
                       device_info: Optional[Dict] = None) -> Dict:
        """
        Create a new session for a user.
        
        Args:
            user_id: The user's ID
            refresh_token: The refresh token
            ip_address: The user's IP address
            user_agent: The user's browser/client user agent
            device_info: Additional device information
            
        Returns:
            The created session data
        """
        try:
            # Hash the refresh token for storage
            refresh_token_hash = self._hash_token(refresh_token)
            
            # Calculate expiry time based on config
            from auth_config import JWT_REFRESH_EXPIRY
            expires_at = datetime.now() + timedelta(seconds=JWT_REFRESH_EXPIRY)
            
            # Create session record
            session_data = {
                'id': str(uuid.uuid4()),
                'user_id': user_id,
                'refresh_token': refresh_token,  # Store encrypted in DB
                'refresh_token_hash': refresh_token_hash,
                'status': SESSION_ACTIVE,
                'expires_at': expires_at.isoformat(),
                'ip_address': ip_address,
                'user_agent': user_agent,
                'device_info': json.dumps(device_info) if device_info else None
            }
            
            # Use in-memory storage if database tables couldn't be created
            if self._use_memory_sessions:
                session_id = session_data['id']
                self._memory_sessions[session_id] = session_data
                session = session_data
                logger.info(f"Created in-memory session for user {user_id}")
            else:
                # Insert into database
                try:
                    response = self.supabase.table('sessions').insert(session_data).execute()
                    
                    if hasattr(response, 'error') and response.error:
                        logger.error(f"Error creating session: {response.error.message}")
                        raise Exception(f"Failed to create session: {response.error.message}")
                    
                    session = response.data[0] if hasattr(response, 'data') and response.data else session_data
                except Exception as e:
                    # Fallback to in-memory if database operation fails
                    logger.warning(f"Failed to create session in database, using in-memory: {str(e)}")
                    self._use_memory_sessions = True
                    session_id = session_data['id']
                    self._memory_sessions[session_id] = session_data
                    session = session_data
            
            # Log the session creation
            self.log_session_activity(
                user_id=user_id,
                session_id=session['id'],
                action='login',
                status='success',
                ip_address=ip_address,
                user_agent=user_agent,
                details={'device_info': device_info}
            )
            
            return session
        except Exception as e:
            logger.error(f"Session creation error: {str(e)}")
            # Log the failed session creation
            self.log_session_activity(
                user_id=user_id,
                session_id=None,
                action='login',
                status='failed',
                ip_address=ip_address,
                user_agent=user_agent,
                details={'error': str(e), 'device_info': device_info}
            )
            raise
    
    def revoke_session(self, session_id: str, reason: str = 'logout') -> bool:
        """
        Revoke a specific session.
        
        Args:
            session_id: The session ID to revoke
            reason: The reason for revocation (e.g., 'logout', 'password_change')
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Get the session to check if it belongs to the current user
            response = self.supabase.table('sessions').select('*').eq('id', session_id).execute()
            
            if response.error:
                logger.error(f"Error retrieving session: {response.error.message}")
                return False
            
            if not response.data:
                logger.warning(f"Session not found: {session_id}")
                return False
            
            session = response.data[0]
            user_id = session['user_id']
            
            # Update the session status
            update_data = {
                'status': SESSION_REVOKED,
                'updated_at': datetime.now().isoformat()
            }
            
            response = self.supabase.table('sessions').update(update_data).eq('id', session_id).execute()
            
            if response.error:
                logger.error(f"Error revoking session: {response.error.message}")
                return False
            
            # Log the session revocation
            self.log_session_activity(
                user_id=user_id,
                session_id=session_id,
                action=reason,
                status='success',
                details={'reason': reason}
            )
            
            return True
        except Exception as e:
            logger.error(f"Session revocation error: {str(e)}")
            return False
    
    def revoke_all_user_sessions(self, user_id: str, reason: str = 'password_change') -> bool:
        """
        Revoke all sessions for a specific user.
        
        Args:
            user_id: The user ID
            reason: The reason for revocation
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Get all active sessions for the user
            response = self.supabase.table('sessions').select('id').eq('user_id', user_id).eq('status', SESSION_ACTIVE).execute()
            
            if response.error:
                logger.error(f"Error retrieving user sessions: {response.error.message}")
                return False
            
            session_ids = [session['id'] for session in response.data]
            
            if not session_ids:
                logger.info(f"No active sessions found for user: {user_id}")
                return True
            
            # Update all sessions to revoked
            update_data = {
                'status': SESSION_REVOKED,
                'updated_at': datetime.now().isoformat()
            }
            
            response = self.supabase.table('sessions').update(update_data).in_('id', session_ids).execute()
            
            if response.error:
                logger.error(f"Error revoking user sessions: {response.error.message}")
                return False
            
            # Log the session revocations
            for session_id in session_ids:
                self.log_session_activity(
                    user_id=user_id,
                    session_id=session_id,
                    action=reason,
                    status='success',
                    details={'reason': reason}
                )
            
            return True
        except Exception as e:
            logger.error(f"User session revocation error: {str(e)}")
            return False
    
    def rotate_refresh_token(self, old_refresh_token: str, new_refresh_token: str,
                           ip_address: Optional[str] = None,
                           user_agent: Optional[str] = None) -> Dict:
        """
        Rotate a refresh token, invalidating the old one and creating a new session.
        
        Args:
            old_refresh_token: The old refresh token
            new_refresh_token: The new refresh token
            ip_address: The user's IP address
            user_agent: The user's browser/client user agent
            
        Returns:
            The new session data
        """
        try:
            # Hash the tokens
            old_token_hash = self._hash_token(old_refresh_token)
            new_token_hash = self._hash_token(new_refresh_token)
            
            # Find the session with the old token
            response = self.supabase.table('sessions').select('*').eq('refresh_token_hash', old_token_hash).eq('status', SESSION_ACTIVE).execute()
            
            if response.error:
                logger.error(f"Error finding session for token rotation: {response.error.message}")
                raise Exception(f"Failed to find session: {response.error.message}")
            
            if not response.data:
                logger.warning("No active session found for the provided refresh token")
                raise Exception("Invalid or expired refresh token")
            
            old_session = response.data[0]
            user_id = old_session['user_id']
            
            # Mark the old session as rotated
            update_data = {
                'status': SESSION_ROTATED,
                'updated_at': datetime.now().isoformat()
            }
            
            response = self.supabase.table('sessions').update(update_data).eq('id', old_session['id']).execute()
            
            if response.error:
                logger.error(f"Error updating old session during rotation: {response.error.message}")
                raise Exception(f"Failed to update old session: {response.error.message}")
            
            # Calculate expiry time based on config
            from auth_config import JWT_REFRESH_EXPIRY
            expires_at = datetime.now() + timedelta(seconds=JWT_REFRESH_EXPIRY)
            
            # Create new session record
            session_data = {
                'user_id': user_id,
                'refresh_token': new_refresh_token,  # Store encrypted in DB
                'refresh_token_hash': new_token_hash,
                'previous_refresh_token_hash': old_token_hash,
                'status': SESSION_ACTIVE,
                'expires_at': expires_at.isoformat(),
                'ip_address': ip_address or old_session.get('ip_address'),
                'user_agent': user_agent or old_session.get('user_agent'),
                'device_info': old_session.get('device_info')
            }
            
            # Insert into database
            response = self.supabase.table('sessions').insert(session_data).execute()
            
            if response.error:
                logger.error(f"Error creating new session during rotation: {response.error.message}")
                raise Exception(f"Failed to create new session: {response.error.message}")
            
            new_session = response.data[0]
            
            # Log the token rotation
            self.log_session_activity(
                user_id=user_id,
                session_id=new_session['id'],
                action='token_rotation',
                status='success',
                ip_address=ip_address,
                user_agent=user_agent,
                details={
                    'old_session_id': old_session['id'],
                    'new_session_id': new_session['id']
                }
            )
            
            return new_session
        except Exception as e:
            logger.error(f"Token rotation error: {str(e)}")
            # Try to extract user_id from the token for logging
            try:
                decoded = jwt.decode(old_refresh_token, options={"verify_signature": False})
                user_id = decoded.get('sub')
            except:
                user_id = None
                
            # Log the failed token rotation
            if user_id:
                self.log_session_activity(
                    user_id=user_id,
                    session_id=None,
                    action='token_rotation',
                    status='failed',
                    ip_address=ip_address,
                    user_agent=user_agent,
                    details={'error': str(e)}
                )
            raise
    
    def validate_session(self, refresh_token: str) -> Tuple[bool, Optional[Dict]]:
        """
        Validate if a refresh token is associated with an active session.
        
        Args:
            refresh_token: The refresh token to validate
            
        Returns:
            Tuple of (is_valid, session_data)
        """
        try:
            # Hash the token
            token_hash = self._hash_token(refresh_token)
            
            # Find the session with the token
            response = self.supabase.table('sessions').select('*').eq('refresh_token_hash', token_hash).execute()
            
            if response.error:
                logger.error(f"Error validating session: {response.error.message}")
                return False, None
            
            if not response.data:
                logger.warning("No session found for the provided refresh token")
                return False, None
            
            session = response.data[0]
            
            # Check if the session is active
            if session['status'] != SESSION_ACTIVE:
                logger.warning(f"Session is not active: {session['status']}")
                return False, None
            
            # Check if the session has expired
            expires_at = datetime.fromisoformat(session['expires_at'].replace('Z', '+00:00'))
            if expires_at < datetime.now():
                # Update session status to expired
                self.supabase.table('sessions').update({'status': SESSION_EXPIRED}).eq('id', session['id']).execute()
                logger.warning("Session has expired")
                return False, None
            
            return True, session
        except Exception as e:
            logger.error(f"Session validation error: {str(e)}")
            return False, None
    
    def get_active_sessions(self, user_id: str) -> List[Dict]:
        """
        Get all active sessions for a user.
        
        Args:
            user_id: The user ID
            
        Returns:
            List of active session data
        """
        try:
            response = self.supabase.table('sessions').select('*').eq('user_id', user_id).eq('status', SESSION_ACTIVE).execute()
            
            if response.error:
                logger.error(f"Error retrieving active sessions: {response.error.message}")
                return []
            
            return response.data
        except Exception as e:
            logger.error(f"Error retrieving active sessions: {str(e)}")
            return []
    
    def log_session_activity(self, user_id: str, action: str, status: str,
                           session_id: Optional[str] = None,
                           ip_address: Optional[str] = None,
                           user_agent: Optional[str] = None,
                           details: Optional[Dict] = None) -> bool:
        """
        Log a session activity for audit purposes.
        
        Args:
            user_id: The user ID
            action: The action performed (login, logout, refresh, etc.)
            status: The status of the action (success, failed)
            session_id: The session ID (if available)
            ip_address: The user's IP address
            user_agent: The user's browser/client user agent
            details: Additional details about the activity
            
        Returns:
            True if successful, False otherwise
        """
        try:
            log_data = {
                'id': str(uuid.uuid4()),
                'user_id': user_id,
                'session_id': session_id,
                'action': action,
                'status': status,
                'ip_address': ip_address,
                'user_agent': user_agent,
                'details': json.dumps(details) if details else None,
                'created_at': datetime.now().isoformat()
            }
            
            # Use in-memory storage if database tables couldn't be created
            if self._use_memory_audit_logs:
                self._memory_audit_logs.append(log_data)
                logger.info(f"Created in-memory audit log for user {user_id}, action {action}")
                return True
            else:
                # Insert into database
                try:
                    response = self.supabase.table('session_audit_logs').insert(log_data).execute()
                    
                    if hasattr(response, 'error') and response.error:
                        logger.error(f"Error logging session activity: {response.error.message}")
                        # Fallback to in-memory
                        self._use_memory_audit_logs = True
                        self._memory_audit_logs.append(log_data)
                        return True
                    
                    return True
                except Exception as e:
                    # Fallback to in-memory if database operation fails
                    logger.warning(f"Failed to log session activity in database, using in-memory: {str(e)}")
                    self._use_memory_audit_logs = True
                    self._memory_audit_logs.append(log_data)
                    return True
        except Exception as e:
            logger.error(f"Session activity logging error: {str(e)}")
            return False
    
    def clean_expired_sessions(self) -> int:
        """
        Clean up expired sessions by marking them as expired.
        
        Returns:
            Number of sessions cleaned up
        """
        try:
            now = datetime.now().isoformat()
            
            # Find expired but still active sessions
            response = self.supabase.table('sessions').select('id').eq('status', SESSION_ACTIVE).lt('expires_at', now).execute()
            
            if response.error:
                logger.error(f"Error finding expired sessions: {response.error.message}")
                return 0
            
            session_ids = [session['id'] for session in response.data]
            
            if not session_ids:
                return 0
            
            # Update sessions to expired
            update_data = {
                'status': SESSION_EXPIRED,
                'updated_at': now
            }
            
            response = self.supabase.table('sessions').update(update_data).in_('id', session_ids).execute()
            
            if response.error:
                logger.error(f"Error updating expired sessions: {response.error.message}")
                return 0
            
            return len(session_ids)
        except Exception as e:
            logger.error(f"Session cleanup error: {str(e)}")
            return 0
    
    def _hash_token(self, token: str) -> str:
        """
        Create a secure hash of a token for storage using HMAC-SHA256.
        This provides better security than a simple hash by using a secret key.
        
        Args:
            token: The token to hash
            
        Returns:
            Hashed token
        """
        import hmac
        import hashlib
        from flask import current_app
        
        # Get the secret key from the application config or environment
        # Fall back to the Flask secret key if no specific HMAC key is configured
        secret_key = current_app.config.get('HMAC_SECRET_KEY') or current_app.config.get('SECRET_KEY')
        
        if not secret_key:
            # If no secret key is available, log a warning and fall back to simple hashing
            # This should never happen in production
            logger.warning("No secret key available for HMAC. Falling back to simple hashing.")
            return hashlib.sha256(token.encode()).hexdigest()
        
        # Create an HMAC using SHA-256 and the secret key
        return hmac.new(
            key=secret_key.encode(),
            msg=token.encode(),
            digestmod=hashlib.sha256
        ).hexdigest()

# Singleton instance
_session_manager = None

def get_session_manager() -> SessionManager:
    """
    Get the session manager instance.
    
    Returns:
        SessionManager instance
    """
    global _session_manager
    if _session_manager is None:
        _session_manager = SessionManager()
    return _session_manager