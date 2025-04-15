"""
CryoProtect Analyzer API - Utilities

This module contains utility functions for the API, including Supabase connection
and authentication helpers.
"""

import os
from functools import wraps
from flask import request, g, current_app, jsonify
from supabase import create_client, Client

def get_supabase_client() -> Client:
    """
    Get or create a Supabase client.
    
    Returns:
        Client: A Supabase client instance
    """
    if not hasattr(g, 'supabase'):
        supabase_url = current_app.config['SUPABASE_URL']
        supabase_key = current_app.config['SUPABASE_KEY']
        
        if not supabase_url or not supabase_key:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in configuration")
        
        g.supabase = create_client(supabase_url, supabase_key)
    
    return g.supabase

def authenticate_user():
    """
    Authenticate a user with Supabase using credentials from configuration.
    
    Returns:
        dict: User data if authentication is successful
    """
    supabase = get_supabase_client()
    supabase_user = current_app.config['SUPABASE_USER']
    supabase_password = current_app.config['SUPABASE_PASSWORD']
    
    if supabase_user and supabase_password:
        try:
            response = supabase.auth.sign_in_with_password({
                "email": supabase_user,
                "password": supabase_password
            })
            
            if response.error:
                current_app.logger.warning(f"Authentication error: {response.error}")
                return None
            
            return response.user
        except Exception as e:
            current_app.logger.warning(f"Authentication error: {str(e)}")
            return None
    
    return None

def token_required(f):
    """
    Decorator to require a valid JWT token for API access.
    
    Args:
        f: The function to decorate
        
    Returns:
        function: The decorated function
    """
    @wraps(f)
    def decorated(*args, **kwargs):
        token = None
        
        # Get token from Authorization header
        if 'Authorization' in request.headers:
            auth_header = request.headers['Authorization']
            if auth_header.startswith('Bearer '):
                token = auth_header.split(' ')[1]
        
        if not token:
            return jsonify({'message': 'Authentication token is missing'}), 401
        
        try:
            # Verify token with Supabase
            supabase = get_supabase_client()
            response = supabase.auth.get_user(token)
            
            if response.error:
                return jsonify({'message': 'Invalid authentication token'}), 401
            
            # Store user info for the route to use
            g.user = response.user
            
        except Exception as e:
            return jsonify({'message': f'Authentication error: {str(e)}'}), 401
        
        return f(*args, **kwargs)
    
    return decorated

def handle_supabase_error(response):
    """
    Handle Supabase error responses.
    
    Args:
        response: Supabase response object
        
    Returns:
        tuple: (error_message, status_code)
    """
    if response.error:
        error_message = response.error.message if hasattr(response.error, 'message') else str(response.error)
        
        if 'not found' in error_message.lower():
            return error_message, 404
        elif 'permission' in error_message.lower() or 'access' in error_message.lower():
            return error_message, 403
        elif 'already exists' in error_message.lower() or 'duplicate' in error_message.lower():
            return error_message, 409
        else:
            return error_message, 400
    
    return None, None

def get_user_id():
    """
    Get the current user ID from the authenticated user.
    
    Returns:
        str: User ID or None if not authenticated
    """
    if hasattr(g, 'user') and g.user:
        return g.user.id
    
    return None