"""
CryoProtect Analyzer - Supabase Client

This module provides functions for creating and managing the Supabase client.
"""

import os
from supabase import create_client as supabase_create_client
from flask import current_app

def create_client():
    """
    Create a Supabase client.
    
    Returns:
        Supabase client
    """
    # Get Supabase URL and key from environment variables or app config
    url = os.environ.get('SUPABASE_URL') or current_app.config.get('SUPABASE_URL')
    key = os.environ.get('SUPABASE_KEY') or current_app.config.get('SUPABASE_KEY')
    
    if not url or not key:
        raise ValueError("Supabase URL and key must be provided in environment variables or app config")
    
    # Create and return the client
    return supabase_create_client(url, key)