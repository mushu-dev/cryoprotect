"""
Database client factory for CryoProtect.

This module provides a factory function for creating the appropriate database
client (Supabase or Convex) based on configuration.
"""

import os
import logging

logger = logging.getLogger(__name__)

def get_db_client(force_convex=None):
    """
    Get the appropriate database client based on configuration.
    
    Args:
        force_convex (bool): If provided, forces use of Convex (True) or Supabase (False).
        
    Returns:
        Object: Either a Convex adapter or Supabase client.
    """
    # Determine whether to use Convex or Supabase
    use_convex_env = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
    should_use_convex = force_convex if force_convex is not None else use_convex_env
    
    if should_use_convex:
        logger.info("Using Convex for database operations")
        from database.convex_adapter import create_client
        return create_client()
    else:
        logger.info("Using Supabase for database operations")
        from supabase import create_client
        
        supabase_url = os.environ.get('SUPABASE_URL', '')
        supabase_key = os.environ.get('SUPABASE_KEY', '')
        
        if not supabase_url or not supabase_key:
            raise ValueError("Supabase URL and key are required when not using Convex")
        
        return create_client(supabase_url, supabase_key)