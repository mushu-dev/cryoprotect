#!/usr/bin/env python3
"""
Simple adapter factory for the database connection system.

This module provides a simplified connection factory that works with the new
direct connection system while maintaining compatibility with code
expecting the adapter pattern.
"""

import os
import logging
import json
from typing import Dict, Any, Optional
from pathlib import Path

# Configure logger
logger = logging.getLogger(__name__)

class ConnectionFactory:
    """
    Simplified connection factory that works with the direct connection system.
    
    This class implements a minimal subset of the ConnectionFactory interface
    to work with the new direct connection system.
    """
    
    def __init__(self, config=None):
        """
        Initialize the connection factory.
        
        Args:
            config: Optional configuration dictionary
        """
        self.config = config or {}
        
        # Initialize the connection pool if needed
        try:
            from database import db
            if not hasattr(db, '_pool') or db._pool is None:
                logger.info("Initializing connection pool from adapter factory")
                db.init_connection_pool(config=self.get_connection_params())
        except Exception as e:
            logger.error(f"Error initializing connection pool: {e}")
    
    def get_connection_params(self) -> Dict[str, Any]:
        """
        Get connection parameters from configuration.
        
        Returns:
            Dict of connection parameters
        """
        # If config is empty, try to load from config file
        if not self.config:
            try:
                config_path = Path('config/db_config.json')
                if config_path.exists():
                    with open(config_path, 'r') as f:
                        return json.load(f)
            except Exception as e:
                logger.error(f"Error loading config file: {e}")
        
        # Extract connection parameters from supabase section if available
        if 'database' in self.config and 'connection' in self.config['database']:
            connection = self.config['database']['connection']
            if 'supabase' in connection and connection.get('mode') == 'supabase':
                return connection['supabase']
        
        # Return the config as is (might be already flattened)
        return self.config
    
    def get_connection(self):
        """
        Get a database connection.
        
        Returns:
            Database connection
        """
        from database import db
        return db.get_connection()

# For backward compatibility
def get_db_connection_info() -> Dict[str, Any]:
    """
    Get information about the current database connection.
    
    Returns:
        Dict with connection information
    """
    from database import db
    
    try:
        conn = db.get_connection()
        if not conn:
            return {'error': 'Could not get database connection'}
        
        try:
            with conn.cursor() as cursor:
                # Get server information
                cursor.execute("SELECT version() as version, current_database() as database, current_user as user")
                server_info = cursor.fetchone()
                
                # Get client encoding
                cursor.execute("SHOW client_encoding")
                encoding = cursor.fetchone()['client_encoding']
                
                return {
                    'server_version': server_info['version'],
                    'database': server_info['database'],
                    'user': server_info['user'],
                    'encoding': encoding,
                }
        finally:
            # Release the connection
            db.release_connection(conn)
            
    except Exception as e:
        logger.error(f"Error getting database connection info: {e}")
        return {'error': str(e)}
