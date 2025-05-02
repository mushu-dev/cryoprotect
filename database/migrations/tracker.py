"""
Migration tracking module.

This module provides functions for tracking applied migrations.
"""

import logging
from typing import Dict, List, Any

# Configure logging
logger = logging.getLogger(__name__)

def get_applied_migrations(conn: Any) -> List[Dict]:
    """
    Get list of applied migrations.

    Args:
        conn: Database connection object

    Returns:
        List of applied migration dictionaries
    """
    try:
        result = conn.table('migrations').select('*').order('version').execute()
        return result.data if result.data else []
    except Exception as e:
        logger.error(f"Error getting applied migrations: {str(e)}")
        return []

def record_migration(
    conn: Any,
    version: str,
    name: str
) -> bool:
    """
    Record a migration as applied.

    Args:
        conn: Database connection object
        version: Migration version
        name: Migration name

    Returns:
        True if successful, False otherwise
    """
    try:
        conn.table('migrations').insert({
            'version': version,
            'name': name,
            'applied_at': 'now()'
        }).execute()
        return True
    except Exception as e:
        logger.error(f"Error recording migration {version}: {str(e)}")
        return False

def remove_migration_record(conn: Any, version: str) -> bool:
    """
    Remove a migration record.

    Args:
        conn: Database connection object
        version: Migration version

    Returns:
        True if successful, False otherwise
    """
    try:
        conn.table('migrations').delete().eq('version', version).execute()
        return True
    except Exception as e:
        logger.error(f"Error removing migration record {version}: {str(e)}")
        return False