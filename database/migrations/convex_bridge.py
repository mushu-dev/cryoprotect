"""
Convex migration bridge module.

This module provides integration between the existing Python-based migration system
and the new Convex TypeScript-based migration system.
"""

import os
import json
import subprocess
import logging
from typing import Any, Dict, List, Optional, Tuple

from database.db_factory import get_db_client
from database.convex_adapter import ConvexAdapter

# Configure logging
logger = logging.getLogger(__name__)

def _is_convex_enabled() -> bool:
    """
    Check if Convex is enabled.
    
    Returns:
        bool: True if Convex is enabled, False otherwise
    """
    use_convex_env = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
    return use_convex_env

def _run_convex_cli(command: str, args: Dict[str, Any] = None) -> Dict:
    """
    Run the Convex CLI for migrations.
    
    Args:
        command: Command to run
        args: Arguments to pass to the command
    
    Returns:
        Dict: Response from the command
    """
    convex_cli_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'convex', 'migrations', 'cli.ts'
    )
    
    # Build the command
    cmd = ['node', '--loader', 'ts-node/esm', convex_cli_path, command]
    
    # Add arguments
    if args:
        for key, value in args.items():
            if isinstance(value, bool):
                if value:
                    cmd.append(f'--{key}')
            else:
                cmd.append(f'--{key}={value}')
    
    logger.info(f"Running Convex CLI command: {' '.join(cmd)}")
    
    # Run the command
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"Error running Convex CLI: {stderr}")
            return {
                'success': False,
                'error': stderr
            }
        
        # Parse the output as JSON if possible
        try:
            if '{' in stdout:
                json_start = stdout.find('{')
                json_str = stdout[json_start:]
                return json.loads(json_str)
            else:
                return {
                    'success': True,
                    'message': stdout
                }
        except json.JSONDecodeError:
            return {
                'success': True,
                'message': stdout
            }
    except Exception as e:
        logger.error(f"Error running Convex CLI: {str(e)}")
        return {
            'success': False,
            'error': str(e)
        }

def get_convex_migration_status() -> Dict:
    """
    Get status of Convex migrations.
    
    Returns:
        Dict: Migration status information
    """
    if not _is_convex_enabled():
        return {
            'success': False,
            'error': 'Convex is not enabled'
        }
    
    return _run_convex_cli('status', {'format': 'json'})

def apply_convex_migrations(
    target_version: Optional[str] = None,
    dry_run: bool = False
) -> Dict:
    """
    Apply Convex migrations.
    
    Args:
        target_version: Target version to migrate to
        dry_run: If True, don't make any changes
    
    Returns:
        Dict: Migration results
    """
    if not _is_convex_enabled():
        return {
            'success': False,
            'error': 'Convex is not enabled'
        }
    
    args = {
        'dry-run': dry_run
    }
    
    if target_version:
        args['target'] = target_version
    
    return _run_convex_cli('apply', args)

def rollback_convex_migrations(
    target_version: Optional[str] = None,
    dry_run: bool = False
) -> Dict:
    """
    Roll back Convex migrations.
    
    Args:
        target_version: Target version to roll back to
        dry_run: If True, don't make any changes
    
    Returns:
        Dict: Migration results
    """
    if not _is_convex_enabled():
        return {
            'success': False,
            'error': 'Convex is not enabled'
        }
    
    args = {
        'dry-run': dry_run
    }
    
    if target_version:
        args['target'] = target_version
    
    return _run_convex_cli('rollback', args)

def create_convex_migration(name: str) -> Dict:
    """
    Create a new Convex migration.
    
    Args:
        name: Migration name
    
    Returns:
        Dict: Result of migration creation
    """
    if not _is_convex_enabled():
        return {
            'success': False,
            'error': 'Convex is not enabled'
        }
    
    return _run_convex_cli('create', {'name': name})

def sync_supabase_to_convex(
    table_name: str,
    limit: int = 1000,
    dry_run: bool = False
) -> Dict:
    """
    Sync data from Supabase to Convex for a given table.
    
    Args:
        table_name: Name of the table to sync
        limit: Maximum number of records to sync
        dry_run: If True, don't make any changes
    
    Returns:
        Dict: Sync results
    """
    if not _is_convex_enabled():
        return {
            'success': False,
            'error': 'Convex is not enabled'
        }
    
    try:
        # Get Supabase client
        supabase = get_db_client(force_convex=False)
        
        # Get Convex client
        convex = get_db_client(force_convex=True)
        
        # Get data from Supabase
        response = supabase.table(table_name).select('*').limit(limit).execute()
        
        if response.error:
            return {
                'success': False,
                'error': f"Error getting data from Supabase: {response.error}"
            }
        
        supabase_data = response.data
        logger.info(f"Got {len(supabase_data)} records from Supabase table '{table_name}'")
        
        if dry_run:
            return {
                'success': True,
                'message': f"Would sync {len(supabase_data)} records from Supabase to Convex (dry run)",
                'count': len(supabase_data)
            }
        
        # Insert data into Convex
        for record in supabase_data:
            response = convex.table(table_name).insert(record).execute()
            
            if response.error:
                logger.error(f"Error inserting record into Convex: {response.error}")
                return {
                    'success': False,
                    'error': f"Error inserting into Convex: {response.error}"
                }
        
        return {
            'success': True,
            'message': f"Synced {len(supabase_data)} records from Supabase to Convex",
            'count': len(supabase_data)
        }
    except Exception as e:
        logger.error(f"Error syncing Supabase to Convex: {str(e)}")
        return {
            'success': False,
            'error': str(e)
        }

def run_universal_migrations(
    target_version: Optional[str] = None,
    dry_run: bool = False,
    environment: str = 'development',
) -> Dict:
    """
    Run migrations on both Supabase and Convex if enabled.
    
    This function is the main entry point for running migrations
    that should apply to both databases.
    
    Args:
        target_version: Target version to migrate to
        dry_run: If True, don't make any changes
        environment: Target environment
    
    Returns:
        Dict: Migration results
    """
    from database.migrations.runner import apply_migrations
    
    results = {
        'supabase': None,
        'convex': None
    }
    
    # Always run Supabase migrations
    conn = get_db_client(force_convex=False)
    
    supabase_result = apply_migrations(
        conn=conn,
        target_version=target_version,
        environment=environment,
        dry_run=dry_run
    )
    
    results['supabase'] = {
        'success': True,
        'applied': supabase_result
    }
    
    # Run Convex migrations if enabled
    if _is_convex_enabled():
        convex_result = apply_convex_migrations(
            target_version=target_version,
            dry_run=dry_run
        )
        
        results['convex'] = convex_result
    
    return results