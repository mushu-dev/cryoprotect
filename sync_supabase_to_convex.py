#!/usr/bin/env python3
"""
Utility script to sync data from Supabase to Convex.

This script provides a command-line interface for syncing data from Supabase
to Convex.
"""

import os
import sys
import argparse
import logging
import json
from typing import List, Dict, Any, Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add the project directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import our migration bridge
from database.migrations.convex_bridge import sync_supabase_to_convex

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Sync data from Supabase to Convex'
    )
    
    parser.add_argument(
        'tables',
        nargs='*',
        help='Tables to sync (default: sync all tables)'
    )
    
    parser.add_argument(
        '--limit',
        type=int,
        default=1000,
        help='Maximum number of records to sync per table (default: 1000)'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help="Don't actually sync data, just show what would happen"
    )
    
    parser.add_argument(
        '--format',
        choices=['text', 'json'],
        default='text',
        help='Output format (default: text)'
    )
    
    return parser.parse_args()

def get_all_tables() -> List[str]:
    """
    Get a list of all tables in the database.
    
    Returns:
        List[str]: List of table names
    """
    from database.db_factory import get_db_client
    
    tables = []
    
    try:
        # Get Supabase client
        supabase = get_db_client(force_convex=False)
        
        # Query the information_schema.tables view
        response = supabase.rpc(
            'exec_sql',
            {
                'query': """
                SELECT table_name
                FROM information_schema.tables
                WHERE table_schema = 'public'
                AND table_type = 'BASE TABLE'
                ORDER BY table_name
                """
            }
        ).execute()
        
        if response.error:
            logger.error(f"Error getting tables: {response.error}")
            return []
        
        if response.data:
            for row in response.data:
                if isinstance(row, dict) and 'table_name' in row:
                    tables.append(row['table_name'])
                elif isinstance(row, list) and len(row) > 0:
                    tables.append(row[0])
        
        return tables
    except Exception as e:
        logger.error(f"Error getting tables: {str(e)}")
        return []

def main():
    """Main entry point."""
    args = parse_args()
    
    # Check if we're using Convex
    if os.environ.get('USE_CONVEX', '').lower() not in ('true', 'yes', '1'):
        logger.error("Convex is not enabled. Please set USE_CONVEX=true.")
        sys.exit(1)
    
    # Get tables to sync
    tables = args.tables
    
    if not tables:
        logger.info("No tables specified, discovering all tables...")
        tables = get_all_tables()
        
        if not tables:
            logger.error("No tables found to sync")
            sys.exit(1)
    
    logger.info(f"Will sync the following tables: {', '.join(tables)}")
    
    if args.dry_run:
        logger.info("DRY RUN: No data will actually be synced")
    
    # Initialize results
    results = {
        'tables': {},
        'total_count': 0,
        'success_count': 0,
        'failure_count': 0
    }
    
    # Sync each table
    for table in tables:
        logger.info(f"Syncing table '{table}'...")
        
        result = sync_supabase_to_convex(
            table_name=table,
            limit=args.limit,
            dry_run=args.dry_run
        )
        
        results['tables'][table] = result
        
        if result.get('success', False):
            results['success_count'] += 1
            count = result.get('count', 0)
            results['total_count'] += count
            
            logger.info(f"Successfully synced {count} records from '{table}'")
        else:
            results['failure_count'] += 1
            error = result.get('error', 'Unknown error')
            
            logger.error(f"Error syncing '{table}': {error}")
    
    # Print summary
    if args.format == 'json':
        print(json.dumps(results, indent=2))
    else:
        print("\nSync Summary:")
        print(f"Synced {results['success_count']} tables with {results['total_count']} records total")
        
        if results['failure_count'] > 0:
            print(f"Failed to sync {results['failure_count']} tables")
            
            for table, result in results['tables'].items():
                if not result.get('success', False):
                    print(f"- {table}: {result.get('error', 'Unknown error')}")

if __name__ == '__main__':
    main()