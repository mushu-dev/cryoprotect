"""
Main entry point for database population operations.

This module provides high-level functions for populating the database
with various types of data.
"""

import logging
import os
from typing import Dict, List, Optional, Union
from database.utils.connection import create_connection
from database.population.molecules import populate_molecules
from database.population.mixtures import populate_mixtures

logger = logging.getLogger(__name__)

def populate_all(
    environment: str = 'development',
    config: Optional[Dict] = None
) -> Dict[str, int]:
    """
    Populate all database tables based on environment.
    
    Args:
        environment: Target environment ('development', 'staging', 'production')
        config: Optional configuration dictionary
        
    Returns:
        Dictionary with table names as keys and count of inserted records as values
    """
    logger.info(f"Populating all tables in {environment} environment")
    
    conn = create_connection(config)
    results = {}
    
    # Populate in dependency order
    results['molecules'] = populate_molecules(conn, environment)
    results['mixtures'] = populate_mixtures(conn, environment)
    
    logger.info(f"Population complete. Results: {results}")
    return results
    
def populate_specific(
    tables: List[str],
    environment: str = 'development',
    config: Optional[Dict] = None
) -> Dict[str, int]:
    """
    Populate specific tables based on environment.
    
    Args:
        tables: List of tables to populate
        environment: Target environment ('development', 'staging', 'production')
        config: Optional configuration dictionary
        
    Returns:
        Dictionary with table names as keys and count of inserted records as values
    """
    logger.info(f"Populating specific tables {tables} in {environment} environment")
    
    conn = create_connection(config)
    results = {}
    
    if 'molecules' in tables:
        results['molecules'] = populate_molecules(conn, environment)
        
    if 'mixtures' in tables:
        results['mixtures'] = populate_mixtures(conn, environment)
    
    logger.info(f"Population complete. Results: {results}")
    return results
    
def populate_from_file(
    file_path: str, 
    table: str,
    environment: str = 'development',
    config: Optional[Dict] = None
) -> int:
    """
    Populate a table from a data file.
    
    Args:
        file_path: Path to the data file (JSON or CSV)
        table: Name of the table to populate
        environment: Target environment ('development', 'staging', 'production')
        config: Optional configuration dictionary
        
    Returns:
        Number of records inserted
    """
    logger.info(f"Populating {table} from file {file_path}")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return 0
        
    conn = create_connection(config)
    
    # Dispatch to appropriate handler based on table name
    if table == 'molecules':
        return populate_molecules(conn, environment, file_path)
    elif table == 'mixtures':
        return populate_mixtures(conn, environment, file_path)
    else:
        logger.error(f"Unsupported table: {table}")
        return 0
    
def main():
    """
    CLI entry point for database population.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Populate CryoProtect database tables')
    parser.add_argument(
        '--env', 
        choices=['development', 'staging', 'production'],
        default='development',
        help='Target environment'
    )
    parser.add_argument(
        '--tables', 
        nargs='+',
        help='Specific tables to populate'
    )
    parser.add_argument(
        '--file',
        help='Path to data file (for single table population)'
    )
    parser.add_argument(
        '--table',
        help='Table to populate from file (used with --file)'
    )
    
    args = parser.parse_args()
    
    if args.file and args.table:
        count = populate_from_file(args.file, args.table, args.env)
        print(f"Populated {args.table} with {count} records from {args.file}")
    elif args.tables:
        results = populate_specific(args.tables, args.env)
        for table, count in results.items():
            print(f"Populated {table} with {count} records")
    else:
        results = populate_all(args.env)
        for table, count in results.items():
            print(f"Populated {table} with {count} records")
    
if __name__ == '__main__':
    main()