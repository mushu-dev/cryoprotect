#\!/usr/bin/env python3
"""
Schema verification module.

This module provides functions for verifying database schema consistency.
"""

import logging
import os
import json
from typing import Any, Dict, List, Optional, Set, Tuple
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

logger = logging.getLogger(__name__)

def _get_expected_tables() -> Dict[str, List[Dict]]:
    """
    Get expected tables and their columns.

    Returns:
        Dictionary mapping table names to column definitions
    """
    # This would ideally load from a schema definition file
    return {
        'molecules': [
            {'name': 'id', 'type': 'varchar', 'nullable': False},
            {'name': 'name', 'type': 'varchar', 'nullable': False},
            {'name': 'properties', 'type': 'jsonb', 'nullable': True},
            {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
            {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
        ],
        'mixtures': [
            {'name': 'id', 'type': 'varchar', 'nullable': False},
            {'name': 'name', 'type': 'varchar', 'nullable': False},
            {'name': 'description', 'type': 'text', 'nullable': True},
            {'name': 'properties', 'type': 'jsonb', 'nullable': True},
            {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
            {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
        ],
        'mixture_components': [
            {'name': 'id', 'type': 'varchar', 'nullable': False},
            {'name': 'mixture_id', 'type': 'varchar', 'nullable': False},
            {'name': 'molecule_id', 'type': 'varchar', 'nullable': False},
            {'name': 'concentration', 'type': 'numeric', 'nullable': True},
            {'name': 'units', 'type': 'varchar', 'nullable': True},
            {'name': 'properties', 'type': 'jsonb', 'nullable': True},
            {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
            {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
        ],
        'molecular_properties': [
            {'name': 'id', 'type': 'varchar', 'nullable': False},
            {'name': 'molecule_id', 'type': 'varchar', 'nullable': False},
            {'name': 'property_name', 'type': 'varchar', 'nullable': False},
            {'name': 'property_value', 'type': 'varchar', 'nullable': True},
            {'name': 'property_type', 'type': 'varchar', 'nullable': True},
            {'name': 'source', 'type': 'varchar', 'nullable': True},
            {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
            {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
        ],
        'calculation_methods': [
            {'name': 'id', 'type': 'varchar', 'nullable': False},
            {'name': 'name', 'type': 'varchar', 'nullable': False},
            {'name': 'description', 'type': 'text', 'nullable': True},
            {'name': 'version', 'type': 'varchar', 'nullable': True},
            {'name': 'parameters', 'type': 'jsonb', 'nullable': True},
            {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
            {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
        ],
        'experiments': [
            {'name': 'id', 'type': 'varchar', 'nullable': False},
            {'name': 'name', 'type': 'varchar', 'nullable': False},
            {'name': 'description', 'type': 'text', 'nullable': True},
            {'name': 'mixture_id', 'type': 'varchar', 'nullable': True},
            {'name': 'protocol', 'type': 'jsonb', 'nullable': True},
            {'name': 'results', 'type': 'jsonb', 'nullable': True},
            {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
            {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
        ],
        'predictions': [
            {'name': 'id', 'type': 'varchar', 'nullable': False},
            {'name': 'name', 'type': 'varchar', 'nullable': False},
            {'name': 'description', 'type': 'text', 'nullable': True},
            {'name': 'mixture_id', 'type': 'varchar', 'nullable': True},
            {'name': 'calculation_method_id', 'type': 'varchar', 'nullable': True},
            {'name': 'results', 'type': 'jsonb', 'nullable': True},
            {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
            {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
        ]
    }

def _get_actual_tables(conn: Any) -> Dict[str, List[Dict]]:
    """
    Get actual tables and their columns from the database.

    Args:
        conn: Database connection or connection-like object with .sql() method

    Returns:
        Dictionary mapping table names to column definitions
    """
    tables = {}

    try:
        # Handle direct connection object
        if hasattr(conn, 'cursor'):
            return _get_actual_tables_direct(conn)
        # Handle MCP-like object with .sql() method
        elif hasattr(conn, 'sql'):
            return _get_actual_tables_mcp(conn)
        else:
            raise ValueError("Unsupported connection type")
    except Exception as e:
        logger.error(f"Error while getting actual tables: {str(e)}")
        raise

def _get_actual_tables_direct(conn) -> Dict[str, List[Dict]]:
    """Get actual tables using direct database connection."""
    tables = {}

    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Get list of tables
        table_query = """
            SELECT tablename FROM pg_tables
            WHERE schemaname='public'
            ORDER BY tablename
        """
        cursor.execute(table_query)
        table_results = cursor.fetchall()

        for table_row in table_results:
            table_name = table_row['tablename']

            # Get columns for this table
            column_query = """
                SELECT
                    column_name,
                    data_type,
                    is_nullable
                FROM information_schema.columns
                WHERE table_schema='public' AND table_name=%s
                ORDER BY ordinal_position
            """
            cursor.execute(column_query, (table_name,))
            column_results = cursor.fetchall()

            columns = []
            for column_row in column_results:
                columns.append({
                    'name': column_row['column_name'],
                    'type': column_row['data_type'],
                    'nullable': column_row['is_nullable'] == 'YES'
                })

            tables[table_name] = columns

    return tables

def _get_actual_tables_mcp(conn) -> Dict[str, List[Dict]]:
    """Get actual tables using MCP connection object."""
    tables = {}

    # Get list of tables
    table_query = """
        SELECT tablename FROM pg_tables
        WHERE schemaname='public'
        ORDER BY tablename
    """
    table_result = conn.sql(table_query).execute()

    for table_row in table_result.data:
        table_name = table_row['tablename']

        # Get columns for this table
        column_query = """
            SELECT
                column_name,
                data_type,
                is_nullable
            FROM information_schema.columns
            WHERE table_schema='public' AND table_name=$1
            ORDER BY ordinal_position
        """
        column_result = conn.sql(column_query, {"params": [table_name]}).execute()

        columns = []
        for column_row in column_result.data:
            columns.append({
                'name': column_row['column_name'],
                'type': column_row['data_type'],
                'nullable': column_row['is_nullable'] == 'YES'
            })

        tables[table_name] = columns

    return tables

def _verify_table_presence(
    expected_tables: Dict[str, List[Dict]],
    actual_tables: Dict[str, List[Dict]]
) -> List[Dict]:
    """
    Verify that all expected tables are present.

    Args:
        expected_tables: Expected tables and columns
        actual_tables: Actual tables and columns

    Returns:
        List of issues found
    """
    issues = []

    # Check for missing tables
    missing_tables = set(expected_tables.keys()) - set(actual_tables.keys())
    for table in missing_tables:
        issues.append({
            'type': 'missing_table',
            'severity': 'error',
            'message': f"Table '{table}' is missing"
        })

    return issues

def _verify_table_columns(
    expected_tables: Dict[str, List[Dict]],
    actual_tables: Dict[str, List[Dict]]
) -> List[Dict]:
    """
    Verify that tables have expected columns.

    Args:
        expected_tables: Expected tables and columns
        actual_tables: Actual tables and columns

    Returns:
        List of issues found
    """
    issues = []

    for table_name, expected_columns in expected_tables.items():
        if table_name not in actual_tables:
            continue  # Skip missing tables (already reported)

        actual_columns = actual_tables[table_name]
        actual_column_names = {col['name'] for col in actual_columns}

        # Check for missing columns
        for expected_column in expected_columns:
            column_name = expected_column['name']

            if column_name not in actual_column_names:
                issues.append({
                    'type': 'missing_column',
                    'severity': 'error',
                    'message': f"Column '{column_name}' is missing from table '{table_name}'"
                })
            else:
                # Check column properties
                actual_column = next(col for col in actual_columns if col['name'] == column_name)

                # Check type
                if expected_column['type'] not in actual_column['type']:
                    issues.append({
                        'type': 'column_type_mismatch',
                        'severity': 'error',
                        'message': (
                            f"Column '{column_name}' in table '{table_name}' "
                            f"has type '{actual_column['type']}', "
                            f"expected '{expected_column['type']}'"
                        )
                    })

                # Check nullability
                if expected_column['nullable'] != actual_column['nullable']:
                    severity = 'warning' if expected_column['nullable'] else 'error'
                    issues.append({
                        'type': 'column_nullability_mismatch',
                        'severity': severity,
                        'message': (
                            f"Column '{column_name}' in table '{table_name}' "
                            f"has nullability={actual_column['nullable']}, "
                            f"expected {expected_column['nullable']}"
                        )
                    })

    return issues

def run_schema_verification(conn: Any, level: str = 'standard') -> Dict:
    """
    Run schema verification checks.

    Args:
        conn: Database connection
        level: Verification level

    Returns:
        Verification results
    """
    logger.info(f"Running schema verification at {level} level")

    issues = []

    try:
        expected_tables = _get_expected_tables()
        actual_tables = _get_actual_tables(conn)

        # Verify table presence
        issues.extend(_verify_table_presence(expected_tables, actual_tables))

        # Verify table columns
        issues.extend(_verify_table_columns(expected_tables, actual_tables))

        # Additional checks for comprehensive level
        if level == 'comprehensive':
            # Additional checks here...
            pass

        success = not any(issue['severity'] == 'error' for issue in issues)

        return {
            'success': success,
            'issues': issues
        }
    except Exception as e:
        logger.error(f"Error during schema verification: {str(e)}")
        return {
            'success': False,
            'issues': [
                {
                    'type': 'verification_error',
                    'severity': 'error',
                    'message': f"Error during schema verification: {str(e)}"
                }
            ]
        }

def get_db_config() -> Dict:
    """
    Get database configuration from environment variables or Supabase URL.
    
    This function first tries to get database configuration from standard
    environment variables (DB_HOST, etc.). If those are not available,
    it attempts to extract the configuration from the Supabase URL.
    
    Returns:
        Dict containing database configuration parameters
    """
    # Try standard environment variables first
    db_host = os.environ.get('DB_HOST')
    db_name = os.environ.get('DB_NAME')
    db_user = os.environ.get('DB_USER')
    db_password = os.environ.get('DB_PASSWORD')
    db_port = os.environ.get('DB_PORT')
    
    # If standard environment variables are set, use them
    if db_host and db_name and db_user and db_password:
        return {
            'host': db_host,
            'dbname': db_name,
            'user': db_user,
            'password': db_password,
            'port': int(db_port) if db_port else 5432
        }
    
    # Otherwise, try to extract from Supabase URL
    supabase_url = os.environ.get('SUPABASE_URL')
    if supabase_url:
        try:
            # Extract project ID from Supabase URL
            # Format: https://[project_id].supabase.co
            project_id = supabase_url.split('//')[1].split('.')[0]
            
            return {
                'host': f"db.{project_id}.supabase.co",
                'dbname': 'postgres',  # Default database name for Supabase
                'user': 'postgres',    # Default user for Supabase
                'password': os.environ.get('SUPABASE_DB_PASSWORD', 'postgres'),
                'port': 5432           # Default PostgreSQL port
            }
        except (IndexError, AttributeError):
            logger.warning(f"Invalid SUPABASE_URL format: {supabase_url}")
    
    # Fall back to default values
    return {
        'host': 'localhost',
        'dbname': 'cryoprotect',
        'user': 'postgres',
        'password': 'postgres',
        'port': 5432
    }

def create_direct_connection():
    """Create a direct database connection."""
    config = get_db_config()
    try:
        conn = psycopg2.connect(
            host=config['host'],
            dbname=config['dbname'],
            user=config['user'],
            password=config['password'],
            port=config['port']
        )
        return conn
    except Exception as e:
        logger.error(f"Failed to connect to database: {str(e)}")
        raise

def main():
    """CLI entry point for schema verification."""
    import argparse
    import json
    from datetime import datetime
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
    )
    
    # Parse arguments
    parser = argparse.ArgumentParser(description='Verify database schema')
    parser.add_argument(
        '--level',
        choices=['basic', 'standard', 'comprehensive'],
        default='standard',
        help='Verification level'
    )
    parser.add_argument(
        '--output',
        help='Output file path for report'
    )
    
    args = parser.parse_args()
    
    try:
        # Load environment variables from .env file
        load_dotenv()
        
        # Create database connection
        conn = create_direct_connection()
        
        # Run verification
        results = run_schema_verification(conn, args.level)
        
        # Add timestamp to results
        results['timestamp'] = datetime.now().isoformat()
        
        # Save results if output file specified
        if args.output:
            os.makedirs(os.path.dirname(args.output), exist_ok=True)
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to {args.output}")
        
        # Display summary
        success = results.get('success', False)
        issues = results.get('issues', [])
        
        if success:
            logger.info("Schema verification PASSED")
        else:
            logger.error("Schema verification FAILED")
        
        for issue in issues:
            log_method = logger.warning if issue['severity'] == 'warning' else logger.error
            log_method(f"{issue['type']}: {issue['message']}")
        
        # Close connection
        conn.close()
        
        # Return exit code
        return 0 if success else 1
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())