"""
Constraint verification module.

This module provides functions for verifying database constraint integrity.
"""

import logging
from typing import Any, Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

def _get_expected_foreign_keys() -> List[Dict]:
    """
    Get expected foreign key constraints.

    Returns:
        List of expected foreign key definitions
    """
    # This would ideally load from a schema definition file
    return [
        {
            'table': 'mixture_components',
            'column': 'mixture_id',
            'references_table': 'mixtures',
            'references_column': 'id'
        },
        {
            'table': 'mixture_components',
            'column': 'molecule_id',
            'references_table': 'molecules',
            'references_column': 'id'
        }
    ]

def _get_actual_foreign_keys(conn: Any) -> List[Dict]:
    """
    Get actual foreign key constraints from the database.

    Args:
        conn: Database connection

    Returns:
        List of actual foreign key definitions
    """
    query = """
        SELECT
            tc.table_name,
            kcu.column_name,
            ccu.table_name AS references_table,
            ccu.column_name AS references_column
        FROM information_schema.table_constraints AS tc
        JOIN information_schema.key_column_usage AS kcu
          ON tc.constraint_name = kcu.constraint_name
        JOIN information_schema.constraint_column_usage AS ccu
          ON ccu.constraint_name = tc.constraint_name
        WHERE tc.constraint_type = 'FOREIGN KEY'
        AND tc.table_schema = 'public'
    """
    
    # Handle direct connection object
    if hasattr(conn, 'cursor'):
        with conn.cursor() as cursor:
            cursor.execute(query)
            result = cursor.fetchall()
            
            foreign_keys = []
            for row in result:
                foreign_keys.append({
                    'table': row[0],  # table_name
                    'column': row[1],  # column_name
                    'references_table': row[2],  # references_table
                    'references_column': row[3]  # references_column
                })
            
            return foreign_keys
    # Handle MCP-like object with .sql() method
    elif hasattr(conn, 'sql'):
        result = conn.sql(query).execute()

        foreign_keys = []
        for row in result.data:
            foreign_keys.append({
                'table': row['table_name'],
                'column': row['column_name'],
                'references_table': row['references_table'],
                'references_column': row['references_column']
            })

        return foreign_keys
    else:
        raise ValueError("Unsupported connection type")

def _verify_foreign_keys(
    expected_foreign_keys: List[Dict],
    actual_foreign_keys: List[Dict]
) -> List[Dict]:
    """
    Verify that expected foreign keys are present.

    Args:
        expected_foreign_keys: Expected foreign key constraints
        actual_foreign_keys: Actual foreign key constraints

    Returns:
        List of issues found
    """
    issues = []

    for expected_fk in expected_foreign_keys:
        found = False
        for actual_fk in actual_foreign_keys:
            if (
                expected_fk['table'] == actual_fk['table']
                and expected_fk['column'] == actual_fk['column']
                and expected_fk['references_table'] == actual_fk['references_table']
                and expected_fk['references_column'] == actual_fk['references_column']
            ):
                found = True
                break

        if not found:
            issues.append({
                'type': 'missing_foreign_key',
                'severity': 'error',
                'message': (
                    f"Foreign key constraint missing: "
                    f"{expected_fk['table']}({expected_fk['column']}) "
                    f"-> {expected_fk['references_table']}({expected_fk['references_column']})"
                )
            })

    return issues

def _verify_integrity_constraints(conn: Any) -> List[Dict]:
    """
    Verify data integrity constraints.

    Args:
        conn: Database connection

    Returns:
        List of issues found
    """
    issues = []

    # Check for orphaned records in mixture_components
    query = """
        SELECT COUNT(*) AS count
        FROM mixture_components
        WHERE mixture_id NOT IN (SELECT id FROM mixtures)
    """
    
    # Handle direct connection object
    if hasattr(conn, 'cursor'):
        with conn.cursor() as cursor:
            cursor.execute(query)
            result = cursor.fetchone()
            orphaned_count = result[0] if result else 0
            
            if orphaned_count > 0:
                issues.append({
                    'type': 'orphaned_records',
                    'severity': 'error',
                    'message': f"Found {orphaned_count} orphaned records in mixture_components (mixture_id)"
                })
            
            # Check for orphaned records in mixture_components (molecule)
            query = """
                SELECT COUNT(*) AS count
                FROM mixture_components
                WHERE molecule_id NOT IN (SELECT id FROM molecules)
            """
            cursor.execute(query)
            result = cursor.fetchone()
            orphaned_count = result[0] if result else 0
            
            if orphaned_count > 0:
                issues.append({
                    'type': 'orphaned_records',
                    'severity': 'error',
                    'message': f"Found {orphaned_count} orphaned records in mixture_components (molecule_id)"
                })
    # Handle MCP-like object with .sql() method
    elif hasattr(conn, 'sql'):
        result = conn.sql(query).execute()
        orphaned_count = result.data[0]['count'] if result.data else 0

        if orphaned_count > 0:
            issues.append({
                'type': 'orphaned_records',
                'severity': 'error',
                'message': f"Found {orphaned_count} orphaned records in mixture_components (mixture_id)"
            })

        # Check for orphaned records in mixture_components (molecule)
        query = """
            SELECT COUNT(*) AS count
            FROM mixture_components
            WHERE molecule_id NOT IN (SELECT id FROM molecules)
        """
        result = conn.sql(query).execute()
        orphaned_count = result.data[0]['count'] if result.data else 0

        if orphaned_count > 0:
            issues.append({
                'type': 'orphaned_records',
                'severity': 'error',
                'message': f"Found {orphaned_count} orphaned records in mixture_components (molecule_id)"
            })
    else:
        raise ValueError("Unsupported connection type")

    return issues

def run_constraint_verification(conn: Any, level: str = 'standard') -> Dict:
    """
    Run constraint verification checks.

    Args:
        conn: Database connection
        level: Verification level

    Returns:
        Verification results
    """
    logger.info(f"Running constraint verification at {level} level")

    issues = []

    try:
        # Verify foreign keys
        expected_foreign_keys = _get_expected_foreign_keys()
        actual_foreign_keys = _get_actual_foreign_keys(conn)
        issues.extend(_verify_foreign_keys(expected_foreign_keys, actual_foreign_keys))

        # Verify integrity constraints (standard and comprehensive levels)
        if level in ('standard', 'comprehensive'):
            issues.extend(_verify_integrity_constraints(conn))

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
        logger.error(f"Error during constraint verification: {str(e)}")
        return {
            'success': False,
            'issues': [
                {
                    'type': 'verification_error',
                    'severity': 'error',
                    'message': f"Error during constraint verification: {str(e)}"
                }
            ]
        }

def main():
    """CLI entry point for constraint verification."""
    import argparse
    import os
    import sys
    import json
    from datetime import datetime
    from dotenv import load_dotenv
    import psycopg2
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
    )
    
    # Parse arguments
    parser = argparse.ArgumentParser(description='Verify database constraints')
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
    
    # Load environment variables from .env file
    load_dotenv()
    
    try:
        # Get database configuration
        db_config = {
            'host': os.environ.get('DB_HOST', 'localhost'),
            'dbname': os.environ.get('DB_NAME', 'cryoprotect'),
            'user': os.environ.get('DB_USER', 'postgres'),
            'password': os.environ.get('DB_PASSWORD', 'postgres'),
            'port': int(os.environ.get('DB_PORT', 5432))
        }
        
        # Create database connection
        conn = psycopg2.connect(
            host=db_config['host'],
            dbname=db_config['dbname'],
            user=db_config['user'],
            password=db_config['password'],
            port=db_config['port']
        )
        
        # Run verification
        results = run_constraint_verification(conn, args.level)
        
        # Save results if output file specified
        if args.output:
            os.makedirs(os.path.dirname(args.output), exist_ok=True)
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to {args.output}")
        
        # Display summary
        success = results.get('success', False)
        if success:
            logger.info("Constraint verification PASSED")
        else:
            logger.error("Constraint verification FAILED")
        
        for issue in results.get('issues', []):
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