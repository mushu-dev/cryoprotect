"""
Data quality verification module.

This module provides functions for verifying data quality.
"""

import logging
from typing import Any, Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

def _check_required_fields(conn: Any) -> List[Dict]:
    """
    Check for missing required fields.

    Args:
        conn: Database connection

    Returns:
        List of issues found
    """
    issues = []

    # Check for missing required fields in molecules
    query = """
        SELECT COUNT(*) AS count
        FROM molecules
        WHERE name IS NULL OR name = ''
    """
    result = conn.sql(query).execute()
    missing_count = result.data[0]['count'] if result.data else 0

    if missing_count > 0:
        issues.append({
            'type': 'missing_required_field',
            'severity': 'error',
            'message': f"Found {missing_count} molecule records with missing name"
        })

    # Check for missing required fields in mixtures
    query = """
        SELECT COUNT(*) AS count
        FROM mixtures
        WHERE name IS NULL OR name = ''
    """
    result = conn.sql(query).execute()
    missing_count = result.data[0]['count'] if result.data else 0

    if missing_count > 0:
        issues.append({
            'type': 'missing_required_field',
            'severity': 'error',
            'message': f"Found {missing_count} mixture records with missing name"
        })

    return issues

def _check_data_consistency(conn: Any) -> List[Dict]:
    """
    Check for data consistency issues.

    Args:
        conn: Database connection

    Returns:
        List of issues found
    """
    issues = []

    # Check for mixtures with no components
    query = """
        SELECT COUNT(*) AS count
        FROM mixtures m
        WHERE NOT EXISTS (
            SELECT 1 FROM mixture_components mc
            WHERE mc.mixture_id = m.id
        )
    """
    result = conn.sql(query).execute()
    empty_count = result.data[0]['count'] if result.data else 0

    if empty_count > 0:
        issues.append({
            'type': 'empty_mixture',
            'severity': 'warning',
            'message': f"Found {empty_count} mixtures with no components"
        })

    # Check for invalid concentration values
    query = """
        SELECT COUNT(*) AS count
        FROM mixture_components
        WHERE concentration < 0 OR concentration > 100
    """
    result = conn.sql(query).execute()
    invalid_count = result.data[0]['count'] if result.data else 0

    if invalid_count > 0:
        issues.append({
            'type': 'invalid_concentration',
            'severity': 'warning',
            'message': f"Found {invalid_count} mixture components with invalid concentration"
        })

    return issues

def _check_duplicate_records(conn: Any) -> List[Dict]:
    """
    Check for duplicate records.

    Args:
        conn: Database connection

    Returns:
        List of issues found
    """
    issues = []

    # Check for duplicate molecule names
    query = """
        SELECT name, COUNT(*) AS count
        FROM molecules
        GROUP BY name
        HAVING COUNT(*) > 1
    """
    result = conn.sql(query).execute()

    if result.data:
        total_duplicates = sum(row['count'] - 1 for row in result.data)
        issues.append({
            'type': 'duplicate_molecule_name',
            'severity': 'warning',
            'message': f"Found {total_duplicates} duplicate molecule names across {len(result.data)} groups"
        })

    # Check for duplicate mixture names
    query = """
        SELECT name, COUNT(*) AS count
        FROM mixtures
        GROUP BY name
        HAVING COUNT(*) > 1
    """
    result = conn.sql(query).execute()

    if result.data:
        total_duplicates = sum(row['count'] - 1 for row in result.data)
        issues.append({
            'type': 'duplicate_mixture_name',
            'severity': 'warning',
            'message': f"Found {total_duplicates} duplicate mixture names across {len(result.data)} groups"
        })

    return issues

def run_data_quality_verification(conn: Any, level: str = 'standard') -> Dict:
    """
    Run data quality verification checks.

    Args:
        conn: Database connection
        level: Verification level

    Returns:
        Verification results
    """
    logger.info(f"Running data quality verification at {level} level")

    issues = []

    try:
        # Basic level checks
        issues.extend(_check_required_fields(conn))

        # Standard and comprehensive level checks
        if level in ('standard', 'comprehensive'):
            issues.extend(_check_data_consistency(conn))

        # Comprehensive level checks
        if level == 'comprehensive':
            issues.extend(_check_duplicate_records(conn))

        success = not any(issue['severity'] == 'error' for issue in issues)

        return {
            'success': success,
            'issues': issues
        }
    except Exception as e:
        logger.error(f"Error during data quality verification: {str(e)}")
        return {
            'success': False,
            'issues': [
                {
                    'type': 'verification_error',
                    'severity': 'error',
                    'message': f"Error during data quality verification: {str(e)}"
                }
            ]
        }