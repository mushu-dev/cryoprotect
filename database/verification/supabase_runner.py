#!/usr/bin/env python3
"""
Database verification runner using Supabase API.

This script runs schema and constraint verification using the Supabase API
instead of a direct database connection.
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client, Client

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
)
logger = logging.getLogger(__name__)

def get_supabase_client() -> Client:
    """
    Create a Supabase client using credentials from environment variables.
    
    Returns:
        Supabase client
    """
    # Load environment variables from .env file
    load_dotenv()
    
    # Get Supabase URL and credentials
    supabase_url = os.environ.get('SUPABASE_URL')
    supabase_key = os.environ.get('SUPABASE_KEY')
    
    if not supabase_url or not supabase_key:
        raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
    
    # Create Supabase client
    return create_client(supabase_url, supabase_key)

def verify_table_existence(supabase: Client) -> dict:
    """
    Verify that all required tables exist.
    
    Args:
        supabase: Supabase client
        
    Returns:
        Verification results
    """
    required_tables = [
        'molecules',
        'molecular_properties',
        'mixture_components',
        'mixtures',
        'calculation_methods',
        'experiments',
        'predictions'
    ]
    
    issues = []
    
    try:
        # Get list of tables using system schema query
        result = supabase.table('pg_tables').select('tablename').eq('schemaname', 'public').execute()
        
        # Extract table names
        actual_tables = [row['tablename'] for row in result.data]
        
        # Check for missing tables
        missing_tables = set(required_tables) - set(actual_tables)
        for table in missing_tables:
            issues.append({
                'type': 'missing_table',
                'severity': 'error',
                'message': f"Table '{table}' is missing"
            })
        
        return {
            'success': len(issues) == 0,
            'issues': issues,
            'tables': {
                'required': required_tables,
                'actual': actual_tables
            }
        }
    
    except Exception as e:
        logger.error(f"Error during table existence verification: {str(e)}")
        return {
            'success': False,
            'issues': [{
                'type': 'verification_error',
                'severity': 'error',
                'message': f"Error during table existence verification: {str(e)}"
            }]
        }

def verify_table_columns(supabase: Client) -> dict:
    """
    Verify that tables have the required columns.
    
    Args:
        supabase: Supabase client
        
    Returns:
        Verification results
    """
    expected_columns = {
        'molecules': ['id', 'name', 'properties', 'created_at', 'updated_at'],
        'mixtures': ['id', 'name', 'description', 'properties', 'created_at', 'updated_at'],
        'mixture_components': ['id', 'mixture_id', 'molecule_id', 'concentration', 'units', 'properties', 'created_at', 'updated_at'],
        'molecular_properties': ['id', 'molecule_id', 'property_name', 'property_value', 'property_type', 'source', 'created_at', 'updated_at'],
        'calculation_methods': ['id', 'name', 'description', 'version', 'parameters', 'created_at', 'updated_at'],
        'experiments': ['id', 'name', 'description', 'mixture_id', 'protocol', 'results', 'created_at', 'updated_at'],
        'predictions': ['id', 'name', 'description', 'mixture_id', 'calculation_method_id', 'results', 'created_at', 'updated_at']
    }
    
    issues = []
    
    try:
        # For each table, get its columns
        for table_name, expected_cols in expected_columns.items():
            # Query information_schema.columns to get column information
            result = supabase.rpc(
                'exec_sql', 
                {
                    'query': f"""
                        SELECT column_name
                        FROM information_schema.columns
                        WHERE table_schema='public' AND table_name='{table_name}'
                        ORDER BY ordinal_position
                    """
                }
            ).execute()
            
            # Check if table exists
            if not result.data or len(result.data) == 0:
                continue  # Skip missing tables (already reported)
            
            # Extract column names
            actual_cols = []
            for row in result.data:
                if isinstance(row, dict) and 'column_name' in row:
                    actual_cols.append(row['column_name'])
                elif isinstance(row, list) and len(row) > 0:
                    actual_cols.append(row[0])
            
            # Check for missing columns
            missing_cols = set(expected_cols) - set(actual_cols)
            for col in missing_cols:
                issues.append({
                    'type': 'missing_column',
                    'severity': 'error',
                    'message': f"Column '{col}' is missing from table '{table_name}'"
                })
        
        return {
            'success': len(issues) == 0,
            'issues': issues
        }
    
    except Exception as e:
        logger.error(f"Error during column verification: {str(e)}")
        return {
            'success': False,
            'issues': [{
                'type': 'verification_error',
                'severity': 'error',
                'message': f"Error during column verification: {str(e)}"
            }]
        }

def verify_foreign_keys(supabase: Client) -> dict:
    """
    Verify that required foreign key constraints exist.
    
    Args:
        supabase: Supabase client
        
    Returns:
        Verification results
    """
    expected_foreign_keys = [
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
    
    issues = []
    
    try:
        # Query information_schema to get foreign key constraints
        result = supabase.rpc(
            'exec_sql', 
            {
                'query': """
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
            }
        ).execute()
        
        # Extract foreign keys
        actual_foreign_keys = []
        for row in result.data:
            if isinstance(row, dict):
                actual_foreign_keys.append({
                    'table': row.get('table_name'),
                    'column': row.get('column_name'),
                    'references_table': row.get('references_table'),
                    'references_column': row.get('references_column')
                })
            elif isinstance(row, list) and len(row) >= 4:
                actual_foreign_keys.append({
                    'table': row[0],
                    'column': row[1],
                    'references_table': row[2],
                    'references_column': row[3]
                })
        
        # Check for missing foreign keys
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
        
        return {
            'success': len(issues) == 0,
            'issues': issues
        }
    
    except Exception as e:
        logger.error(f"Error during foreign key verification: {str(e)}")
        return {
            'success': False,
            'issues': [{
                'type': 'verification_error',
                'severity': 'error',
                'message': f"Error during foreign key verification: {str(e)}"
            }]
        }

def verify_integrity_constraints(supabase: Client) -> dict:
    """
    Verify data integrity constraints.
    
    Args:
        supabase: Supabase client
        
    Returns:
        Verification results
    """
    issues = []
    
    try:
        # Check for orphaned records in mixture_components (mixture_id)
        result = supabase.rpc(
            'exec_sql', 
            {
                'query': """
                    SELECT COUNT(*) AS count
                    FROM mixture_components
                    WHERE mixture_id NOT IN (SELECT id FROM mixtures)
                """
            }
        ).execute()
        
        # Extract count
        orphaned_count = 0
        if result.data and len(result.data) > 0:
            if isinstance(result.data[0], dict) and 'count' in result.data[0]:
                orphaned_count = result.data[0]['count']
            elif isinstance(result.data[0], list) and len(result.data[0]) > 0:
                orphaned_count = result.data[0][0]
        
        if orphaned_count > 0:
            issues.append({
                'type': 'orphaned_records',
                'severity': 'error',
                'message': f"Found {orphaned_count} orphaned records in mixture_components (mixture_id)"
            })
        
        # Check for orphaned records in mixture_components (molecule_id)
        result = supabase.rpc(
            'exec_sql', 
            {
                'query': """
                    SELECT COUNT(*) AS count
                    FROM mixture_components
                    WHERE molecule_id NOT IN (SELECT id FROM molecules)
                """
            }
        ).execute()
        
        # Extract count
        orphaned_count = 0
        if result.data and len(result.data) > 0:
            if isinstance(result.data[0], dict) and 'count' in result.data[0]:
                orphaned_count = result.data[0]['count']
            elif isinstance(result.data[0], list) and len(result.data[0]) > 0:
                orphaned_count = result.data[0][0]
        
        if orphaned_count > 0:
            issues.append({
                'type': 'orphaned_records',
                'severity': 'error',
                'message': f"Found {orphaned_count} orphaned records in mixture_components (molecule_id)"
            })
        
        return {
            'success': len(issues) == 0,
            'issues': issues
        }
    
    except Exception as e:
        logger.error(f"Error during integrity constraint verification: {str(e)}")
        return {
            'success': False,
            'issues': [{
                'type': 'verification_error',
                'severity': 'error',
                'message': f"Error during integrity constraint verification: {str(e)}"
            }]
        }

def verify_json_properties(supabase: Client) -> dict:
    """
    Verify JSON property columns in the molecules table.
    
    Args:
        supabase: Supabase client
        
    Returns:
        Verification results
    """
    try:
        # Query to check JSON property columns
        result = supabase.rpc(
            'exec_sql', 
            {
                'query': """
                    SELECT 
                        count(*) as total_count,
                        count(properties->'pubchem') as pubchem_count,
                        count(properties->'chembl') as chembl_count,
                        count(properties->'rdkit') as rdkit_count
                    FROM molecules;
                """
            }
        ).execute()
        
        # Extract counts
        counts = {}
        if result.data and len(result.data) > 0:
            if isinstance(result.data[0], dict):
                counts = result.data[0]
            elif isinstance(result.data[0], list) and len(result.data[0]) >= 4:
                counts = {
                    'total_count': result.data[0][0],
                    'pubchem_count': result.data[0][1],
                    'chembl_count': result.data[0][2],
                    'rdkit_count': result.data[0][3]
                }
        
        # Check if any molecules exist
        if counts.get('total_count', 0) == 0:
            return {
                'success': False,
                'issues': [{
                    'type': 'no_molecules',
                    'severity': 'warning',
                    'message': 'No molecules found in the database'
                }],
                'counts': counts
            }
        
        # Check if property columns exist
        issues = []
        if counts.get('pubchem_count', 0) == 0:
            issues.append({
                'type': 'missing_property',
                'severity': 'warning',
                'message': 'No molecules with PubChem properties found'
            })
        
        if counts.get('chembl_count', 0) == 0:
            issues.append({
                'type': 'missing_property',
                'severity': 'warning',
                'message': 'No molecules with ChEMBL properties found'
            })
        
        if counts.get('rdkit_count', 0) == 0:
            issues.append({
                'type': 'missing_property',
                'severity': 'warning',
                'message': 'No molecules with RDKit properties found'
            })
        
        return {
            'success': len(issues) == 0,
            'issues': issues,
            'counts': counts
        }
    
    except Exception as e:
        logger.error(f"Error during JSON property verification: {str(e)}")
        return {
            'success': False,
            'issues': [{
                'type': 'verification_error',
                'severity': 'error',
                'message': f"Error during JSON property verification: {str(e)}"
            }]
        }

def run_verification(level='standard', output_path=None):
    """
    Run database verification using Supabase API.
    
    Args:
        level: Verification level (basic, standard, comprehensive)
        output_path: Path to save the verification report
        
    Returns:
        Combined verification results
    """
    try:
        # Create Supabase client
        supabase = get_supabase_client()
        
        # Run table existence verification
        logger.info(f"Running table existence verification at {level} level")
        table_results = verify_table_existence(supabase)
        
        # Run column verification
        logger.info(f"Running column verification at {level} level")
        column_results = verify_table_columns(supabase)
        
        # Run foreign key verification
        logger.info(f"Running foreign key verification at {level} level")
        foreign_key_results = verify_foreign_keys(supabase)
        
        # Run integrity constraint verification (standard and comprehensive levels)
        integrity_results = {'success': True, 'issues': []}
        if level in ('standard', 'comprehensive'):
            logger.info(f"Running integrity constraint verification at {level} level")
            integrity_results = verify_integrity_constraints(supabase)
        
        # Run JSON property verification
        logger.info("Verifying JSON property columns")
        json_property_results = verify_json_properties(supabase)
        
        # Combine all issues
        all_issues = (
            table_results.get('issues', []) +
            column_results.get('issues', []) +
            foreign_key_results.get('issues', []) +
            integrity_results.get('issues', []) +
            json_property_results.get('issues', [])
        )
        
        # Determine overall success
        success = not any(issue.get('severity') == 'error' for issue in all_issues)
        
        # Combine results
        combined_results = {
            'timestamp': datetime.now().isoformat(),
            'level': level,
            'success': success,
            'issues': all_issues,
            'details': {
                'table_verification': table_results,
                'column_verification': column_results,
                'foreign_key_verification': foreign_key_results,
                'integrity_verification': integrity_results,
                'json_property_verification': json_property_results
            }
        }
        
        # Save results if output path specified
        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w') as f:
                json.dump(combined_results, f, indent=2)
            logger.info(f"Results saved to {output_path}")
        
        # Display summary
        if success:
            logger.info("Database verification PASSED")
        else:
            logger.error("Database verification FAILED")
        
        for issue in all_issues:
            log_method = logger.warning if issue.get('severity') == 'warning' else logger.error
            log_method(f"{issue.get('type')}: {issue.get('message')}")
        
        return combined_results
    
    except Exception as e:
        logger.error(f"Error during verification: {str(e)}")
        error_results = {
            'timestamp': datetime.now().isoformat(),
            'level': level,
            'success': False,
            'error': str(e),
            'issues': [{
                'type': 'verification_error',
                'severity': 'error',
                'message': f"Error during verification: {str(e)}"
            }]
        }
        
        # Save error results if output path specified
        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w') as f:
                json.dump(error_results, f, indent=2)
            logger.info(f"Error results saved to {output_path}")
        
        return error_results

def main():
    """CLI entry point for database verification."""
    parser = argparse.ArgumentParser(description='Verify database schema and constraints using Supabase API')
    parser.add_argument(
        '--level',
        choices=['basic', 'standard', 'comprehensive'],
        default='standard',
        help='Verification level'
    )
    parser.add_argument(
        '--output',
        default='reports/schema_validation_report.json',
        help='Output file path for report'
    )
    
    args = parser.parse_args()
    
    try:
        # Run verification
        results = run_verification(args.level, args.output)
        
        # Return exit code
        return 0 if results.get('success', False) else 1
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())