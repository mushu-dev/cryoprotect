#!/usr/bin/env python3
"""
Simple database verification script.

This script performs basic verification of the database schema
using the Supabase API and generates a simple report.
"""

import os
import sys
import json
import logging
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
)
logger = logging.getLogger(__name__)

def get_supabase_client():
    """Create a Supabase client using credentials from environment variables."""
    # Load environment variables from .env file
    load_dotenv()
    
    # Get Supabase URL and credentials
    supabase_url = os.environ.get('SUPABASE_URL')
    supabase_key = os.environ.get('SUPABASE_KEY')
    
    if not supabase_url or not supabase_key:
        raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
    
    # Create Supabase client
    return create_client(supabase_url, supabase_key)

def check_table_exists(supabase, table_name):
    """
    Check if a table exists by attempting to query it.
    
    Args:
        supabase: Supabase client
        table_name: Name of the table to check
        
    Returns:
        True if the table exists, False otherwise
    """
    try:
        # Try to select a single row from the table
        result = supabase.table(table_name).select('*').limit(1).execute()
        # If we get here, the table exists
        return True
    except Exception as e:
        # Check if the error message indicates the table doesn't exist
        error_str = str(e)
        if "relation" in error_str and "does not exist" in error_str:
            return False
        # For other errors, log and return False
        logger.error(f"Error checking if table {table_name} exists: {error_str}")
        return False

def verify_database_schema():
    """
    Verify the database schema by checking if required tables exist.
    
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
    
    try:
        # Create Supabase client
        supabase = get_supabase_client()
        
        # Check each required table
        table_status = {}
        for table in required_tables:
            exists = check_table_exists(supabase, table)
            table_status[table] = exists
            if exists:
                logger.info(f"Table '{table}' exists")
            else:
                logger.warning(f"Table '{table}' does not exist")
        
        # Check if all required tables exist
        all_tables_exist = all(table_status.values())
        
        # Create verification results
        results = {
            'timestamp': datetime.now().isoformat(),
            'success': all_tables_exist,
            'table_status': table_status,
            'message': "All required tables exist" if all_tables_exist else "Some required tables are missing"
        }
        
        return results
    
    except Exception as e:
        logger.error(f"Error during database schema verification: {str(e)}")
        return {
            'timestamp': datetime.now().isoformat(),
            'success': False,
            'error': str(e),
            'message': f"Error during database schema verification: {str(e)}"
        }

def verify_json_properties(supabase):
    """
    Verify JSON property columns in the molecules table.
    
    Args:
        supabase: Supabase client
        
    Returns:
        Verification results
    """
    try:
        # Check if molecules table exists
        if not check_table_exists(supabase, 'molecules'):
            return {
                'success': False,
                'message': "Molecules table does not exist"
            }
        
        # Try to select a row with properties column
        try:
            result = supabase.table('molecules').select('properties').limit(1).execute()
            # If we get here, the properties column exists
            return {
                'success': True,
                'message': "Properties column exists in molecules table"
            }
        except Exception as e:
            error_str = str(e)
            if "column" in error_str and "does not exist" in error_str:
                return {
                    'success': False,
                    'message': "Properties column does not exist in molecules table"
                }
            # For other errors, log and return error
            logger.error(f"Error checking properties column: {error_str}")
            return {
                'success': False,
                'error': error_str,
                'message': f"Error checking properties column: {error_str}"
            }
    
    except Exception as e:
        logger.error(f"Error during JSON property verification: {str(e)}")
        return {
            'success': False,
            'error': str(e),
            'message': f"Error during JSON property verification: {str(e)}"
        }

def main():
    """CLI entry point for simple database verification."""
    try:
        # Create Supabase client
        supabase = get_supabase_client()
        
        # Verify database schema
        logger.info("Verifying database schema...")
        schema_results = verify_database_schema()
        
        # Verify JSON properties
        logger.info("Verifying JSON properties...")
        json_results = verify_json_properties(supabase)
        
        # Combine results
        results = {
            'timestamp': datetime.now().isoformat(),
            'schema_verification': schema_results,
            'json_property_verification': json_results,
            'success': schema_results.get('success', False) and json_results.get('success', False)
        }
        
        # Save results to file
        output_path = 'reports/schema_validation_report.json'
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Results saved to {output_path}")
        
        # Display summary
        if results.get('success', False):
            logger.info("Database verification PASSED")
        else:
            logger.error("Database verification FAILED")
            if not schema_results.get('success', False):
                logger.error(f"Schema verification: {schema_results.get('message', 'Failed')}")
            if not json_results.get('success', False):
                logger.error(f"JSON property verification: {json_results.get('message', 'Failed')}")
        
        # Return exit code
        return 0 if results.get('success', False) else 1
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())