#!/usr/bin/env python3
"""
Database verification runner.

This script runs both schema and constraint verification and generates a combined report.
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
from dotenv import load_dotenv
import psycopg2
from psycopg2.extras import RealDictCursor

# Import verification modules
from .schema import run_schema_verification
from .constraints import run_constraint_verification

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
)
logger = logging.getLogger(__name__)

def get_supabase_db_config():
    """
    Get database configuration, preferring session pooler if DB_CONNECTION_MODE=pooler.
    """
    load_dotenv()
    db_mode = os.environ.get('DB_CONNECTION_MODE', 'direct').lower()
    if db_mode == 'pooler':
        return {
            'host': os.environ.get('POOLER_DB_HOST'),
            'dbname': os.environ.get('POOLER_DB_NAME', 'postgres'),
            'user': os.environ.get('POOLER_DB_USER'),
            'password': os.environ.get('POOLER_DB_PASSWORD'),
            'port': int(os.environ.get('POOLER_DB_PORT', 5432))
        }
    else:
        # Fallback to direct connection
        supabase_url = os.environ.get('SUPABASE_URL')
        if not supabase_url:
            raise ValueError("SUPABASE_URL must be set in .env file")
        try:
            project_id = supabase_url.split('//')[1].split('.')[0]
        except (IndexError, AttributeError):
            raise ValueError(f"Invalid SUPABASE_URL format: {supabase_url}")
        return {
            'host': f"db.{project_id}.supabase.co",
            'dbname': 'postgres',
            'user': 'postgres',
            'password': os.environ.get('SUPABASE_DB_PASSWORD', 'postgres'),
            'port': 5432
        }

def create_direct_connection():
    """Create a direct database connection using Supabase credentials."""
    config = get_supabase_db_config()
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

def run_verification(level='standard', output_path=None):
    """
    Run both schema and constraint verification.
    
    Args:
        level: Verification level (basic, standard, comprehensive)
        output_path: Path to save the verification report
        
    Returns:
        Combined verification results
    """
    try:
        # Create database connection
        conn = create_direct_connection()
        
        # Run schema verification
        logger.info(f"Running schema verification at {level} level")
        schema_results = run_schema_verification(conn, level)
        
        # Run constraint verification
        logger.info(f"Running constraint verification at {level} level")
        constraint_results = run_constraint_verification(conn, level)
        
        # Run JSON property column verification
        logger.info("Verifying JSON property columns")
        json_property_results = verify_json_properties(conn)
        
        # Combine results
        combined_results = {
            'timestamp': datetime.now().isoformat(),
            'level': level,
            'schema_verification': schema_results,
            'constraint_verification': constraint_results,
            'json_property_verification': json_property_results,
            'success': (
                schema_results.get('success', False) and 
                constraint_results.get('success', False) and
                json_property_results.get('success', False)
            )
        }
        
        # Save results if output path specified
        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w') as f:
                json.dump(combined_results, f, indent=2)
            logger.info(f"Results saved to {output_path}")
        
        # Display summary
        success = combined_results.get('success', False)
        if success:
            logger.info("Database verification PASSED")
        else:
            logger.error("Database verification FAILED")
        
        # Close connection
        conn.close()
        
        return combined_results
    
    except Exception as e:
        logger.error(f"Error during verification: {str(e)}")
        return {
            'timestamp': datetime.now().isoformat(),
            'level': level,
            'success': False,
            'error': str(e)
        }

def verify_json_properties(conn):
    """
    Verify JSON property columns in the molecules table.
    
    Args:
        conn: Database connection
        
    Returns:
        Verification results
    """
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            query = """
                SELECT 
                    count(*) as total_count,
                    count(properties->'pubchem') as pubchem_count,
                    count(properties->'chembl') as chembl_count,
                    count(properties->'rdkit') as rdkit_count
                FROM molecules;
            """
            cursor.execute(query)
            result = cursor.fetchone()
            
            # Check if any molecules exist
            if result['total_count'] == 0:
                return {
                    'success': False,
                    'issues': [{
                        'type': 'no_molecules',
                        'severity': 'warning',
                        'message': 'No molecules found in the database'
                    }]
                }
            
            # Check if property columns exist
            issues = []
            if result['pubchem_count'] == 0:
                issues.append({
                    'type': 'missing_property',
                    'severity': 'warning',
                    'message': 'No molecules with PubChem properties found'
                })
            
            if result['chembl_count'] == 0:
                issues.append({
                    'type': 'missing_property',
                    'severity': 'warning',
                    'message': 'No molecules with ChEMBL properties found'
                })
            
            if result['rdkit_count'] == 0:
                issues.append({
                    'type': 'missing_property',
                    'severity': 'warning',
                    'message': 'No molecules with RDKit properties found'
                })
            
            return {
                'success': len(issues) == 0,
                'issues': issues,
                'counts': result
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

def main():
    """CLI entry point for database verification."""
    parser = argparse.ArgumentParser(description='Verify database schema and constraints')
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