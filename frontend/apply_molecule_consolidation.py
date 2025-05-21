#!/usr/bin/env python3
"""
Script to consolidate duplicate molecules in the CryoProtect database.
This script applies the SQL migration to clean up duplicate molecules
and ensure data consistency.
"""

import os
import sys
import argparse
import json
from datetime import datetime
import psycopg2
from psycopg2.extras import RealDictCursor

# Environment variables for database connection
DB_HOST = os.environ.get('DB_HOST')
DB_PORT = os.environ.get('DB_PORT', '5432')
DB_NAME = os.environ.get('DB_NAME')
DB_USER = os.environ.get('DB_USER')
DB_PASSWORD = os.environ.get('DB_PASSWORD')
DB_URL = os.environ.get('DATABASE_URL')

def create_db_connection():
    """Create a connection to the database."""
    try:
        if DB_URL:
            conn = psycopg2.connect(DB_URL)
        else:
            conn = psycopg2.connect(
                host=DB_HOST,
                port=DB_PORT,
                dbname=DB_NAME,
                user=DB_USER,
                password=DB_PASSWORD
            )
        return conn
    except Exception as e:
        print(f"Error connecting to the database: {e}")
        sys.exit(1)

def execute_migration(conn, sql_file):
    """Execute the SQL migration file."""
    try:
        with open(sql_file, 'r') as f:
            sql = f.read()
            
        with conn.cursor() as cur:
            cur.execute(sql)
            conn.commit()
            print(f"Successfully executed migration: {sql_file}")
            return True
    except Exception as e:
        conn.rollback()
        print(f"Error executing migration: {e}")
        return False

def get_duplicate_stats(conn):
    """Get statistics on duplicate molecules before and after consolidation."""
    results = {}
    
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            # Get stats before consolidation
            cur.execute("""
                WITH duplicate_check AS (
                  SELECT 
                    smiles,
                    COUNT(*) as count
                  FROM 
                    consolidated_molecules
                  WHERE 
                    smiles IS NOT NULL
                  GROUP BY 
                    smiles
                  HAVING 
                    COUNT(*) > 1
                )
                SELECT 
                  COUNT(*) as total_duplicates,
                  SUM(count) as total_duplicate_records,
                  MAX(count) as max_duplicates_for_one_molecule
                FROM 
                  duplicate_check;
            """)
            results['before'] = cur.fetchone()
            
            # Get stats on data completeness
            cur.execute("""
                SELECT 
                  COUNT(*) as total_molecules,
                  SUM(CASE WHEN smiles IS NULL THEN 1 ELSE 0 END) as missing_smiles,
                  SUM(CASE WHEN formula IS NULL THEN 1 ELSE 0 END) as missing_formula,
                  SUM(CASE WHEN molecular_weight IS NULL THEN 1 ELSE 0 END) as missing_weight,
                  SUM(CASE WHEN pubchem_cid IS NULL THEN 1 ELSE 0 END) as missing_pubchem_cid,
                  SUM(CASE WHEN is_consolidated IS NULL THEN 1 ELSE 0 END) as missing_consolidation_status
                FROM 
                  consolidated_molecules;
            """)
            results['data_completeness'] = cur.fetchone()
            
            # Generate a timestamp for the report
            results['timestamp'] = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            return results
    except Exception as e:
        print(f"Error getting stats: {e}")
        return None

def save_report(stats, filename_prefix="duplicate_consolidation_report"):
    """Save the statistics report to a file."""
    timestamp = stats.get('timestamp', datetime.now().strftime("%Y%m%d_%H%M%S"))
    filename = f"{filename_prefix}_{timestamp}.json"
    
    try:
        with open(filename, 'w') as f:
            json.dump(stats, f, indent=2)
        print(f"Report saved to {filename}")
    except Exception as e:
        print(f"Error saving report: {e}")

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Consolidate duplicate molecules in the CryoProtect database.")
    parser.add_argument("--no-execute", action="store_true", help="Check for duplicates but don't execute the migration")
    parser.add_argument("--report-only", action="store_true", help="Generate a report without executing the migration")
    args = parser.parse_args()
    
    conn = create_db_connection()
    
    # Get stats before consolidation
    before_stats = get_duplicate_stats(conn)
    
    if args.report_only:
        save_report(before_stats)
        print("Report generated successfully.")
        return
    
    # Execute the migration if not in no-execute mode
    if not args.no_execute:
        sql_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), 
            'migrations', 
            '025_consolidate_duplicate_molecules.sql'
        )
        success = execute_migration(conn, sql_file)
        
        if not success:
            print("Migration failed. Exiting.")
            return
        
        # Get stats after consolidation
        after_stats = get_duplicate_stats(conn)
        if after_stats:
            stats = {
                'before': before_stats.get('before', {}),
                'data_completeness_before': before_stats.get('data_completeness', {}),
                'data_completeness_after': after_stats.get('data_completeness', {}),
                'timestamp': datetime.now().strftime("%Y%m%d_%H%M%S")
            }
            save_report(stats)
    else:
        print("Checking for duplicates without executing the migration...")
        if before_stats:
            print(f"Found {before_stats['before']['total_duplicates']} molecules with duplicates, "
                  f"totaling {before_stats['before']['total_duplicate_records']} records.")
            print(f"The molecule with the most duplicates has {before_stats['before']['max_duplicates_for_one_molecule']} records.")
            save_report(before_stats, filename_prefix="duplicate_check_report")
    
    conn.close()

if __name__ == "__main__":
    main()