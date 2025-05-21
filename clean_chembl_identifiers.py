#!/usr/bin/env python3
"""
Clean ChEMBL source identifiers in molecules table.

This script removes the "vUnknown ID:" prefix from ChEMBL sources
and standardizes the data format.
"""

import os
import sys
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
import re
from datetime import datetime

def connect_to_db(db_params):
    """Connect to the database and return a connection."""
    try:
        conn = psycopg2.connect(**db_params)
        return conn
    except psycopg2.Error as e:
        print(f"Database connection error: {e}")
        sys.exit(1)

def get_chembl_source_patterns(conn):
    """Get unique patterns in ChEMBL source identifiers."""
    patterns = {}
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT DISTINCT data_source
            FROM molecules
            WHERE data_source LIKE 'ChEMBL%'
            ORDER BY data_source
        """)
        
        for row in cursor.fetchall():
            source = row['data_source']
            # Extract pattern (e.g., "vUnknown ID:", "ID:", etc.)
            match = re.search(r'ChEMBL\s+([^:]+):', source)
            pattern = match.group(1) if match else "None"
            
            if pattern not in patterns:
                patterns[pattern] = []
            
            patterns[pattern].append(source)
    
    return patterns

def count_chembl_sources(conn):
    """Count ChEMBL sources by pattern."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT
                CASE
                    WHEN data_source LIKE 'ChEMBL vUnknown ID:%' THEN 'vUnknown ID:'
                    WHEN data_source LIKE 'ChEMBL ID:%' THEN 'ID:'
                    ELSE 'Other'
                END AS pattern,
                COUNT(*) as count
            FROM molecules
            WHERE data_source LIKE 'ChEMBL%'
            GROUP BY pattern
            ORDER BY count DESC
        """)
        
        return cursor.fetchall()

def clean_chembl_identifiers(conn, dry_run=False):
    """Clean up ChEMBL source identifiers."""
    # Count before changes
    before_counts = count_chembl_sources(conn)
    print("\nCurrent ChEMBL source patterns:")
    for row in before_counts:
        print(f"  {row['pattern']}: {row['count']} records")
    
    # Get pattern examples
    patterns = get_chembl_source_patterns(conn)
    print("\nExample sources for each pattern:")
    for pattern, examples in patterns.items():
        print(f"  Pattern '{pattern}':")
        for i, example in enumerate(examples[:3]):  # Show up to 3 examples
            print(f"    - {example}")
        if len(examples) > 3:
            print(f"    - ... and {len(examples) - 3} more")
    
    # Confirm before proceeding
    if not dry_run:
        confirmation = input("\nProceed with cleaning ChEMBL identifiers? (y/n): ")
        if confirmation.lower() != 'y':
            print("Operation cancelled.")
            return
    
    # Update the data sources
    with conn.cursor() as cursor:
        # Fix "vUnknown ID:" pattern
        cursor.execute("""
            UPDATE molecules
            SET 
                data_source = REPLACE(data_source, 'ChEMBL vUnknown ID:', 'ChEMBL ID:'),
                updated_at = NOW()
            WHERE data_source LIKE 'ChEMBL vUnknown ID:%'
            RETURNING COUNT(*)
        """)
        unknown_id_count = cursor.fetchone()[0]
        
        # Standardize spacing in "ID:" pattern
        cursor.execute("""
            UPDATE molecules
            SET 
                data_source = REGEXP_REPLACE(data_source, 'ChEMBL\\s+ID:', 'ChEMBL ID:'),
                updated_at = NOW()
            WHERE data_source LIKE 'ChEMBL%ID:%'
               AND data_source <> REGEXP_REPLACE(data_source, 'ChEMBL\\s+ID:', 'ChEMBL ID:')
            RETURNING COUNT(*)
        """)
        spacing_count = cursor.fetchone()[0]
        
        if not dry_run:
            conn.commit()
            print(f"\nUpdated {unknown_id_count} records with 'vUnknown ID:' pattern")
            print(f"Standardized spacing in {spacing_count} records")
        else:
            conn.rollback()
            print(f"\n[DRY RUN] Would update {unknown_id_count} records with 'vUnknown ID:' pattern")
            print(f"[DRY RUN] Would standardize spacing in {spacing_count} records")
    
    # Count after changes (only meaningful if not dry run)
    if not dry_run:
        after_counts = count_chembl_sources(conn)
        print("\nUpdated ChEMBL source patterns:")
        for row in after_counts:
            print(f"  {row['pattern']}: {row['count']} records")

def create_audit_log(conn, dry_run=False):
    """Create an audit log entry for the ChEMBL identifier cleaning."""
    if dry_run:
        return
    
    with conn.cursor() as cursor:
        cursor.execute("""
            INSERT INTO scientific_data_audit (
                operation_type,
                table_name,
                operation_details,
                performed_by,
                performed_at,
                affected_rows,
                operation_status,
                operation_metadata
            ) VALUES (
                'UPDATE',
                'molecules',
                'Clean ChEMBL source identifiers',
                'database_remediation_script',
                NOW(),
                (SELECT COUNT(*) FROM molecules WHERE data_source LIKE 'ChEMBL ID:%'),
                'COMPLETED',
                jsonb_build_object(
                    'operation', 'clean_chembl_identifiers',
                    'timestamp', NOW()::text,
                    'details', 'Removed vUnknown ID: prefix and standardized ChEMBL identifiers'
                )
            )
            RETURNING id
        """)
        
        audit_id = cursor.fetchone()[0]
        conn.commit()
        
        print(f"\nCreated audit log entry with ID: {audit_id}")

def main():
    parser = argparse.ArgumentParser(description="Clean ChEMBL source identifiers")
    parser.add_argument("--host", default=os.environ.get("DB_HOST", "aws-0-us-east-1.pooler.supabase.com"), 
                        help="Database host")
    parser.add_argument("--port", default=os.environ.get("DB_PORT", "5432"), 
                        help="Database port")
    parser.add_argument("--dbname", default=os.environ.get("DB_NAME", "postgres"), 
                        help="Database name")
    parser.add_argument("--user", default=os.environ.get("DB_USER", "postgres.tsdlmynydfuypiugmkev"), 
                        help="Database user")
    parser.add_argument("--password", default=os.environ.get("DB_PASSWORD"), 
                        help="Database password")
    parser.add_argument("--dry-run", action="store_true", 
                        help="Show what would be done without making changes")
    
    args = parser.parse_args()
    
    # Prepare database connection parameters
    db_params = {
        'host': args.host,
        'port': args.port,
        'dbname': args.dbname,
        'user': args.user,
        'password': args.password,
        'sslmode': 'require'
    }
    
    # Ensure password is provided
    if not db_params['password']:
        print("Error: Database password is required. Set DB_PASSWORD environment variable or use --password.")
        sys.exit(1)
    
    print(f"Starting ChEMBL identifier cleaning process at {datetime.now().isoformat()}")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    
    # Connect to the database
    conn = connect_to_db(db_params)
    
    try:
        # Clean ChEMBL identifiers
        clean_chembl_identifiers(conn, args.dry_run)
        
        # Create audit log entry
        create_audit_log(conn, args.dry_run)
        
    except Exception as e:
        print(f"Error during ChEMBL identifier cleaning: {e}")
        conn.rollback()
        sys.exit(1)
    
    finally:
        # Close the connection
        if conn:
            conn.close()
    
    print(f"ChEMBL identifier cleaning process completed at {datetime.now().isoformat()}")

if __name__ == "__main__":
    main()