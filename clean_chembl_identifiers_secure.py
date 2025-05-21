#!/usr/bin/env python3
"""
Clean ChEMBL source identifiers in molecules table.

This script removes the "vUnknown ID:" prefix from ChEMBL sources
and standardizes the data format.

This version uses environment variables for database connection,
improving security by not passing passwords on the command line.
"""

import os
import sys
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
import re
from datetime import datetime
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

def connect_to_db():
    """Connect to the database using environment variables."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    
    # Ensure required parameters are present
    if not all([db_params['host'], db_params['user'], db_params['password']]):
        print("Error: Missing required database connection parameters in environment variables.")
        print("Make sure you have a valid .env file with SUPABASE_DB_* variables.")
        sys.exit(1)
    
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

def clean_chembl_identifiers(conn, dry_run=False, auto_confirm=False):
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
    if not dry_run and not auto_confirm:
        confirmation = input("\nProceed with cleaning ChEMBL identifiers? (y/n): ")
        if confirmation.lower() != 'y':
            print("Operation cancelled.")
            return
    elif not dry_run and auto_confirm:
        print("\nAuto-confirming ChEMBL identifier cleaning (--yes flag provided)")
    
    # Update the data sources
    with conn.cursor() as cursor:
        # Count how many records will be affected by the "vUnknown ID:" pattern fix
        cursor.execute("""
            SELECT COUNT(*)
            FROM molecules
            WHERE data_source LIKE 'ChEMBL vUnknown ID:%'
        """)
        unknown_id_count = cursor.fetchone()[0]

        # Fix "vUnknown ID:" pattern if not in dry run mode
        if not dry_run:
            cursor.execute("""
                UPDATE molecules
                SET
                    data_source = REPLACE(data_source, 'ChEMBL vUnknown ID:', 'ChEMBL ID:'),
                    updated_at = NOW()
                WHERE data_source LIKE 'ChEMBL vUnknown ID:%'
            """)

        # Count how many records will be affected by spacing standardization
        cursor.execute("""
            SELECT COUNT(*)
            FROM molecules
            WHERE data_source LIKE 'ChEMBL%ID:%'
               AND data_source <> REGEXP_REPLACE(data_source, 'ChEMBL\\s+ID:', 'ChEMBL ID:')
        """)
        spacing_count = cursor.fetchone()[0]

        # Standardize spacing in "ID:" pattern if not in dry run mode
        if not dry_run:
            cursor.execute("""
                UPDATE molecules
                SET
                    data_source = REGEXP_REPLACE(data_source, 'ChEMBL\\s+ID:', 'ChEMBL ID:'),
                    updated_at = NOW()
                WHERE data_source LIKE 'ChEMBL%ID:%'
                   AND data_source <> REGEXP_REPLACE(data_source, 'ChEMBL\\s+ID:', 'ChEMBL ID:')
            """)
        
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
        print("\n[DRY RUN] Would create audit log entry")
        return

    try:
        with conn.cursor() as cursor:
            cursor.execute("""
                INSERT INTO scientific_data_audit (
                    table_name,
                    operation,
                    user_id,
                    timestamp,
                    old_data,
                    new_data,
                    application_context
                ) VALUES (
                    'molecules',
                    'UPDATE',
                    NULL,
                    NOW(),
                    jsonb_build_object(
                        'pattern', 'vUnknown ID:'
                    ),
                    jsonb_build_object(
                        'pattern', 'ID:'
                    ),
                    'database_remediation_script'
                )
                RETURNING id
            """)

            audit_id = cursor.fetchone()[0]
            conn.commit()

            print(f"\nCreated audit log entry with ID: {audit_id}")
    except Exception as e:
        print(f"Error creating audit log: {e}")
        print("The data changes were applied, but the audit log could not be created.")
        print("This might be due to incompatibility with the existing audit table structure.")

def main():
    parser = argparse.ArgumentParser(description="Clean ChEMBL source identifiers")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be done without making changes")
    parser.add_argument("--yes", action="store_true",
                        help="Automatically confirm changes without prompting")

    args = parser.parse_args()

    print(f"Starting ChEMBL identifier cleaning process at {datetime.now().isoformat()}")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    if args.yes:
        print("AUTO-CONFIRM MODE: Changes will be applied without confirmation")

    # Connect to the database using environment variables
    conn = connect_to_db()

    try:
        # Clean ChEMBL identifiers
        clean_chembl_identifiers(conn, args.dry_run, args.yes)

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