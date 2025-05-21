#!/usr/bin/env python3
"""
Remove redundant cid field from molecules table.

This script:
1. Verifies that pubchem_cid has all the data from cid
2. Creates a backup of the cid values just in case
3. Removes the redundant cid column from molecules table
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

def connect_to_db():
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    return psycopg2.connect(**db_params)

def verify_cid_fields(conn):
    """Verify that pubchem_cid has all the data from cid."""
    with conn.cursor() as cursor:
        cursor.execute("""
            SELECT COUNT(*) 
            FROM molecules 
            WHERE pubchem_cid IS NULL AND cid IS NOT NULL
        """)
        missing_count = cursor.fetchone()[0]
        
        cursor.execute("""
            SELECT COUNT(*) 
            FROM molecules 
            WHERE pubchem_cid != cid AND pubchem_cid IS NOT NULL AND cid IS NOT NULL
        """)
        diff_count = cursor.fetchone()[0]
    
    return missing_count == 0 and diff_count == 0

def backup_cid_values(conn):
    """Create a backup of cid values."""
    backup_data = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, pubchem_cid, cid
            FROM molecules
            WHERE cid IS NOT NULL
        """)
        backup_data = cursor.fetchall()
    
    # Save backup to file
    with open("cid_backup.json", "w") as f:
        json.dump(backup_data, f, indent=2, default=str)
    
    return len(backup_data)

def remove_cid_column(conn):
    """Remove the redundant cid column."""
    try:
        with conn.cursor() as cursor:
            # First check if materialized view exists
            cursor.execute("""
                SELECT 1
                FROM pg_matviews
                WHERE matviewname = 'public_molecules_summary'
            """)

            if cursor.fetchone():
                print("Found dependent materialized view 'public_molecules_summary'")

                # Drop the materialized view
                cursor.execute("DROP MATERIALIZED VIEW IF EXISTS public_molecules_summary")
                print("Dropped materialized view")

                # Now drop the column
                cursor.execute("ALTER TABLE molecules DROP COLUMN IF EXISTS cid")
                print("Dropped cid column")

                # Recreate the materialized view without cid
                cursor.execute("""
                    CREATE MATERIALIZED VIEW public_molecules_summary AS
                    SELECT
                        m.id,
                        m.name,
                        m.molecular_formula,
                        m.smiles,
                        m.pubchem_cid,
                        m.pubchem_link,
                        m.is_public,
                        count(mp.id) AS property_count
                    FROM
                        molecules m
                        LEFT JOIN molecular_properties mp ON m.id = mp.molecule_id
                    WHERE
                        m.is_public = true
                    GROUP BY
                        m.id
                """)
                print("Recreated materialized view with pubchem_cid instead of cid")
            else:
                # No materialized view, just drop the column
                cursor.execute("ALTER TABLE molecules DROP COLUMN IF EXISTS cid")
                print("Dropped cid column")

            return True
    except Exception as e:
        print(f"Error removing cid column: {e}")
        return False

def main():
    """Remove redundant cid field from molecules table."""
    print("Removing redundant cid field from molecules table...")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Verify cid fields
        print("Verifying cid fields...")
        if not verify_cid_fields(conn):
            print("Error: pubchem_cid does not have all the data from cid.")
            print("Run fix_redundant_cid_fields.py first.")
            return
        
        # Create backup
        print("Creating backup of cid values...")
        backup_count = backup_cid_values(conn)
        print(f"Backed up {backup_count} cid values to cid_backup.json")
        
        # Remove column
        print("Removing cid column...")
        if remove_cid_column(conn):
            print("Successfully removed cid column")
        else:
            print("Failed to remove cid column")
            return
        
        # Commit changes
        conn.commit()
        print("Column removal committed to database")
        
        # Verify removal
        with conn.cursor() as cursor:
            cursor.execute("""
                SELECT column_name
                FROM information_schema.columns
                WHERE table_name = 'molecules' AND column_name = 'cid'
            """)
            if cursor.fetchone() is None:
                print("Verified: cid column has been removed")
            else:
                print("Warning: cid column still exists for some reason")
        
        print("\nRedundant cid field removal completed successfully")
    
    except Exception as e:
        conn.rollback()
        print(f"Error removing cid field: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()