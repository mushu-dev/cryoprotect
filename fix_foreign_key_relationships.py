"""
fix_foreign_key_relationships.py

Phase 2: Essential Maintenance Utilities
----------------------------------------
This script checks and fixes foreign key relationships for:
- lab_verifications.verifier (should be UUID referencing auth.users)
- shares.item_id and shared_resources.resource_id (checks for orphaned references)

Usage:
    python fix_foreign_key_relationships.py

Requires:
    - psycopg2 (install with pip if needed)
    - Database credentials set via environment variables or hardcoded below
"""

import os
import sys
import uuid
import psycopg2
from psycopg2.extras import RealDictCursor

# Database connection settings (edit as needed)
DB_HOST = os.getenv("DB_HOST", "localhost")
DB_PORT = os.getenv("DB_PORT", "5432")
DB_NAME = os.getenv("DB_NAME", "cryoprotect")
DB_USER = os.getenv("DB_USER", "postgres")
DB_PASSWORD = os.getenv("DB_PASSWORD", "password")

def connect():
    return psycopg2.connect(
        host=DB_HOST,
        port=DB_PORT,
        dbname=DB_NAME,
        user=DB_USER,
        password=DB_PASSWORD
    )

def is_uuid(val):
    try:
        uuid.UUID(str(val))
        return True
    except Exception:
        return False

def fix_lab_verifications_verifier(conn):
    print("Checking lab_verifications.verifier for non-UUID values...")
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute("SELECT id, verifier FROM public.lab_verifications;")
        rows = cur.fetchall()
        non_uuid = [r for r in rows if not is_uuid(r['verifier'])]
        if non_uuid:
            print(f"Found {len(non_uuid)} non-UUID verifier values. Manual migration required before running migration 017.")
            for r in non_uuid:
                print(f"  id={r['id']} verifier={r['verifier']}")
        else:
            print("All verifier values are valid UUIDs.")

def check_orphaned_shares(conn):
    print("Checking for orphaned shares.item_id references...")
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        # Check for molecules
        cur.execute("""
            SELECT s.id, s.item_id
            FROM public.shares s
            WHERE s.data_type = 'molecules'
              AND NOT EXISTS (SELECT 1 FROM public.molecules m WHERE m.id = s.item_id);
        """)
        orphaned = cur.fetchall()
        if orphaned:
            print(f"Found {len(orphaned)} orphaned shares for molecules:")
            for r in orphaned:
                print(f"  share_id={r['id']} item_id={r['item_id']}")
        else:
            print("No orphaned shares for molecules.")

def check_orphaned_shared_resources(conn):
    print("Checking for orphaned shared_resources.resource_id references...")
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        # Check for projects
        cur.execute("""
            SELECT sr.id, sr.resource_id
            FROM public.shared_resources sr
            WHERE sr.resource_type = 'project'
              AND NOT EXISTS (SELECT 1 FROM public.project p WHERE p.id = sr.resource_id);
        """)
        orphaned = cur.fetchall()
        if orphaned:
            print(f"Found {len(orphaned)} orphaned shared_resources for projects:")
            for r in orphaned:
                print(f"  shared_resource_id={r['id']} resource_id={r['resource_id']}")
        else:
            print("No orphaned shared_resources for projects.")

def main():
    try:
        conn = connect()
        fix_lab_verifications_verifier(conn)
        check_orphaned_shares(conn)
        check_orphaned_shared_resources(conn)
        conn.close()
        print("Foreign key relationship checks complete.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()