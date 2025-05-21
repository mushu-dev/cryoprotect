#!/usr/bin/env python3
"""
Fix redundant pubchem_cid/cid fields in molecules table.

This script identifies and consolidates redundant PubChem CID fields in the molecules table:
- Ensures pubchem_cid has the correct values
- Updates any records where pubchem_cid is NULL but cid has a value
- Documents the modifications for audit purposes
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

def analyze_cid_fields(conn):
    """Analyze the pubchem_cid and cid fields in molecules table."""
    stats = {
        "total_molecules": 0,
        "with_pubchem_cid": 0,
        "with_cid": 0,
        "with_both": 0,
        "with_matching_values": 0,
        "with_different_values": 0,
        "with_cid_only": 0,
        "with_pubchem_cid_only": 0
    }
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Get total count
        cursor.execute("SELECT COUNT(*) FROM molecules")
        stats["total_molecules"] = cursor.fetchone()["count"]
        
        # Count molecules with pubchem_cid
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NOT NULL")
        stats["with_pubchem_cid"] = cursor.fetchone()["count"]
        
        # Count molecules with cid
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE cid IS NOT NULL")
        stats["with_cid"] = cursor.fetchone()["count"]
        
        # Count molecules with both fields
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NOT NULL AND cid IS NOT NULL")
        stats["with_both"] = cursor.fetchone()["count"]
        
        # Count molecules with matching values
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid = cid AND pubchem_cid IS NOT NULL AND cid IS NOT NULL")
        stats["with_matching_values"] = cursor.fetchone()["count"]
        
        # Count molecules with different values
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid != cid AND pubchem_cid IS NOT NULL AND cid IS NOT NULL")
        stats["with_different_values"] = cursor.fetchone()["count"]
        
        # Count molecules with cid only
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NULL AND cid IS NOT NULL")
        stats["with_cid_only"] = cursor.fetchone()["count"]
        
        # Count molecules with pubchem_cid only
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NOT NULL AND cid IS NULL")
        stats["with_pubchem_cid_only"] = cursor.fetchone()["count"]
        
        # Get sample records with different values (if any)
        if stats["with_different_values"] > 0:
            cursor.execute("""
                SELECT id, name, pubchem_cid, cid, pubchem_link
                FROM molecules
                WHERE pubchem_cid != cid AND pubchem_cid IS NOT NULL AND cid IS NOT NULL
                LIMIT 5
            """)
            stats["different_value_examples"] = cursor.fetchall()
        
        # Get sample records with cid only
        if stats["with_cid_only"] > 0:
            cursor.execute("""
                SELECT id, name, pubchem_cid, cid, pubchem_link
                FROM molecules
                WHERE pubchem_cid IS NULL AND cid IS NOT NULL
                LIMIT 5
            """)
            stats["cid_only_examples"] = cursor.fetchall()
    
    return stats

def consolidate_cid_fields(conn, stats):
    """Consolidate pubchem_cid and cid fields."""
    updates = {
        "updated_records": 0,
        "details": []
    }

    with conn.cursor() as cursor:
        # Update records where pubchem_cid is NULL but cid has a value
        if stats["with_cid_only"] > 0:
            cursor.execute("""
                UPDATE molecules
                SET pubchem_cid = cid,
                    updated_at = NOW()
                WHERE pubchem_cid IS NULL AND cid IS NOT NULL
            """)
            updates["updated_records"] += cursor.rowcount
            print(f"Updated {cursor.rowcount} records where pubchem_cid was NULL")

            # Get details of updated records
            with conn.cursor(cursor_factory=RealDictCursor) as detail_cursor:
                detail_cursor.execute("""
                    SELECT id, name, pubchem_cid, cid, pubchem_link
                    FROM molecules
                    WHERE pubchem_cid = cid AND pubchem_cid IS NOT NULL
                    ORDER BY name
                    LIMIT 5
                """)
                updates["cid_to_pubchem_cid_examples"] = detail_cursor.fetchall()

        # Add modification records
        cursor.execute("""
            UPDATE molecules
            SET modification_history = COALESCE(modification_history, '[]'::jsonb) || jsonb_build_object(
                'timestamp', NOW(),
                'action', 'consolidate_cid_fields',
                'details', 'Consolidated pubchem_cid and cid fields'
            )::jsonb
            WHERE pubchem_cid IS NOT NULL
        """)

    return updates

def main():
    """Fix redundant pubchem_cid/cid fields in molecules table."""
    print("Analyzing and fixing redundant PubChem CID fields...")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Analyze current state
        stats = analyze_cid_fields(conn)
        print("\nCurrent State Analysis:")
        print(f"Total molecules: {stats['total_molecules']}")
        print(f"With pubchem_cid: {stats['with_pubchem_cid']}")
        print(f"With cid: {stats['with_cid']}")
        print(f"With both fields: {stats['with_both']}")
        print(f"With matching values: {stats['with_matching_values']}")
        print(f"With different values: {stats['with_different_values']}")
        print(f"With cid only: {stats['with_cid_only']}")
        print(f"With pubchem_cid only: {stats['with_pubchem_cid_only']}")
        
        # Show examples if there are records with different values
        if stats.get("different_value_examples"):
            print("\nExamples of records with different values:")
            for example in stats["different_value_examples"]:
                print(f"  Name: {example['name']}, pubchem_cid: {example['pubchem_cid']}, cid: {example['cid']}")
        
        # Show examples of records with cid only
        if stats.get("cid_only_examples"):
            print("\nExamples of records with cid only:")
            for example in stats["cid_only_examples"]:
                print(f"  Name: {example['name']}, cid: {example['cid']}")
        
        # Consolidate fields
        print("\nConsolidating fields...")
        updates = consolidate_cid_fields(conn, stats)
        
        # Commit changes
        conn.commit()
        
        # Reanalyze after changes
        print("\nReanalyzing after changes...")
        new_stats = analyze_cid_fields(conn)
        print("\nUpdated State Analysis:")
        print(f"Total molecules: {new_stats['total_molecules']}")
        print(f"With pubchem_cid: {new_stats['with_pubchem_cid']} (Before: {stats['with_pubchem_cid']})")
        print(f"With cid: {new_stats['with_cid']} (Before: {stats['with_cid']})")
        print(f"With both fields: {new_stats['with_both']} (Before: {stats['with_both']})")
        print(f"With cid only: {new_stats['with_cid_only']} (Before: {stats['with_cid_only']})")
        print(f"With pubchem_cid only: {new_stats['with_pubchem_cid_only']} (Before: {stats['with_pubchem_cid_only']})")
        
        # Save results to file
        results = {
            "before": stats,
            "updates": updates,
            "after": new_stats
        }
        with open("cid_field_consolidation_report.json", "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        print("\nCID field consolidation completed successfully")
        print("Report saved to 'cid_field_consolidation_report.json'")
    
    except Exception as e:
        conn.rollback()
        print(f"Error fixing CID fields: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()