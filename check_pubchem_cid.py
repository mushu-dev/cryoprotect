#!/usr/bin/env python3
"""
Check for issues with pubchem_cid field in molecules table.

This script:
1. Identifies molecules with missing pubchem_cid values
2. Checks if there are any inconsistencies in pubchem_cid and pubchem_link values
3. Provides a summary of the analysis
"""

import os
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

def analyze_pubchem_cid_field(conn):
    """Analyze the pubchem_cid field in molecules table."""
    stats = {
        "total_molecules": 0,
        "with_pubchem_cid": 0,
        "without_pubchem_cid": 0,
        "with_pubchem_link": 0,
        "with_pubchem_cid_but_no_link": 0,
        "with_pubchem_link_but_no_cid": 0,
        "with_inconsistent_values": 0
    }
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Get total count
        cursor.execute("SELECT COUNT(*) FROM molecules")
        stats["total_molecules"] = cursor.fetchone()["count"]
        
        # Count molecules with pubchem_cid
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NOT NULL")
        stats["with_pubchem_cid"] = cursor.fetchone()["count"]
        
        # Count molecules without pubchem_cid
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NULL")
        stats["without_pubchem_cid"] = cursor.fetchone()["count"]
        
        # Count molecules with pubchem_link
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_link IS NOT NULL")
        stats["with_pubchem_link"] = cursor.fetchone()["count"]
        
        # Count molecules with pubchem_cid but no pubchem_link
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NOT NULL AND pubchem_link IS NULL")
        stats["with_pubchem_cid_but_no_link"] = cursor.fetchone()["count"]
        
        # Count molecules with pubchem_link but no pubchem_cid
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NULL AND pubchem_link IS NOT NULL")
        stats["with_pubchem_link_but_no_cid"] = cursor.fetchone()["count"]
        
        # Count molecules with inconsistent values (pubchem_link doesn't contain pubchem_cid)
        cursor.execute("""
            SELECT COUNT(*) FROM molecules 
            WHERE pubchem_cid IS NOT NULL AND pubchem_link IS NOT NULL 
            AND pubchem_link NOT LIKE '%' || pubchem_cid || '%'
        """)
        stats["with_inconsistent_values"] = cursor.fetchone()["count"]
        
        # Get sample records with missing pubchem_cid values (if any)
        if stats["without_pubchem_cid"] > 0:
            cursor.execute("""
                SELECT id, name, pubchem_cid, pubchem_link, molecular_formula
                FROM molecules
                WHERE pubchem_cid IS NULL
                LIMIT 5
            """)
            stats["missing_cid_examples"] = cursor.fetchall()
        
        # Get sample records with inconsistent values (if any)
        if stats["with_inconsistent_values"] > 0:
            cursor.execute("""
                SELECT id, name, pubchem_cid, pubchem_link
                FROM molecules
                WHERE pubchem_cid IS NOT NULL AND pubchem_link IS NOT NULL 
                AND pubchem_link NOT LIKE '%' || pubchem_cid || '%'
                LIMIT 5
            """)
            stats["inconsistent_examples"] = cursor.fetchall()
    
    return stats

def fix_pubchem_links(conn, stats):
    """Fix inconsistent pubchem_link values based on pubchem_cid."""
    updates = {
        "updated_records": 0,
        "details": []
    }
    
    with conn.cursor() as cursor:
        # Update inconsistent pubchem_link values
        if stats["with_inconsistent_values"] > 0:
            cursor.execute("""
                UPDATE molecules
                SET pubchem_link = 'https://pubchem.ncbi.nlm.nih.gov/compound/' || pubchem_cid,
                    updated_at = NOW()
                WHERE pubchem_cid IS NOT NULL AND pubchem_link IS NOT NULL 
                AND pubchem_link NOT LIKE '%' || pubchem_cid || '%'
            """)
            updates["updated_records"] += cursor.rowcount
            print(f"Updated {cursor.rowcount} records with inconsistent pubchem_link values")
        
        # Add pubchem_link where missing but pubchem_cid exists
        if stats["with_pubchem_cid_but_no_link"] > 0:
            cursor.execute("""
                UPDATE molecules
                SET pubchem_link = 'https://pubchem.ncbi.nlm.nih.gov/compound/' || pubchem_cid,
                    updated_at = NOW()
                WHERE pubchem_cid IS NOT NULL AND pubchem_link IS NULL
            """)
            updates["updated_link_added"] = cursor.rowcount
            print(f"Added pubchem_link for {cursor.rowcount} records with missing links")
        
        # Add modification records
        cursor.execute("""
            UPDATE molecules
            SET modification_history = COALESCE(modification_history, '[]'::jsonb) || jsonb_build_object(
                'timestamp', NOW(),
                'action', 'fix_pubchem_links',
                'details', 'Fixed inconsistent pubchem_link values'
            )::jsonb
            WHERE pubchem_cid IS NOT NULL
        """)
    
    return updates

def main():
    """Check and fix pubchem_cid and pubchem_link fields in molecules table."""
    print("Analyzing pubchem_cid and pubchem_link fields...")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Analyze current state
        stats = analyze_pubchem_cid_field(conn)
        print("\nCurrent State Analysis:")
        print(f"Total molecules: {stats['total_molecules']}")
        print(f"With pubchem_cid: {stats['with_pubchem_cid']}")
        print(f"Without pubchem_cid: {stats['without_pubchem_cid']}")
        print(f"With pubchem_link: {stats['with_pubchem_link']}")
        print(f"With pubchem_cid but no link: {stats['with_pubchem_cid_but_no_link']}")
        print(f"With pubchem_link but no cid: {stats['with_pubchem_link_but_no_cid']}")
        print(f"With inconsistent values: {stats['with_inconsistent_values']}")
        
        # Show examples if there are records with missing pubchem_cid
        if stats.get("missing_cid_examples"):
            print("\nExamples of records with missing pubchem_cid:")
            for example in stats["missing_cid_examples"]:
                print(f"  ID: {example['id']}, Name: {example['name']}, Formula: {example['molecular_formula']}")
        
        # Show examples of records with inconsistent values
        if stats.get("inconsistent_examples"):
            print("\nExamples of records with inconsistent values:")
            for example in stats["inconsistent_examples"]:
                print(f"  Name: {example['name']}, pubchem_cid: {example['pubchem_cid']}, pubchem_link: {example['pubchem_link']}")
        
        # Ask if user wants to fix inconsistencies
        if stats["with_inconsistent_values"] > 0 or stats["with_pubchem_cid_but_no_link"] > 0:
            user_input = input("\nWould you like to fix inconsistent pubchem_link values? (y/n): ")
            if user_input.lower() == 'y':
                # Fix inconsistencies
                print("\nFixing inconsistencies...")
                updates = fix_pubchem_links(conn, stats)
                
                # Commit changes
                conn.commit()
                
                # Reanalyze after changes
                print("\nReanalyzing after changes...")
                new_stats = analyze_pubchem_cid_field(conn)
                print("\nUpdated State Analysis:")
                print(f"Total molecules: {new_stats['total_molecules']}")
                print(f"With pubchem_cid: {new_stats['with_pubchem_cid']}")
                print(f"With pubchem_link: {new_stats['with_pubchem_link']} (Before: {stats['with_pubchem_link']})")
                print(f"With pubchem_cid but no link: {new_stats['with_pubchem_cid_but_no_link']} (Before: {stats['with_pubchem_cid_but_no_link']})")
                print(f"With inconsistent values: {new_stats['with_inconsistent_values']} (Before: {stats['with_inconsistent_values']})")
                
                # Save results to file
                results = {
                    "before": stats,
                    "updates": updates,
                    "after": new_stats
                }
                with open("pubchem_cid_fix_report.json", "w") as f:
                    json.dump(results, f, indent=2, default=str)
                
                print("\nPubChem CID/link fixes completed successfully")
                print("Report saved to 'pubchem_cid_fix_report.json'")
            else:
                print("No changes were made.")
        else:
            print("\nNo inconsistencies found that need fixing.")
    
    except Exception as e:
        conn.rollback()
        print(f"Error analyzing PubChem CID fields: {e}")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()