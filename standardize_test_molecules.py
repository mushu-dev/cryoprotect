#!/usr/bin/env python3
"""
Standardize test molecules in the database.

This script:
1. Identifies molecules that appear to be test entries
2. Updates their names with a TEST_ prefix if not already present
3. Adds a metadata flag to clearly indicate test data
4. Generates a report of changes made
"""

import os
import sys
import json
import argparse
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

def identify_test_molecules(conn):
    """Identify molecules that appear to be test entries."""
    test_molecules = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, molecular_formula, pubchem_cid, 
                   smiles, properties, is_public
            FROM molecules
            WHERE 
                name ILIKE '%test%' OR 
                name ILIKE '%dummy%' OR 
                name ILIKE '%example%' OR
                (properties IS NOT NULL AND properties ? 'is_test')
            ORDER BY name
        """)
        
        test_molecules = cursor.fetchall()
    
    return test_molecules

def update_test_molecule(conn, molecule_id, old_name, dry_run=False):
    """Update a test molecule with standardized naming and metadata."""
    # Determine new name with TEST_ prefix if not already present
    new_name = old_name
    if not old_name.startswith("TEST_"):
        new_name = f"TEST_{old_name}"
    
    # Update properties to include test flag
    try:
        if not dry_run:
            with conn.cursor() as cursor:
                # Update name and add test flag to properties
                cursor.execute("""
                    UPDATE molecules
                    SET name = %s,
                        properties = COALESCE(properties, '{}'::jsonb) || '{"is_test": true}'::jsonb,
                        is_public = false,
                        updated_at = NOW(),
                        modification_history = COALESCE(modification_history, '[]'::jsonb) || jsonb_build_object(
                            'timestamp', NOW(),
                            'action', 'standardize_test_molecules',
                            'details', 'Standardized test molecule name and metadata',
                            'old_name', %s,
                            'new_name', %s
                        )::jsonb
                    WHERE id = %s
                """, (new_name, old_name, new_name, molecule_id))
                
                return cursor.rowcount, new_name
        else:
            # Return what would be updated in dry run mode
            return 1, new_name
    except Exception as e:
        print(f"  Error updating test molecule {old_name}: {e}")
        return 0, old_name

def main():
    """Standardize test molecules in the database."""
    parser = argparse.ArgumentParser(description="Standardize test molecules in the database")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be updated without making changes")
    parser.add_argument("--auto-commit", action="store_true", help="Automatically commit changes without confirmation")
    args = parser.parse_args()
    
    print("Standardizing test molecules in the database...")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Identify test molecules
        test_molecules = identify_test_molecules(conn)
        print(f"Found {len(test_molecules)} potential test molecules")
        
        if not test_molecules:
            print("No test molecules to update.")
            return
        
        # Process each test molecule
        updated_count = 0
        needs_update_count = 0
        results = []
        
        for molecule in test_molecules:
            molecule_id = molecule["id"]
            name = molecule["name"]
            
            # Check if the molecule already has the TEST_ prefix
            if name.startswith("TEST_"):
                print(f"  {name} already has TEST_ prefix, checking metadata...")
                
                # Check if properties already has is_test flag
                has_test_flag = molecule.get("properties") and molecule["properties"].get("is_test")
                
                if has_test_flag:
                    print(f"  {name} already has complete test metadata, skipping")
                    results.append({
                        "molecule_id": molecule_id,
                        "name": name,
                        "status": "already_standardized"
                    })
                else:
                    print(f"  {name} needs metadata update")
                    needs_update_count += 1
                    updated, new_name = update_test_molecule(conn, molecule_id, name, args.dry_run)
                    
                    if updated:
                        updated_count += 1
                        results.append({
                            "molecule_id": molecule_id,
                            "name": name,
                            "new_name": new_name,
                            "status": "metadata_updated"
                        })
            else:
                print(f"  Standardizing: {name}")
                needs_update_count += 1
                updated, new_name = update_test_molecule(conn, molecule_id, name, args.dry_run)
                
                if updated:
                    updated_count += 1
                    print(f"  Updated name from '{name}' to '{new_name}'")
                    results.append({
                        "molecule_id": molecule_id,
                        "name": name,
                        "new_name": new_name,
                        "status": "updated"
                    })
                else:
                    print(f"  Failed to update {name}")
                    results.append({
                        "molecule_id": molecule_id,
                        "name": name,
                        "status": "update_failed"
                    })
        
        # Summary
        print("\nSummary:")
        print(f"Total test molecules: {len(test_molecules)}")
        print(f"Molecules needing updates: {needs_update_count}")
        print(f"Molecules updated: {updated_count}")
        
        # Save report
        with open("test_molecule_standardization_report.json", "w") as f:
            json.dump({
                "total": len(test_molecules),
                "needed_updates": needs_update_count,
                "updated": updated_count,
                "results": results,
                "dry_run": args.dry_run
            }, f, indent=2)
        
        print(f"Report saved to 'test_molecule_standardization_report.json'")
        
        # Commit changes if not dry run
        if not args.dry_run and updated_count > 0:
            if args.auto_commit:
                conn.commit()
                print("Changes automatically committed to the database")
            else:
                confirm = input("\nCommit these changes to the database? (y/n): ")
                if confirm.lower() == 'y':
                    conn.commit()
                    print("Changes committed to the database")
                else:
                    conn.rollback()
                    print("Changes rolled back, no updates were made")
    
    except Exception as e:
        conn.rollback()
        print(f"Error standardizing test molecules: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()