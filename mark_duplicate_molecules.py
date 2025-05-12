#!/usr/bin/env python3
"""
Mark duplicate molecules in the database.

This script:
1. Identifies likely duplicate molecules based on name and formula
2. Sets a 'duplicate_group' property in their properties field
3. Generates a report of the duplicate groups
"""

import os
import sys
import json
import uuid
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

def find_duplicate_names(conn):
    """Find molecules with duplicate names."""
    duplicates = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Find names that appear multiple times
        cursor.execute("""
            SELECT lower(name) as name, COUNT(*) as count
            FROM molecules
            WHERE name IS NOT NULL AND name NOT LIKE 'TEST\\_%'
            GROUP BY lower(name)
            HAVING COUNT(*) > 1
            ORDER BY count DESC
        """)
        
        duplicate_names = cursor.fetchall()
        
        # For each duplicate name, get the molecules
        for dup in duplicate_names:
            name = dup["name"]
            cursor.execute("""
                SELECT id, name, molecular_formula, pubchem_cid, smiles, 
                       data_source, created_at, is_public, properties
                FROM molecules
                WHERE lower(name) = %s
                ORDER BY pubchem_cid NULLS LAST, created_at
            """, (name,))
            
            molecules = cursor.fetchall()
            
            # Generate a unique group ID
            group_id = str(uuid.uuid4())
            
            duplicates.append({
                "group_id": group_id,
                "name": name,
                "count": dup["count"],
                "molecules": molecules,
                "duplicate_type": "name"
            })
    
    return duplicates

def find_duplicate_formulas(conn):
    """Find molecules with duplicate molecular formulas."""
    duplicates = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Find formulas that appear multiple times
        cursor.execute("""
            SELECT molecular_formula, COUNT(*) as count
            FROM molecules
            WHERE molecular_formula IS NOT NULL 
            AND name NOT LIKE 'TEST\\_%'
            GROUP BY molecular_formula
            HAVING COUNT(*) > 1
            ORDER BY count DESC
        """)
        
        duplicate_formulas = cursor.fetchall()
        
        # For each duplicate formula, get the molecules
        for dup in duplicate_formulas:
            formula = dup["molecular_formula"]
            cursor.execute("""
                SELECT id, name, molecular_formula, pubchem_cid, smiles, 
                       data_source, created_at, is_public, properties
                FROM molecules
                WHERE molecular_formula = %s
                ORDER BY pubchem_cid NULLS LAST, created_at
            """, (formula,))
            
            molecules = cursor.fetchall()
            
            # Generate a unique group ID
            group_id = str(uuid.uuid4())
            
            duplicates.append({
                "group_id": group_id,
                "formula": formula,
                "count": dup["count"],
                "molecules": molecules,
                "duplicate_type": "formula"
            })
    
    return duplicates

def mark_molecule_as_duplicate(conn, molecule_id, group_id, duplicate_type, group_value, dry_run=False):
    """Mark a molecule as part of a duplicate group."""
    try:
        if not dry_run:
            with conn.cursor() as cursor:
                cursor.execute("""
                    UPDATE molecules
                    SET properties = COALESCE(properties, '{}'::jsonb) || %s::jsonb,
                        updated_at = NOW(),
                        modification_history = COALESCE(modification_history, '[]'::jsonb) || jsonb_build_object(
                            'timestamp', NOW(),
                            'action', 'mark_duplicate_molecules',
                            'details', %s
                        )::jsonb
                    WHERE id = %s
                """, (
                    json.dumps({
                        "duplicate_group": group_id,
                        "duplicate_type": duplicate_type,
                        "duplicate_value": group_value
                    }),
                    f"Marked as duplicate in group {group_id} by {duplicate_type}",
                    molecule_id
                ))
                
                return cursor.rowcount
        else:
            # Return what would be updated in dry run mode
            return 1
    except Exception as e:
        print(f"  Error marking molecule as duplicate: {e}")
        return 0

def mark_molecules_in_group(conn, group, dry_run=False):
    """Mark all molecules in a duplicate group."""
    updated_count = 0
    group_id = group["group_id"]
    duplicate_type = group["duplicate_type"]
    group_value = group.get("name") if duplicate_type == "name" else group.get("formula")
    
    for molecule in group["molecules"]:
        molecule_id = molecule["id"]
        name = molecule["name"]
        
        if not dry_run:
            print(f"  Marking {name} (ID: {molecule_id}) as part of group {group_id}")
        else:
            print(f"  Would mark {name} (ID: {molecule_id}) as part of group {group_id}")
        
        updated = mark_molecule_as_duplicate(conn, molecule_id, group_id, duplicate_type, group_value, dry_run)
        if updated:
            updated_count += 1
    
    return updated_count

def main():
    """Mark duplicate molecules in the database."""
    import argparse
    parser = argparse.ArgumentParser(description="Mark duplicate molecules in the database")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be updated without making changes")
    parser.add_argument("--auto-commit", action="store_true", help="Automatically commit changes without confirmation")
    args = parser.parse_args()
    
    print("Marking duplicate molecules in the database...")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Find duplicate names
        name_duplicates = find_duplicate_names(conn)
        print(f"Found {len(name_duplicates)} groups of molecules with duplicate names")
        
        # Find duplicate formulas
        formula_duplicates = find_duplicate_formulas(conn)
        print(f"Found {len(formula_duplicates)} groups of molecules with duplicate formulas")
        
        if not name_duplicates and not formula_duplicates:
            print("No duplicate molecules to mark.")
            return
        
        # Mark duplicates
        name_updated_count = 0
        formula_updated_count = 0
        
        print("\nProcessing name duplicates...")
        for group in name_duplicates:
            name = group["name"]
            count = group["count"]
            print(f"\nDuplicate name group: {name} ({count} molecules)")
            
            updated = mark_molecules_in_group(conn, group, args.dry_run)
            name_updated_count += updated
        
        print("\nProcessing formula duplicates...")
        for group in formula_duplicates:
            formula = group["formula"]
            count = group["count"]
            print(f"\nDuplicate formula group: {formula} ({count} molecules)")
            
            updated = mark_molecules_in_group(conn, group, args.dry_run)
            formula_updated_count += updated
        
        # Save results to file
        results = {
            "name_duplicates": {
                "groups": len(name_duplicates),
                "total_molecules": sum(g["count"] for g in name_duplicates),
                "updated_molecules": name_updated_count,
                "groups_data": name_duplicates
            },
            "formula_duplicates": {
                "groups": len(formula_duplicates),
                "total_molecules": sum(g["count"] for g in formula_duplicates),
                "updated_molecules": formula_updated_count,
                "groups_data": formula_duplicates
            },
            "dry_run": args.dry_run
        }
        
        with open("duplicate_molecules_report.json", "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        print("\nSummary:")
        print(f"Name duplicate groups: {len(name_duplicates)}")
        print(f"Formula duplicate groups: {len(formula_duplicates)}")
        print(f"Total molecules marked as name duplicates: {name_updated_count}")
        print(f"Total molecules marked as formula duplicates: {formula_updated_count}")
        
        print(f"Report saved to 'duplicate_molecules_report.json'")
        
        # Commit changes if not dry run
        if not args.dry_run:
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
        print(f"Error marking duplicate molecules: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()