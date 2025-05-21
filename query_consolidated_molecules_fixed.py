#!/usr/bin/env python3
"""
Query consolidated molecules after consolidation.

This script demonstrates how to:
1. Query molecules and handle both primary and secondary molecules
2. Look up primary molecules for any secondary molecules
3. Show relationships between primary and secondary molecules
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Connect to the database
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

def get_primary_molecule_id(conn, molecule_id):
    """
    Get the primary molecule ID for any molecule.
    Returns the same ID if it's already a primary molecule.
    """
    with conn.cursor() as cursor:
        cursor.execute("""
            SELECT 
                CASE 
                    WHEN properties->>'consolidated_to' IS NOT NULL THEN 
                        (properties->>'consolidated_to')::uuid
                    ELSE 
                        id
                END AS primary_id
            FROM molecules
            WHERE id = %s
        """, (molecule_id,))
        
        result = cursor.fetchone()
        return result[0] if result else None

def is_primary_molecule(conn, molecule_id):
    """
    Check if a molecule is a primary molecule (not consolidated).
    """
    with conn.cursor() as cursor:
        cursor.execute("""
            SELECT 
                properties->>'consolidated_to' IS NULL AS is_primary
            FROM molecules
            WHERE id = %s
        """, (molecule_id,))
        
        result = cursor.fetchone()
        return result[0] if result else False

def get_molecule_with_primary_info(conn, molecule_id):
    """
    Get a molecule with information about its primary molecule if it's consolidated.
    """
    primary_id = get_primary_molecule_id(conn, molecule_id)
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Get the original molecule
        cursor.execute("""
            SELECT * FROM molecules WHERE id = %s
        """, (molecule_id,))
        
        molecule = cursor.fetchone()
        if not molecule:
            return None
        
        # Add consolidated info
        molecule['is_consolidated'] = primary_id != molecule_id
        
        # If consolidated, get primary molecule info
        if molecule['is_consolidated']:
            cursor.execute("""
                SELECT * FROM molecules WHERE id = %s
            """, (primary_id,))
            
            primary_molecule = cursor.fetchone()
            if primary_molecule:
                molecule['primary_molecule'] = primary_molecule
                molecule['primary_id'] = primary_id
        
        return molecule

def find_molecules_by_name(conn, name_pattern, include_secondary=False):
    """
    Find molecules by name pattern, optionally including or excluding secondary molecules.
    """
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        if include_secondary:
            # Include all molecules (both primary and secondary)
            cursor.execute("""
                SELECT * FROM molecules
                WHERE name ILIKE %s
            """, (f'%{name_pattern}%',))
        else:
            # Only primary molecules (those without consolidated_to property)
            cursor.execute("""
                SELECT * FROM molecules
                WHERE name ILIKE %s
                AND (properties->>'consolidated_to' IS NULL)
            """, (f'%{name_pattern}%',))
        
        molecules = cursor.fetchall()
        
        # Enhance with consolidated information
        for molecule in molecules:
            primary_id = get_primary_molecule_id(conn, molecule['id'])
            molecule['is_consolidated'] = primary_id != molecule['id']
            molecule['primary_id'] = primary_id
            
            if molecule['is_consolidated']:
                cursor.execute("""
                    SELECT name FROM molecules WHERE id = %s
                """, (primary_id,))
                primary_name = cursor.fetchone()
                if primary_name:
                    molecule['primary_name'] = primary_name['name']
        
        return molecules

def get_secondary_molecules(conn, primary_id):
    """
    Get all secondary molecules that point to a specific primary molecule.
    """
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT *
            FROM molecules
            WHERE properties->>'consolidated_to' = %s
        """, (str(primary_id),))
        
        return cursor.fetchall()

def main():
    """
    Demonstrate how to use the consolidated molecules.
    """
    conn = connect_to_db()
    
    try:
        print("\nDemonstrating consolidated molecule queries:")
        
        # 1. Find molecules by name (including secondary)
        print("\n1. Finding all molecules containing 'glycol' (including secondary):")
        molecules = find_molecules_by_name(conn, 'glycol', include_secondary=True)
        for mol in molecules:
            if mol.get('is_consolidated'):
                print(f"  {mol['name']} (ID: {mol['id']}) -> Consolidated to {mol.get('primary_name', 'Unknown')} (ID: {mol['primary_id']})")
            else:
                print(f"  {mol['name']} (ID: {mol['id']}) - Primary molecule")
        
        # 2. Find primary molecules only
        print("\n2. Finding only primary molecules containing 'glycol':")
        primary_molecules = find_molecules_by_name(conn, 'glycol', include_secondary=False)
        for mol in primary_molecules:
            print(f"  {mol['name']} (ID: {mol['id']})")
            
            # Get any secondary molecules
            secondary_mols = get_secondary_molecules(conn, mol['id'])
            if secondary_mols:
                print(f"    Secondary molecules pointing to this primary:")
                for sec_mol in secondary_mols:
                    print(f"      {sec_mol['name']} (ID: {sec_mol['id']})")
        
        # 3. Get details for a specific molecule (handle both primary and secondary)
        if molecules:
            test_mol = molecules[0]
            print(f"\n3. Looking up details for molecule: {test_mol['name']} (ID: {test_mol['id']}):")
            molecule_details = get_molecule_with_primary_info(conn, test_mol['id'])
            
            if molecule_details['is_consolidated']:
                print(f"  This is a SECONDARY molecule consolidated to:")
                print(f"  {molecule_details['primary_molecule']['name']} (ID: {molecule_details['primary_id']})")
            else:
                print(f"  This is a PRIMARY molecule")
                
                # Check if it has any secondary molecules
                secondary_mols = get_secondary_molecules(conn, test_mol['id'])
                if secondary_mols:
                    print(f"  Secondary molecules:")
                    for sec_mol in secondary_mols:
                        print(f"    {sec_mol['name']} (ID: {sec_mol['id']})")
        
    except Exception as e:
        print(f"Error: {e}")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()