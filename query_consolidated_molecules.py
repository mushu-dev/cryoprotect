#!/usr/bin/env python3
"""
Query consolidated molecules using the consolidated_molecules view.

This script provides utility functions to:
1. Query molecules through the consolidated_molecules view
2. Get primary molecule information for any molecule ID
3. Demonstrate how to use the consolidated molecule system in applications
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

def setup_consolidated_view(conn):
    """
    Set up the consolidated_molecules view if it doesn't exist.
    """
    with conn.cursor() as cursor:
        # Check if the view already exists
        cursor.execute("""
            SELECT EXISTS (
                SELECT 1 FROM information_schema.views
                WHERE table_schema = 'public'
                AND table_name = 'consolidated_molecules'
            )
        """)
        view_exists = cursor.fetchone()[0]
        
        if not view_exists:
            print("Creating consolidated_molecules view...")
            
            # Create helper function to get primary molecule ID
            cursor.execute("""
                CREATE OR REPLACE FUNCTION get_primary_molecule_id(molecule_id UUID)
                RETURNS UUID AS $$
                DECLARE
                    primary_id UUID;
                BEGIN
                    -- Check if the molecule has a consolidated_to property
                    SELECT (properties->>'consolidated_to')::UUID INTO primary_id
                    FROM molecules
                    WHERE id = molecule_id;
                    
                    -- If not, it's already a primary molecule
                    IF primary_id IS NULL THEN
                        RETURN molecule_id;
                    ELSE
                        RETURN primary_id;
                    END IF;
                END;
                $$ LANGUAGE plpgsql;
            """)
            
            # Create helper function to check if a molecule is primary
            cursor.execute("""
                CREATE OR REPLACE FUNCTION is_primary_molecule(molecule_id UUID)
                RETURNS BOOLEAN AS $$
                BEGIN
                    -- Check if the molecule has no consolidated_to property
                    RETURN NOT EXISTS (
                        SELECT 1
                        FROM molecules
                        WHERE id = molecule_id
                        AND properties ? 'consolidated_to'
                    );
                END;
                $$ LANGUAGE plpgsql;
            """)
            
            # Check if consolidated_molecules view needs to be refreshed
            cursor.execute("""
                DROP VIEW IF EXISTS consolidated_molecules;
                DROP VIEW IF EXISTS primary_molecules;
            """)

            # Create consolidated_molecules view
            cursor.execute("""
                CREATE VIEW consolidated_molecules AS
                WITH consolidated_mapping AS (
                    SELECT
                        m.id,
                        CASE
                            WHEN m.properties->>'consolidated_to' IS NOT NULL THEN
                                m.properties->>'consolidated_to'
                            ELSE
                                m.id::text
                        END AS primary_id,
                        m.properties->>'consolidated_to' IS NOT NULL AS is_consolidated
                    FROM
                        molecules m
                )
                SELECT
                    m.id,
                    m.name,
                    m.pubchem_cid,
                    m.smiles,
                    m.molecular_formula,
                    m.data_source,
                    m.created_at,
                    m.updated_at,
                    m.properties,
                    cm.is_consolidated,
                    CASE WHEN cm.is_consolidated THEN pm.name ELSE m.name END AS consolidated_name,
                    CASE WHEN cm.is_consolidated THEN pm.id ELSE m.id END AS primary_molecule_id
                FROM
                    molecules m
                JOIN
                    consolidated_mapping cm ON m.id = cm.id
                LEFT JOIN
                    molecules pm ON pm.id = cm.primary_id::uuid;
            """)
            
            # Create primary_molecules view for querying just primary molecules
            cursor.execute("""
                CREATE OR REPLACE VIEW primary_molecules AS
                SELECT 
                    m.*
                FROM 
                    molecules m
                WHERE 
                    NOT (m.properties ? 'consolidated_to');
            """)
            
            conn.commit()
            print("Views created successfully.")
        else:
            print("Consolidated molecules view already exists.")

def get_primary_molecule(conn, molecule_id):
    """
    Get the primary molecule for any molecule ID.
    Returns the primary molecule itself if it's already primary.
    """
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT get_primary_molecule_id(%s) as primary_id
        """, (molecule_id,))
        
        result = cursor.fetchone()
        primary_id = result['primary_id']
        
        # Get the primary molecule details
        cursor.execute("""
            SELECT * FROM molecules WHERE id = %s
        """, (primary_id,))
        
        return cursor.fetchone()

def get_consolidated_view_molecule(conn, molecule_id):
    """
    Get a molecule from the consolidated_molecules view.
    This shows how redirected properties work.
    """
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT * FROM consolidated_molecules WHERE id = %s
        """, (molecule_id,))
        
        return cursor.fetchone()

def get_secondary_molecules(conn, primary_id):
    """
    Get all secondary molecules that redirect to a primary molecule.
    """
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT * FROM molecules
            WHERE properties->>'consolidated_to' = %s
        """, (primary_id,))
        
        return cursor.fetchall()

def find_molecules_by_name(conn, name_pattern, include_secondary=False):
    """
    Find molecules by name pattern, optionally including secondary molecules.
    """
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        if include_secondary:
            cursor.execute("""
                SELECT * FROM consolidated_molecules
                WHERE name ILIKE %s OR consolidated_name ILIKE %s
            """, (f'%{name_pattern}%', f'%{name_pattern}%'))
        else:
            cursor.execute("""
                SELECT * FROM consolidated_molecules
                WHERE (name ILIKE %s OR consolidated_name ILIKE %s)
                AND NOT is_consolidated
            """, (f'%{name_pattern}%', f'%{name_pattern}%'))
        
        return cursor.fetchall()

def main():
    """
    Demonstrate how to use the consolidated molecules system.
    """
    conn = connect_to_db()
    
    try:
        # Setup consolidated view if it doesn't exist
        setup_consolidated_view(conn)
        
        print("\nDemonstrating consolidated molecule system:")
        
        # 1. Find molecules by name
        print("\n1. Finding molecules by name 'glycol':")
        molecules = find_molecules_by_name(conn, 'glycol', include_secondary=True)
        for mol in molecules:
            if mol['is_consolidated']:
                print(f"  {mol['name']} (ID: {mol['id']}) -> Consolidated to {mol['consolidated_name']} (ID: {mol['primary_molecule_id']})")
            else:
                print(f"  {mol['name']} (ID: {mol['id']}) - Primary molecule")
        
        # 2. Get primary molecule for any ID
        if molecules and any(mol['is_consolidated'] for mol in molecules):
            secondary_mol = next(mol for mol in molecules if mol['is_consolidated'])
            secondary_id = secondary_mol['id']
            
            print(f"\n2. Getting primary molecule for secondary molecule {secondary_mol['name']} (ID: {secondary_id}):")
            primary_mol = get_primary_molecule(conn, secondary_id)
            print(f"  Primary molecule: {primary_mol['name']} (ID: {primary_mol['id']})")
            
            # 3. Get view version
            print(f"\n3. Getting consolidated view entry for {secondary_mol['name']} (ID: {secondary_id}):")
            view_mol = get_consolidated_view_molecule(conn, secondary_id)
            print(f"  Original name: {view_mol['name']}")
            print(f"  Consolidated name: {view_mol['consolidated_name']}")
            print(f"  Is consolidated: {view_mol['is_consolidated']}")
            print(f"  Primary ID: {view_mol['primary_molecule_id']}")
            
            # 4. Get all secondary molecules for a primary
            print(f"\n4. Finding all secondary molecules for primary {primary_mol['name']} (ID: {primary_mol['id']}):")
            secondary_mols = get_secondary_molecules(conn, primary_mol['id'])
            for sec_mol in secondary_mols:
                print(f"  {sec_mol['name']} (ID: {sec_mol['id']})")
        
        # 5. Find only primary molecules
        print("\n5. Finding primary molecules with 'glycerol':")
        primary_only = find_molecules_by_name(conn, 'glycerol', include_secondary=False)
        for mol in primary_only:
            print(f"  {mol['name']} (ID: {mol['id']})")
            
    except Exception as e:
        print(f"Error: {e}")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()