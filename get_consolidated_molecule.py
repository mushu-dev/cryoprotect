#!/usr/bin/env python3
"""
Helper utility to get the primary molecule for any molecule ID.

This script:
1. Takes a molecule ID as input
2. Checks if it's a consolidated molecule
3. Returns the primary molecule ID if it's consolidated
4. Otherwise returns the original ID
"""

import os
import sys
import json
import uuid
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

def get_consolidated_molecule(conn, molecule_id):
    """
    Get the primary molecule ID for a given molecule ID.
    
    Args:
        conn: Database connection
        molecule_id: UUID of the molecule to check
        
    Returns:
        Dictionary with molecule details including primary ID
    """
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # First, check if the molecule exists
            cursor.execute("""
                SELECT id, name, pubchem_cid, molecular_formula, properties,
                       created_at, updated_at
                FROM molecules
                WHERE id = %s
            """, (molecule_id,))
            
            molecule = cursor.fetchone()
            if not molecule:
                return {
                    "error": f"Molecule with ID {molecule_id} not found",
                    "id": molecule_id,
                    "primary_id": None,
                    "is_consolidated": False
                }
            
            # Check if it's a consolidated molecule
            is_consolidated = (molecule['properties'] and 
                               'consolidated_to' in molecule['properties'])
            
            if is_consolidated:
                primary_id = molecule['properties']['consolidated_to']
                
                # Get primary molecule details
                cursor.execute("""
                    SELECT id, name, pubchem_cid, molecular_formula, properties,
                           created_at, updated_at
                    FROM molecules
                    WHERE id = %s
                """, (primary_id,))
                
                primary_molecule = cursor.fetchone()
                if primary_molecule:
                    return {
                        "id": molecule_id,
                        "name": molecule['name'],
                        "primary_id": primary_id,
                        "primary_name": primary_molecule['name'],
                        "is_consolidated": True,
                        "pubchem_cid": molecule['pubchem_cid'],
                        "primary_pubchem_cid": primary_molecule['pubchem_cid'],
                        "molecular_formula": molecule['molecular_formula']
                    }
                else:
                    # This is unusual - the primary doesn't exist
                    return {
                        "error": f"Primary molecule with ID {primary_id} not found",
                        "id": molecule_id,
                        "primary_id": primary_id,
                        "is_consolidated": True,
                        "name": molecule['name']
                    }
            else:
                # Not consolidated - this is the primary
                return {
                    "id": molecule_id,
                    "name": molecule['name'],
                    "primary_id": molecule_id,  # Same as id
                    "is_consolidated": False,
                    "pubchem_cid": molecule['pubchem_cid'],
                    "molecular_formula": molecule['molecular_formula']
                }
    except Exception as e:
        return {
            "error": f"Error retrieving molecule: {str(e)}",
            "id": molecule_id,
            "primary_id": None,
            "is_consolidated": False
        }

def main():
    """Command-line interface for retrieving primary molecule information."""
    parser = argparse.ArgumentParser(description="Get consolidated molecule information.")
    parser.add_argument("molecule_id", help="UUID of the molecule to check")
    parser.add_argument("--json", action="store_true", help="Output in JSON format")
    args = parser.parse_args()
    
    try:
        molecule_id = str(uuid.UUID(args.molecule_id))
    except ValueError:
        print(f"Error: Invalid UUID format: {args.molecule_id}")
        sys.exit(1)
    
    conn = connect_to_db()
    try:
        result = get_consolidated_molecule(conn, molecule_id)
        
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            print(f"Molecule ID: {result['id']}")
            print(f"Name: {result.get('name', 'Unknown')}")
            
            if 'error' in result:
                print(f"Error: {result['error']}")
            else:
                if result['is_consolidated']:
                    print(f"Status: Consolidated to {result['primary_id']}")
                    print(f"Primary Name: {result['primary_name']}")
                    print(f"PubChem CID: {result['pubchem_cid']} → {result['primary_pubchem_cid']}")
                else:
                    print("Status: Primary molecule (not consolidated)")
                    print(f"PubChem CID: {result['pubchem_cid']}")
                
                print(f"Molecular Formula: {result.get('molecular_formula', 'Unknown')}")
    finally:
        conn.close()

if __name__ == "__main__":
    main()