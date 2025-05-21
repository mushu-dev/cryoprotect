#!/usr/bin/env python3
"""
Example of how to use the consolidated molecules view in applications.

This script demonstrates:
1. Querying consolidated molecules view
2. Using primary_molecules view for duplicate-free queries
3. Using the get_primary_molecule_id function
4. Implementing consolidation handling in application code
"""

import os
import sys
import json
import uuid
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
from tabulate import tabulate

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

def display_consolidated_molecules(conn):
    """
    Display consolidated molecules using the consolidated_molecules view.
    """
    print("\n=== Consolidated Molecules ===")
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Query consolidated molecules
        cursor.execute("""
            SELECT id, name, molecule_status, primary_molecule_id, primary_molecule_name,
                   pubchem_cid, primary_pubchem_cid
            FROM consolidated_molecules
            WHERE is_consolidated = TRUE
            ORDER BY name
        """)
        
        consolidated = cursor.fetchall()
        
        # Format the results for display
        table_data = []
        for molecule in consolidated:
            table_data.append([
                molecule['name'],
                molecule['id'],
                molecule['pubchem_cid'],
                molecule['primary_molecule_name'],
                molecule['primary_molecule_id'],
                molecule['primary_pubchem_cid']
            ])
        
        # Display results
        headers = ["Name", "ID", "PubChem CID", "Primary Name", "Primary ID", "Primary PubChem CID"]
        print(tabulate(table_data, headers=headers, tablefmt="grid"))
        
        print(f"\nTotal consolidated molecules: {len(consolidated)}")

def display_primary_molecules(conn):
    """
    Display primary molecules used in consolidation.
    """
    print("\n=== Primary Molecules ===")
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Query primary molecules
        cursor.execute("""
            SELECT id, name, molecule_status, pubchem_cid
            FROM consolidated_molecules
            WHERE molecule_status = 'Primary'
            ORDER BY name
        """)
        
        primaries = cursor.fetchall()
        
        # Format the results for display
        table_data = []
        for molecule in primaries:
            table_data.append([
                molecule['name'],
                molecule['id'],
                molecule['pubchem_cid'],
                molecule['molecule_status']
            ])
        
        # Display results
        headers = ["Name", "ID", "PubChem CID", "Status"]
        print(tabulate(table_data, headers=headers, tablefmt="grid"))
        
        print(f"\nTotal primary molecules: {len(primaries)}")

def demonstrate_function_usage(conn):
    """
    Demonstrate usage of the consolidated molecule functions.
    """
    print("\n=== Consolidated Molecule Function Examples ===")
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Get a consolidated molecule to use as an example
        cursor.execute("""
            SELECT id, name FROM consolidated_molecules
            WHERE is_consolidated = TRUE
            LIMIT 1
        """)
        
        consolidated = cursor.fetchone()
        
        if consolidated:
            # Example 1: Get primary molecule ID using function
            consolidated_id = consolidated['id']
            cursor.execute("""
                SELECT get_primary_molecule_id(%s) as primary_id
            """, (consolidated_id,))
            
            result = cursor.fetchone()
            primary_id = result['primary_id']
            
            print(f"Example 1: Get primary ID for {consolidated['name']} ({consolidated_id})")
            print(f"  get_primary_molecule_id('{consolidated_id}') = {primary_id}")
            
            # Example 2: Check if a molecule is primary
            cursor.execute("""
                SELECT is_primary_molecule(%s) as is_primary
            """, (consolidated_id,))
            
            result = cursor.fetchone()
            is_primary_consolidated = result['is_primary']
            
            cursor.execute("""
                SELECT is_primary_molecule(%s) as is_primary
            """, (primary_id,))
            
            result = cursor.fetchone()
            is_primary_primary = result['is_primary']
            
            print(f"\nExample 2: Check if molecules are primary")
            print(f"  is_primary_molecule('{consolidated_id}') = {is_primary_consolidated}")
            print(f"  is_primary_molecule('{primary_id}') = {is_primary_primary}")
            
            # Example 3: Get all consolidated molecules for a primary
            cursor.execute("""
                SELECT id, name
                FROM consolidated_molecules
                WHERE primary_molecule_id = %s AND is_consolidated = TRUE
            """, (primary_id,))
            
            consolidated_list = cursor.fetchall()
            
            print(f"\nExample 3: Get all consolidated molecules for primary ID {primary_id}")
            for mol in consolidated_list:
                print(f"  {mol['name']} ({mol['id']})")
        else:
            print("No consolidated molecules found to demonstrate functions")

def demonstrate_application_integration(conn):
    """
    Demonstrate how to handle consolidated molecules in application code.
    """
    print("\n=== Application Integration Example ===")
    print("Scenario: User requests molecular properties for a specific molecule")
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Get a consolidated molecule to use as an example
        cursor.execute("""
            SELECT id, name, primary_molecule_id, primary_molecule_name,
                  is_consolidated
            FROM consolidated_molecules
            WHERE is_consolidated = TRUE
            LIMIT 1
        """)
        
        molecule = cursor.fetchone()
        
        if not molecule:
            print("No consolidated molecules found for the example")
            return
        
        # Simulate user request for this molecule ID
        requested_id = molecule['id']
        requested_name = molecule['name']
        
        print(f"User requests properties for: {requested_name} ({requested_id})")
        
        # Application code pattern:
        # 1. Use the get_primary_molecule_id function to get the primary ID
        cursor.execute("""
            SELECT get_primary_molecule_id(%s) as primary_id
        """, (requested_id,))
        
        result = cursor.fetchone()
        primary_id = result['primary_id']
        
        # 2. Check if this molecule is consolidated
        cursor.execute("""
            SELECT is_consolidated, primary_molecule_name
            FROM consolidated_molecules
            WHERE id = %s
        """, (requested_id,))
        
        result = cursor.fetchone()
        is_consolidated = result['is_consolidated']
        primary_name = result['primary_molecule_name'] if is_consolidated else None
        
        # 3. Get molecular properties using the primary ID
        cursor.execute("""
            SELECT mp.property_type, mp.property_value, mp.numeric_value, 
                   mp.unit, mp.source
            FROM molecular_properties mp
            WHERE mp.molecule_id = %s
            LIMIT 5
        """, (primary_id,))
        
        properties = cursor.fetchall()
        
        # 4. Display the results with consolidated info
        print("\nApplication Response:")
        
        if is_consolidated:
            print(f"NOTE: This molecule has been consolidated into {primary_name}")
            print(f"Showing properties from the primary molecule ({primary_id})")
        
        print("\nMolecular Properties:")
        table_data = []
        for prop in properties:
            value = prop['property_value'] or prop['numeric_value']
            table_data.append([
                prop['property_type'],
                value,
                prop['unit'],
                prop['source']
            ])
        
        headers = ["Property Type", "Value", "Unit", "Source"]
        print(tabulate(table_data, headers=headers, tablefmt="simple"))

def main():
    """Example script to demonstrate the consolidated molecules view and functions."""
    # Connect to database
    conn = connect_to_db()
    
    try:
        # Display consolidated molecules
        display_consolidated_molecules(conn)
        
        # Display primary molecules
        display_primary_molecules(conn)
        
        # Demonstrate function usage
        demonstrate_function_usage(conn)
        
        # Demonstrate application integration
        demonstrate_application_integration(conn)
        
    except Exception as e:
        print(f"Error in example script: {e}")
    finally:
        conn.close()

if __name__ == "__main__":
    main()