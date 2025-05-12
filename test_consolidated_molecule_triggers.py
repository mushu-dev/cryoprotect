#!/usr/bin/env python3
"""
Test script to verify that consolidated molecule triggers are working correctly.
This script verifies:
1. Automatic redirection to primary molecules
2. Protection against modifying secondary molecules
3. Prevention of circular references
4. Protection against deleting primary molecules with secondaries
"""

import os
import sys
import uuid
import psycopg2
import psycopg2.extras
from datetime import datetime
import json

# Import database connection utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from database.utils import get_db_connection

def run_test():
    """Run tests to verify the consolidated molecule triggers work correctly."""
    conn = get_db_connection()
    conn.autocommit = False
    
    try:
        print("Starting consolidated molecule trigger tests")
        print("=" * 80)
        
        # Create test molecules
        primary_id = create_test_molecule(conn, "Test Primary Molecule")
        secondary_id = create_test_molecule(conn, "Test Secondary Molecule")
        
        print(f"Created test primary molecule: {primary_id}")
        print(f"Created test secondary molecule: {secondary_id}")
        
        # Set up consolidation relationship
        print("\nEstablishing consolidation relationship...")
        with conn.cursor() as cursor:
            cursor.execute("""
                UPDATE molecules 
                SET primary_molecule_id = %s, is_secondary = TRUE
                WHERE id = %s
            """, (primary_id, secondary_id))
            conn.commit()
        print("✓ Successfully established consolidation relationship")
        
        # Test 1: Test automatic redirection to primary molecule
        print("\nTest 1: Automatic redirection to primary molecules")
        print("-" * 80)
        
        # Create a property for the secondary molecule - should be redirected
        prop_id = create_test_property(conn, secondary_id, "Test Property", "Test Value")
        
        # Verify the property was created for the primary instead
        with conn.cursor() as cursor:
            cursor.execute("""
                SELECT molecule_id FROM molecular_properties
                WHERE id = %s
            """, (prop_id,))
            result = cursor.fetchone()
            actual_molecule_id = result[0]
            
            if str(actual_molecule_id) == str(primary_id):
                print("✓ Property was correctly redirected to primary molecule")
            else:
                print(f"✗ Property was not redirected. Expected: {primary_id}, Got: {actual_molecule_id}")
        
        # Test 2: Protection against modifying secondary molecules
        print("\nTest 2: Protection against modifying secondary molecules")
        print("-" * 80)
        
        try:
            with conn.cursor() as cursor:
                cursor.execute("""
                    UPDATE molecules 
                    SET name = 'Modified Secondary Name'
                    WHERE id = %s
                """, (secondary_id,))
                conn.commit()
                print("✗ Was able to modify a secondary molecule - trigger not working!")
        except psycopg2.Error as e:
            conn.rollback()
            print(f"✓ Correctly prevented modification of secondary molecule: {e}")
        
        # Test 3: Prevention of circular references
        print("\nTest 3: Prevention of circular references")
        print("-" * 80)
        
        try:
            with conn.cursor() as cursor:
                cursor.execute("""
                    UPDATE molecules 
                    SET primary_molecule_id = %s
                    WHERE id = %s
                """, (secondary_id, primary_id))
                conn.commit()
                print("✗ Was able to create circular reference - trigger not working!")
        except psycopg2.Error as e:
            conn.rollback()
            print(f"✓ Correctly prevented circular reference: {e}")
        
        # Test 4: Protection against deleting primary with secondaries
        print("\nTest 4: Protection against deleting primary with secondaries")
        print("-" * 80)
        
        try:
            with conn.cursor() as cursor:
                cursor.execute("""
                    DELETE FROM molecules 
                    WHERE id = %s
                """, (primary_id,))
                conn.commit()
                print("✗ Was able to delete primary molecule with secondaries - trigger not working!")
        except psycopg2.Error as e:
            conn.rollback()
            print(f"✓ Correctly prevented deletion of primary molecule: {e}")
        
        # Test 5: Verify get_primary_molecule_id function works
        print("\nTest 5: Verify get_primary_molecule_id function")
        print("-" * 80)
        
        with conn.cursor() as cursor:
            # Test with secondary molecule
            cursor.execute("SELECT get_primary_molecule_id(%s)", (secondary_id,))
            result = cursor.fetchone()
            retrieved_primary_id = result[0]
            
            if str(retrieved_primary_id) == str(primary_id):
                print(f"✓ get_primary_molecule_id correctly returns primary for secondary molecule")
            else:
                print(f"✗ get_primary_molecule_id failed. Expected: {primary_id}, Got: {retrieved_primary_id}")
                
            # Test with primary molecule
            cursor.execute("SELECT get_primary_molecule_id(%s)", (primary_id,))
            result = cursor.fetchone()
            retrieved_id = result[0]
            
            if str(retrieved_id) == str(primary_id):
                print(f"✓ get_primary_molecule_id correctly returns same ID for primary molecule")
            else:
                print(f"✗ get_primary_molecule_id failed. Expected: {primary_id}, Got: {retrieved_id}")
        
        print("\nTest Summary")
        print("=" * 80)
        print("All tests completed. Check the results above for any failures.")
        
    except Exception as e:
        conn.rollback()
        print(f"Error during testing: {e}")
    finally:
        # Clean up test data
        print("\nCleaning up test data...")
        cleanup_test_data(conn, primary_id, secondary_id)
        
        # Close connection
        conn.close()

def create_test_molecule(conn, name):
    """Create a test molecule and return its ID."""
    with conn.cursor() as cursor:
        molecule_id = uuid.uuid4()
        cursor.execute("""
            INSERT INTO molecules (
                id, name, pubchem_cid, molecular_formula, smiles, 
                molecular_weight, canonical_smiles, created_at, 
                updated_at, properties
            ) VALUES (
                %s, %s, %s, %s, %s, 
                %s, %s, %s, %s, %s
            ) RETURNING id
        """, (
            molecule_id, name, 12345, 'C10H20O', 'CCCCCCCCCC(=O)',
            100.0, 'CCCCCCCCCC(=O)', datetime.now(), datetime.now(),
            json.dumps({'test': 'data'})
        ))
        conn.commit()
        return molecule_id

def create_test_property(conn, molecule_id, name, value):
    """Create a test molecular property and return its ID."""
    with conn.cursor() as cursor:
        property_id = uuid.uuid4()
        cursor.execute("""
            INSERT INTO molecular_properties (
                id, molecule_id, property_name, property_value,
                created_at, updated_at
            ) VALUES (
                %s, %s, %s, %s, %s, %s
            ) RETURNING id
        """, (
            property_id, molecule_id, name, value,
            datetime.now(), datetime.now()
        ))
        conn.commit()
        return property_id

def cleanup_test_data(conn, primary_id, secondary_id):
    """Clean up test data created during the test."""
    try:
        with conn.cursor() as cursor:
            # First remove the consolidation relationship
            cursor.execute("""
                UPDATE molecules 
                SET primary_molecule_id = NULL, is_secondary = FALSE
                WHERE id = %s
            """, (secondary_id,))
            
            # Delete properties
            cursor.execute("DELETE FROM molecular_properties WHERE molecule_id IN (%s, %s)",
                         (primary_id, secondary_id))
            
            # Now delete the test molecules
            cursor.execute("DELETE FROM molecules WHERE id IN (%s, %s)",
                         (primary_id, secondary_id))
            
            conn.commit()
    except Exception as e:
        conn.rollback()
        print(f"Error during cleanup: {e}")

if __name__ == "__main__":
    run_test()