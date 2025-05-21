#!/usr/bin/env python3
"""
Check molecules with missing molecular formulas in the CryoProtect database
"""
import os
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Connect to the database
def connect_to_db():
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    
    return psycopg2.connect(**db_params)

def delete_test_molecule(conn, molecule_id):
    """Delete a test molecule with the given ID"""
    with conn.cursor() as cursor:
        # Delete from molecular_properties first (foreign key constraint)
        cursor.execute("DELETE FROM molecular_properties WHERE molecule_id = %s", (molecule_id,))
        
        # Delete from molecules
        cursor.execute("DELETE FROM molecules WHERE id = %s", (molecule_id,))
        
        conn.commit()
        print(f"Deleted molecule with ID: {molecule_id}")

def main():
    conn = connect_to_db()
    
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Find molecules with missing formulas but valid SMILES
            cursor.execute("""
                SELECT id, name, smiles, data_source, pubchem_cid, properties 
                FROM molecules 
                WHERE molecular_formula IS NULL 
                  AND smiles IS NOT NULL
            """)
            
            results = cursor.fetchall()
            print("Molecules with missing formulas but valid SMILES:")
            for row in results:
                print(f"ID: {row['id']}")
                print(f"  Name: {row['name']}")
                print(f"  SMILES: {row['smiles']}")
                print(f"  Source: {row['data_source']}")
                print(f"  PubChem CID: {row['pubchem_cid']}")
                print(f"  Properties: {row['properties']}")
                print()
            
            # Find all molecules with missing formulas
            cursor.execute("""
                SELECT id, name, smiles, data_source, pubchem_cid, properties
                FROM molecules 
                WHERE molecular_formula IS NULL
            """)
            
            results = cursor.fetchall()
            print("\nAll molecules with missing formulas:")
            for row in results:
                print(f"ID: {row['id']}")
                print(f"  Name: {row['name']}")
                print(f"  SMILES: {row['smiles']}")
                print(f"  Source: {row['data_source']}")
                print(f"  PubChem CID: {row['pubchem_cid']}")
                print(f"  Properties: {row['properties']}")
                print()
            
            # Automatically delete test molecules
            print("Automatically deleting test molecules with missing formulas:")
            for row in results:
                is_test = False
                
                # Check if it's a test molecule based on name, source, or properties
                if row['data_source'] and 'test' in row['data_source'].lower():
                    is_test = True
                elif row['name'] and 'test' in row['name'].lower():
                    is_test = True
                elif row['properties'] and isinstance(row['properties'], dict) and row['properties'].get('is_test'):
                    is_test = True
                
                if is_test:
                    print(f"Deleting test molecule: {row['name']}")
                    delete_test_molecule(conn, row['id'])
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()