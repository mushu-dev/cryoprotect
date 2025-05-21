#!/usr/bin/env python3
"""
Check the status of molecular properties in the CryoProtect database
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

def main():
    conn = connect_to_db()
    
    try:
        with conn.cursor() as cursor:
            # Get total count of molecular properties
            cursor.execute("SELECT COUNT(*) FROM molecular_properties")
            total_count = cursor.fetchone()[0]
            print(f"Total molecular properties: {total_count}")
            
            # Get property counts by type
            cursor.execute("""
                SELECT property_type, COUNT(*) 
                FROM molecular_properties 
                GROUP BY property_type 
                ORDER BY COUNT(*) DESC 
                LIMIT 10
            """)
            
            print("\nProperties by type:")
            for row in cursor.fetchall():
                print(f"  {row[0]}: {row[1]}")
            
            # Get counts by molecule
            cursor.execute("""
                SELECT COUNT(DISTINCT molecule_id) 
                FROM molecular_properties
            """)
            distinct_molecules = cursor.fetchone()[0]
            print(f"\nDistinct molecules with properties: {distinct_molecules}")
            
            # Get total molecules count
            cursor.execute("SELECT COUNT(*) FROM molecules")
            total_molecules = cursor.fetchone()[0]
            print(f"Total molecules: {total_molecules}")
            
            # Calculate coverage
            coverage = (distinct_molecules / total_molecules) * 100
            print(f"Property coverage: {coverage:.2f}%")
            
            # Check for molecules with missing common properties
            cursor.execute("""
                SELECT COUNT(*) 
                FROM molecules m
                WHERE NOT EXISTS (
                    SELECT 1 FROM molecular_properties mp 
                    WHERE mp.molecule_id = m.id 
                    AND mp.property_type = 'molecular_weight'
                )
            """)
            molecules_missing_mw = cursor.fetchone()[0]
            print(f"\nMolecules missing molecular weight: {molecules_missing_mw}")
            
            cursor.execute("""
                SELECT COUNT(*) 
                FROM molecules m
                WHERE NOT EXISTS (
                    SELECT 1 FROM molecular_properties mp 
                    WHERE mp.molecule_id = m.id 
                    AND mp.property_type = 'logp'
                )
            """)
            molecules_missing_logp = cursor.fetchone()[0]
            print(f"Molecules missing LogP: {molecules_missing_logp}")

    finally:
        conn.close()

if __name__ == "__main__":
    main()