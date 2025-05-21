#!/usr/bin/env python3
"""
Check the status of the property completion process in the CryoProtect database
"""
import os
import psycopg2
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
            # Check molecules with completed properties
            for property_type in ['molecular_weight', 'logp', 'tpsa', 'h_donors', 'h_acceptors', 'rotatable_bonds']:
                cursor.execute(f"""
                    SELECT COUNT(*) 
                    FROM molecular_properties 
                    WHERE property_type = '{property_type}'
                """)
                count = cursor.fetchone()[0]
                print(f"Molecules with {property_type}: {count}")
            
            # Check most recently added properties
            cursor.execute("""
                SELECT property_type, molecule_id, created_at
                FROM molecular_properties
                ORDER BY created_at DESC
                LIMIT 5
            """)
            
            print("\nMost recently added properties:")
            for row in cursor.fetchall():
                print(f"  {row[0]} for molecule {row[1]} at {row[2]}")
            
            # Check which molecules still need property completion
            cursor.execute("""
                SELECT COUNT(*) 
                FROM molecules m
                WHERE m.smiles IS NOT NULL
                AND NOT EXISTS (
                    SELECT 1 FROM molecular_properties mp 
                    WHERE mp.molecule_id = m.id 
                    AND mp.property_type = 'molecular_weight'
                )
            """)
            remaining = cursor.fetchone()[0]
            print(f"\nMolecules still needing property completion: {remaining}")

    finally:
        conn.close()

if __name__ == "__main__":
    main()