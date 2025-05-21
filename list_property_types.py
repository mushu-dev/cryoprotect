#!/usr/bin/env python3
"""
List distinct property types in the CryoProtect database
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
            # List all distinct property types
            cursor.execute("SELECT DISTINCT property_type FROM molecular_properties")
            
            print("Distinct property types:")
            for row in cursor.fetchall():
                print(f"  {row[0]}")
            
            # Check case-sensitivity issues
            cursor.execute("""
                SELECT property_type, COUNT(*) 
                FROM molecular_properties 
                WHERE LOWER(property_type) LIKE '%weight%' 
                   OR LOWER(property_type) LIKE '%mass%'
                GROUP BY property_type
            """)
            
            print("\nWeight-related property types:")
            for row in cursor.fetchall():
                print(f"  {row[0]}: {row[1]}")
            
            cursor.execute("""
                SELECT property_type, COUNT(*) 
                FROM molecular_properties 
                WHERE LOWER(property_type) LIKE '%logp%'
                GROUP BY property_type
            """)
            
            print("\nLogP-related property types:")
            for row in cursor.fetchall():
                print(f"  {row[0]}: {row[1]}")

    finally:
        conn.close()

if __name__ == "__main__":
    main()