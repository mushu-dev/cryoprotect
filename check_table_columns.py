#!/usr/bin/env python3
"""
Check table columns in the database.
"""

import os
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

def get_table_columns(conn, table_name):
    """Get columns for a table."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT column_name, data_type, is_nullable
            FROM information_schema.columns
            WHERE table_schema = 'public' AND table_name = %s
            ORDER BY ordinal_position
        """, (table_name,))
        return cursor.fetchall()

def main():
    # Tables to check
    tables = [
        'molecules',
        'molecular_properties',
        'mixture_components', 
        'mixtures',
        'experiment_properties',
        'predictions'
    ]
    
    # Connect to database
    print("Connecting to database...")
    conn = connect_to_db()
    
    try:
        for table in tables:
            print(f"\nColumns for table '{table}':")
            columns = get_table_columns(conn, table)
            if columns:
                for col in columns:
                    nullable = "NULL" if col['is_nullable'] == 'YES' else "NOT NULL"
                    print(f"  {col['column_name']} ({col['data_type']}, {nullable})")
            else:
                print(f"  Table {table} not found or no columns available")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()