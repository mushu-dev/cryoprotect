#!/usr/bin/env python3
"""
Check the structure of the mixture_components table.
"""

import os
import sys
import psycopg2
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

def main():
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }

    # Ensure required parameters are present
    if not all([db_params['host'], db_params['user'], db_params['password']]):
        print("Error: Missing required database connection parameters in environment variables.")
        print("Make sure you have a valid .env file with SUPABASE_DB_* variables.")
        sys.exit(1)

    conn = None
    try:
        # Connect to the database
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()

        # Check table structure: mixture_components
        cursor.execute("""
            SELECT column_name, data_type, is_nullable
            FROM information_schema.columns
            WHERE table_name = 'mixture_components'
            ORDER BY ordinal_position
        """)

        columns = cursor.fetchall()

        print("Table structure for 'mixture_components':")
        for column in columns:
            name, data_type, is_nullable = column
            print(f"  {name}: {data_type} {'' if is_nullable == 'YES' else '(NOT NULL)'}")

        # Get a sample row
        cursor.execute("""
            SELECT *
            FROM mixture_components
            LIMIT 1
        """)

        row = cursor.fetchone()
        if row:
            print("\nSample row:")
            for i, column in enumerate(columns):
                print(f"  {column[0]}: {row[i]}")
        else:
            print("\nNo rows found in the mixture_components table.")

        # Check table structure: mixtures
        cursor.execute("""
            SELECT column_name, data_type, is_nullable
            FROM information_schema.columns
            WHERE table_name = 'mixtures'
            ORDER BY ordinal_position
        """)

        columns = cursor.fetchall()

        print("\nTable structure for 'mixtures':")
        for column in columns:
            name, data_type, is_nullable = column
            print(f"  {name}: {data_type} {'' if is_nullable == 'YES' else '(NOT NULL)'}")

        # Get a sample row
        cursor.execute("""
            SELECT *
            FROM mixtures
            LIMIT 1
        """)

        row = cursor.fetchone()
        if row:
            print("\nSample row:")
            for i, column in enumerate(columns):
                print(f"  {column[0]}: {row[i]}")
        else:
            print("\nNo rows found in the mixtures table.")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    main()