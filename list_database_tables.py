#!/usr/bin/env python3
"""
List all tables in the database to understand the correct schema.
"""

import os
import sys
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

def get_tables(conn):
    """Get all tables in the public schema."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT table_name, table_type
            FROM information_schema.tables
            WHERE table_schema = 'public'
            ORDER BY table_name
        """)
        return cursor.fetchall()

def get_table_row_count(conn, table_name):
    """Get the row count for a table."""
    try:
        with conn.cursor() as cursor:
            cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
            return cursor.fetchone()[0]
    except Exception as e:
        return f"Error counting rows: {str(e)}"

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
    # Connect to database
    print("Connecting to database...")
    try:
        conn = connect_to_db()
    except Exception as e:
        print(f"Error connecting to database: {str(e)}")
        sys.exit(1)

    try:
        # Get all tables
        tables = get_tables(conn)
        print(f"Found {len(tables)} tables in the public schema:")
        
        # Print table information
        for table in tables:
            table_name = table['table_name']
            table_type = table['table_type']
            row_count = get_table_row_count(conn, table_name)
            
            print(f"\n- {table_name} ({table_type}, {row_count} rows)")
            
            # Get and print column information
            columns = get_table_columns(conn, table_name)
            if columns:
                for col in columns:
                    nullable = "NULL" if col['is_nullable'] == 'YES' else "NOT NULL"
                    print(f"  - {col['column_name']} ({col['data_type']}, {nullable})")
            else:
                print(f"  No columns available for {table_name}")
    
    except Exception as e:
        print(f"Error listing tables: {str(e)}")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()