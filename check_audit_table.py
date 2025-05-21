#!/usr/bin/env python3
"""
Check if the scientific_data_audit table exists and create it if needed.
"""

import os
import sys
import psycopg2
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

def main():
    # Get database connection parameters from environment
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
        
        # Check if the table exists
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = 'public'
                AND table_name = 'scientific_data_audit'
            )
        """)
        table_exists = cursor.fetchone()[0]
        
        if table_exists:
            print("Table 'scientific_data_audit' already exists.")
            
            # Check table structure
            cursor.execute("""
                SELECT column_name, data_type
                FROM information_schema.columns
                WHERE table_name = 'scientific_data_audit'
                ORDER BY ordinal_position
            """)
            columns = cursor.fetchall()
            
            print("\nTable structure:")
            for column in columns:
                print(f"  {column[0]}: {column[1]}")
                
            # Count records
            cursor.execute("SELECT COUNT(*) FROM scientific_data_audit")
            count = cursor.fetchone()[0]
            print(f"\nCurrent record count: {count}")
            
            sys.exit(0)
        
        print("Table 'scientific_data_audit' does not exist. Creating it now...")
        
        # Create the table
        cursor.execute("""
            CREATE TABLE scientific_data_audit (
                id SERIAL PRIMARY KEY,
                operation_type VARCHAR(50) NOT NULL,
                table_name VARCHAR(100) NOT NULL,
                operation_details TEXT,
                performed_by VARCHAR(100),
                performed_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
                affected_rows INTEGER,
                operation_status VARCHAR(50) DEFAULT 'COMPLETED',
                operation_metadata JSONB
            )
        """)
        conn.commit()
        
        print("Table 'scientific_data_audit' created successfully.")
        sys.exit(0)
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    main()