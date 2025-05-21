#!/usr/bin/env python3
"""
Script to check database schema for the molecules table.
"""

import os
import sys
import psycopg2
from psycopg2.extras import RealDictCursor

# Database connection parameters
DB_PARAMS = {
    'host': 'aws-0-us-east-1.pooler.supabase.com',
    'port': 5432,
    'dbname': 'postgres',
    'user': 'postgres.tsdlmynydfuypiugmkev',
    'password': 'LDHt$rkaM&Gmf3X@LQ37',
    'sslmode': 'require'
}

def get_db_connection():
    """Get database connection using psycopg2."""
    return psycopg2.connect(**DB_PARAMS)

def check_table_schema(table_name):
    """Check schema of a specific table."""
    conn = get_db_connection()
    cursor = conn.cursor(cursor_factory=RealDictCursor)
    
    # Get column information
    cursor.execute("""
        SELECT column_name, data_type, is_nullable
        FROM information_schema.columns
        WHERE table_name = %s
        ORDER BY ordinal_position
    """, (table_name,))
    
    columns = cursor.fetchall()
    
    print(f"Schema for table '{table_name}':")
    for col in columns:
        print(f"  {col['column_name']} ({col['data_type']}, nullable: {col['is_nullable']})")
    
    cursor.close()
    conn.close()

def count_rows(table_name):
    """Count rows in a specific table."""
    conn = get_db_connection()
    cursor = conn.cursor()
    
    cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
    count = cursor.fetchone()[0]
    
    print(f"Table '{table_name}' has {count} rows.")
    
    cursor.close()
    conn.close()

if __name__ == "__main__":
    try:
        # Check molecules table schema
        check_table_schema("molecules")
        count_rows("molecules")
        
        # Check molecular_properties table schema
        check_table_schema("molecular_properties")
        count_rows("molecular_properties")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)