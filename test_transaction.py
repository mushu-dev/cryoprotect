#!/usr/bin/env python3
"""
Test transactions and RLS on database.
"""

import os
from dotenv import load_dotenv
from database import db

def main():
    # Load environment variables
    load_dotenv()
    
    # Initialize database connection
    config = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
        'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }
    
    print(f"Initializing database connection to {config['host']}:{config['port']}/{config['database']}")
    db.init_connection_pool(config=config)
    
    # Test transaction commit
    test_id = None
    print("Testing transaction commit:")
    with db.transaction() as cursor:
        cursor.execute("""
            INSERT INTO molecules 
            (chembl_id, name, pubchem_cid, created_at, updated_at) 
            VALUES (%s, %s, %s, NOW(), NOW()) RETURNING id
        """, ('TEST_TRANS', 'Transaction Test', 12345))
        result = cursor.fetchone()
        print("Insert result:", result)
        test_id = result['id']
        
        # Verify within same transaction
        cursor.execute("SELECT id FROM molecules WHERE chembl_id = %s", ('TEST_TRANS',))
        print("Immediate verify in same transaction:", cursor.fetchone())
    
    # Verify after transaction commit
    result = db.execute_query("SELECT id FROM molecules WHERE chembl_id = %s", ('TEST_TRANS',))
    print("Verify after transaction commit:", result)
    
    # Test if we can find by UUID directly 
    if test_id:
        result = db.execute_query("SELECT id FROM molecules WHERE id = %s", (test_id,))
        print("Verify by UUID:", result)
    
    # Test RLS
    print("\nTesting ROW LEVEL SECURITY:")
    try:
        rls_enabled = db.execute_query("SHOW enable_row_level_security;")
        print("enable_row_level_security:", rls_enabled)
    except Exception as e:
        print("Error checking enable_row_level_security:", str(e))
    
    try:
        check_tables = db.execute_query("""
            SELECT 
                table_name, table_schema, 
                relpersistence as persistence, 
                relrowsecurity as rls_enabled
            FROM pg_class c
            JOIN pg_namespace n ON c.relnamespace = n.oid
            JOIN information_schema.tables t ON t.table_name = c.relname AND t.table_schema = n.nspname
            WHERE c.relkind = 'r' AND t.table_schema = 'public' AND t.table_name = 'molecules'
        """)
        print("RLS settings for molecules table:", check_tables)
    except Exception as e:
        print("Error checking tables RLS settings:", str(e))
    
    # Check RLS policies 
    try:
        policies = db.execute_query("""
            SELECT 
                schemaname, tablename, 
                policyname, cmd, 
                roles, qual
            FROM pg_policies
            WHERE tablename = 'molecules'
        """)
        print("\nRLS Policies for molecules table:")
        for policy in policies:
            print(f"  {policy['policyname']}: {policy['cmd']} {policy['qual']}")
    except Exception as e:
        print("Error checking policies:", str(e))
    
    # Clean up test data 
    print("\nCleaning up test data...")
    db.execute_query("DELETE FROM molecules WHERE chembl_id LIKE 'TEST%'")
    
    print("Test complete.")
    db.close_all_connections()

if __name__ == "__main__":
    main()