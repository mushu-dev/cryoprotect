#!/usr/bin/env python3
"""
Apply the consolidated molecules view to the database.

This script:
1. Creates a view for consolidated molecules
2. Creates a view for primary molecules
3. Creates utility functions for working with consolidated molecules
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

def apply_sql_from_file(conn, filename):
    """
    Apply SQL from a file to the database.
    
    Args:
        conn: Database connection
        filename: Path to the SQL file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        with open(filename, 'r') as f:
            sql = f.read()
        
        with conn.cursor() as cursor:
            cursor.execute(sql)
        
        conn.commit()
        return True
    except Exception as e:
        print(f"Error applying SQL: {e}")
        conn.rollback()
        return False

def test_consolidated_molecules_view(conn):
    """
    Test the consolidated molecules view.
    
    Args:
        conn: Database connection
        
    Returns:
        True if the view exists and works, False otherwise
    """
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("""
                SELECT COUNT(*) as total_molecules,
                       SUM(CASE WHEN is_consolidated THEN 1 ELSE 0 END) as consolidated_count,
                       SUM(CASE WHEN molecule_status = 'Primary' THEN 1 ELSE 0 END) as primary_count,
                       SUM(CASE WHEN molecule_status = 'Unique' THEN 1 ELSE 0 END) as unique_count
                FROM consolidated_molecules
            """)
            
            result = cursor.fetchone()
            
            if result:
                print("Consolidated Molecules View Test:")
                print(f"  Total molecules: {result['total_molecules']}")
                print(f"  Consolidated molecules: {result['consolidated_count']}")
                print(f"  Primary molecules: {result['primary_count']}")
                print(f"  Unique molecules: {result['unique_count']}")
                
                # Test the primary molecules view
                cursor.execute("SELECT COUNT(*) as primary_total FROM primary_molecules")
                primary_result = cursor.fetchone()
                
                if primary_result:
                    print(f"  Primary molecules view count: {primary_result['primary_total']}")
                    
                    # Verify counts match
                    if primary_result['primary_total'] == (result['primary_count'] + result['unique_count']):
                        print("  Count verification: SUCCESS")
                    else:
                        print("  Count verification: FAILED - counts don't match")
                        return False
                    
                # Test the get_primary_molecule_id function
                cursor.execute("""
                    SELECT id, name FROM consolidated_molecules 
                    WHERE is_consolidated = TRUE 
                    LIMIT 1
                """)
                consolidated = cursor.fetchone()
                
                if consolidated:
                    consolidated_id = consolidated['id']
                    cursor.execute("""
                        SELECT get_primary_molecule_id(%s) as primary_id
                    """, (consolidated_id,))
                    
                    func_result = cursor.fetchone()
                    
                    if func_result:
                        primary_id = func_result['primary_id']
                        print(f"  Testing get_primary_molecule_id with {consolidated['name']} ({consolidated_id})")
                        print(f"    Primary ID: {primary_id}")
                        
                        # Verify primary ID matches
                        cursor.execute("""
                            SELECT primary_molecule_id FROM consolidated_molecules 
                            WHERE id = %s
                        """, (consolidated_id,))
                        
                        check_result = cursor.fetchone()
                        if check_result and check_result['primary_molecule_id'] == primary_id:
                            print("    Function verification: SUCCESS")
                        else:
                            print("    Function verification: FAILED - IDs don't match")
                            return False
                
                return True
            else:
                print("No results returned from consolidated molecules view")
                return False
    except Exception as e:
        print(f"Error testing view: {e}")
        return False

def main():
    """Apply the consolidated molecules view to the database."""
    # Connect to database
    print("Connecting to database...")
    conn = connect_to_db()
    
    try:
        # Apply SQL from file
        sql_file = "create_consolidated_molecules_view.sql"
        print(f"Applying SQL from {sql_file}...")
        
        if apply_sql_from_file(conn, sql_file):
            print("SQL applied successfully")
            
            # Test the view
            print("\nTesting consolidated molecules view...")
            if test_consolidated_molecules_view(conn):
                print("\nConsolidated molecules view is working correctly!")
                print("\nYou can now use the following in your queries:")
                print("- consolidated_molecules view: All molecules with consolidation info")
                print("- primary_molecules view: Only primary and unique molecules")
                print("- get_primary_molecule_id(uuid): Function to get primary ID")
                print("- is_primary_molecule(uuid): Function to check if a molecule is primary")
            else:
                print("\nFailed to verify consolidated molecules view")
                sys.exit(1)
        else:
            print("Failed to apply SQL")
            sys.exit(1)
    finally:
        conn.close()

if __name__ == "__main__":
    main()