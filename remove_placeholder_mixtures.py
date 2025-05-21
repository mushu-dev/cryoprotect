#!/usr/bin/env python3
"""
Remove placeholder mixtures with no components and no scientific value.
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
    
    # List of mixtures to keep even if they have no components
    KEEP_MIXTURES = [
        "VS55 Vitrification Solution",
        "M22 Vitrification Solution",
        "DP6 Vitrification Solution",
        "EAFS10/10 Solution",
        "Glycerol/Trehalose Solution",
        "Betaine/Glycerol Solution",
        "DMSO/Proline Solution",
        "Hydroxyectoine/DMSO Solution"
    ]
    
    conn = None
    try:
        # Connect to the database
        conn = psycopg2.connect(**db_params)
        conn.autocommit = False  # Start a transaction
        cursor = conn.cursor()
        
        # Find placeholders (mixtures with no components)
        cursor.execute("""
            SELECT m.id, m.name
            FROM mixtures m
            LEFT JOIN mixture_components mc ON m.id = mc.mixture_id
            WHERE mc.id IS NULL
            ORDER BY m.name
        """)
        
        placeholder_mixtures = cursor.fetchall()
        print(f"Found {len(placeholder_mixtures)} placeholders:")
        
        removed_count = 0
        for mixture_id, mixture_name in placeholder_mixtures:
            # Skip important mixtures
            if mixture_name in KEEP_MIXTURES:
                print(f"  Skipping '{mixture_name}' (on keep list)")
                continue
            
            # Check for references
            cursor.execute("""
                SELECT 'experiments' as table_name, COUNT(*) as ref_count 
                FROM experiments 
                WHERE mixture_id = %s
                UNION ALL
                SELECT 'predictions' as table_name, COUNT(*) as ref_count 
                FROM predictions 
                WHERE mixture_id = %s
            """, (mixture_id, mixture_id))
            
            references = cursor.fetchall()
            
            has_references = False
            for table_name, ref_count in references:
                if ref_count > 0:
                    print(f"  Cannot remove '{mixture_name}': Referenced by {ref_count} rows in {table_name}")
                    has_references = True
                    break
            
            if has_references:
                continue
            
            # Remove the mixture
            cursor.execute("""
                DELETE FROM mixtures
                WHERE id = %s
            """, (mixture_id,))
            
            print(f"  Removed '{mixture_name}'")
            removed_count += 1
            
        # Commit all changes
        conn.commit()
        print(f"Removed {removed_count} placeholder mixtures")
        
    except Exception as e:
        print(f"Error: {e}")
        if conn:
            conn.rollback()
        sys.exit(1)
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    main()