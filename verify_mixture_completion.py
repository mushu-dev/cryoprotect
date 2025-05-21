#!/usr/bin/env python3
"""
Verify that mixtures have been completed with components.
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
        
        # Get mixture component counts
        cursor.execute("""
            SELECT 
                m.name AS mixture_name,
                COUNT(mc.id) AS component_count
            FROM 
                mixtures m
            LEFT JOIN 
                mixture_components mc ON m.id = mc.mixture_id
            GROUP BY 
                m.name
            ORDER BY 
                component_count DESC, m.name
        """)
        
        mixtures = cursor.fetchall()
        
        print("Mixture component counts:")
        for mixture_name, component_count in mixtures:
            status = "✓" if component_count > 0 else "✗"
            print(f"  {status} {mixture_name}: {component_count} components")
        
        # Get detailed mixtures and components
        cursor.execute("""
            SELECT 
                m.name AS mixture_name,
                mol.name AS molecule_name,
                mc.concentration,
                mc.concentration_unit
            FROM 
                mixtures m
            JOIN 
                mixture_components mc ON m.id = mc.mixture_id
            JOIN 
                molecules mol ON mc.molecule_id = mol.id
            ORDER BY 
                m.name, mol.name
        """)
        
        components = cursor.fetchall()
        
        print("\nDetailed mixture compositions:")
        current_mixture = None
        for mixture_name, molecule_name, concentration, unit in components:
            if mixture_name != current_mixture:
                if current_mixture is not None:
                    print()
                print(f"  {mixture_name}:")
                current_mixture = mixture_name
            
            print(f"    - {molecule_name}: {concentration} {unit}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    main()