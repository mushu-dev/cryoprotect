#!/usr/bin/env python3
"""
Check molecule names to confirm standardization.
"""

import os
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Connect to database
db_params = {
    'host': os.getenv('SUPABASE_DB_HOST'),
    'port': os.getenv('SUPABASE_DB_PORT', '5432'),
    'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
    'user': os.getenv('SUPABASE_DB_USER'),
    'password': os.getenv('SUPABASE_DB_PASSWORD'),
    'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
}

with psycopg2.connect(**db_params) as conn:
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Check for specific molecule names
        names_to_check = [
            "Dimethyl sulfoxide",  # Was "DMSO"
            "Ethylene glycol",     # Was "ethane-1,2-diol"
            "Glycerol",            # Was "propane-1,2,3-triol"
            "Propylene glycol",    # Was "propane-1,2-diol"
            "Trehalose",           # Was "Trehalose" (capitalization only)
            "Methanol"             # Was "methanol" (capitalization only)
        ]
        
        print("Checking for standardized molecule names:")
        print("------------------------------------------")
        
        for name in names_to_check:
            cursor.execute("""
                SELECT id, name, molecular_formula, pubchem_cid, properties 
                FROM molecules 
                WHERE name = %s
            """, (name,))
            
            rows = cursor.fetchall()
            
            print(f"\nName: {name}")
            print(f"Found {len(rows)} molecules")
            
            for row in rows:
                props_str = ""
                if row['properties'] and 'original_name' in row['properties']:
                    props_str = f" (Original name: {row['properties']['original_name']})"
                
                print(f"  ID: {row['id']}, Formula: {row['molecular_formula']}, PubChem CID: {row['pubchem_cid']}{props_str}")
                
        # Check for None names that should have been updated
        print("\n\nChecking if any 'None' names remain:")
        print("------------------------------------------")
        
        cursor.execute("""
            SELECT COUNT(*) as count
            FROM molecules
            WHERE name = 'None'
        """)
        
        none_count = cursor.fetchone()['count']
        print(f"Molecules with 'None' as name: {none_count}")
        
        if none_count > 0:
            print("Examples:")
            cursor.execute("""
                SELECT id, name, molecular_formula, pubchem_cid
                FROM molecules
                WHERE name = 'None'
                LIMIT 5
            """)
            
            for row in cursor.fetchall():
                print(f"  ID: {row['id']}, Formula: {row['molecular_formula']}, PubChem CID: {row['pubchem_cid']}")