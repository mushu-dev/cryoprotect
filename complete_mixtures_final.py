#!/usr/bin/env python3
"""
Simplified script to complete mixtures by adding components.
Focuses on just adding components without additional features.
"""

import os
import sys
import uuid
import psycopg2
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Simplified mixture components
MIXTURES = {
    "Betaine/Glycerol Solution": [
        {"name": "glycerol", "concentration": 10, "unit": "%w/v"},
        {"name": "betaine", "concentration": 5, "unit": "%w/v"}
    ],
    "DMSO/Proline Solution": [
        {"name": "dimethyl sulfoxide", "concentration": 10, "unit": "%w/v"},
        {"name": "proline", "concentration": 5, "unit": "%w/v"}
    ],
    "DP6 Vitrification Solution": [
        {"name": "dimethyl sulfoxide", "concentration": 6, "unit": "%w/v"},
        {"name": "propylene glycol", "concentration": 6, "unit": "%w/v"}
    ],
    "EAFS10/10 Solution": [
        {"name": "ethylene glycol", "concentration": 10, "unit": "%w/v"},
        {"name": "acetamide", "concentration": 10, "unit": "%w/v"},
        {"name": "formamide", "concentration": 10, "unit": "%w/v"},
        {"name": "dimethyl sulfoxide", "concentration": 10, "unit": "%w/v"}
    ],
    "Glycerol/Trehalose Solution": [
        {"name": "glycerol", "concentration": 15, "unit": "%w/v"},
        {"name": "trehalose", "concentration": 15, "unit": "%w/v"}
    ],
    "Hydroxyectoine/DMSO Solution": [
        {"name": "dimethyl sulfoxide", "concentration": 7.5, "unit": "%w/v"},
        {"name": "hydroxyectoine", "concentration": 5, "unit": "%w/v"}
    ],
    "M22 Vitrification Solution": [
        {"name": "dimethyl sulfoxide", "concentration": 22.3, "unit": "%w/v"},
        {"name": "formamide", "concentration": 12.8, "unit": "%w/v"},
        {"name": "ethylene glycol", "concentration": 16.8, "unit": "%w/v"}
    ]
}

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
        conn.autocommit = False  # Start a transaction
        cursor = conn.cursor()
        
        # Process each mixture
        for mixture_name, components in MIXTURES.items():
            print(f"Processing {mixture_name}...")
            
            # Get mixture ID
            cursor.execute("""
                SELECT id FROM mixtures WHERE name = %s
            """, (mixture_name,))
            
            mixture_row = cursor.fetchone()
            if not mixture_row:
                print(f"  Mixture '{mixture_name}' not found in database")
                continue
                
            mixture_id = mixture_row[0]
            
            # Check if mixture already has components
            cursor.execute("""
                SELECT COUNT(*) FROM mixture_components WHERE mixture_id = %s
            """, (mixture_id,))
            
            component_count = cursor.fetchone()[0]
            if component_count > 0:
                print(f"  Mixture '{mixture_name}' already has {component_count} components - skipping")
                continue
            
            # Add each component
            for component in components:
                # Find molecule
                cursor.execute("""
                    SELECT id FROM molecules WHERE LOWER(name) LIKE %s
                """, (f"%{component['name'].lower()}%",))
                
                molecule_row = cursor.fetchone()
                if not molecule_row:
                    print(f"  Could not find molecule '{component['name']}' - skipping")
                    continue
                    
                molecule_id = molecule_row[0]
                
                # Add component
                component_id = str(uuid.uuid4())
                cursor.execute("""
                    INSERT INTO mixture_components (
                        id, mixture_id, molecule_id, concentration, concentration_unit, 
                        role, created_at, updated_at
                    ) VALUES (
                        %s, %s, %s, %s, %s, 
                        %s, NOW(), NOW()
                    )
                """, (
                    component_id,
                    mixture_id,
                    molecule_id,
                    component["concentration"],
                    component["unit"],
                    "cryoprotectant"
                ))
                
                print(f"  Added component '{component['name']}' at {component['concentration']} {component['unit']}")
            
        # Commit all changes
        conn.commit()
        print("All changes committed successfully")
            
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