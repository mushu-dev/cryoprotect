#!/usr/bin/env python3
"""
Standardize concentration units in mixture components.

This script identifies mixture components with non-standard units (mol/L)
and converts them to the standard unit (%w/v) for consistent comparison.
"""

import os
import sys
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Conversion factor for specific known compounds (mol/L to %w/v)
# Formula: %w/v = (mol/L * molecular_weight) / 10
MOLECULAR_WEIGHTS = {
    # Common cryoprotectants
    "DMSO": 78.13,         # Dimethyl sulfoxide
    "Glycerol": 92.09,      # Glycerol
    "Ethylene glycol": 62.07, # Ethylene glycol
    "Propylene glycol": 76.10, # Propylene glycol
    "Sucrose": 342.30,      # Sucrose
    "Trehalose": 342.30,    # Trehalose
    "Mannitol": 182.17,     # Mannitol
    "Dextran": 1000.0,      # Approximate - varies by size
    # Components used in specific mixtures like VS55
    "Formamide": 45.04,     # Formamide
    "2,3-Butanediol": 90.12  # 2,3-Butanediol
}

# Default molecular weight for unknown compounds
DEFAULT_MOLECULAR_WEIGHT = 100.0

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

def get_molecule_name(conn, molecule_id):
    """Get the name of a molecule by its ID."""
    with conn.cursor() as cursor:
        cursor.execute("SELECT name FROM molecules WHERE id = %s", (molecule_id,))
        result = cursor.fetchone()
        return result[0] if result else "Unknown"

def get_non_standard_components(conn):
    """Find mixture components using non-standard units (mol/L)."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT mc.id, mc.mixture_id, m.name as mixture_name, 
                   mc.molecule_id, mc.concentration, mc.concentration_unit, mc.role
            FROM mixture_components mc
            JOIN mixtures m ON mc.mixture_id = m.id
            WHERE mc.concentration_unit != '%w/v'
            ORDER BY m.name
        """)
        return cursor.fetchall()

def convert_molar_to_percent(concentration, molecular_weight):
    """Convert concentration from mol/L to %w/v."""
    # Convert to float to avoid decimal/float type issues
    concentration_float = float(concentration)
    # Formula: %w/v = (mol/L * molecular_weight) / 10
    return (concentration_float * molecular_weight) / 10.0

def update_component_unit(conn, component_id, new_concentration, new_unit):
    """Update the concentration and unit of a mixture component."""
    with conn.cursor() as cursor:
        cursor.execute("""
            UPDATE mixture_components 
            SET concentration = %s, concentration_unit = %s, updated_at = NOW() 
            WHERE id = %s
        """, (new_concentration, new_unit, component_id))
        return cursor.rowcount

def main():
    """Standardize concentration units in mixture components."""
    print("Standardizing mixture component units...")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False  # Use transaction
    
    try:
        # Get components with non-standard units
        components = get_non_standard_components(conn)
        print(f"Found {len(components)} components with non-standard units")
        
        if not components:
            print("All mixture components already have standardized units.")
            return
        
        standardized_count = 0
        
        # Process each component
        for component in components:
            # Get molecule name for molecular weight lookup
            molecule_name = get_molecule_name(conn, component['molecule_id'])
            print(f"\nProcessing component: {molecule_name} in mixture {component['mixture_name']}")
            print(f"  Current: {component['concentration']} {component['concentration_unit']}")
            
            # Skip components that aren't in mol/L
            if component['concentration_unit'] != 'mol/L':
                print(f"  Skipping: Unit '{component['concentration_unit']}' not convertible to %w/v")
                continue
            
            # Get molecular weight
            molecular_weight = MOLECULAR_WEIGHTS.get(molecule_name, DEFAULT_MOLECULAR_WEIGHT)
            if molecule_name not in MOLECULAR_WEIGHTS:
                print(f"  Warning: Unknown molecular weight for '{molecule_name}', using default {DEFAULT_MOLECULAR_WEIGHT}")
            
            # Convert concentration
            new_concentration = convert_molar_to_percent(component['concentration'], molecular_weight)
            new_unit = '%w/v'
            
            # Update component
            updated = update_component_unit(conn, component['id'], new_concentration, new_unit)
            if updated:
                standardized_count += 1
                print(f"  Updated to: {new_concentration:.2f} {new_unit}")
            else:
                print(f"  Failed to update component {component['id']}")
        
        # Commit changes
        conn.commit()
        print(f"\nStandardized {standardized_count} of {len(components)} components")
    
    except Exception as e:
        conn.rollback()
        print(f"Error standardizing units: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()