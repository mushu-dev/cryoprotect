#!/usr/bin/env python3
"""
Create property types for critical molecular properties.

This script creates the necessary property types for essential molecular properties:
- logP (lipophilicity)
- h_bond_donors (hydrogen bond donors)
- h_bond_acceptors (hydrogen bond acceptors)
"""

import os
import sys
import uuid
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Critical property types to create
PROPERTY_TYPES = [
    {
        "name": "logP",
        "description": "Octanol-water partition coefficient (lipophilicity)",
        "data_type": "float",
        "units": "log units",
        "min_value": "-10.0",
        "max_value": "10.0"
    },
    {
        "name": "h_bond_donors",
        "description": "Number of hydrogen bond donors in the molecule",
        "data_type": "integer",
        "units": "count",
        "min_value": "0",
        "max_value": "100"
    },
    {
        "name": "h_bond_acceptors",
        "description": "Number of hydrogen bond acceptors in the molecule",
        "data_type": "integer",
        "units": "count",
        "min_value": "0",
        "max_value": "100"
    },
    {
        "name": "molecular_weight",
        "description": "Molecular weight of the compound",
        "data_type": "float",
        "units": "g/mol",
        "min_value": "1.0",
        "max_value": "10000.0"
    },
    {
        "name": "rotatable_bonds",
        "description": "Number of rotatable bonds in the molecule",
        "data_type": "integer",
        "units": "count",
        "min_value": "0",
        "max_value": "100"
    }
]

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

def check_existing_property_types(conn):
    """Check which property types already exist."""
    existing = set()
    
    with conn.cursor() as cursor:
        cursor.execute("""
            SELECT name FROM property_types
        """)
        for row in cursor.fetchall():
            existing.add(row[0])
    
    return existing

def create_property_type(conn, property_type):
    """Create a new property type."""
    with conn.cursor() as cursor:
        cursor.execute("""
            INSERT INTO property_types (
                id, name, description, data_type, units, min_value, max_value,
                created_at, updated_at
            ) VALUES (
                %s, %s, %s, %s, %s, %s, %s, NOW(), NOW()
            )
        """, (
            str(uuid.uuid4()),
            property_type["name"],
            property_type["description"],
            property_type["data_type"],
            property_type["units"],
            property_type["min_value"],
            property_type["max_value"]
        ))
        return cursor.rowcount

def main():
    """Create property types for critical molecular properties."""
    print("Creating property types for critical molecular properties...")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Check existing property types
        existing_types = check_existing_property_types(conn)
        print(f"Found {len(existing_types)} existing property types")
        
        # Create new property types
        created_count = 0
        skipped_count = 0
        
        for prop_type in PROPERTY_TYPES:
            if prop_type["name"] in existing_types:
                print(f"Skipping existing property type: {prop_type['name']}")
                skipped_count += 1
                continue
            
            created = create_property_type(conn, prop_type)
            if created:
                created_count += 1
                print(f"Created property type: {prop_type['name']}")
            else:
                print(f"Failed to create property type: {prop_type['name']}")
        
        # Commit changes
        conn.commit()
        print(f"\nSummary: Created {created_count} property types, skipped {skipped_count} existing types")
    
    except Exception as e:
        conn.rollback()
        print(f"Error creating property types: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()