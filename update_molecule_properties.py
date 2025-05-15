#!/usr/bin/env python3
"""
Update Molecule Properties with RDKit Microservice

This script updates the molecular properties of existing molecules in the database
using the new RDKit microservice deployed on Fly.io.
"""

import os
import sys
import logging
import psycopg2
from psycopg2.extras import RealDictCursor
import requests
import json
from datetime import datetime
import uuid
from urllib.parse import urlparse
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"rdkit_update_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Load environment variables
load_dotenv()

# RDKit service configuration
RDKIT_SERVICE_URL = os.environ.get('RDKIT_SERVICE_URL', 'https://cryoprotect-rdkit.fly.dev')
RDKIT_API_KEY = os.environ.get('RDKIT_API_KEY', '')

def get_db_connection():
    """Connect to the database using environment variables."""
    # First try Supabase connection params
    if os.environ.get('SUPABASE_DB_HOST'):
        db_params = {
            'host': os.environ.get('SUPABASE_DB_HOST'),
            'port': os.environ.get('SUPABASE_DB_PORT', '5432'),
            'dbname': os.environ.get('SUPABASE_DB_NAME', 'postgres'),
            'user': os.environ.get('SUPABASE_DB_USER'),
            'password': os.environ.get('SUPABASE_DB_PASSWORD')
        }
    # Otherwise try DATABASE_URL (Heroku style)
    elif os.environ.get('DATABASE_URL'):
        db_url = os.environ.get('DATABASE_URL')
        url = urlparse(db_url)
        db_params = {
            'host': url.hostname,
            'port': url.port or 5432,
            'dbname': url.path[1:],
            'user': url.username,
            'password': url.password
        }
    else:
        logger.error("No database connection parameters found")
        return None
    
    try:
        conn = psycopg2.connect(**db_params)
        logger.info("Connected to database")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def get_property_types(conn):
    """Get mapping of property type names to IDs."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT id, name FROM property_types")
            property_types = {row['name']: row['id'] for row in cursor.fetchall()}
            return property_types
    except Exception as e:
        logger.error(f"Error getting property types: {e}")
        return {}

def create_property_type(conn, name, description, data_type='numeric', unit=''):
    """Create a new property type."""
    try:
        with conn.cursor() as cursor:
            property_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO property_types (id, name, description, data_type, unit, created_at)
                VALUES (%s, %s, %s, %s, %s, %s)
                RETURNING id
            """, (
                property_id,
                name,
                description,
                data_type,
                unit,
                datetime.now()
            ))
            new_id = cursor.fetchone()[0]
            logger.info(f"Created new property type: {name} (ID: {new_id})")
            return new_id
    except Exception as e:
        logger.error(f"Error creating property type {name}: {e}")
        return None

def ensure_property_types(conn):
    """Ensure all needed property types exist in the database."""
    property_types = get_property_types(conn)
    
    # Define the RDKit properties we need
    rdkit_properties = {
        'logP': 'The octanol-water partition coefficient',
        'TPSA': 'Topological Polar Surface Area',
        'numHDonors': 'Number of hydrogen bond donors',
        'numHAcceptors': 'Number of hydrogen bond acceptors',
        'molecularWeight': 'Molecular weight calculated from RDKit',
        'numRotatableBonds': 'Number of rotatable bonds',
        'aromaticRings': 'Number of aromatic rings',
        'heavyAtomCount': 'Number of heavy (non-hydrogen) atoms',
    }
    
    # Create any missing property types
    for name, description in rdkit_properties.items():
        if name not in property_types:
            new_id = create_property_type(conn, name, description)
            if new_id:
                property_types[name] = new_id
    
    return property_types

def get_molecules(conn):
    """Get all molecules from the database."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT id, name, smiles FROM molecules")
            return cursor.fetchall()
    except Exception as e:
        logger.error(f"Error getting molecules: {e}")
        return []

def calculate_properties(smiles):
    """Calculate molecular properties using the RDKit microservice."""
    try:
        headers = {'Content-Type': 'application/json'}
        if RDKIT_API_KEY:
            headers['X-API-Key'] = RDKIT_API_KEY
            
        response = requests.post(
            f"{RDKIT_SERVICE_URL.rstrip('/')}/api/calculate-properties",
            headers=headers,
            json={'molecule_data': smiles},
            timeout=60
        )
        
        if response.status_code != 200:
            logger.error(f"Error from RDKit service: {response.status_code} {response.text}")
            return None
            
        data = response.json()
        if data.get('status') != 'success':
            logger.error(f"RDKit service error: {data}")
            return None
            
        return data.get('data', {})
    except Exception as e:
        logger.error(f"Error calling RDKit service: {e}")
        return None

def map_rdkit_to_properties(rdkit_data):
    """Map RDKit API response to property values."""
    if not rdkit_data:
        return {}
        
    properties = {
        'logP': rdkit_data.get('logp'),
        'TPSA': rdkit_data.get('tpsa'),
        'heavyAtomCount': rdkit_data.get('molecular_properties', {}).get('heavy_atom_count'),
        'molecularWeight': rdkit_data.get('molecular_properties', {}).get('molecular_weight'),
        'numRotatableBonds': rdkit_data.get('molecular_properties', {}).get('rotatable_bond_count'),
        'aromaticRings': rdkit_data.get('molecular_properties', {}).get('aromatic_ring_count'),
    }
    
    # Get hydrogen bonding data
    h_bonding = rdkit_data.get('hydrogen_bonding', {})
    properties['numHDonors'] = h_bonding.get('donors')
    properties['numHAcceptors'] = h_bonding.get('acceptors')
    
    # Filter out None values
    return {k: v for k, v in properties.items() if v is not None}

def update_molecular_properties(conn, molecule_id, properties, property_type_ids):
    """Update molecular properties for a molecule."""
    updated_count = 0
    
    try:
        with conn.cursor() as cursor:
            for prop_name, value in properties.items():
                if prop_name not in property_type_ids:
                    logger.warning(f"Property type '{prop_name}' not found in database")
                    continue
                    
                property_type_id = property_type_ids[prop_name]
                
                # Check if property already exists
                cursor.execute("""
                    SELECT id FROM molecular_properties 
                    WHERE molecule_id = %s AND property_type_id = %s
                """, (molecule_id, property_type_id))
                
                existing = cursor.fetchone()
                
                if existing:
                    # Update existing property
                    cursor.execute("""
                        UPDATE molecular_properties
                        SET numeric_value = %s, updated_at = %s
                        WHERE id = %s
                    """, (value, datetime.now(), existing[0]))
                else:
                    # Create new property
                    cursor.execute("""
                        INSERT INTO molecular_properties
                        (id, molecule_id, property_type_id, numeric_value, created_at, updated_at)
                        VALUES (%s, %s, %s, %s, %s, %s)
                    """, (
                        str(uuid.uuid4()),
                        molecule_id,
                        property_type_id,
                        value,
                        datetime.now(),
                        datetime.now()
                    ))
                
                updated_count += 1
    except Exception as e:
        logger.error(f"Error updating properties for molecule {molecule_id}: {e}")
        raise
        
    return updated_count

def main():
    """Main function."""
    conn = get_db_connection()
    if not conn:
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Set autocommit to False to use transactions
        conn.autocommit = False
        
        # Ensure property types exist
        property_type_ids = ensure_property_types(conn)
        if not property_type_ids:
            logger.error("Failed to ensure property types")
            return 1
            
        # Get all molecules
        molecules = get_molecules(conn)
        logger.info(f"Found {len(molecules)} molecules in database")
        
        total_updated = 0
        
        # Process each molecule
        for i, molecule in enumerate(molecules):
            molecule_id = molecule['id']
            smiles = molecule['smiles']
            name = molecule['name']
            
            logger.info(f"Processing molecule {i+1}/{len(molecules)}: {name}")
            
            if not smiles:
                logger.warning(f"Molecule {name} has no SMILES string, skipping")
                continue
                
            # Calculate properties using RDKit service
            rdkit_data = calculate_properties(smiles)
            if not rdkit_data:
                logger.warning(f"Failed to calculate properties for {name}, skipping")
                continue
                
            # Map RDKit data to property values
            properties = map_rdkit_to_properties(rdkit_data)
            
            # Update properties in database
            updated = update_molecular_properties(conn, molecule_id, properties, property_type_ids)
            logger.info(f"Updated {updated} properties for {name}")
            
            total_updated += updated
        
        conn.commit()
        logger.info(f"Successfully updated {total_updated} properties for {len(molecules)} molecules")
        return 0
    except Exception as e:
        logger.exception(f"Error updating molecular properties: {e}")
        conn.rollback()
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())