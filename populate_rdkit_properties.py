#\!/usr/bin/env python3
"""
Populate RDKit Properties

This script populates essential RDKit properties for all molecules.
"""

import os
import sys
import logging
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import json
import uuid
from datetime import datetime
import math

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"rdkit_properties_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Load environment variables
load_dotenv()

def connect_to_db():
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }
    
    try:
        conn = psycopg2.connect(**db_params)
        logger.info("Connected to database")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def get_property_type_mapping(conn):
    """Get a mapping of old property types to standardized ones."""
    property_mapping = {
        'LogP': 'logP',
        'Molecular Weight': 'molecularWeight',
        'Hydrogen Bond Donor Count': 'numHDonors',
        'Hydrogen Bond Acceptor Count': 'numHAcceptors',
        'Rotatable Bond Count': 'numRotatableBonds',
        'Aromatic Ring Count': 'aromaticRings',
        'Topological Polar Surface Area': 'TPSA'
    }
    
    property_id_mapping = {}
    
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Get IDs for both old and new property types
            for old_name, new_name in property_mapping.items():
                # Check old name
                cursor.execute("SELECT id FROM property_types WHERE name = %s", (old_name,))
                old_result = cursor.fetchone()
                
                # Check new name
                cursor.execute("SELECT id FROM property_types WHERE name = %s", (new_name,))
                new_result = cursor.fetchone()
                
                if old_result and new_result:
                    property_id_mapping[old_result['id']] = new_result['id']
                elif old_result:
                    logger.warning(f"Property type '{old_name}' exists but '{new_name}' doesn't")
                elif new_result:
                    logger.warning(f"Property type '{new_name}' exists but '{old_name}' doesn't")
                else:
                    logger.warning(f"Neither '{old_name}' nor '{new_name}' property types exist")
            
        return property_id_mapping
    except Exception as e:
        logger.error(f"Error getting property type mapping: {e}")
        return {}

def migrate_property_values(conn, property_id_mapping):
    """Migrate property values from old property types to new ones."""
    try:
        success_count = 0
        
        for old_id, new_id in property_id_mapping.items():
            with conn.cursor() as cursor:
                # Count properties with old ID
                cursor.execute("SELECT COUNT(*) FROM molecular_properties WHERE property_type_id = %s", (old_id,))
                count = cursor.fetchone()[0]
                
                logger.info(f"Found {count} properties with old ID {old_id}")
                
                if count > 0:
                    # Update properties to use new ID
                    cursor.execute("""
                        INSERT INTO molecular_properties 
                        (id, molecule_id, property_type_id, numeric_value, text_value, created_at, updated_at)
                        SELECT 
                            uuid_generate_v4(), molecule_id, %s, numeric_value, text_value, created_at, NOW()
                        FROM molecular_properties
                        WHERE property_type_id = %s
                        AND molecule_id NOT IN (
                            SELECT molecule_id FROM molecular_properties WHERE property_type_id = %s
                        )
                    """, (new_id, old_id, new_id))
                    
                    inserted = cursor.rowcount
                    logger.info(f"Inserted {inserted} properties with new ID {new_id}")
                    
                    if inserted > 0:
                        success_count += 1
        
        return success_count
    except Exception as e:
        logger.error(f"Error migrating property values: {e}")
        return 0

def create_missing_property_types(conn):
    """Create any missing standardized property types."""
    standard_properties = {
        'logP': ('The octanol-water partition coefficient', 'numeric'),
        'molecularWeight': ('The molecular weight of the compound', 'numeric'),
        'numHDonors': ('The number of hydrogen bond donors', 'numeric'),
        'numHAcceptors': ('The number of hydrogen bond acceptors', 'numeric'),
        'numRotatableBonds': ('The number of rotatable bonds', 'numeric'),
        'aromaticRings': ('The number of aromatic rings', 'numeric'),
        'TPSA': ('Topological Polar Surface Area', 'numeric')
    }
    
    created_count = 0
    
    try:
        with conn.cursor() as cursor:
            for name, (description, data_type) in standard_properties.items():
                # Check if property type exists
                cursor.execute("SELECT id FROM property_types WHERE name = %s", (name,))
                result = cursor.fetchone()
                
                if not result:
                    # Create property type
                    logger.info(f"Creating property type: {name}")
                    property_id = str(uuid.uuid4())
                    cursor.execute("""
                        INSERT INTO property_types (id, name, description, data_type, created_at)
                        VALUES (%s, %s, %s, %s, %s)
                    """, (
                        property_id,
                        name,
                        description,
                        data_type,
                        datetime.now()
                    ))
                    created_count += 1
        
        return created_count
    except Exception as e:
        logger.error(f"Error creating missing property types: {e}")
        return 0

def copy_logp_to_lowercase(conn):
    """Copy LogP property values to logP for all molecules."""
    try:
        with conn.cursor() as cursor:
            # Get property type IDs
            cursor.execute("SELECT id FROM property_types WHERE name = 'LogP'")
            old_id = cursor.fetchone()
            
            cursor.execute("SELECT id FROM property_types WHERE name = 'logP'")
            new_id = cursor.fetchone()
            
            if not old_id or not new_id:
                logger.error("Could not find LogP and/or logP property types")
                return 0
            
            old_id = old_id[0]
            new_id = new_id[0]
            
            # Count existing LogP values
            cursor.execute("SELECT COUNT(*) FROM molecular_properties WHERE property_type_id = %s", (old_id,))
            count = cursor.fetchone()[0]
            
            logger.info(f"Found {count} LogP values")
            
            # Copy LogP values to logP for molecules that don't have logP yet
            cursor.execute("""
                INSERT INTO molecular_properties 
                (id, molecule_id, property_type_id, numeric_value, text_value, created_at, updated_at)
                SELECT 
                    uuid_generate_v4(), molecule_id, %s, numeric_value, text_value, created_at, NOW()
                FROM molecular_properties
                WHERE property_type_id = %s
                AND molecule_id NOT IN (
                    SELECT molecule_id FROM molecular_properties WHERE property_type_id = %s
                )
            """, (new_id, old_id, new_id))
            
            copied = cursor.rowcount
            logger.info(f"Copied {copied} LogP values to logP")
            
            return copied
    except Exception as e:
        logger.error(f"Error copying LogP to logP: {e}")
        return 0

def copy_property_values(conn, source_name, target_name):
    """Copy property values from one property type to another."""
    try:
        with conn.cursor() as cursor:
            # Get property type IDs
            cursor.execute("SELECT id FROM property_types WHERE name = %s", (source_name,))
            source_id = cursor.fetchone()
            
            cursor.execute("SELECT id FROM property_types WHERE name = %s", (target_name,))
            target_id = cursor.fetchone()
            
            if not source_id or not target_id:
                logger.error(f"Could not find {source_name} and/or {target_name} property types")
                return 0
            
            source_id = source_id[0]
            target_id = target_id[0]
            
            # Count existing source values
            cursor.execute("SELECT COUNT(*) FROM molecular_properties WHERE property_type_id = %s", (source_id,))
            count = cursor.fetchone()[0]
            
            logger.info(f"Found {count} {source_name} values")
            
            # Copy values for molecules that don't have target yet
            cursor.execute("""
                INSERT INTO molecular_properties 
                (id, molecule_id, property_type_id, numeric_value, text_value, created_at, updated_at)
                SELECT 
                    uuid_generate_v4(), molecule_id, %s, numeric_value, text_value, created_at, NOW()
                FROM molecular_properties
                WHERE property_type_id = %s
                AND molecule_id NOT IN (
                    SELECT molecule_id FROM molecular_properties WHERE property_type_id = %s
                )
            """, (target_id, source_id, target_id))
            
            copied = cursor.rowcount
            logger.info(f"Copied {copied} {source_name} values to {target_name}")
            
            return copied
    except Exception as e:
        logger.error(f"Error copying {source_name} to {target_name}: {e}")
        return 0

def main():
    """Main function."""
    conn = connect_to_db()
    if not conn:
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Set autocommit to False to use transactions
        conn.autocommit = False
        
        # Create missing property types
        created = create_missing_property_types(conn)
        logger.info(f"Created {created} missing property types")
        
        # Special case: Copy LogP to logP
        copied_logp = copy_property_values(conn, 'LogP', 'logP')
        
        # Try to find and apply other property mappings
        property_mapping = [
            ('Molecular Weight', 'molecularWeight'),
            ('Hydrogen Bond Donor Count', 'numHDonors'),
            ('Hydrogen Bond Acceptor Count', 'numHAcceptors'),
            ('Rotatable Bond Count', 'numRotatableBonds'),
            ('Aromatic Ring Count', 'aromaticRings'),
            ('Topological Polar Surface Area', 'TPSA')
        ]
        
        for source, target in property_mapping:
            copy_property_values(conn, source, target)
        
        conn.commit()
        logger.info("RDKit properties populated successfully")
        return 0
    except Exception as e:
        logger.exception(f"Error populating RDKit properties: {e}")
        conn.rollback()
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())
