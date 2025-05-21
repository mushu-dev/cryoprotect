#\!/usr/bin/env python3
"""
Merge Duplicate Properties

This script merges duplicate property types and copies values from h_bond_donors/acceptors
to the standardized numHDonors/numHAcceptors.
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

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"merge_properties_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
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

def find_duplicate_property_types(conn):
    """Find duplicate property types (case insensitive)."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("""
                SELECT LOWER(name) as lower_name, COUNT(*) as count
                FROM property_types
                GROUP BY LOWER(name)
                HAVING COUNT(*) > 1
            """)
            
            duplicate_names = cursor.fetchall()
            
            if not duplicate_names:
                logger.info("No duplicate property types found")
                return []
            
            logger.info(f"Found {len(duplicate_names)} duplicate property type names")
            
            duplicates = []
            
            for row in duplicate_names:
                lower_name = row['lower_name']
                
                cursor.execute("""
                    SELECT id, name, 
                           (SELECT COUNT(*) FROM molecular_properties WHERE property_type_id = property_types.id) as usage_count
                    FROM property_types
                    WHERE LOWER(name) = %s
                    ORDER BY usage_count DESC
                """, (lower_name,))
                
                property_types = cursor.fetchall()
                
                if len(property_types) > 1:
                    primary = property_types[0]
                    others = property_types[1:]
                    
                    duplicates.append({
                        'name': lower_name,
                        'primary': primary,
                        'duplicates': others
                    })
                    
                    logger.info(f"Property name '{lower_name}' has {len(property_types)} entries")
                    logger.info(f"  Primary: {primary['name']} (ID: {primary['id']}) - Used {primary['usage_count']} times")
                    
                    for pt in others:
                        logger.info(f"  Duplicate: {pt['name']} (ID: {pt['id']}) - Used {pt['usage_count']} times")
            
            return duplicates
    except Exception as e:
        logger.error(f"Error finding duplicate property types: {e}")
        return []

def merge_duplicate_property_types(conn, duplicates):
    """Merge duplicate property types."""
    if not duplicates:
        return True
    
    try:
        for duplicate in duplicates:
            primary = duplicate['primary']
            
            for dup in duplicate['duplicates']:
                logger.info(f"Merging '{dup['name']}' into '{primary['name']}'")
                
                with conn.cursor() as cursor:
                    # Move molecular_properties to primary
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
                    """, (primary['id'], dup['id'], primary['id']))
                    
                    moved = cursor.rowcount
                    logger.info(f"Moved {moved} property values to primary property type")
                    
                    # Move predictions to primary
                    cursor.execute("""
                        INSERT INTO predictions 
                        (id, molecule_id, property_type_id, numeric_value, confidence, created_at, updated_at)
                        SELECT 
                            uuid_generate_v4(), molecule_id, %s, numeric_value, confidence, created_at, NOW()
                        FROM predictions
                        WHERE property_type_id = %s
                        AND molecule_id NOT IN (
                            SELECT molecule_id FROM predictions WHERE property_type_id = %s
                        )
                    """, (primary['id'], dup['id'], primary['id']))
                    
                    moved_predictions = cursor.rowcount
                    logger.info(f"Moved {moved_predictions} predictions to primary property type")
                    
                    # Delete the duplicate property type's data
                    cursor.execute("DELETE FROM molecular_properties WHERE property_type_id = %s", (dup['id'],))
                    deleted_props = cursor.rowcount
                    
                    cursor.execute("DELETE FROM predictions WHERE property_type_id = %s", (dup['id'],))
                    deleted_preds = cursor.rowcount
                    
                    # Delete the duplicate property type
                    cursor.execute("DELETE FROM property_types WHERE id = %s", (dup['id'],))
                    
                    logger.info(f"Deleted duplicate property type and {deleted_props + deleted_preds} related records")
        
        return True
    except Exception as e:
        logger.error(f"Error merging duplicate property types: {e}")
        return False

def copy_hydrogen_bond_properties(conn):
    """Copy h_bond_donors/acceptors to numHDonors/numHAcceptors."""
    property_mappings = [
        ('h_bond_donors', 'numHDonors'),
        ('h_bond_acceptors', 'numHAcceptors'),
        ('rtb', 'numRotatableBonds'),
        ('hydrogenBondDonorCount', 'numHDonors'),
        ('hydrogenBondAcceptorCount', 'numHAcceptors'),
        ('rotatableBondCount', 'numRotatableBonds'),
        ('aromaticRingCount', 'aromaticRings')
    ]
    
    try:
        for source, target in property_mappings:
            with conn.cursor(cursor_factory=RealDictCursor) as cursor:
                # Get property type IDs
                cursor.execute("SELECT id FROM property_types WHERE name = %s", (source,))
                source_result = cursor.fetchone()
                
                cursor.execute("SELECT id FROM property_types WHERE name = %s", (target,))
                target_result = cursor.fetchone()
                
                if not source_result or not target_result:
                    logger.warning(f"Could not find either '{source}' or '{target}' property type")
                    continue
                
                source_id = source_result['id']
                target_id = target_result['id']
                
                # Copy values
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
                logger.info(f"Copied {copied} values from '{source}' to '{target}'")
        
        return True
    except Exception as e:
        logger.error(f"Error copying hydrogen bond properties: {e}")
        return False

def main():
    """Main function."""
    conn = connect_to_db()
    if not conn:
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Set autocommit to False to use transactions
        conn.autocommit = False
        
        # Find duplicate property types
        duplicates = find_duplicate_property_types(conn)
        
        # Merge duplicate property types
        if duplicates and merge_duplicate_property_types(conn, duplicates):
            logger.info("Successfully merged duplicate property types")
        
        # Copy hydrogen bond properties
        if copy_hydrogen_bond_properties(conn):
            logger.info("Successfully copied hydrogen bond properties")
        
        conn.commit()
        logger.info("Property merge completed successfully")
        return 0
    except Exception as e:
        logger.exception(f"Error merging properties: {e}")
        conn.rollback()
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())
