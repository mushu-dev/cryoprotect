#\!/usr/bin/env python3
"""
Calculate Cryoprotection Scores

This script calculates cryoprotection scores for all molecules based on
their molecular properties. The score is a composite of multiple factors:
- HBD/HBA ratio
- LogP
- Molecular weight
- Number of rotatable bonds
- TPSA
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
        logging.FileHandler(f"cryoprotection_scores_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
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

def get_property_type_ids(conn):
    """Get IDs for all required property types."""
    try:
        property_names = [
            'numHDonors', 
            'numHAcceptors', 
            'logP', 
            'molecularWeight', 
            'numRotatableBonds', 
            'TPSA',
            'cryoprotectionScore'
        ]
        
        property_ids = {}
        
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            for name in property_names:
                cursor.execute("""
                    SELECT id FROM property_types WHERE name = %s
                """, (name,))
                
                result = cursor.fetchone()
                if result:
                    property_ids[name] = result['id']
                else:
                    # Create the property type if it doesn't exist
                    logger.info(f"Creating property type: {name}")
                    property_id = str(uuid.uuid4())
                    cursor.execute("""
                        INSERT INTO property_types (id, name, description, units, data_type, created_at)
                        VALUES (%s, %s, %s, %s, %s, %s)
                    """, (
                        property_id,
                        name,
                        f"The {name} property",
                        None,
                        'numeric',
                        datetime.now()
                    ))
                    property_ids[name] = property_id
        
        return property_ids
    except Exception as e:
        logger.error(f"Error getting property type ids: {e}")
        return {}

def get_molecules_without_scores(conn, property_ids):
    """Get all active molecules without cryoprotection scores."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("""
                SELECT m.id, m.name 
                FROM consolidated_molecules m
                WHERE m.molecule_status != 'Consolidated'
                AND NOT EXISTS (
                    SELECT 1 FROM molecular_properties mp 
                    WHERE mp.molecule_id = m.id 
                    AND mp.property_type_id = %s
                )
            """, (property_ids['cryoprotectionScore'],))
            
            return cursor.fetchall()
    except Exception as e:
        logger.error(f"Error getting molecules without scores: {e}")
        return []

def get_molecule_properties(conn, molecule_id, property_ids):
    """Get all relevant properties for a molecule."""
    try:
        properties = {}
        
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Get properties excluding cryoprotection score
            for name, type_id in property_ids.items():
                if name != 'cryoprotectionScore':
                    cursor.execute("""
                        SELECT mp.numeric_value
                        FROM molecular_properties mp
                        WHERE mp.molecule_id = %s AND mp.property_type_id = %s
                    """, (molecule_id, type_id))
                    
                    result = cursor.fetchone()
                    if result:
                        properties[name] = result['numeric_value']
                    else:
                        properties[name] = None
        
        return properties
    except Exception as e:
        logger.error(f"Error getting properties for molecule {molecule_id}: {e}")
        return {}

def calculate_score(properties):
    """Calculate cryoprotection score based on molecular properties."""
    # Check if we have enough properties to calculate a score
    required_properties = ['numHDonors', 'numHAcceptors', 'logP', 'molecularWeight', 'numRotatableBonds']
    for prop in required_properties:
        if prop not in properties or properties[prop] is None:
            logger.warning(f"Missing required property: {prop}")
            return None
    
    try:
        # Weight factors (totaling 1.0)
        weights = {
            'hb_ratio': 0.30,        # Hydrogen bond ratio
            'logp': 0.20,            # LogP (balanced hydrophobicity)
            'mol_weight': 0.20,      # Molecular weight
            'flexibility': 0.20,     # Rotatable bonds
            'polar_surface': 0.10    # TPSA
        }
        
        # Calculate component scores (0-100 scale)
        hb_donor = properties['numHDonors']
        hb_acceptor = properties['numHAcceptors']
        hb_donor_acceptor_ratio = hb_donor / max(hb_acceptor, 1)
        
        # Optimal H-bond donor/acceptor ratio is around 0.5-0.7
        # Score peaks at 0.6 ratio and falls off in both directions
        hb_ratio_score = 100 * math.exp(-5 * ((hb_donor_acceptor_ratio - 0.6) ** 2))
        
        # LogP score - optimal range around -1 to 1 (more water-like)
        # Score peaks at 0 (balanced) and falls off in both directions
        logp = properties['logP']
        logp_score = 100 * math.exp(-0.5 * (logp ** 2))
        
        # Molecular weight - prefer smaller molecules (100-300 Da range)
        mol_weight = properties['molecularWeight']
        if mol_weight <= 300:
            mol_weight_score = 100 * (1 - (mol_weight - 100) / 300)
        else:
            mol_weight_score = max(0, 100 * (1 - (mol_weight - 300) / 700))
        
        # Flexibility score - prefer moderately flexible molecules (2-5 rotatable bonds)
        rotatable_bonds = properties['numRotatableBonds']
        if rotatable_bonds <= 5:
            flexibility_score = 100 * (rotatable_bonds / 5)
        else:
            flexibility_score = max(0, 100 * (2 - rotatable_bonds / 10))
        
        # Polar surface area - prefer moderate values (40-90 Å²)
        tpsa = properties.get('TPSA', 0)
        if tpsa is None:
            tpsa = 0
            
        if tpsa <= 90:
            polar_surface_score = 100 * (tpsa / 90)
        else:
            polar_surface_score = max(0, 100 * (2 - tpsa / 180))
        
        # Calculate weighted total score
        total_score = (
            weights['hb_ratio'] * hb_ratio_score +
            weights['logp'] * logp_score +
            weights['mol_weight'] * mol_weight_score +
            weights['flexibility'] * flexibility_score +
            weights['polar_surface'] * polar_surface_score
        )
        
        # Scale to 0-1 range
        return total_score / 100
    
    except Exception as e:
        logger.error(f"Error calculating score: {e}")
        return None

def store_cryoprotection_score(conn, molecule_id, score, property_ids):
    """Store the calculated cryoprotection score in the database."""
    try:
        if score is None:
            return False
            
        with conn.cursor() as cursor:
            # Check if score already exists
            cursor.execute("""
                SELECT id FROM molecular_properties
                WHERE molecule_id = %s AND property_type_id = %s
            """, (molecule_id, property_ids['cryoprotectionScore']))
            
            existing = cursor.fetchone()
            
            if existing:
                # Update existing score
                cursor.execute("""
                    UPDATE molecular_properties
                    SET numeric_value = %s, updated_at = %s
                    WHERE id = %s
                """, (score, datetime.now(), existing[0]))
            else:
                # Insert new score
                cursor.execute("""
                    INSERT INTO molecular_properties
                    (id, molecule_id, property_type_id, numeric_value, created_at, updated_at)
                    VALUES (%s, %s, %s, %s, %s, %s)
                """, (
                    str(uuid.uuid4()),
                    molecule_id,
                    property_ids['cryoprotectionScore'],
                    score,
                    datetime.now(),
                    datetime.now()
                ))
                
            return True
    except Exception as e:
        logger.error(f"Error storing score for molecule {molecule_id}: {e}")
        return False

def calculate_all_scores(conn, property_ids, batch_size=100):
    """Calculate cryoprotection scores for all molecules missing them."""
    molecules = get_molecules_without_scores(conn, property_ids)
    total_molecules = len(molecules)
    
    logger.info(f"Found {total_molecules} molecules without cryoprotection scores")
    
    if not total_molecules:
        return True
    
    success_count = 0
    errors = []
    
    # Process in batches to avoid long transactions
    for i in range(0, total_molecules, batch_size):
        batch = molecules[i:i+batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{math.ceil(total_molecules/batch_size)} ({len(batch)} molecules)")
        
        # Start a new transaction for each batch
        conn.autocommit = False
        
        try:
            for molecule in batch:
                properties = get_molecule_properties(conn, molecule['id'], property_ids)
                score = calculate_score(properties)
                
                if score is not None:
                    if store_cryoprotection_score(conn, molecule['id'], score, property_ids):
                        success_count += 1
                    else:
                        errors.append({
                            'molecule_id': molecule['id'],
                            'name': molecule['name'],
                            'error': 'Failed to store score'
                        })
                else:
                    errors.append({
                        'molecule_id': molecule['id'],
                        'name': molecule['name'],
                        'error': 'Failed to calculate score'
                    })
            
            # Commit the batch
            conn.commit()
            logger.info(f"Committed batch with {len(batch)} molecules")
            
        except Exception as e:
            logger.error(f"Error processing batch: {e}")
            conn.rollback()
            
            # Mark all molecules in batch as errors
            for molecule in batch:
                if not any(e['molecule_id'] == molecule['id'] for e in errors):
                    errors.append({
                        'molecule_id': molecule['id'],
                        'name': molecule['name'],
                        'error': str(e)
                    })
    
    coverage_percentage = (success_count / total_molecules) * 100 if total_molecules > 0 else 0
    logger.info(f"Calculated scores for {success_count} of {total_molecules} molecules ({coverage_percentage:.2f}%)")
    
    if errors:
        error_file = f"cryoprotection_score_errors_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(error_file, 'w') as f:
            json.dump(errors, f, indent=2)
        logger.info(f"Saved {len(errors)} errors to {error_file}")
    
    return coverage_percentage >= 90

def main():
    """Main function."""
    conn = connect_to_db()
    if not conn:
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Get property type IDs
        property_ids = get_property_type_ids(conn)
        if not property_ids or len(property_ids) < 7:
            logger.error("Failed to get or create all required property types")
            return 1
        
        # Calculate scores for all molecules
        success = calculate_all_scores(conn, property_ids)
        
        if success:
            logger.info("Cryoprotection scores calculated successfully!")
            return 0
        else:
            logger.error("Failed to calculate scores for at least 90% of molecules")
            return 1
    except Exception as e:
        logger.exception(f"Error calculating cryoprotection scores: {e}")
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())
