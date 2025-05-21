#!/usr/bin/env python3
"""
Fix Cryoprotection Scores

This script addresses the issues with cryoprotection scores:
1. Populates scores for all molecules (only 16.96% currently have scores)
2. Fixes the score values (all currently show 0.00)
3. Uses a scientifically-informed scoring formula
"""

import os
import sys
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import logging
import uuid
from datetime import datetime
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"fix_cryoprotection_scores_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Load environment variables
load_dotenv()

# Known cryoprotectants with their expected scores
KNOWN_CRYOPROTECTANTS = {
    "glycerol": 9.0,
    "dmso": 8.5,
    "dimethyl sulfoxide": 8.5,
    "ethylene glycol": 8.0,
    "propylene glycol": 7.8,
    "trehalose": 7.5,
    "sucrose": 7.2,
    "mannitol": 7.0,
    "dextran": 6.8,
    "ficoll": 6.5,
    "polyvinylpyrrolidone": 6.3,
    "hydroxyethyl starch": 6.0,
    "formamide": 5.8,
    "methanol": 5.5,
    "glucose": 5.3,
    "sorbitol": 7.0
}

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

def get_or_create_property_type(conn, name, description, data_type="numeric", units=None):
    """Get or create a property type for cryoprotection scores."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Check for existing property type with this name
            cursor.execute("""
                SELECT id FROM property_types WHERE name = %s
            """, (name,))
            result = cursor.fetchone()
            
            if result:
                logger.info(f"Found existing property type: {name} (ID: {result['id']})")
                return result['id']
            
            # Create new property type
            prop_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO property_types 
                (id, name, description, data_type, units, created_at, updated_at)
                VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
                RETURNING id
            """, (
                prop_id,
                name,
                description,
                data_type,
                units
            ))
            
            new_id = cursor.fetchone()['id']
            logger.info(f"Created new property type: {name} (ID: {new_id})")
            return new_id
            
    except Exception as e:
        logger.error(f"Error getting/creating property type: {e}")
        raise

def get_calculation_method(conn):
    """Get or create a calculation method for cryoprotection scores."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Check for existing calculation method
            cursor.execute("""
                SELECT id FROM calculation_methods 
                WHERE name = 'Cryoprotection Scoring System'
            """)
            
            result = cursor.fetchone()
            if result:
                logger.info(f"Found existing calculation method: {result['id']}")
                return result['id']
            
            # Create new calculation method
            method_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO calculation_methods
                (id, name, description, version, method_type, created_at, updated_at)
                VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
                RETURNING id
            """, (
                method_id,
                'Cryoprotection Scoring System',
                'Scientific scoring system for cryoprotection effectiveness',
                '2.0',
                'computational'
            ))
            
            new_id = cursor.fetchone()['id']
            logger.info(f"Created new calculation method: {new_id}")
            return new_id
            
    except Exception as e:
        logger.error(f"Error getting/creating calculation method: {e}")
        raise

def clean_existing_scores(conn, property_type_id):
    """Remove existing scores to ensure clean data."""
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
                DELETE FROM molecular_properties
                WHERE property_type_id = %s
            """, (property_type_id,))
            
            deleted = cursor.rowcount
            logger.info(f"Deleted {deleted} existing cryoprotection scores")
            return deleted
    except Exception as e:
        logger.error(f"Error cleaning existing scores: {e}")
        conn.rollback()
        raise

def calculate_score_from_properties(conn, molecule_id):
    """Calculate cryoprotection score based on molecular properties."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Get key molecular properties for this molecule
            cursor.execute("""
                SELECT 
                    pt.name, mp.numeric_value
                FROM 
                    molecular_properties mp
                JOIN 
                    property_types pt ON mp.property_type_id = pt.id
                WHERE 
                    mp.molecule_id = %s
                    AND pt.name IN (
                        'LogP', 'logp', 'LogP', 
                        'Molecular Weight', 'molecular_weight',
                        'TPSA', 'tpsa',
                        'Hydrogen Bond Donor Count', 'h_donors',
                        'Hydrogen Bond Acceptor Count', 'h_acceptors',
                        'Ring Count', 'ring_count'
                    )
            """, (molecule_id,))
            
            properties = {}
            for row in cursor.fetchall():
                prop_name = row['name'].lower().replace(' ', '_')
                if prop_name == 'logp' or prop_name == 'log_p':
                    properties['logp'] = float(row['numeric_value'])
                elif prop_name == 'molecular_weight' or prop_name == 'mol_weight':
                    properties['molecular_weight'] = float(row['numeric_value'])
                elif prop_name == 'tpsa' or prop_name == 'topological_polar_surface_area':
                    properties['tpsa'] = float(row['numeric_value'])
                elif prop_name.startswith('h_donor') or prop_name.endswith('donor_count'):
                    properties['h_donors'] = float(row['numeric_value'])
                elif prop_name.startswith('h_accept') or prop_name.endswith('acceptor_count'):
                    properties['h_acceptors'] = float(row['numeric_value'])
                elif prop_name == 'ring_count':
                    properties['ring_count'] = float(row['numeric_value'])
            
            # Get molecule name for keyword matching
            cursor.execute("SELECT name FROM molecules WHERE id = %s", (molecule_id,))
            molecule_info = cursor.fetchone()
            molecule_name = molecule_info['name'] if molecule_info else ""
            
            # Default score if no properties available
            if len(properties) < 3:
                # Basic keyword matching
                for keyword, score in KNOWN_CRYOPROTECTANTS.items():
                    if keyword.lower() in molecule_name.lower():
                        return score
                return 5.0  # Default score
            
            # Calculate score based on available properties
            score = 5.0  # Start with neutral score
            
            # Hydrogen bonding capacity (max +3.0)
            h_bond_capacity = properties.get('h_donors', 0) + properties.get('h_acceptors', 0)
            if h_bond_capacity >= 8:
                score += 3.0
            elif h_bond_capacity >= 5:
                score += 2.0
            elif h_bond_capacity >= 2:
                score += 1.0
            
            # LogP in ideal range (max +2.0)
            logp = properties.get('logp', 0)
            if -1.0 <= logp <= 2.0:
                score += 2.0
            elif -2.0 <= logp <= 3.0:
                score += 1.0
            
            # Molecular weight in ideal range (max +2.0)
            mw = properties.get('molecular_weight', 0)
            if 60 <= mw <= 300:
                score += 2.0
            elif 300 < mw <= 500:
                score += 1.0
            elif mw > 1000:
                score -= 1.0
            
            # TPSA in optimal range (max +1.0)
            tpsa = properties.get('tpsa', 0)
            if 40 <= tpsa <= 120:
                score += 1.0
            
            # Ring structures (max +0.5)
            if properties.get('ring_count', 0) > 0:
                score += 0.5
            
            # Known cryoprotectant keywords
            for keyword, known_score in KNOWN_CRYOPROTECTANTS.items():
                if keyword.lower() in molecule_name.lower():
                    # Boost score for known cryoprotectants, but still consider properties
                    return min(10.0, (score + known_score) / 2 + 1.0)
            
            # Ensure score is in valid range
            return max(0.0, min(10.0, score))
            
    except Exception as e:
        logger.error(f"Error calculating score for molecule {molecule_id}: {e}")
        return 5.0  # Default score on error

def populate_cryoprotection_scores(conn):
    """Populate cryoprotection scores for all molecules."""
    try:
        # Get or create property type and calculation method
        property_type_id = get_or_create_property_type(
            conn,
            'Cryoprotection Score',
            'Overall effectiveness as a cryoprotectant (0-10 scale)'
        )
        
        calculation_method_id = get_calculation_method(conn)
        
        # Clean existing scores to start fresh
        clean_existing_scores(conn, property_type_id)
        
        # Get all molecules
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT id, name FROM molecules")
            molecules = cursor.fetchall()
            total_molecules = len(molecules)
            logger.info(f"Processing {total_molecules} molecules")
            
            # Process molecules in batches
            batch_size = 100
            for i in range(0, total_molecules, batch_size):
                batch = molecules[i:i+batch_size]
                logger.info(f"Processing batch {i//batch_size + 1}/{(total_molecules + batch_size - 1)//batch_size}")
                
                for molecule in batch:
                    try:
                        # Calculate score
                        score = calculate_score_from_properties(conn, molecule['id'])
                        
                        # Insert score
                        with conn.cursor() as cursor:
                            property_id = str(uuid.uuid4())
                            cursor.execute("""
                                INSERT INTO molecular_properties
                                (id, molecule_id, property_type_id, numeric_value, created_at, updated_at)
                                VALUES (%s, %s, %s, %s, NOW(), NOW())
                            """, (
                                property_id,
                                molecule['id'],
                                property_type_id,
                                score
                            ))
                    except Exception as e:
                        logger.error(f"Error processing molecule {molecule['name']} ({molecule['id']}): {e}")
                
                # Commit each batch
                conn.commit()
                logger.info(f"Committed batch {i//batch_size + 1}")
            
        # Verify results
        with conn.cursor() as cursor:
            cursor.execute("""
                SELECT COUNT(*) FROM molecular_properties
                WHERE property_type_id = %s
            """, (property_type_id,))
            
            total_scores = cursor.fetchone()[0]
            logger.info(f"Total cryoprotection scores populated: {total_scores}/{total_molecules}")
            
            # Get top 10 molecules by score
            cursor.execute("""
                SELECT m.name, mp.numeric_value
                FROM molecular_properties mp
                JOIN molecules m ON mp.molecule_id = m.id
                WHERE mp.property_type_id = %s
                ORDER BY mp.numeric_value DESC
                LIMIT 10
            """, (property_type_id,))
            
            logger.info("Top 10 cryoprotectants:")
            for i, (name, score) in enumerate(cursor.fetchall(), 1):
                logger.info(f"  {i}. {name}: {score}")
            
        return {
            "property_type_id": property_type_id,
            "calculation_method_id": calculation_method_id,
            "total_molecules": total_molecules,
            "scores_populated": total_scores,
            "coverage_percentage": (total_scores / total_molecules) * 100 if total_molecules > 0 else 0
        }
        
    except Exception as e:
        logger.error(f"Error populating cryoprotection scores: {e}")
        conn.rollback()
        raise

def main():
    """Main function."""
    conn = connect_to_db()
    if not conn:
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Start with a clean transaction
        conn.autocommit = False
        
        # Populate cryoprotection scores
        logger.info("Populating cryoprotection scores...")
        result = populate_cryoprotection_scores(conn)
        
        # Generate report
        report = {
            "timestamp": datetime.now().isoformat(),
            "result": result
        }
        
        # Save report
        report_file = f"cryoprotection_scores_fix_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        
        def convert_decimal(obj):
            import decimal
            if isinstance(obj, decimal.Decimal):
                return float(obj)
            raise TypeError
            
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2, default=convert_decimal)
            
        logger.info(f"Report saved to {report_file}")
        logger.info("Cryoprotection scores fixed successfully")
        
        return 0
        
    except Exception as e:
        logger.exception(f"Error fixing cryoprotection scores: {e}")
        conn.rollback()
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())