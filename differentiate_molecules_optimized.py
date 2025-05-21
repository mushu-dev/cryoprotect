#!/usr/bin/env python3
"""
Optimized differentiate script for large molecule groups.

This script:
1. Processes the large 'DIFFERENTIATE' group in batches
2. Uses a more efficient approach for adding differentiation properties
3. Handles potential database timeouts and connection issues
"""

import os
import sys
import json
import uuid
import logging
import psycopg2
from psycopg2.extras import RealDictCursor, Json
from datetime import datetime
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Constants
PLAN_FILE = "duplicate_consolidation_plan.json"
DIFFERENTIATION_LOG_FILE = "differentiation_log.json"
BATCH_SIZE = 20  # Process this many molecules at once

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('differentiation_large_group.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

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

def load_consolidation_plan():
    """Load the consolidation plan."""
    try:
        with open(PLAN_FILE, 'r') as f:
            plan = json.load(f)
        return plan
    except Exception as e:
        logger.error(f"Error loading consolidation plan: {e}")
        sys.exit(1)

def process_large_differentiate_group(group_id):
    """
    Process the large differentiate group specifically.
    
    Args:
        group_id: ID of the large group to process
        
    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Processing large differentiate group {group_id}")
    
    # Load the plan and identify the large group
    plan = load_consolidation_plan()
    large_group = None
    
    for group in plan.get('differentiate', []):
        if group['group_id'] == group_id:
            large_group = group
            break
    
    if not large_group:
        logger.error(f"Group {group_id} not found in consolidation plan")
        return False
    
    logger.info(f"Found group with {large_group['molecule_count']} molecules")
    
    # Create differentiation log
    differentiation_log = {
        'group_id': group_id,
        'differentiation_type': 'DIFFERENTIATE',
        'molecules': [],
        'timestamp': datetime.now().isoformat(),
        'changes': []
    }
    
    # Process molecules in batches
    total_molecules = large_group['molecules_detail']
    batches = [total_molecules[i:i + BATCH_SIZE] for i in range(0, len(total_molecules), BATCH_SIZE)]
    
    logger.info(f"Processing molecules in {len(batches)} batches of up to {BATCH_SIZE} molecules each")
    
    total_processed = 0
    total_success = 0
    
    for batch_idx, batch in enumerate(batches):
        logger.info(f"Processing batch {batch_idx + 1} of {len(batches)} with {len(batch)} molecules")
        
        # Connect to database (fresh connection for each batch)
        try:
            conn = connect_to_db()
            conn.autocommit = False  # Use transaction for each batch
            
            with conn.cursor(cursor_factory=RealDictCursor) as cursor:
                for mol in batch:
                    mol_id = mol['id']
                    mol_name = mol['name'] or "Unnamed molecule"
                    total_processed += 1
                    
                    logger.info(f"Processing molecule {total_processed} of {large_group['molecule_count']}: {mol_name} (ID: {mol_id})")
                    
                    # Add molecule to the log
                    differentiation_log['molecules'].append({
                        'id': mol_id,
                        'name': mol_name,
                        'pubchem_cid': mol.get('pubchem_cid'),
                        'smiles': mol.get('smiles'),
                        'molecular_formula': mol.get('molecular_formula')
                    })
                    
                    # Create differentiation properties
                    diff_properties = {
                        'differentiation_group': group_id,
                        'has_structural_differences': True
                    }
                    
                    # Add structural information if available
                    if mol.get('smiles'):
                        diff_properties['smiles_structure'] = mol.get('smiles')
                        
                    if mol.get('molecular_formula'):
                        diff_properties['formula'] = mol.get('molecular_formula')
                    
                    # Update molecule properties
                    try:
                        cursor.execute("""
                            UPDATE molecules
                            SET properties = jsonb_set(
                                COALESCE(properties, '{}'::jsonb),
                                '{differentiation}',
                                %s::jsonb,
                                true
                            )
                            WHERE id = %s
                            RETURNING id, name
                        """, (Json(diff_properties), mol_id))
                        
                        updated_mol = cursor.fetchone()
                        if updated_mol:
                            change = {
                                'type': 'property_update',
                                'molecule_id': mol_id,
                                'field': 'properties.differentiation',
                                'value': diff_properties
                            }
                            differentiation_log['changes'].append(change)
                            total_success += 1
                            logger.info(f"Updated molecule {mol_id}")
                    except Exception as e:
                        logger.error(f"Error updating molecule {mol_id}: {e}")
                        # Continue with next molecule
                
                # Commit batch
                conn.commit()
                logger.info(f"Committed batch {batch_idx + 1}, processed {len(batch)} molecules")
        
        except Exception as e:
            conn.rollback()
            logger.error(f"Error processing batch {batch_idx + 1}: {e}")
        
        finally:
            conn.close()
    
    # Save the complete differentiation log
    append_to_differentiation_log(differentiation_log)
    
    logger.info(f"Large group processing complete. Processed {total_processed} molecules, successfully updated {total_success}.")
    return total_success > 0

def append_to_differentiation_log(log_entry):
    """Append a differentiation log entry to the log file."""
    try:
        # Create or load existing log
        if os.path.exists(DIFFERENTIATION_LOG_FILE):
            with open(DIFFERENTIATION_LOG_FILE, 'r') as f:
                log = json.load(f)
        else:
            log = []
        
        # Append new entry
        log.append(log_entry)
        
        # Write back
        with open(DIFFERENTIATION_LOG_FILE, 'w') as f:
            json.dump(log, f, indent=2, default=str)
            
    except Exception as e:
        logger.error(f"Error appending to differentiation log: {e}")

def main():
    """
    Process the large differentiate group.
    """
    import argparse
    parser = argparse.ArgumentParser(description="Optimized script to differentiate large molecule groups")
    parser.add_argument("--group", default="b278ce51-6575-46f2-ae4e-79369475cc8d", 
                      help="Group ID of the large molecule group (default: b278ce51-6575-46f2-ae4e-79369475cc8d)")
    args = parser.parse_args()
    
    logger.info(f"Starting optimized differentiation for group {args.group}")
    
    success = process_large_differentiate_group(args.group)
    
    if success:
        logger.info(f"Successfully processed large group {args.group}")
    else:
        logger.error(f"Failed to process large group {args.group}")
        sys.exit(1)

if __name__ == "__main__":
    main()