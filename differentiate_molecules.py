#!/usr/bin/env python3
"""
Differentiate molecules with similar names but different structures.

This script:
1. Handles the 'DIFFERENTIATE' groups from the duplicate analysis
2. Adds properties to clarify the differences between similar molecules
3. Logs all changes for traceability
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

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('differentiation.log'),
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

def process_differentiate_group(conn, group, dry_run=False):
    """
    Process a differentiate group.
    These are molecules that appeared similar but have different structures.
    We'll add descriptive properties to clearly differentiate them.
    
    Args:
        conn: Database connection
        group: Group information from the consolidation plan
        dry_run: If True, don't make any changes
        
    Returns:
        True if successful, False otherwise
    """
    group_id = group['group_id']
    logger.info(f"Processing differentiate group {group_id}")
    
    # Create differentiation log
    differentiation_log = {
        'group_id': group_id,
        'differentiation_type': 'DIFFERENTIATE',
        'molecules': [],
        'timestamp': datetime.now().isoformat(),
        'changes': []
    }
    
    if dry_run:
        logger.info(f"DRY RUN: Would process {len(group['molecules_detail'])} molecules")
        for mol in group['molecules_detail']:
            logger.info(f"  Would differentiate: {mol['name']} (ID: {mol['id']})")
        return True
    
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Process each molecule in the group
            for mol in group['molecules_detail']:
                mol_id = mol['id']
                mol_name = mol['name']
                logger.info(f"Processing molecule: {mol_name} (ID: {mol_id})")
                
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
                
                # Find what makes this molecule unique in the group
                unique_properties = find_unique_properties(mol, group['molecules_detail'])
                if unique_properties:
                    diff_properties['differentiating_properties'] = unique_properties
                
                # Update molecule properties
                cursor.execute("""
                    UPDATE molecules
                    SET properties = jsonb_set(
                        COALESCE(properties, '{}'::jsonb),
                        '{differentiation}',
                        %s::jsonb,
                        true
                    )
                    WHERE id = %s
                    RETURNING id, name, properties
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
                    logger.info(f"Updated molecule {mol_id} with differentiation properties")
            
            # Commit all changes
            conn.commit()
            logger.info(f"Successfully differentiated group {group_id}")
            
            # Save the differentiation log
            append_to_differentiation_log(differentiation_log)
            return True
            
    except Exception as e:
        conn.rollback()
        logger.error(f"Error differentiating group {group_id}: {e}")
        return False

def find_unique_properties(molecule, all_molecules):
    """Find properties that make this molecule unique in the group."""
    unique_properties = {}
    
    # Skip if no relationships data
    if 'relationships' not in molecule:
        return unique_properties
    
    # Check properties
    if 'properties' in molecule['relationships']:
        props = molecule['relationships']['properties']
        
        # Find unique property values
        for prop in props:
            prop_type = prop.get('property_type')
            prop_value = prop.get('property_value')
            
            if not prop_type or not prop_value:
                continue
                
            # Check if this property with this value is unique
            is_unique = True
            for other_mol in all_molecules:
                if other_mol['id'] == molecule['id']:
                    continue
                    
                if 'relationships' not in other_mol:
                    continue
                    
                if 'properties' not in other_mol['relationships']:
                    continue
                
                for other_prop in other_mol['relationships']['properties']:
                    if (other_prop.get('property_type') == prop_type and 
                        other_prop.get('property_value') == prop_value):
                        is_unique = False
                        break
                
                if not is_unique:
                    break
            
            if is_unique:
                unique_properties[prop_type] = prop_value
    
    return unique_properties

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
    Differentiate molecules with similar names but different structures.
    """
    import argparse
    parser = argparse.ArgumentParser(description="Differentiate similar molecules with different structures")
    parser.add_argument("--group", help="Process only a specific group ID")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without making changes")
    args = parser.parse_args()
    
    # Load consolidation plan
    plan = load_consolidation_plan()
    differentiate_groups = plan.get('differentiate', [])
    
    # Connect to database
    logger.info("Connecting to database...")
    conn = connect_to_db()
    
    try:
        total_processed = 0
        total_success = 0
        
        # Filter to a specific group if requested
        if args.group:
            differentiate_groups = [g for g in differentiate_groups if g['group_id'] == args.group]
            if not differentiate_groups:
                logger.warning(f"Group {args.group} not found or not a differentiate group")
                sys.exit(1)
        
        logger.info(f"Processing {len(differentiate_groups)} differentiate groups")
        
        for group in differentiate_groups:
            group_id = group['group_id']
            logger.info(f"Differentiating group {group_id} with {group['molecule_count']} molecules")
            total_processed += 1
            
            success = process_differentiate_group(conn, group, dry_run=args.dry_run)
            
            if success:
                total_success += 1
        
        # Summary
        logger.info(f"Differentiation {'dry run' if args.dry_run else 'complete'}. " +
                    f"Processed {total_processed} groups, " + 
                    f"{'would have' if args.dry_run else ''} successfully differentiated {total_success} groups.")
        
    except Exception as e:
        logger.error(f"Error during differentiation: {e}")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()