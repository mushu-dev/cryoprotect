#!/usr/bin/env python3
"""
Consolidate duplicate molecules based on the analysis and plan.

This script:
1. Reads the duplicate groups analysis and consolidation plan
2. Implements different consolidation strategies based on the group type
3. Updates molecule records in the database
4. Preserves relationships to ensure data integrity
5. Logs all changes for traceability
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
ANALYSIS_FILE = "duplicate_groups_analysis.json"
PLAN_FILE = "duplicate_consolidation_plan.json"
CONSOLIDATION_LOG_FILE = "duplicate_consolidation_log.json"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('consolidation.log'),
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

def load_analysis_and_plan():
    """Load the duplicate groups analysis and consolidation plan."""
    try:
        with open(ANALYSIS_FILE, 'r') as f:
            analysis = json.load(f)
        
        with open(PLAN_FILE, 'r') as f:
            plan = json.load(f)
        
        return analysis, plan
    except Exception as e:
        logger.error(f"Error loading analysis and plan: {e}")
        sys.exit(1)

def get_property_type_id(cursor, property_name):
    """Get property type ID by name."""
    if not property_name:
        return None
        
    cursor.execute("""
        SELECT id FROM property_types WHERE name = %s
    """, (property_name,))
    result = cursor.fetchone()
    
    if result:
        return result['id']
    return None

def consolidate_selective_merge(conn, group, analysis):
    """
    Implement selective merge consolidation strategy.
    
    This strategy is used when molecules have multiple PubChem CIDs.
    We select one molecule as primary and update others to reference it.
    """
    logger.info(f"Implementing SELECTIVE_MERGE for group {group['group_id']}")
    
    primary_id = group['primary_candidate']
    if not primary_id:
        logger.warning(f"No primary candidate for group {group['group_id']}. Skipping.")
        return False
    
    # Get the primary molecule's details
    primary_molecule = None
    for mol in group['molecules_detail']:
        if mol['id'] == primary_id:
            primary_molecule = mol
            break
    
    if not primary_molecule:
        logger.warning(f"Could not find primary molecule {primary_id} for group {group['group_id']}. Skipping.")
        return False
    
    logger.info(f"Primary molecule selected: {primary_molecule['name']} (ID: {primary_id})")
    
    # Process each secondary molecule
    secondary_molecules = [mol for mol in group['molecules_detail'] if mol['id'] != primary_id]
    
    # Create consolidation log
    consolidation_log = {
        'group_id': group['group_id'],
        'consolidation_type': 'SELECTIVE_MERGE',
        'primary_molecule': {
            'id': primary_molecule['id'],
            'name': primary_molecule['name'],
            'pubchem_cid': primary_molecule.get('pubchem_cid')
        },
        'secondary_molecules': [],
        'timestamp': datetime.now().isoformat(),
        'changes': []
    }
    
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            for sec_mol in secondary_molecules:
                logger.info(f"Processing secondary molecule: {sec_mol['name']} (ID: {sec_mol['id']})")
                
                # Record this in the consolidation log
                consolidation_log['secondary_molecules'].append({
                    'id': sec_mol['id'],
                    'name': sec_mol['name'],
                    'pubchem_cid': sec_mol.get('pubchem_cid')
                })
                
                # 1. Update secondary molecule to indicate it's a duplicate
                cursor.execute("""
                    UPDATE molecules
                    SET properties = jsonb_set(
                        COALESCE(properties, '{}'::jsonb),
                        '{consolidated_to}',
                        %s::jsonb,
                        true
                    )
                    WHERE id = %s
                    RETURNING id, name, properties
                """, (Json(primary_id), sec_mol['id']))
                
                updated_mol = cursor.fetchone()
                if updated_mol:
                    change = {
                        'type': 'property_update',
                        'molecule_id': sec_mol['id'],
                        'field': 'properties.consolidated_to',
                        'value': primary_id
                    }
                    consolidation_log['changes'].append(change)
                    logger.info(f"Updated secondary molecule {sec_mol['id']} properties")
                
                # 2. Find properties of the secondary molecule
                cursor.execute("""
                    SELECT *
                    FROM molecular_properties
                    WHERE molecule_id = %s
                """, (sec_mol['id'],))
                
                sec_properties = cursor.fetchall()
                
                # Process each property
                for prop in sec_properties:
                    property_type = prop.get('property_type')
                    property_type_id = prop.get('property_type_id')
                    
                    # Determine if primary already has this property
                    exists_query = "SELECT 1 FROM molecular_properties WHERE molecule_id = %s AND "
                    
                    if property_type_id:
                        exists_query += "property_type_id = %s"
                        exists_params = (primary_id, property_type_id)
                    elif property_type:
                        exists_query += "property_type = %s"
                        exists_params = (primary_id, property_type)
                    else:
                        logger.warning(f"Property has neither type_id nor type name, skipping")
                        continue
                    
                    cursor.execute(exists_query, exists_params)
                    
                    # If property doesn't exist on primary, copy it
                    if not cursor.fetchone():
                        # Create new property for primary
                        new_prop_id = str(uuid.uuid4())
                        
                        # Ensure we have a property_type_id
                        if not property_type_id and property_type:
                            property_type_id = get_property_type_id(cursor, property_type)
                            if not property_type_id:
                                logger.warning(f"Could not find property type ID for {property_type}, skipping")
                                continue
                        
                        # Build the insert query
                        insert_query = """
                            INSERT INTO molecular_properties
                            (id, molecule_id, data_source, created_at, updated_at
                        """
                        
                        # Add optional columns if they have values
                        params = [new_prop_id, primary_id, f"Consolidated from {sec_mol['id']}", 
                                  datetime.now(), datetime.now()]
                        
                        # Add property_type_id if available
                        if property_type_id:
                            insert_query += ", property_type_id"
                            params.append(property_type_id)
                        
                        # Add other fields if available
                        fields = [
                            ('property_type', prop.get('property_type')),
                            ('property_value', prop.get('property_value')),
                            ('numeric_value', prop.get('numeric_value')),
                            ('text_value', prop.get('text_value')),
                            ('boolean_value', prop.get('boolean_value')),
                            ('unit', prop.get('unit')),
                            ('source', prop.get('source'))
                        ]
                        
                        for field_name, field_value in fields:
                            if field_value is not None:
                                insert_query += f", {field_name}"
                                params.append(field_value)
                        
                        # Complete the query
                        placeholders = ", ".join(["%s"] * len(params))
                        insert_query += ") VALUES (" + placeholders + ") RETURNING id"
                        
                        # Execute the insert
                        cursor.execute(insert_query, params)
                        
                        new_prop = cursor.fetchone()
                        if new_prop:
                            change = {
                                'type': 'property_migration',
                                'from_molecule_id': sec_mol['id'],
                                'to_molecule_id': primary_id,
                                'property_type': property_type,
                                'new_property_id': new_prop['id']
                            }
                            consolidation_log['changes'].append(change)
                            logger.info(f"Migrated property {property_type} from {sec_mol['id']} to {primary_id}")
            
            # Commit all changes
            conn.commit()
            logger.info(f"Successfully consolidated group {group['group_id']}")
            
            # Save consolidation log
            append_to_consolidation_log(consolidation_log)
            return True
    
    except Exception as e:
        conn.rollback()
        logger.error(f"Error consolidating group {group['group_id']}: {e}")
        return False

def consolidate_safe_merge(conn, group, analysis):
    """
    Implement safe merge consolidation strategy.
    
    This strategy is used when none of the duplicate molecules have relationships.
    We can safely merge all molecules into the primary one.
    """
    logger.info(f"Implementing SAFE_MERGE for group {group['group_id']}")
    
    # Implementation similar to selective_merge but can be more aggressive
    # since we don't need to worry about preserving relationships
    # Not implemented in this version since there are no safe merge groups
    return False

def consolidate_primary_selection(conn, group, analysis):
    """
    Implement primary selection consolidation strategy.
    
    This strategy is used when only one molecule has relationships.
    We designate it as primary and update others to reference it.
    """
    logger.info(f"Implementing PRIMARY_SELECTION for group {group['group_id']}")
    
    # Implementation similar to selective_merge
    # Not implemented in this version since there are no primary selection groups
    return False

def consolidate_complex_merge(conn, group, analysis):
    """
    Implement complex merge consolidation strategy.
    
    This strategy is used when multiple molecules have relationships.
    We need to carefully analyze and merge to preserve all relationships.
    """
    logger.info(f"Implementing COMPLEX_MERGE for group {group['group_id']}")
    
    # This is more complex and requires careful relationship migration
    # Not implemented in this version
    return False

def consolidate_differentiate(conn, group, analysis):
    """
    Implement differentiate consolidation strategy.
    
    This strategy is used when molecules have different chemical structures.
    We update the molecules to more clearly distinguish them.
    """
    logger.info(f"Implementing DIFFERENTIATE for group {group['group_id']}")
    
    # This involves updating names or adding distinguishing properties
    # Not implemented in this version
    return False

def append_to_consolidation_log(log_entry):
    """Append a consolidation log entry to the log file."""
    try:
        # Create or load existing log
        if os.path.exists(CONSOLIDATION_LOG_FILE):
            with open(CONSOLIDATION_LOG_FILE, 'r') as f:
                log = json.load(f)
        else:
            log = []
        
        # Append new entry
        log.append(log_entry)
        
        # Write back
        with open(CONSOLIDATION_LOG_FILE, 'w') as f:
            json.dump(log, f, indent=2, default=str)
            
    except Exception as e:
        logger.error(f"Error appending to consolidation log: {e}")

def main():
    """
    Implement the consolidation plan for duplicate molecules.
    
    By default, only processes selective merge groups in this version.
    """
    import argparse
    parser = argparse.ArgumentParser(description="Consolidate duplicate molecules")
    parser.add_argument("--strategy", choices=["selective", "safe", "primary", "complex", "differentiate", "all"],
                       default="selective", help="Consolidation strategy to implement")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without making changes")
    parser.add_argument("--group", help="Process only a specific group ID")
    args = parser.parse_args()
    
    # Load analysis and plan
    analysis, plan = load_analysis_and_plan()
    
    # Connect to database
    logger.info("Connecting to database...")
    conn = connect_to_db()
    
    try:
        total_processed = 0
        total_success = 0
        
        # Process groups based on selected strategy
        strategies = []
        
        if args.strategy == "all":
            strategies = ["selective", "safe", "primary", "complex", "differentiate"]
        else:
            strategies = [args.strategy]
        
        for strategy in strategies:
            groups = []
            
            if strategy == "selective":
                groups = plan["selective_merge"]
                logger.info(f"Processing {len(groups)} SELECTIVE_MERGE groups")
            elif strategy == "safe":
                groups = plan["safe_merge"]
                logger.info(f"Processing {len(groups)} SAFE_MERGE groups")
            elif strategy == "primary":
                groups = plan["primary_selection"]
                logger.info(f"Processing {len(groups)} PRIMARY_SELECTION groups")
            elif strategy == "complex":
                groups = plan["complex_merge"]
                logger.info(f"Processing {len(groups)} COMPLEX_MERGE groups")
            elif strategy == "differentiate":
                groups = plan["differentiate"]
                logger.info(f"Processing {len(groups)} DIFFERENTIATE groups")
            
            # If specific group ID is requested, filter to just that group
            if args.group:
                groups = [g for g in groups if g['group_id'] == args.group]
                if not groups:
                    logger.warning(f"Group {args.group} not found or not part of the {strategy} strategy")
                    continue
            
            for group in groups:
                group_id = group['group_id']
                group_analysis = analysis[group_id]
                
                logger.info(f"Processing group {group_id} with {group['molecule_count']} molecules")
                total_processed += 1
                
                # Skip if dry run
                if args.dry_run:
                    logger.info(f"DRY RUN: Would process group {group_id}")
                    continue
                
                # Apply the appropriate consolidation strategy
                success = False
                
                if strategy == "selective":
                    success = consolidate_selective_merge(conn, group, group_analysis)
                elif strategy == "safe":
                    success = consolidate_safe_merge(conn, group, group_analysis)
                elif strategy == "primary":
                    success = consolidate_primary_selection(conn, group, group_analysis)
                elif strategy == "complex":
                    success = consolidate_complex_merge(conn, group, group_analysis)
                elif strategy == "differentiate":
                    success = consolidate_differentiate(conn, group, group_analysis)
                
                if success:
                    total_success += 1
        
        # Summary
        logger.info(f"Consolidation complete. Processed {total_processed} groups, successfully consolidated {total_success} groups.")
        
    except Exception as e:
        logger.error(f"Error during consolidation: {e}")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()