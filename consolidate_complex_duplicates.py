#!/usr/bin/env python3
"""
Consolidate complex duplicate molecule groups.

This script:
1. Handles the 'COMPLEX_MERGE' groups from the duplicate analysis
2. Preserves all relationships while consolidating duplicates
3. Migrates properties, mixture components, and predictions to the primary molecule
4. Logs all changes for traceability
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
CONSOLIDATION_LOG_FILE = "complex_consolidation_log.json"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('complex_consolidation.log'),
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

def consolidate_complex_group(conn, group, analysis, dry_run=False):
    """
    Consolidate a complex merge group.
    
    Args:
        conn: Database connection
        group: Group information from the consolidation plan
        analysis: The full analysis data for reference
        dry_run: If True, don't make any changes
        
    Returns:
        True if successful, False otherwise
    """
    group_id = group['group_id']
    logger.info(f"Processing complex merge group {group_id}")
    
    primary_id = group['primary_candidate']
    if not primary_id:
        logger.warning(f"No primary candidate for group {group_id}. Skipping.")
        return False
    
    # Find the primary molecule details
    primary_molecule = None
    for mol in group['molecules_detail']:
        if mol['id'] == primary_id:
            primary_molecule = mol
            break
    
    if not primary_molecule:
        logger.warning(f"Could not find primary molecule {primary_id} for group {group_id}. Skipping.")
        return False
    
    logger.info(f"Primary molecule selected: {primary_molecule['name']} (ID: {primary_id})")
    
    # Process each secondary molecule
    secondary_molecules = [mol for mol in group['molecules_detail'] if mol['id'] != primary_id]
    
    # Create consolidation log
    consolidation_log = {
        'group_id': group_id,
        'consolidation_type': 'COMPLEX_MERGE',
        'primary_molecule': {
            'id': primary_molecule['id'],
            'name': primary_molecule['name'],
            'pubchem_cid': primary_molecule.get('pubchem_cid')
        },
        'secondary_molecules': [],
        'timestamp': datetime.now().isoformat(),
        'changes': []
    }
    
    if dry_run:
        logger.info(f"DRY RUN: Would process {len(secondary_molecules)} secondary molecules")
        for sec_mol in secondary_molecules:
            logger.info(f"  Would consolidate: {sec_mol['name']} (ID: {sec_mol['id']})")
        return True
    
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
                
                # 2. Migrate properties from secondary to primary
                migrate_molecular_properties(cursor, sec_mol['id'], primary_id, consolidation_log)
                
                # 3. Migrate mixture components
                migrate_mixture_components(cursor, sec_mol['id'], primary_id, consolidation_log)
                
                # 4. Migrate predictions
                migrate_predictions(cursor, sec_mol['id'], primary_id, consolidation_log)
            
            # Commit all changes
            conn.commit()
            logger.info(f"Successfully consolidated complex group {group_id}")
            
            # Save the consolidation log
            append_to_consolidation_log(consolidation_log)
            return True
            
    except Exception as e:
        conn.rollback()
        logger.error(f"Error consolidating complex group {group_id}: {e}")
        return False

def migrate_molecular_properties(cursor, source_id, target_id, consolidation_log):
    """
    Migrate molecular properties from source to target molecule.
    
    Only migrates properties that don't already exist on the target.
    """
    # Get all properties from source molecule
    cursor.execute("""
        SELECT * FROM molecular_properties
        WHERE molecule_id = %s
    """, (source_id,))
    
    source_properties = cursor.fetchall()
    logger.info(f"Found {len(source_properties)} properties to migrate")
    
    for prop in source_properties:
        # Skip properties that don't have either property_type_id or property_type
        if not prop.get('property_type_id') and not prop.get('property_type'):
            logger.warning(f"Skipping property with no type information: {prop.get('id')}")
            continue
        
        # Check if property already exists on target
        if prop.get('property_type_id'):
            cursor.execute("""
                SELECT 1 FROM molecular_properties
                WHERE molecule_id = %s AND property_type_id = %s
            """, (target_id, prop['property_type_id']))
        elif prop.get('property_type'):
            cursor.execute("""
                SELECT 1 FROM molecular_properties
                WHERE molecule_id = %s AND property_type = %s
            """, (target_id, prop['property_type']))
        
        if not cursor.fetchone():
            # Property doesn't exist on target, migrate it
            new_prop_id = str(uuid.uuid4())
            
            # Build column list and values for insert
            columns = ["id", "molecule_id", "data_source", "created_at", "updated_at"]
            values = [new_prop_id, target_id, f"Consolidated from {source_id}", datetime.now(), datetime.now()]
            
            # Add other fields if they have values
            fields = [
                ('property_type_id', prop.get('property_type_id')),
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
                    columns.append(field_name)
                    values.append(field_value)
            
            # Create SQL query
            column_str = ", ".join(columns)
            placeholder_str = ", ".join(["%s"] * len(values))
            
            sql = f"""
                INSERT INTO molecular_properties ({column_str})
                VALUES ({placeholder_str})
                RETURNING id
            """
            
            # Execute the insert
            cursor.execute(sql, values)
            new_prop = cursor.fetchone()
            
            if new_prop:
                change = {
                    'type': 'property_migration',
                    'from_molecule_id': source_id,
                    'to_molecule_id': target_id,
                    'property_type': prop.get('property_type'),
                    'new_property_id': new_prop['id']
                }
                consolidation_log['changes'].append(change)
                logger.info(f"Migrated property {prop.get('property_type') or prop.get('property_type_id')} from {source_id} to {target_id}")

def migrate_mixture_components(cursor, source_id, target_id, consolidation_log):
    """
    Migrate mixture components from source to target molecule.
    """
    # Get all mixture components for the source molecule
    cursor.execute("""
        SELECT * FROM mixture_components
        WHERE molecule_id = %s
    """, (source_id,))
    
    mixture_components = cursor.fetchall()
    logger.info(f"Found {len(mixture_components)} mixture components to migrate")
    
    for component in mixture_components:
        mixture_id = component['mixture_id']
        
        # Check if target molecule already exists in this mixture
        cursor.execute("""
            SELECT 1 FROM mixture_components
            WHERE mixture_id = %s AND molecule_id = %s
        """, (mixture_id, target_id))
        
        if cursor.fetchone():
            logger.info(f"Target molecule already exists in mixture {mixture_id}, skipping")
            continue
        
        # Migrate the mixture component to the target molecule
        new_component_id = str(uuid.uuid4())

        # Build insert query with properties handling
        properties_json = component.get('properties')
        if properties_json:
            # Convert to string if it's a dict
            if isinstance(properties_json, dict):
                properties_json = json.dumps(properties_json)

            # Add consolidated_from info
            properties_json = f'{{"consolidated_from": "{source_id}", {properties_json[1:]}}'
        else:
            properties_json = f'{{"consolidated_from": "{source_id}"}}'

        cursor.execute("""
            INSERT INTO mixture_components (
                id, mixture_id, molecule_id, concentration, concentration_unit,
                role, created_by, created_at, updated_at, properties, units
            )
            VALUES (
                %s, %s, %s, %s, %s,
                %s, %s, NOW(), NOW(), %s, %s
            )
            RETURNING id
        """, (
            new_component_id, mixture_id, target_id,
            component.get('concentration'), component.get('concentration_unit'),
            component.get('role'), component.get('created_by'),
            properties_json, component.get('units')
        ))
        
        new_component = cursor.fetchone()
        
        if new_component:
            change = {
                'type': 'mixture_component_migration',
                'from_molecule_id': source_id,
                'to_molecule_id': target_id,
                'mixture_id': mixture_id,
                'new_component_id': new_component['id']
            }
            consolidation_log['changes'].append(change)
            logger.info(f"Migrated mixture component for mixture {mixture_id} from {source_id} to {target_id}")

def migrate_predictions(cursor, source_id, target_id, consolidation_log):
    """
    Migrate predictions from source to target molecule.
    """
    # Get all predictions for the source molecule
    cursor.execute("""
        SELECT * FROM predictions
        WHERE molecule_id = %s
    """, (source_id,))
    
    predictions = cursor.fetchall()
    logger.info(f"Found {len(predictions)} predictions to migrate")
    
    for prediction in predictions:
        property_type_id = prediction.get('property_type_id')
        
        # Check if a prediction for this property already exists for the target
        cursor.execute("""
            SELECT 1 FROM predictions
            WHERE molecule_id = %s AND property_type_id = %s
        """, (target_id, property_type_id))
        
        if cursor.fetchone():
            logger.info(f"Target molecule already has prediction for property {property_type_id}, skipping")
            continue
        
        # Migrate the prediction to the target molecule
        new_prediction_id = str(uuid.uuid4())

        # Handle data source
        data_source = prediction.get('data_source')
        if not data_source:
            data_source = f"Consolidated from {source_id}"

        # Handle modification history
        mod_history = prediction.get('modification_history')
        if mod_history:
            # Convert to string if it's a dict
            if isinstance(mod_history, dict):
                mod_history = json.dumps(mod_history)

            # Add consolidated_from info
            mod_history = f'{{"consolidated_from": "{source_id}", {mod_history[1:]}}'
        else:
            mod_history = f'{{"consolidated_from": "{source_id}"}}'

        # Handle results
        results = prediction.get('results')
        if results and isinstance(results, dict):
            results = json.dumps(results)

        # Build insert query
        cursor.execute("""
            INSERT INTO predictions (
                id, molecule_id, mixture_id, property_type_id, calculation_method_id,
                numeric_value, text_value, boolean_value, confidence,
                created_at, updated_at, created_by, data_source, version,
                modification_history, name, description, results
            )
            VALUES (
                %s, %s, %s, %s, %s,
                %s, %s, %s, %s,
                NOW(), NOW(), %s, %s, %s,
                %s, %s, %s, %s
            )
            RETURNING id
        """, (
            new_prediction_id, target_id, prediction.get('mixture_id'),
            prediction.get('property_type_id'), prediction.get('calculation_method_id'),
            prediction.get('numeric_value'), prediction.get('text_value'),
            prediction.get('boolean_value'), prediction.get('confidence'),
            prediction.get('created_by'), data_source, prediction.get('version'),
            mod_history, prediction.get('name'), prediction.get('description'),
            results
        ))
        
        new_pred = cursor.fetchone()
        
        if new_pred:
            change = {
                'type': 'prediction_migration',
                'from_molecule_id': source_id,
                'to_molecule_id': target_id,
                'property_type_id': property_type_id,
                'new_prediction_id': new_pred['id']
            }
            consolidation_log['changes'].append(change)
            logger.info(f"Migrated prediction for property {property_type_id} from {source_id} to {target_id}")

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
    Consolidate complex duplicate groups.
    """
    import argparse
    parser = argparse.ArgumentParser(description="Consolidate complex duplicate molecules")
    parser.add_argument("--group", help="Process only a specific group ID")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without making changes")
    args = parser.parse_args()
    
    # Load analysis and plan
    analysis, plan = load_analysis_and_plan()
    complex_groups = plan.get('complex_merge', [])
    
    # Connect to database
    logger.info("Connecting to database...")
    conn = connect_to_db()
    
    try:
        total_processed = 0
        total_success = 0
        
        # Filter to a specific group if requested
        if args.group:
            complex_groups = [g for g in complex_groups if g['group_id'] == args.group]
            if not complex_groups:
                logger.warning(f"Group {args.group} not found or not a complex merge group")
                sys.exit(1)
        
        logger.info(f"Processing {len(complex_groups)} complex merge groups")
        
        for group in complex_groups:
            group_id = group['group_id']
            logger.info(f"Consolidating complex group {group_id} with {group['molecule_count']} molecules")
            total_processed += 1
            
            success = consolidate_complex_group(conn, group, analysis, dry_run=args.dry_run)
            
            if success:
                total_success += 1
        
        # Summary
        logger.info(f"Consolidation {'dry run' if args.dry_run else 'complete'}. " +
                    f"Processed {total_processed} groups, " + 
                    f"{'would have' if args.dry_run else ''} successfully consolidated {total_success} groups.")
        
    except Exception as e:
        logger.error(f"Error during consolidation: {e}")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()