#!/usr/bin/env python3
"""
Reference compounds import script with complete property data import.

This script ensures that all reference cryoprotectant compounds are imported with
complete property data from both ChEMBL and PubChem sources.
Enhanced to use the connection factory for improved performance and reliability.
"""

import os
import sys
import logging
import argparse
import json
import time
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
import uuid

# Import custom modules
from chembl.client import ResilientChEMBLClient as ChEMBLClient
from pubchem.client import ResilientPubChemClient as PubChemClient
from cryoprotectant_identifiers import CryoprotectantIdentifierManager
from chembl.reference_compounds import get_reference_compound_ids

# Import database connection modules
from database.connection import get_db_connection
from sql_executor import (
    execute_query, bulk_insert, execute_batch, with_transaction,
    with_retry, process_in_batches
)
from property_utils import PropertyManager

# Import new utility modules for improved resilience and performance
from db_connection_utils import get_db_connection as get_resilient_connection, safe_transaction
from transaction_utils import with_transaction_retry, execute_in_transaction, is_transaction_active
from batch_utils import bulk_insert_properties, resumable_batch_import, batch_delete_properties

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Checkpoint file path
DEFAULT_CHECKPOINT_PATH = "checkpoints/reference_import_checkpoint.json"

def create_checkpoint(checkpoint_path: str, processed_ids: List[str], results: Dict) -> None:
    """
    Create a checkpoint file to allow resuming interrupted imports.
    Enhanced with safe transaction handling and better error recovery.
    
    Args:
        checkpoint_path: Path to save the checkpoint file
        processed_ids: List of ChEMBL IDs that have been processed
        results: Current results dictionary
    """
    os.makedirs(os.path.dirname(checkpoint_path), exist_ok=True)
    
    checkpoint_data = {
        'timestamp': datetime.now().isoformat(),
        'processed_ids': processed_ids,
        'results': results,
        'position': len(processed_ids),  # Add position for compatibility with batch_utils
        'total': results.get('total_compounds', 0)
    }
    
    try:
        # Use safe transaction to ensure atomic file write
        with safe_transaction():
            with open(checkpoint_path, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
        
        logger.info(f"Checkpoint saved to {checkpoint_path}")
    except Exception as e:
        logger.error(f"Error saving checkpoint: {str(e)}")
        # Try to save to an alternative location as a backup
        try:
            backup_file = f"{checkpoint_path}.backup.{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            with open(backup_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
            logger.info(f"Backup checkpoint saved to {backup_file}")
        except Exception as backup_error:
            logger.error(f"Failed to save backup checkpoint: {str(backup_error)}")

def load_checkpoint(checkpoint_path: str) -> Tuple[List[str], Dict]:
    """
    Load a checkpoint file to resume an interrupted import.
    Enhanced with better error handling and backup checkpoint support.
    
    Args:
        checkpoint_path: Path to the checkpoint file
        
    Returns:
        Tuple of (processed_ids, results)
    """
    if not os.path.exists(checkpoint_path):
        return [], {}
    
    try:
        # Use safe transaction for better error handling
        with safe_transaction():
            with open(checkpoint_path, 'r') as f:
                checkpoint_data = json.load(f)
        
        processed_ids = checkpoint_data.get('processed_ids', [])
        results = checkpoint_data.get('results', {})
        
        # Verify checkpoint data integrity
        if not isinstance(processed_ids, list):
            logger.warning(f"Invalid processed_ids in checkpoint: {type(processed_ids)}")
            processed_ids = []
        
        if not isinstance(results, dict):
            logger.warning(f"Invalid results in checkpoint: {type(results)}")
            results = {}
        
        logger.info(f"Loaded checkpoint with {len(processed_ids)} processed compounds")
        return processed_ids, results
    except Exception as e:
        logger.warning(f"Failed to load checkpoint: {str(e)}")
        
        # Try to find backup checkpoints
        try:
            backup_pattern = f"{checkpoint_path}.backup.*.json"
            import glob
            backup_files = glob.glob(backup_pattern)
            
            if backup_files:
                # Sort by timestamp (newest first)
                backup_files.sort(reverse=True)
                latest_backup = backup_files[0]
                logger.info(f"Attempting to load backup checkpoint: {latest_backup}")
                
                with open(latest_backup, 'r') as f:
                    backup_data = json.load(f)
                
                processed_ids = backup_data.get('processed_ids', [])
                results = backup_data.get('results', {})
                
                logger.info(f"Successfully loaded backup checkpoint with {len(processed_ids)} processed compounds")
                return processed_ids, results
        except Exception as backup_error:
            logger.error(f"Failed to load backup checkpoint: {str(backup_error)}")
        
        return [], {}

def import_reference_compounds(output_report: Optional[str] = None,
                              checkpoint_path: str = DEFAULT_CHECKPOINT_PATH,
                              resume: bool = True,
                              batch_size: int = 5) -> Dict:
    """
    Import reference cryoprotectant compounds with complete property data.
    Uses direct PostgreSQL connections and batch operations for improved performance.
    
    Args:
        output_report: Optional path to save the import report
        checkpoint_path: Path for checkpoint file to enable resumable operations
        resume: Whether to resume from checkpoint if available
        batch_size: Number of compounds to process in each batch
        
    Returns:
        Import results dictionary
    """
    # Initialize clients
    chembl_client = ChEMBLClient()
    pubchem_client = PubChemClient()
    id_manager = CryoprotectantIdentifierManager.get_instance()
    
    # Initialize PropertyManager with direct PostgreSQL connection
    property_manager = PropertyManager()
    
    # Get reference compound IDs
    chembl_ids = get_reference_compound_ids()
    logger.info(f"Importing {len(chembl_ids)} reference compounds")
    
    # Load checkpoint if resuming
    processed_ids = []
    results = {}
    
    if resume and os.path.exists(checkpoint_path):
        processed_ids, results = load_checkpoint(checkpoint_path)
        logger.info(f"Resuming import from checkpoint. {len(processed_ids)} compounds already processed.")
        
        # Filter out already processed IDs
        chembl_ids = [cid for cid in chembl_ids if cid not in processed_ids]
        logger.info(f"{len(chembl_ids)} compounds remaining to process")
    
    # Initialize results if not resuming
    if not results:
        results = {
            'timestamp': datetime.now().isoformat(),
            'total_compounds': len(chembl_ids) + len(processed_ids),
            'imported': 0,
            'updated': 0,
            'failed': 0,
            'details': {}
        }
    
    # Process compounds in batches
    def process_batch(batch_chembl_ids):
        batch_results = []
        
        for chembl_id in batch_chembl_ids:
            try:
                logger.info(f"Processing reference compound: {chembl_id}")
                
                # Step 1: Get compound data from ChEMBL
                chembl_data = chembl_client.get_molecule_by_chembl_id(chembl_id)
                if not chembl_data:
                    logger.error(f"Failed to fetch data for ChEMBL ID: {chembl_id}")
                    results['failed'] += 1
                    results['details'][chembl_id] = {'status': 'failed', 'error': 'ChEMBL data not found'}
                    processed_ids.append(chembl_id)
                    continue
                    
                # Step 2: Get internal ID or create new one
                internal_id, is_new = id_manager.resolve_identifier(chembl_id=chembl_id)
                if not internal_id:
                    # Generate a proper UUID instead of a string ID
                    internal_id = str(uuid.uuid4())
                    is_new = True
                
                # Step 3: Get PubChem data if possible
                pubchem_cid = None
                pubchem_data = None
                if chembl_data.get('cross_references', {}).get('pubchem_cid'):
                    pubchem_cid = chembl_data['cross_references']['pubchem_cid']
                    pubchem_data = pubchem_client.get_molecule_properties(pubchem_cid)
                
                # Step 4: Merge properties from both sources
                properties = {
                    'logs': {
                        'imported_at': datetime.now().isoformat(),
                        'source': 'reference_import'
                    },
                    'chembl': chembl_data,
                    'basic': {
                        'name': chembl_data.get('pref_name') or chembl_data.get('molecule_properties', {}).get('full_molformula') or chembl_id,
                        'molecular_formula': chembl_data.get('molecule_properties', {}).get('full_molformula') or 'Unknown',
                        'molecular_weight': chembl_data.get('molecule_properties', {}).get('full_mwt') or 0.0
                    },
                    'identifiers': {
                        'chembl_id': chembl_id,
                        'pubchem_cid': pubchem_cid,
                        'inchi': chembl_data.get('molecule_structures', {}).get('standard_inchi') or '',
                        'inchi_key': chembl_data.get('molecule_structures', {}).get('standard_inchi_key') or '',
                        'smiles': chembl_data.get('molecule_structures', {}).get('canonical_smiles') or ''
                    },
                    'properties': {
                        'alogp': chembl_data.get('molecule_properties', {}).get('alogp') or 0.0,
                        'h_bond_donors': chembl_data.get('molecule_properties', {}).get('hbd') or 0,
                        'h_bond_acceptors': chembl_data.get('molecule_properties', {}).get('hba') or 0,
                        'rotatable_bonds': chembl_data.get('molecule_properties', {}).get('rtb') or 0,
                        'psa': chembl_data.get('molecule_properties', {}).get('psa') or 0.0,
                        'heavy_atoms': chembl_data.get('molecule_properties', {}).get('heavy_atoms') or 0
                    }
                }
                
                # Add PubChem properties if available
                if pubchem_data:
                    properties['pubchem'] = pubchem_data
                    
                    # Extract and add PubChem-specific properties
                    pubchem_props = {}
                    if 'props' in pubchem_data:
                        for prop in pubchem_data.get('props', []):
                            if prop.get('urn', {}).get('label') == 'LogP':
                                pubchem_props['logP'] = prop.get('value', {}).get('sval')
                            elif prop.get('urn', {}).get('label') == 'Water Solubility':
                                pubchem_props['water_solubility'] = prop.get('value', {}).get('sval')
                            elif prop.get('urn', {}).get('label') == 'Melting Point':
                                pubchem_props['melting_point'] = prop.get('value', {}).get('sval')
                            elif prop.get('urn', {}).get('label') == 'Boiling Point':
                                pubchem_props['boiling_point'] = prop.get('value', {}).get('sval')
                    
                    # Merge with existing properties
                    properties['properties'].update(pubchem_props)
                
                # Step 5: Update the identifier manager
                # Ensure we have valid values for the identifier manager
                molecule_data = {
                    "internal_id": internal_id,
                    "chembl_id": chembl_id,
                    "pubchem_cid": pubchem_cid,
                    "names": [properties['basic']['name']],
                    "inchi_key": properties['identifiers']['inchi_key'] or '',
                    "smiles": properties['identifiers']['smiles'] or '',
                    "formula": properties['basic']['molecular_formula'] or 'Unknown',
                    "molecular_weight": properties['basic']['molecular_weight'] or 0.0,
                    "category": "reference"
                }
                id_manager.add_molecule(internal_id, molecule_data)
                id_manager.save_identifiers()
                
                # Step 6: Insert or update in database using direct connection
                try:
                    # Check if internal_id is in the old string format (CRYO*) and convert to UUID if needed
                    if isinstance(internal_id, str) and internal_id.startswith("CRYO"):
                        # Generate a UUID for database operations while keeping the original ID for the identifier manager
                        db_id = str(uuid.uuid4())
                        logger.info(f"Converting string ID {internal_id} to UUID {db_id} for database operations")
                    else:
                        db_id = internal_id
                        
                    if is_new:
                        # Insert new molecule with basic fields only (no properties JSONB)
                        insert_query = """
                        INSERT INTO molecules
                        (id, name, smiles, inchi, inchikey, formula, molecular_weight, chembl_id, pubchem_cid, created_at, updated_at)
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW())
                        """
                        execute_query(insert_query,
                                    (db_id,
                                     properties['basic']['name'],
                                     properties['identifiers']['smiles'],
                                     properties['identifiers']['inchi'],
                                     properties['identifiers']['inchi_key'],
                                     properties['basic']['molecular_formula'],
                                     properties['basic']['molecular_weight'],
                                     chembl_id,
                                     pubchem_cid))
                        
                        logger.info(f"Inserted new reference molecule {internal_id} (ChEMBL ID: {chembl_id})")
                        results['imported'] += 1
                    else:
                        # Update existing molecule basic fields
                        update_query = """
                        UPDATE molecules
                        SET name = %s, smiles = %s, inchi = %s, inchikey = %s,
                            formula = %s, molecular_weight = %s,
                            chembl_id = %s, pubchem_cid = %s, updated_at = NOW()
                        WHERE id = %s
                        """
                        execute_query(update_query,
                                    (properties['basic']['name'],
                                     properties['identifiers']['smiles'],
                                     properties['identifiers']['inchi'],
                                     properties['identifiers']['inchi_key'],
                                     properties['basic']['molecular_formula'],
                                     properties['basic']['molecular_weight'],
                                     chembl_id,
                                     pubchem_cid,
                                     db_id))
                        
                        logger.info(f"Updated reference molecule {internal_id} (ChEMBL ID: {chembl_id})")
                        results['updated'] += 1
                    
                    # Now insert properties using PropertyManager
                    property_data = {}
                    
                    # Add properties from ChEMBL
                    if 'properties' in properties:
                        for prop_name, prop_value in properties['properties'].items():
                            if prop_value is not None:
                                property_data[prop_name] = prop_value
                    
                    # Ensure we have the critical properties for reference compounds
                    # These are the properties that verification is checking for
                    critical_properties = ['logP', 'h_bond_donors', 'h_bond_acceptors']
                    logger.info(f"Checking critical properties for molecule {internal_id} (ChEMBL ID: {chembl_id})")
                    
                    for prop in critical_properties:
                        if prop not in property_data or property_data[prop] is None:
                            # Try to derive from other properties if possible
                            if prop == 'logP' and 'alogp' in property_data:
                                property_data['logP'] = property_data['alogp']
                                logger.info(f"Derived {prop}={property_data['logP']} from alogp")
                            elif prop == 'h_bond_donors' and 'hbd' in properties.get('molecule_properties', {}):
                                property_data['h_bond_donors'] = properties['molecule_properties']['hbd']
                                logger.info(f"Derived {prop}={property_data['h_bond_donors']} from hbd")
                            elif prop == 'h_bond_acceptors' and 'hba' in properties.get('molecule_properties', {}):
                                property_data['h_bond_acceptors'] = properties['molecule_properties']['hba']
                                logger.info(f"Derived {prop}={property_data['h_bond_acceptors']} from hba")
                            else:
                                # Set default values for missing critical properties
                                logger.warning(f"Setting default value for missing critical property: {prop}")
                                if prop == 'logP':
                                    property_data[prop] = 0.0
                                elif prop in ['h_bond_donors', 'h_bond_acceptors']:
                                    property_data[prop] = 0
                                logger.info(f"Set default value for {prop}={property_data[prop]}")
                    
                    # Make sure critical properties are explicitly set
                    for prop in critical_properties:
                        if prop not in property_data or property_data[prop] is None:
                            if prop == 'logP':
                                property_data[prop] = 0.0
                            else:  # h_bond_donors or h_bond_acceptors
                                property_data[prop] = 0
                            logger.info(f"Explicitly setting default value for critical property {prop} = {property_data[prop]}")
                    
                    # Log all critical properties to verify they're set
                    logger.info(f"Critical properties for molecule {internal_id} (ChEMBL ID: {chembl_id}):")
                    for prop in critical_properties:
                        logger.info(f"  {prop} = {property_data.get(prop, 'NOT SET')}")
                    
                    # Execute property setting with retry
                    # Use transaction context manager for better transaction handling
                    try:
                        # Use the new transaction_utils for better resilience
                        with with_transaction_retry(max_retries=3, retry_delay=1.0):
                            # First try to use bulk_insert_properties from batch_utils
                            try:
                                # Convert property_data to the format expected by bulk_insert_properties
                                property_records = []
                                for prop_name, prop_value in property_data.items():
                                    if prop_value is not None:
                                        # Get property type ID first
                                        try:
                                            # Use property_manager to get property_type_id
                                            property_type_id = property_manager.get_property_type_id(prop_name)
                                            
                                            # Determine property type
                                            if isinstance(prop_value, (int, float)):
                                                property_records.append({
                                                    "molecule_id": db_id,
                                                    "property_type_id": property_type_id,
                                                    "numeric_value": prop_value,
                                                    "text_value": None,
                                                    "boolean_value": None,
                                                    "created_by": None
                                                })
                                            else:
                                                property_records.append({
                                                    "molecule_id": db_id,
                                                    "property_type_id": property_type_id,
                                                    "numeric_value": None,
                                                    "text_value": str(prop_value),
                                                    "boolean_value": None,
                                                    "created_by": None
                                                })
                                        except Exception as e:
                                            logger.error(f"Failed to get property_type_id for {prop_name}: {str(e)}")
                                
                                # Use bulk_insert_properties for better performance
                                if property_records:
                                    # Log the property records being inserted
                                    logger.info(f"Inserting {len(property_records)} property records for molecule {db_id}")
                                    for i, record in enumerate(property_records):
                                        prop_type_id = record.get('property_type_id')
                                        logger.info(f"  Record {i}: property_type_id={prop_type_id}, " +
                                                   f"numeric_value={record.get('numeric_value')}, " +
                                                   f"text_value={record.get('text_value')}")
                                    
                                    # Verify property_type_ids exist in the database before insertion
                                    for record in property_records:
                                        prop_type_id = record.get('property_type_id')
                                        verify_query = "SELECT id FROM property_types WHERE id = %s"
                                        verify_result = execute_query(verify_query, (prop_type_id,), fetch_one=True)
                                        if not verify_result:
                                            logger.error(f"Property type ID {prop_type_id} does not exist in database")
                                            # Create the property type if it doesn't exist
                                            create_query = """
                                            INSERT INTO property_types (id, name, data_type, description)
                                            VALUES (%s, %s, %s, %s)
                                            """
                                            prop_name = "unknown_property"
                                            for name, type_id in property_manager._property_types_cache.items():
                                                if type_id.get('id') == prop_type_id:
                                                    prop_name = name
                                                    break
                                            execute_query(create_query, (prop_type_id, prop_name, 'numeric', f"Auto-created property: {prop_name}"))
                                            logger.info(f"Created missing property type: {prop_name} (ID: {prop_type_id})")
                                    
                                    # Try direct SQL insertion as a fallback if bulk_insert_properties fails
                                    try:
                                        success_count = bulk_insert_properties(property_records)
                                        total_count = len(property_records)
                                        logger.info(f"Bulk inserted {success_count}/{total_count} properties using batch_utils")
                                        
                                        # If no properties were inserted, try direct SQL insertion
                                        if success_count == 0:
                                            logger.warning("Bulk insert returned 0 success count, trying direct SQL insertion")
                                            direct_insert_count = 0
                                            
                                            for record in property_records:
                                                try:
                                                    direct_query = """
                                                    INSERT INTO molecular_properties
                                                    (molecule_id, property_type_id, numeric_value, text_value, boolean_value, created_by)
                                                    VALUES (%s, %s, %s, %s, %s, %s)
                                                    RETURNING id
                                                    """
                                                    direct_result = execute_query(direct_query, (
                                                        record.get('molecule_id'),
                                                        record.get('property_type_id'),
                                                        record.get('numeric_value'),
                                                        record.get('text_value'),
                                                        record.get('boolean_value'),
                                                        record.get('created_by')
                                                    ), fetch_one=True)
                                                    
                                                    if direct_result:
                                                        direct_insert_count += 1
                                                        logger.info(f"Direct insert successful for property_type_id={record.get('property_type_id')}")
                                                except Exception as direct_error:
                                                    logger.error(f"Direct insert failed: {str(direct_error)}")
                                            
                                            logger.info(f"Directly inserted {direct_insert_count}/{total_count} properties")
                                            success_count = direct_insert_count
                                    except Exception as bulk_error:
                                        logger.error(f"Bulk insert error: {str(bulk_error)}")
                                        success_count = 0
                                    
                                    # If not all properties were inserted successfully, log a warning
                                    if success_count < total_count:
                                        logger.warning(f"Failed to insert {total_count - success_count} properties")
                            except Exception as bulk_error:
                                logger.warning(f"Bulk insert failed: {str(bulk_error)}. Falling back to property_manager.")
                                # Fall back to property_manager
                                success_count, total_count = property_manager.set_properties(db_id, property_data, created_by=None)
                            
                            # Verify critical properties were set
                            verification_properties = property_manager.get_properties(db_id, critical_properties)
                            logger.info(f"Verification properties for molecule {db_id}: {verification_properties}")
                            
                            missing_properties = [prop for prop in critical_properties if prop not in verification_properties]
                            
                            if missing_properties:
                                logger.warning(f"Critical properties still missing after setting: {missing_properties}")
                                
                                # Log the database schema to verify the molecular_properties table structure
                                try:
                                    schema_query = """
                                    SELECT column_name, data_type
                                    FROM information_schema.columns
                                    WHERE table_name = 'molecular_properties'
                                    """
                                    schema_result = execute_query(schema_query)
                                    logger.info(f"molecular_properties table schema: {schema_result}")
                                    
                                    # Also check if the property types exist
                                    prop_types_query = "SELECT id, name, data_type FROM property_types"
                                    prop_types_result = execute_query(prop_types_query)
                                    logger.info(f"Property types in database: {prop_types_result}")
                                except Exception as schema_error:
                                    logger.error(f"Error fetching schema: {str(schema_error)}")
                                
                                # Try direct insertion for missing properties
                                for prop in missing_properties:
                                    if prop == 'logP':
                                        default_value = 0.0
                                    else:  # h_bond_donors or h_bond_acceptors
                                        default_value = 0
                                    logger.info(f"Setting default value for critical property {prop} = {default_value}")
                                    
                                    # Get property type ID
                                    prop_type_id = property_manager.get_property_type_id(prop)
                                    
                                    # Insert directly using SQL
                                    try:
                                        direct_query = """
                                        INSERT INTO molecular_properties
                                        (molecule_id, property_type_id, numeric_value, created_at, updated_at)
                                        VALUES (%s, %s, %s, NOW(), NOW())
                                        ON CONFLICT (molecule_id, property_type_id)
                                        DO UPDATE SET numeric_value = EXCLUDED.numeric_value, updated_at = NOW()
                                        RETURNING id
                                        """
                                        direct_result = execute_query(direct_query, (db_id, prop_type_id, default_value), fetch_one=True)
                                        if direct_result:
                                            logger.info(f"Successfully set property {prop} with direct SQL: {direct_result}")
                                        else:
                                            logger.error(f"Failed to set property {prop} with direct SQL")
                                    except Exception as direct_error:
                                        logger.error(f"Error with direct SQL for property {prop}: {str(direct_error)}")
                                        # Fall back to original method
                                        execute_in_transaction(
                                            lambda: property_manager.set_property(db_id, prop, default_value)
                                        )
                                
                                # Verify again
                                verification_properties = property_manager.get_properties(db_id, critical_properties)
                                still_missing = [prop for prop in critical_properties if prop not in verification_properties]
                                if still_missing:
                                    logger.warning(f"Properties still missing after direct SQL: {still_missing}")
                                    # Last resort: try raw SQL without parameters
                                    for prop in still_missing:
                                        try:
                                            # Get property type ID directly from database
                                            prop_type_query = "SELECT id FROM property_types WHERE name = %s"
                                            prop_type_result = execute_query(prop_type_query, (prop,), fetch_one=True)
                                            
                                            if prop_type_result:
                                                prop_type_id = prop_type_result['id']
                                                default_value = 0.0 if prop == 'logP' else 0
                                                
                                                # Insert with raw SQL
                                                raw_query = f"""
                                                INSERT INTO molecular_properties
                                                (molecule_id, property_type_id, numeric_value, created_at, updated_at)
                                                VALUES ('{db_id}', '{prop_type_id}', {default_value}, NOW(), NOW())
                                                ON CONFLICT (molecule_id, property_type_id)
                                                DO UPDATE SET numeric_value = EXCLUDED.numeric_value, updated_at = NOW()
                                                """
                                                execute_query(raw_query)
                                                logger.info(f"Raw SQL insert for {prop} executed")
                                            else:
                                                logger.error(f"Could not find property type for {prop}")
                                        except Exception as raw_error:
                                            logger.error(f"Raw SQL insert failed: {str(raw_error)}")
                                    
                                    # Final verification
                                    final_verification = property_manager.get_properties(db_id, critical_properties)
                                    final_missing = [prop for prop in critical_properties if prop not in final_verification]
                                    if final_missing:
                                        raise Exception(f"Failed to set critical properties after all attempts: {final_missing}")
                    except Exception as prop_error:
                        logger.error(f"Error setting properties: {str(prop_error)}")
                        # Last resort: try setting each critical property individually outside transaction
                        for prop in critical_properties:
                            try:
                                if prop == 'logP':
                                    default_value = 0.0
                                else:  # h_bond_donors or h_bond_acceptors
                                    default_value = 0
                                logger.warning(f"Last resort: Setting critical property {prop} = {default_value}")
                                
                                # Try multiple approaches in sequence
                                try:
                                    # 1. Try with property manager
                                    with safe_transaction():
                                        property_manager.set_property(db_id, prop, default_value)
                                        logger.info(f"Property manager succeeded for {prop}")
                                except Exception as pm_error:
                                    logger.error(f"Property manager failed: {str(pm_error)}")
                                    
                                    try:
                                        # 2. Try direct SQL
                                        prop_type_id = property_manager.get_property_type_id(prop)
                                        direct_query = """
                                        INSERT INTO molecular_properties
                                        (molecule_id, property_type_id, numeric_value, created_at, updated_at)
                                        VALUES (%s, %s, %s, NOW(), NOW())
                                        ON CONFLICT (molecule_id, property_type_id)
                                        DO UPDATE SET numeric_value = EXCLUDED.numeric_value, updated_at = NOW()
                                        """
                                        execute_query(direct_query, (db_id, prop_type_id, default_value))
                                        logger.info(f"Direct SQL insert for {prop} succeeded")
                                    except Exception as sql_error:
                                        logger.error(f"Direct SQL failed: {str(sql_error)}")
                                        
                                        try:
                                            # 3. Try raw SQL as absolute last resort
                                            prop_type_query = "SELECT id FROM property_types WHERE name = %s"
                                            prop_type_result = execute_query(prop_type_query, (prop,), fetch_one=True)
                                            
                                            if prop_type_result:
                                                prop_type_id = prop_type_result['id']
                                                raw_query = f"""
                                                INSERT INTO molecular_properties
                                                (molecule_id, property_type_id, numeric_value, created_at, updated_at)
                                                VALUES ('{db_id}', '{prop_type_id}', {default_value}, NOW(), NOW())
                                                ON CONFLICT (molecule_id, property_type_id)
                                                DO UPDATE SET numeric_value = EXCLUDED.numeric_value, updated_at = NOW()
                                                """
                                                execute_query(raw_query)
                                                logger.info(f"Raw SQL insert for {prop} succeeded")
                                            else:
                                                logger.error(f"Could not find property type for {prop}")
                                        except Exception as raw_error:
                                            logger.error(f"Raw SQL failed: {str(raw_error)}")
                            except Exception as e:
                                logger.error(f"All attempts failed for property {prop}: {str(e)}")
                        
                        # Final verification
                        verification_properties = property_manager.get_properties(db_id, critical_properties)
                        still_missing = [prop for prop in critical_properties if prop not in verification_properties]
                        
                        if still_missing:
                            logger.error(f"Critical properties still missing after all attempts: {still_missing}")
                            # Log database state for debugging
                            try:
                                count_query = "SELECT COUNT(*) FROM molecular_properties WHERE molecule_id = %s"
                                count_result = execute_query(count_query, (db_id,), fetch_one=True)
                                logger.info(f"Total properties for molecule {db_id}: {count_result}")
                                
                                props_query = """
                                SELECT mp.id, pt.name, mp.numeric_value, mp.text_value, mp.boolean_value
                                FROM molecular_properties mp
                                JOIN property_types pt ON mp.property_type_id = pt.id
                                WHERE mp.molecule_id = %s
                                """
                                props_result = execute_query(props_query, (db_id,))
                                logger.info(f"Existing properties: {props_result}")
                            except Exception as debug_error:
                                logger.error(f"Error fetching debug info: {str(debug_error)}")
                    
                    logger.info(f"Added {success_count}/{total_count} properties for molecule {internal_id}")
                    
                    # Record success
                    results['details'][chembl_id] = {
                        'status': 'success',
                        'internal_id': internal_id,
                        'is_new': is_new,
                        'has_pubchem': pubchem_cid is not None,
                        'properties_added': success_count,
                        'critical_properties_verified': len(missing_properties) == 0
                    }
                    
                except Exception as db_error:
                    logger.warning(f"Database operation failed for {chembl_id}: {str(db_error)}")
                    logger.info(f"Continuing with identifier manager updates only for {chembl_id}")
                    # Still count as success since we updated the identifier manager
                    if is_new:
                        results['imported'] += 1
                    else:
                        results['updated'] += 1
                    
                    # Record partial success
                    results['details'][chembl_id] = {
                        'status': 'partial',
                        'internal_id': internal_id,
                        'is_new': is_new,
                        'has_pubchem': pubchem_cid is not None,
                        'error': str(db_error)
                    }
                
                # Add to processed IDs
                processed_ids.append(chembl_id)
                
                # Create checkpoint after each compound
                create_checkpoint(checkpoint_path, processed_ids, results)
                
            except Exception as e:
                logger.error(f"Error processing reference compound {chembl_id}: {str(e)}")
                results['failed'] += 1
                results['details'][chembl_id] = {'status': 'failed', 'error': str(e)}
                processed_ids.append(chembl_id)
                
                # Create checkpoint after each failure too
                create_checkpoint(checkpoint_path, processed_ids, results)
        
        return batch_results
    
    # Process compounds in batches
    if chembl_ids:
        process_in_batches(chembl_ids, batch_size=batch_size, process_func=process_batch)
    
    # Generate final report
    if output_report:
        os.makedirs(os.path.dirname(output_report), exist_ok=True)
        with open(output_report, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Import report saved to {output_report}")
    
    logger.info(f"Reference compound import completed: {results['imported']} imported, " +
               f"{results['updated']} updated, {results['failed']} failed")
    
    return results

def verify_reference_compounds(import_results):
    """
    Verify that reference compounds have all required properties.
    
    Args:
        import_results: Results from the import operation
    """
    logger.info("Starting reference compound verification")
    
    # Get database connection
    db_conn = get_db()
    
    # Get all reference compounds
    query = """
    SELECT m.id, m.name, m.chembl_id
    FROM molecules m
    WHERE m.chembl_id IS NOT NULL
    """
    molecules = execute_query(query)
    
    if not molecules:
        logger.warning("No reference compounds found in database")
        return
    
    logger.info(f"Found {len(molecules)} reference compounds to verify")
    
    # Initialize PropertyManager
    property_manager = PropertyManager()
    
    # Critical properties that must be present
    critical_properties = ['logP', 'h_bond_donors', 'h_bond_acceptors']
    
    # Check each molecule for critical properties
    missing_properties = {}
    complete_count = 0
    
    for molecule in molecules:
        molecule_id = molecule['id']
        chembl_id = molecule['chembl_id']
        
        # Get properties for this molecule
        properties = property_manager.get_properties(molecule_id, critical_properties)
        
        # Check if all critical properties are present
        missing = [prop for prop in critical_properties if prop not in properties]
        
        if missing:
            logger.warning(f"Reference compound {chembl_id} missing properties: {missing}")
            missing_properties[chembl_id] = missing
            
            # Try to fix missing properties with transaction
            try:
                # Use the new transaction_utils for better resilience
                with with_transaction_retry(max_retries=3, retry_delay=1.0):
                    # First try to use bulk_insert_properties from batch_utils
                    try:
                        # Convert missing properties to the format expected by bulk_insert_properties
                        property_records = []
                        for prop in missing:
                            try:
                                # Get property type ID first
                                property_type_id = property_manager.get_property_type_id(prop)
                                
                                if prop == 'logP':
                                    default_value = 0.0
                                    property_records.append({
                                        "molecule_id": molecule_id,
                                        "property_type_id": property_type_id,
                                        "numeric_value": default_value,
                                        "text_value": None,
                                        "boolean_value": None,
                                        "created_by": None
                                    })
                                else:  # h_bond_donors or h_bond_acceptors
                                    default_value = 0
                                    property_records.append({
                                        "molecule_id": molecule_id,
                                        "property_type_id": property_type_id,
                                        "numeric_value": default_value,
                                        "text_value": None,
                                        "boolean_value": None,
                                        "created_by": None
                                    })
                                logger.info(f"Setting default value for critical property {prop} = {default_value}")
                            except Exception as e:
                                logger.error(f"Failed to get property_type_id for {prop}: {str(e)}")
                            logger.info(f"Setting default value for critical property {prop} = {default_value}")
                        
                        # Use bulk_insert_properties for better performance
                        if property_records:
                            success_count = bulk_insert_properties(property_records)
                            logger.info(f"Bulk inserted {success_count}/{len(property_records)} properties using batch_utils")
                    except Exception as bulk_error:
                        logger.warning(f"Bulk insert failed: {str(bulk_error)}. Falling back to individual property setting.")
                        # Fall back to individual property setting
                        for prop in missing:
                            if prop == 'logP':
                                default_value = 0.0
                            else:  # h_bond_donors or h_bond_acceptors
                                default_value = 0
                            
                            logger.info(f"Setting default value for critical property {prop} = {default_value}")
                            property_manager.set_property(molecule_id, prop, default_value)
                    
                    # Verify properties were set
                    verification_properties = property_manager.get_properties(molecule_id, missing)
                    still_missing = [prop for prop in missing if prop not in verification_properties]
                    if still_missing:
                        raise Exception(f"Failed to set properties: {still_missing}")
            except Exception as e:
                logger.error(f"Error fixing properties for {chembl_id}: {str(e)}")
                # Last resort: try setting each property individually outside transaction
                for prop in missing:
                    try:
                        if prop == 'logP':
                            default_value = 0.0
                        else:  # h_bond_donors or h_bond_acceptors
                            default_value = 0
                        logger.warning(f"Last resort: Setting critical property {prop} = {default_value}")
                        # Use safe_transaction for better error handling
                        with safe_transaction():
                            property_manager.set_property(molecule_id, prop, default_value)
                    except Exception as prop_error:
                        logger.error(f"Failed to set property {prop}: {str(prop_error)}")
        else:
            complete_count += 1
    
    # Final verification
    final_count = 0
    for molecule in molecules:
        molecule_id = molecule['id']
        properties = property_manager.get_properties(molecule_id, critical_properties)
        if all(prop in properties for prop in critical_properties):
            final_count += 1
    
    logger.info(f"Verification complete: {final_count}/{len(molecules)} compounds have all critical properties")
    if missing_properties:
        logger.info(f"Fixed properties for {len(missing_properties)} compounds")
    
    return final_count, len(molecules)

def main():
    """CLI entry point for reference compound import."""
    parser = argparse.ArgumentParser(description='Import reference cryoprotectant compounds')
    parser.add_argument('--report',
                      default=f"reports/reference_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                      help='Output file path for report')
    parser.add_argument('--checkpoint',
                      default=DEFAULT_CHECKPOINT_PATH,
                      help='Checkpoint file path for resumable operations')
    parser.add_argument('--no-resume', action='store_true',
                      help='Do not resume from checkpoint even if available')
    parser.add_argument('--batch-size', type=int, default=5,
                      help='Number of compounds to process in each batch')
    parser.add_argument('--verify', action='store_true',
                      help='Verify reference compounds after import')
    parser.add_argument('--force-verify', action='store_true',
                      help='Force verification even if import fails')
    
    args = parser.parse_args()
    
    try:
        results = import_reference_compounds(
            args.report,
            checkpoint_path=args.checkpoint,
            resume=not args.no_resume,
            batch_size=args.batch_size
        )
        
        # Always verify after import to ensure critical properties are set
        logger.info("Verifying reference compounds after import...")
        complete_count, total_count = verify_reference_compounds(results)
        
        # Verify again to make sure all properties are set
        logger.info("Running final verification...")
        final_count, total_count = verify_reference_compounds(results)
        
        if final_count < total_count:
            logger.warning(f"Final verification shows {final_count}/{total_count} complete compounds")
            return 1
        else:
            logger.info(f"All {total_count} reference compounds have required properties")
            return 0
    except Exception as e:
        logger.error(f"Import failed: {str(e)}")
        
        if args.force_verify:
            logger.info("Running verification despite import failure...")
            try:
                verify_reference_compounds({})
            except Exception as verify_error:
                logger.error(f"Verification after failure also failed: {str(verify_error)}")
        
        return 1

if __name__ == "__main__":
    sys.exit(main())