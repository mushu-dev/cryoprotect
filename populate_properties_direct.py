"""
Direct Molecular Properties Population Script

This script populates the molecular_properties table using direct PostgreSQL connection
for optimal performance, especially with large datasets.

Usage:
    python populate_properties_direct.py [--batch-size 1000] [--limit 5000]
"""

import os
import sys
import json
import time
import uuid
import logging
import argparse
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from datetime import datetime
from dotenv import load_dotenv

# Import direct connection modules
from postgres_direct import PostgresDirectConnection
from sql_executor import execute_query, bulk_insert, transaction, close_all

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/direct_properties_population.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Ensure log directory exists
Path("logs").mkdir(exist_ok=True)

# Checkpoint file for resuming imports
CHECKPOINT_FILE = "checkpoints/properties_import_direct.json"

def ensure_checkpoint_dir():
    """Ensure the checkpoint directory exists."""
    Path("checkpoints").mkdir(exist_ok=True)

def save_checkpoint(progress: Dict[str, Any]):
    """Save progress checkpoint to file."""
    ensure_checkpoint_dir()
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(progress, f)
    logger.info(f"Checkpoint saved: {progress}")

def load_checkpoint() -> Dict[str, Any]:
    """Load progress checkpoint from file."""
    ensure_checkpoint_dir()
    try:
        if os.path.exists(CHECKPOINT_FILE):
            with open(CHECKPOINT_FILE, 'r') as f:
                checkpoint = json.load(f)
            logger.info(f"Checkpoint loaded: {checkpoint}")
            return checkpoint
        else:
            logger.info("No checkpoint file found, starting fresh")
            return {'completed_molecules': 0, 'completed_properties': 0, 'last_molecule_id': None, 'timestamp': None}
    except Exception as e:
        logger.warning(f"Error loading checkpoint: {e}. Starting fresh.")
        return {'completed_molecules': 0, 'completed_properties': 0, 'last_molecule_id': None, 'timestamp': None}

def verify_database_schema() -> bool:
    """Verify that required tables exist with the expected schema."""
    try:
        # Check if the molecules table exists
        result = execute_query("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = 'public' AND table_name = 'molecules'
            ) as exists
        """)
        
        if not result or not result[0]['exists']:
            logger.error("The 'molecules' table does not exist")
            return False
        
        # Check if the molecular_properties table exists
        result = execute_query("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = 'public' AND table_name = 'molecular_properties'
            ) as exists
        """)
        
        if not result or not result[0]['exists']:
            logger.error("The 'molecular_properties' table does not exist")
            return False
        
        # Check if the required columns exist in molecular_properties
        result = execute_query("""
            SELECT column_name 
            FROM information_schema.columns 
            WHERE table_schema = 'public' AND table_name = 'molecular_properties'
        """)
        
        columns = [row['column_name'] for row in result]
        required_columns = ['id', 'molecule_id', 'property_name', 'property_value', 'property_type']
        
        for column in required_columns:
            if column not in columns:
                logger.error(f"Required column '{column}' is missing from 'molecular_properties' table")
                return False
        
        logger.info("Database schema verified successfully")
        return True
    except Exception as e:
        logger.error(f"Error verifying database schema: {e}")
        return False

def get_molecules(offset: int, limit: int) -> List[Dict[str, Any]]:
    """
    Get a batch of molecules from the database.
    
    Args:
        offset: Number of molecules to skip
        limit: Maximum number of molecules to return
        
    Returns:
        List of molecule dictionaries
    """
    query = """
        SELECT id, name
        FROM molecules
        ORDER BY id
        OFFSET %s LIMIT %s
    """
    
    try:
        result = execute_query(query, (offset, limit))
        return result
    except Exception as e:
        logger.error(f"Error fetching molecules: {e}")
        return []

def get_molecule_count() -> int:
    """Get the total number of molecules in the database."""
    try:
        result = execute_query("SELECT COUNT(*) as count FROM molecules")
        return result[0]['count'] if result else 0
    except Exception as e:
        logger.error(f"Error getting molecule count: {e}")
        return 0

def load_property_data(file_path: str) -> Dict[str, Dict[str, Any]]:
    """
    Load property data from a JSON file.
    
    Returns:
        Dictionary mapping molecule name to properties
    """
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        property_data = {}
        for entry in data:
            name = entry.get('name')
            if name:
                property_data[name] = entry
        
        logger.info(f"Loaded properties for {len(property_data)} molecules from {file_path}")
        return property_data
    except Exception as e:
        logger.error(f"Error loading property data from {file_path}: {e}")
        return {}

def prepare_property_batch(molecules: List[Dict[str, Any]],
                          property_data: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Prepare a batch of molecular properties for insertion.
    
    Args:
        molecules: List of molecule dictionaries with id and name
        property_data: Dictionary mapping molecule name to properties
        
    Returns:
        List of property dictionaries ready for insertion
    """
    property_records = []
    processed_molecules = 0
    
    for molecule in molecules:
        molecule_id = molecule['id']
        molecule_name = molecule['name']
        
        # Look up properties for this molecule
        molecule_props = property_data.get(molecule_name, {})
        if not molecule_props:
            # Skip molecules without property data
            continue
        
        # Add each property as a separate record
        for prop_name, prop_value in molecule_props.items():
            # Skip name as it's not a property
            if prop_name == 'name':
                continue
            
            # Skip empty or None values
            if prop_value is None or prop_value == '':
                continue
            
            # Determine property type based on value
            if isinstance(prop_value, (int, float)):
                prop_type = 'numeric'
            elif isinstance(prop_value, bool):
                prop_type = 'boolean'
            else:
                prop_type = 'text'
            
            # Create property record
            property_record = {
                'id': str(uuid.uuid4()),
                'molecule_id': molecule_id,
                'property_name': prop_name,
                'property_value': str(prop_value),
                'property_type': prop_type,
                'created_at': datetime.now().isoformat(),
                'updated_at': datetime.now().isoformat()
            }
            
            property_records.append(property_record)
        
        processed_molecules += 1
    
    return property_records, processed_molecules

def populate_properties(property_file: str, batch_size: int = 1000, molecule_limit: Optional[int] = None):
    """
    Populate the molecular_properties table using direct connection.
    
    Args:
        property_file: Path to the JSON file containing property data
        batch_size: Number of molecules to process in each batch
        molecule_limit: Optional limit on the total number of molecules to process
    """
    # Verify database schema
    if not verify_database_schema():
        logger.error("Database schema verification failed. Aborting.")
        return
    
    # Load property data
    property_data = load_property_data(property_file)
    if not property_data:
        logger.error("No property data loaded. Aborting.")
        return
    
    # Get molecule count
    total_molecules = get_molecule_count()
    if total_molecules == 0:
        logger.error("No molecules found in database. Aborting.")
        return
    
    # Apply limit if specified
    if molecule_limit is not None:
        total_molecules = min(total_molecules, molecule_limit)
    
    logger.info(f"Total molecules to process: {total_molecules}")
    
    # Load checkpoint
    checkpoint = load_checkpoint()
    molecule_offset = checkpoint['completed_molecules']
    total_properties_inserted = checkpoint['completed_properties']
    
    logger.info(f"Starting property population from molecule offset {molecule_offset} with {total_properties_inserted} properties already inserted")
    
    start_time = time.time()
    
    try:
        # Process molecules in batches
        while molecule_offset < total_molecules:
            batch_start_time = time.time()
            
            # Get batch of molecules
            molecules = get_molecules(molecule_offset, batch_size)
            if not molecules:
                logger.warning(f"No molecules found at offset {molecule_offset}. Stopping.")
                break
            
            logger.info(f"Processing batch: molecules {molecule_offset+1}-{molecule_offset+len(molecules)} of {total_molecules}")
            
            # Prepare property records
            property_records, processed_molecules = prepare_property_batch(molecules, property_data)
            
            if property_records:
                logger.info(f"Inserting {len(property_records)} property records for {processed_molecules} molecules")
                
                try:
                    # Insert property records
                    bulk_insert('molecular_properties', property_records)
                    
                    # Update progress
                    molecule_offset += len(molecules)
                    total_properties_inserted += len(property_records)
                    
                    checkpoint = {
                        'completed_molecules': molecule_offset,
                        'completed_properties': total_properties_inserted,
                        'last_molecule_id': molecules[-1]['id'] if molecules else None,
                        'timestamp': datetime.now().isoformat()
                    }
                    save_checkpoint(checkpoint)
                    
                    batch_time = time.time() - batch_start_time
                    logger.info(f"Batch completed in {batch_time:.2f} seconds "
                              f"({len(property_records) / batch_time:.2f} properties/sec)")
                    
                except Exception as e:
                    logger.error(f"Error inserting property batch: {e}")
                    # Save checkpoint at the last successful batch
                    save_checkpoint(checkpoint)
                    raise
            else:
                logger.info(f"No properties found for molecules in this batch. Skipping.")
                # Update progress even though no properties were inserted
                molecule_offset += len(molecules)
                checkpoint = {
                    'completed_molecules': molecule_offset,
                    'completed_properties': total_properties_inserted,
                    'last_molecule_id': molecules[-1]['id'] if molecules else None,
                    'timestamp': datetime.now().isoformat()
                }
                save_checkpoint(checkpoint)
        
        total_time = time.time() - start_time
        logger.info(f"Population completed: {total_properties_inserted} properties inserted for {molecule_offset} molecules in {total_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"Error during property population: {e}")
        raise
    finally:
        # Close all connections
        close_all()

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Populate molecular_properties table using direct connection")
    parser.add_argument("--file", default="data/cryoprotectant_master_list.json", 
                        help="Path to JSON file with property data")
    parser.add_argument("--batch-size", type=int, default=1000, 
                        help="Number of molecules to process in each batch")
    parser.add_argument("--limit", type=int, default=None, 
                        help="Limit the total number of molecules to process")
    
    args = parser.parse_args()
    
    try:
        # Initialize database connection
        db = PostgresDirectConnection.get_instance()
        if not db.test_connection():
            logger.error("Database connection failed. Aborting.")
            return 1
        
        # Run population
        populate_properties(args.file, args.batch_size, args.limit)
        return 0
        
    except Exception as e:
        logger.error(f"Error in main: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())