"""
Direct Molecule Population Script

This script populates the molecules table using direct PostgreSQL connection
for optimal performance, especially with large datasets.

Usage:
    python populate_molecules_direct.py [--batch-size 500] [--limit 1000]
"""

import os
import sys
import json
import time
import uuid
import logging
import argparse
from pathlib import Path
from typing import List, Dict, Any, Optional
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
        logging.FileHandler("logs/direct_population.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Ensure log directory exists
Path("logs").mkdir(exist_ok=True)

# Checkpoint file for resuming imports
CHECKPOINT_FILE = "checkpoints/molecule_import_direct.json"

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
            return {'completed': 0, 'last_id': None, 'timestamp': None}
    except Exception as e:
        logger.warning(f"Error loading checkpoint: {e}. Starting fresh.")
        return {'completed': 0, 'last_id': None, 'timestamp': None}

def verify_database_schema():
    """Verify that the molecules table exists with the expected schema."""
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
        
        # Check if the required columns exist
        result = execute_query("""
            SELECT column_name 
            FROM information_schema.columns 
            WHERE table_schema = 'public' AND table_name = 'molecules'
        """)
        
        columns = [row['column_name'] for row in result]
        required_columns = ['id', 'name', 'formula', 'molecular_weight', 'smiles']
        
        for column in required_columns:
            if column not in columns:
                logger.error(f"Required column '{column}' is missing from 'molecules' table")
                return False
        
        logger.info("Database schema verified successfully")
        return True
    except Exception as e:
        logger.error(f"Error verifying database schema: {e}")
        return False

def load_molecules_from_file(file_path: str) -> List[Dict[str, Any]]:
    """Load molecule data from a JSON file."""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        logger.info(f"Loaded {len(data)} molecules from {file_path}")
        return data
    except Exception as e:
        logger.error(f"Error loading molecules from {file_path}: {e}")
        return []

def prepare_molecule_batch(molecules: List[Dict[str, Any]], 
                          start_index: int, 
                          batch_size: int) -> List[Dict[str, Any]]:
    """Prepare a batch of molecules for insertion."""
    end_index = min(start_index + batch_size, len(molecules))
    batch = molecules[start_index:end_index]
    
    prepared_data = []
    for mol in batch:
        # Generate a UUID if not present
        if 'id' not in mol:
            mol['id'] = str(uuid.uuid4())
        
        # Add created_at and updated_at timestamps if not present
        if 'created_at' not in mol:
            mol['created_at'] = datetime.now().isoformat()
        if 'updated_at' not in mol:
            mol['updated_at'] = datetime.now().isoformat()
        
        # Add owner_id if it exists as a column
        if 'owner_id' not in mol:
            mol['owner_id'] = os.getenv('DEFAULT_OWNER_ID', 'system')
        
        prepared_data.append(mol)
    
    return prepared_data

def populate_molecules(file_path: str, batch_size: int = 500, limit: Optional[int] = None):
    """
    Populate the molecules table using direct connection for better performance.
    
    Args:
        file_path: Path to the JSON file containing molecule data
        batch_size: Number of molecules to insert in each batch
        limit: Optional limit on the total number of molecules to insert
    """
    # Verify database schema
    if not verify_database_schema():
        logger.error("Database schema verification failed. Aborting.")
        return
    
    # Load molecules from file
    molecules = load_molecules_from_file(file_path)
    if not molecules:
        logger.error("No molecules loaded. Aborting.")
        return
    
    # Apply limit if specified
    if limit is not None:
        molecules = molecules[:limit]
        logger.info(f"Limited to {limit} molecules")
    
    # Load checkpoint
    checkpoint = load_checkpoint()
    start_index = checkpoint['completed']
    
    logger.info(f"Starting molecule population from index {start_index} of {len(molecules)}")
    start_time = time.time()
    total_inserted = 0
    
    try:
        # Process in batches
        for i in range(start_index, len(molecules), batch_size):
            batch_start_time = time.time()
            
            # Prepare batch
            batch_molecules = prepare_molecule_batch(molecules, i, batch_size)
            batch_size_actual = len(batch_molecules)
            
            logger.info(f"Processing batch {i//batch_size + 1}/{(len(molecules) + batch_size - 1)//batch_size}: "
                      f"{batch_size_actual} molecules")
            
            try:
                # Insert batch
                bulk_insert('molecules', batch_molecules)
                
                # Update progress
                total_inserted += batch_size_actual
                checkpoint = {
                    'completed': i + batch_size_actual,
                    'last_id': batch_molecules[-1]['id'] if batch_molecules else None,
                    'timestamp': datetime.now().isoformat()
                }
                save_checkpoint(checkpoint)
                
                batch_time = time.time() - batch_start_time
                logger.info(f"Batch completed in {batch_time:.2f} seconds "
                          f"({batch_size_actual / batch_time:.2f} molecules/sec)")
                
            except Exception as e:
                logger.error(f"Error inserting batch: {e}")
                # Save checkpoint at the last successful batch
                save_checkpoint(checkpoint)
                raise
        
        total_time = time.time() - start_time
        logger.info(f"Population completed: {total_inserted} molecules inserted in {total_time:.2f} seconds "
                  f"({total_inserted / total_time:.2f} molecules/sec)")
        
    except Exception as e:
        logger.error(f"Error during population: {e}")
        raise
    finally:
        # Close all connections
        close_all()

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Populate molecules table using direct connection")
    parser.add_argument("--file", default="data/cryoprotectant_master_list.json", 
                        help="Path to JSON file with molecule data")
    parser.add_argument("--batch-size", type=int, default=500, 
                        help="Number of molecules to insert in each batch")
    parser.add_argument("--limit", type=int, default=None, 
                        help="Limit the total number of molecules to insert")
    
    args = parser.parse_args()
    
    try:
        # Initialize database connection
        db = PostgresDirectConnection.get_instance()
        if not db.test_connection():
            logger.error("Database connection failed. Aborting.")
            return 1
        
        # Run population
        populate_molecules(args.file, args.batch_size, args.limit)
        return 0
        
    except Exception as e:
        logger.error(f"Error in main: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())