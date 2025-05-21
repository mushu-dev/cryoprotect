#!/usr/bin/env python3
"""
Batch processor for molecules with None names.

This script identifies molecules with None names and updates them with
appropriate names based on various strategies:
1. Structural information (SMILES)
2. Data from external sources (PubChem, ChEMBL)
3. Standardized naming conventions

Usage:
    python batch_process_none_names.py [--batch-size=100] [--dry-run]
"""

import argparse
import sys
import logging
import time
import json
from typing import List, Dict, Any, Optional, Tuple
from datetime import datetime
import uuid

# Database connection
from database.adapter import get_database_connection
from database.utils import format_molecule_name

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"none_name_batch_process_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Constants
DEFAULT_BATCH_SIZE = 100
MAX_BATCH_SIZE = 1000
SLEEP_BETWEEN_BATCHES = 1  # seconds

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Batch process molecules with None names")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Number of molecules to process in each batch (default: {DEFAULT_BATCH_SIZE}, max: {MAX_BATCH_SIZE})"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run in dry-run mode without making any changes"
    )
    parser.add_argument(
        "--checkpoint-file",
        type=str,
        default=f"none_name_process_checkpoint_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
        help="File to store checkpoint information"
    )
    parser.add_argument(
        "--load-checkpoint",
        type=str,
        help="Load a previous checkpoint file to resume processing"
    )
    return parser.parse_args()

def get_molecules_with_none_names(db, batch_size: int, last_id: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    Get a batch of molecules with None names.
    
    Args:
        db: Database connection
        batch_size: Number of molecules to retrieve
        last_id: Last molecule ID processed (for pagination)
        
    Returns:
        List of molecules with None names
    """
    query = """
    SELECT id, smiles, pubchem_cid, chembl_id, formula, molecular_weight
    FROM molecules
    WHERE name IS NULL
    """
    
    params = []
    if last_id:
        query += " AND id > %s"
        params.append(last_id)
    
    query += " ORDER BY id LIMIT %s"
    params.append(batch_size)
    
    result = db.execute(query, params)
    return result.fetchall()

def get_name_from_smiles(smiles: str) -> Optional[str]:
    """
    Generate a name from SMILES using RDKit if available.
    
    Args:
        smiles: SMILES string for the molecule
        
    Returns:
        Generated name or None if not possible
    """
    if not smiles:
        return None
        
    try:
        # Try to import RDKit
        from rdkit import Chem
        from rdkit.Chem import AllChem, Draw
        
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
            
        # Try to get name from molecule properties
        try:
            return mol.GetProp("_Name")
        except:
            pass
            
        # Generate name from SMARTS pattern recognition
        try:
            from rdkit.Chem.Draw import SimilarityMaps
            return f"Compound-{Chem.MolToSmiles(mol, isomericSmiles=False)[:15]}"
        except:
            pass
            
        # Default to a simplified SMILES
        return f"Compound-{Chem.MolToSmiles(mol, isomericSmiles=False)[:15]}"
    except ImportError:
        logger.warning("RDKit not available, falling back to basic name generation")
        
    # Very basic name generation without RDKit
    if len(smiles) <= 20:
        return f"Compound-{smiles}"
    else:
        return f"Compound-{smiles[:20]}"

def get_name_from_external_id(pubchem_cid: Optional[str], chembl_id: Optional[str]) -> Optional[str]:
    """
    Get name from external database IDs.
    
    Args:
        pubchem_cid: PubChem CID
        chembl_id: ChEMBL ID
        
    Returns:
        Name from external database or None
    """
    if pubchem_cid:
        return f"PubChem-{pubchem_cid}"
    elif chembl_id:
        return f"ChEMBL-{chembl_id}"
    return None

def get_name_from_formula(formula: Optional[str], molecular_weight: Optional[float]) -> Optional[str]:
    """
    Generate a name from chemical formula and molecular weight.
    
    Args:
        formula: Chemical formula
        molecular_weight: Molecular weight
        
    Returns:
        Generated name or None
    """
    if formula:
        if molecular_weight:
            return f"{formula} (MW: {molecular_weight:.2f})"
        return f"Formula-{formula}"
    return None

def generate_name_for_molecule(molecule: Dict[str, Any]) -> str:
    """
    Generate a name for a molecule using various strategies.
    
    Args:
        molecule: Molecule data
        
    Returns:
        Generated name
    """
    # Try to get name from SMILES
    name = get_name_from_smiles(molecule.get('smiles'))
    if name:
        return name
        
    # Try to get name from external IDs
    name = get_name_from_external_id(molecule.get('pubchem_cid'), molecule.get('chembl_id'))
    if name:
        return name
        
    # Try to get name from formula and molecular weight
    name = get_name_from_formula(molecule.get('formula'), molecule.get('molecular_weight'))
    if name:
        return name
        
    # Last resort: generate a unique name
    return f"Molecule-{str(uuid.uuid4())[:8]}"

def update_molecule_name(db, molecule_id: str, name: str, dry_run: bool = False) -> bool:
    """
    Update the name of a molecule in the database.
    
    Args:
        db: Database connection
        molecule_id: Molecule ID
        name: New name
        dry_run: If True, don't actually update the database
        
    Returns:
        True if successful, False otherwise
    """
    if dry_run:
        logger.info(f"Would update molecule {molecule_id} with name: {name}")
        return True
        
    try:
        # Format the name according to standardization rules
        formatted_name = format_molecule_name(name)
        
        # Update the molecule name
        result = db.execute(
            "UPDATE molecules SET name = %s WHERE id = %s RETURNING id",
            (formatted_name, molecule_id)
        )
        
        return result.fetchone() is not None
    except Exception as e:
        logger.error(f"Error updating molecule {molecule_id}: {str(e)}")
        return False

def save_checkpoint(checkpoint_file: str, last_id: str, processed_count: int, success_count: int) -> None:
    """
    Save a checkpoint to resume processing later.
    
    Args:
        checkpoint_file: Path to checkpoint file
        last_id: Last molecule ID processed
        processed_count: Total number of molecules processed
        success_count: Number of molecules successfully updated
    """
    checkpoint_data = {
        'last_id': last_id,
        'processed_count': processed_count,
        'success_count': success_count,
        'timestamp': datetime.now().isoformat()
    }
    
    try:
        with open(checkpoint_file, 'w') as f:
            json.dump(checkpoint_data, f, indent=2)
        logger.info(f"Checkpoint saved to {checkpoint_file}")
    except Exception as e:
        logger.error(f"Error saving checkpoint: {str(e)}")

def load_checkpoint(checkpoint_file: str) -> Tuple[Optional[str], int, int]:
    """
    Load a checkpoint to resume processing.
    
    Args:
        checkpoint_file: Path to checkpoint file
        
    Returns:
        Tuple of (last_id, processed_count, success_count)
    """
    try:
        with open(checkpoint_file, 'r') as f:
            checkpoint_data = json.load(f)
            
        last_id = checkpoint_data.get('last_id')
        processed_count = checkpoint_data.get('processed_count', 0)
        success_count = checkpoint_data.get('success_count', 0)
        
        logger.info(f"Loaded checkpoint: last_id={last_id}, processed={processed_count}, success={success_count}")
        return last_id, processed_count, success_count
    except Exception as e:
        logger.error(f"Error loading checkpoint: {str(e)}")
        return None, 0, 0

def process_molecules(db, batch_size: int, dry_run: bool, checkpoint_file: str, resume_from: Optional[str] = None) -> None:
    """
    Process molecules with None names in batches.
    
    Args:
        db: Database connection
        batch_size: Number of molecules to process in each batch
        dry_run: If True, don't actually update the database
        checkpoint_file: Path to checkpoint file
        resume_from: Last molecule ID to resume from
    """
    last_id = resume_from
    total_processed = 0
    total_success = 0
    
    if resume_from:
        _, processed_count, success_count = load_checkpoint(checkpoint_file)
        total_processed = processed_count
        total_success = success_count
    
    start_time = time.time()
    
    while True:
        # Get a batch of molecules
        molecules = get_molecules_with_none_names(db, batch_size, last_id)
        
        if not molecules:
            logger.info("No more molecules with None names found")
            break
            
        batch_size = len(molecules)
        logger.info(f"Processing batch of {batch_size} molecules")
        
        batch_success = 0
        
        for molecule in molecules:
            molecule_id = molecule['id']
            
            # Generate a name for the molecule
            name = generate_name_for_molecule(molecule)
            
            # Update the molecule name
            success = update_molecule_name(db, molecule_id, name, dry_run)
            
            if success:
                batch_success += 1
                total_success += 1
                
            # Update the last_id for pagination
            last_id = molecule_id
            
        total_processed += batch_size
        
        # Log batch results
        logger.info(f"Batch completed: {batch_success}/{batch_size} successful")
        logger.info(f"Progress: {total_processed} processed, {total_success} successful")
        
        # Save checkpoint
        save_checkpoint(checkpoint_file, last_id, total_processed, total_success)
        
        # Sleep between batches to avoid database overload
        time.sleep(SLEEP_BETWEEN_BATCHES)
    
    # Log final results
    elapsed_time = time.time() - start_time
    logger.info(f"Processing completed in {elapsed_time:.2f} seconds")
    logger.info(f"Total: {total_processed} processed, {total_success} successful")

def count_molecules_with_none_names(db) -> int:
    """
    Count the number of molecules with None names.
    
    Args:
        db: Database connection
        
    Returns:
        Number of molecules with None names
    """
    result = db.execute("SELECT COUNT(*) FROM molecules WHERE name IS NULL")
    return result.fetchone()[0]

def main():
    """Main entry point."""
    args = parse_args()
    
    # Validate batch size
    batch_size = min(args.batch_size, MAX_BATCH_SIZE)
    if batch_size <= 0:
        batch_size = DEFAULT_BATCH_SIZE
        
    logger.info(f"Starting batch processor with batch size: {batch_size}")
    logger.info(f"Dry run: {args.dry_run}")
    
    # Get database connection
    db = get_database_connection()
    
    # Count molecules with None names
    count = count_molecules_with_none_names(db)
    logger.info(f"Found {count} molecules with None names")
    
    if count == 0:
        logger.info("No molecules with None names found")
        return 0
        
    # Load checkpoint if specified
    resume_from = None
    if args.load_checkpoint:
        resume_from, _, _ = load_checkpoint(args.load_checkpoint)
        checkpoint_file = args.load_checkpoint
    else:
        checkpoint_file = args.checkpoint_file
    
    # Process molecules
    process_molecules(db, batch_size, args.dry_run, checkpoint_file, resume_from)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())