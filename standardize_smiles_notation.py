#!/usr/bin/env python3
"""
This script standardizes SMILES notation for all molecules in the database
using the rdkit library to ensure consistent representation.
"""

import os
import sys
import logging
import psycopg2
import psycopg2.extras
from dotenv import load_dotenv
from datetime import datetime
import json
import argparse

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Using mock implementation.")
    # Mock implementation for environments without RDKit
    class MockChem:
        @staticmethod
        def MolFromSmiles(smiles):
            if not smiles or not isinstance(smiles, str):
                return None
            return {'smiles': smiles}
            
        @staticmethod
        def MolToSmiles(mol, **kwargs):
            if not mol:
                return None
            return mol.get('smiles', '')
            
    class MockAllChem:
        pass
        
    Chem = MockChem()
    AllChem = MockAllChem()

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("smiles_standardization.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Database connection parameters
DB_HOST = os.getenv("SUPABASE_DB_HOST")
DB_PORT = os.getenv("SUPABASE_DB_PORT")
DB_NAME = os.getenv("SUPABASE_DB_NAME")
DB_USER = os.getenv("SUPABASE_DB_USER")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD")

def get_db_connection():
    """Get a direct database connection using psycopg2."""
    try:
        conn = psycopg2.connect(
            host=DB_HOST,
            port=DB_PORT,
            dbname=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD
        )
        logger.info("Connected to database using direct PostgreSQL connection")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def standardize_smiles(smiles):
    """
    Standardize SMILES notation using RDKit.
    
    Args:
        smiles (str): Original SMILES string
        
    Returns:
        str: Standardized SMILES string or None if invalid
    """
    if not smiles:
        return None
        
    try:
        # Parse SMILES string to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None
            
        # Generate canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        
        if not canonical_smiles:
            logger.warning(f"Failed to generate canonical SMILES for: {smiles}")
            return None
            
        return canonical_smiles
        
    except Exception as e:
        logger.error(f"Error standardizing SMILES {smiles}: {e}")
        return None

def get_molecules_with_non_standard_smiles(batch_size=500, offset=0):
    """
    Get molecules with SMILES that need standardization.
    Process in batches to avoid memory issues.
    
    Args:
        batch_size: Number of molecules to process in each batch
        offset: Starting offset for pagination
        
    Returns:
        List of molecule records
    """
    conn = get_db_connection()
    if not conn:
        return []
        
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    molecules = []
    
    try:
        # Get molecules with non-null SMILES
        cursor.execute("""
            SELECT id, name, smiles
            FROM molecules
            WHERE smiles IS NOT NULL
            ORDER BY id
            LIMIT %s OFFSET %s
        """, (batch_size, offset))
        
        molecules = cursor.fetchall()
        logger.info(f"Retrieved {len(molecules)} molecules (offset {offset})")
        
    except Exception as e:
        logger.error(f"Error retrieving molecules: {e}")
        conn.rollback()
    
    cursor.close()
    conn.close()
    return molecules

def update_molecule_smiles(molecule_id, original_smiles, standardized_smiles, dry_run=False):
    """
    Update a molecule's SMILES with the standardized version.
    
    Args:
        molecule_id: The molecule ID
        original_smiles: Original SMILES string
        standardized_smiles: Standardized SMILES string
        dry_run: If True, just log changes without making them
        
    Returns:
        bool: True if successful, False otherwise
    """
    if original_smiles == standardized_smiles:
        return True  # No change needed
        
    conn = get_db_connection()
    if not conn:
        return False
        
    cursor = conn.cursor()
    success = False
    
    try:
        if dry_run:
            logger.info(f"DRY RUN: Would update molecule {molecule_id} SMILES from '{original_smiles}' to '{standardized_smiles}'")
            success = True
        else:
            # Start a transaction
            conn.autocommit = False
            
            # Update the molecule record
            cursor.execute("""
                UPDATE molecules
                SET smiles = %s, updated_at = NOW()
                WHERE id = %s
            """, (standardized_smiles, molecule_id))
            
            # Create an audit record
            cursor.execute("""
                INSERT INTO scientific_data_audit 
                (table_name, record_id, operation, old_value, new_value, timestamp)
                VALUES 
                ('molecules', %s, 'standardize_smiles', %s, %s, %s)
            """, (
                str(molecule_id),
                json.dumps({"smiles": original_smiles}),
                json.dumps({"smiles": standardized_smiles}),
                datetime.now()
            ))
            
            conn.commit()
            logger.info(f"Updated molecule {molecule_id} SMILES: '{original_smiles}' → '{standardized_smiles}'")
            success = True
            
    except Exception as e:
        conn.rollback()
        logger.error(f"Error updating molecule {molecule_id}: {e}")
    
    cursor.close()
    conn.close()
    return success

def process_batch(molecules, dry_run=False):
    """
    Process a batch of molecules for SMILES standardization.
    
    Args:
        molecules: List of molecule records
        dry_run: If True, don't make actual changes
        
    Returns:
        dict: Statistics about the processing
    """
    stats = {
        "processed": len(molecules),
        "standardized": 0,
        "unchanged": 0,
        "errors": 0,
        "invalid": 0
    }
    
    for molecule in molecules:
        molecule_id = molecule['id']
        original_smiles = molecule['smiles']
        
        # Standardize the SMILES
        standardized_smiles = standardize_smiles(original_smiles)
        
        if standardized_smiles is None:
            logger.warning(f"Invalid SMILES for molecule {molecule_id}: {original_smiles}")
            stats["invalid"] += 1
            continue
            
        if original_smiles == standardized_smiles:
            logger.debug(f"SMILES already standardized for molecule {molecule_id}")
            stats["unchanged"] += 1
            continue
            
        # Update the molecule with standardized SMILES
        if update_molecule_smiles(molecule_id, original_smiles, standardized_smiles, dry_run):
            stats["standardized"] += 1
        else:
            stats["errors"] += 1
            
    return stats

def main():
    """Main function to run SMILES standardization."""
    parser = argparse.ArgumentParser(description="Standardize SMILES notation for molecules.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing them")
    parser.add_argument("--batch-size", type=int, default=500, help="Number of molecules to process in each batch")
    parser.add_argument("--max-molecules", type=int, default=None, help="Maximum number of molecules to process (default: all)")
    args = parser.parse_args()
    
    if not RDKIT_AVAILABLE and not args.dry_run:
        logger.error("RDKit is required for SMILES standardization. Please install RDKit and try again.")
        logger.error("You can run with --dry-run to test the workflow without RDKit.")
        return False
    
    logger.info(f"Starting SMILES standardization {'in dry-run mode' if args.dry_run else 'in execution mode'}")
    
    total_stats = {
        "processed": 0,
        "standardized": 0,
        "unchanged": 0,
        "errors": 0,
        "invalid": 0
    }
    
    offset = 0
    batch_size = args.batch_size
    
    # Process molecules in batches
    while True:
        molecules = get_molecules_with_non_standard_smiles(batch_size, offset)
        
        if not molecules:
            logger.info("No more molecules to process")
            break
            
        batch_stats = process_batch(molecules, args.dry_run)
        
        # Update total statistics
        for key in total_stats:
            total_stats[key] += batch_stats[key]
            
        logger.info(f"Batch complete: {batch_stats['processed']} processed, "
                   f"{batch_stats['standardized']} standardized, "
                   f"{batch_stats['unchanged']} unchanged, "
                   f"{batch_stats['errors']} errors, "
                   f"{batch_stats['invalid']} invalid")
        
        offset += batch_size
        
        # Check if we've reached the maximum number of molecules to process
        if args.max_molecules and total_stats["processed"] >= args.max_molecules:
            logger.info(f"Reached maximum of {args.max_molecules} molecules")
            break
    
    # Print summary
    print("\n" + "=" * 60)
    print("SMILES Standardization Summary")
    print("=" * 60)
    print(f"Total molecules processed: {total_stats['processed']}")
    print(f"Standardized: {total_stats['standardized']}")
    print(f"Already standardized: {total_stats['unchanged']}")
    print(f"Invalid SMILES: {total_stats['invalid']}")
    print(f"Errors during processing: {total_stats['errors']}")
    print("=" * 60)
    
    logger.info("SMILES standardization complete")
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)