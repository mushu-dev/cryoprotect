#!/usr/bin/env python3
"""
Enhance PubChem molecules with missing properties.

This script identifies PubChem molecules with missing properties and retrieves
the missing data from the PubChem API. Uses the connection factory for
improved performance and reliability.
"""

import os
import sys
import logging
import argparse
import json
import time
import requests
from typing import Dict, List, Any, Optional, Tuple, Union
from datetime import datetime
import concurrent.futures
import uuid

# Import database connection utilities
from database.connection import get_db_connection, close_all_db_connections
from sql_executor import (
    execute_query, bulk_insert, execute_batch, with_retry,
    process_in_batches, get_db, with_transaction
)
from property_utils import PropertyManager

# Import PubChem utilities
from pubchem.simple_client import PubChemClient
from pubchem.simple_rate_limiter import RateLimiter

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
BATCH_SIZE = 10  # Number of compounds to process in parallel
MAX_RETRIES = 3  # Maximum number of retries for failed operations
DRY_RUN = False  # Set to True to run without database updates
CHECKPOINT_FILE = "checkpoints/pubchem_enhancement.json"  # Checkpoint file for resumable operations
PROPERTY_MAPPING = {
    'XLogP': 'logP',
    'HBondDonorCount': 'h_bond_donors',
    'HBondAcceptorCount': 'h_bond_acceptors',
    'RotatableBondCount': 'rotatable_bonds',
    'HeavyAtomCount': 'heavy_atoms',
    'TPSA': 'tpsa',
    'MolecularWeight': 'molecular_weight',
    'MolecularFormula': 'molecular_formula'
}

class PubChemEnhancer:
    """
    Class for enhancing PubChem molecules with missing properties.
    Implements batch processing, checkpointing, and resilient API operations.
    """
    
    def __init__(self, batch_size: int = BATCH_SIZE, dry_run: bool = DRY_RUN, 
                checkpoint_file: str = CHECKPOINT_FILE):
        """
        Initialize the PubChem enhancer.
        
        Args:
            batch_size: Number of compounds to process in parallel
            dry_run: Whether to run without database updates
            checkpoint_file: Path to the checkpoint file for resumable operations
        """
        self.batch_size = batch_size
        self.dry_run = dry_run
        self.checkpoint_file = checkpoint_file
        
        # Initialize clients and utilities
        self.pubchem_client = PubChemClient()
        self.rate_limiter = RateLimiter(requests_per_minute=5)
        self.property_manager = PropertyManager()
        
        # Initialize results tracking
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'total_molecules': 0,
            'molecules_enhanced': 0,
            'molecules_skipped': 0,
            'molecules_failed': 0,
            'details': {}
        }
        
        # Load checkpoint if exists
        self.processed_molecules = self._load_checkpoint()
        
    def _load_checkpoint(self) -> set:
        """
        Load checkpoint data from file.
        
        Returns:
            Set of molecule IDs that have already been processed
        """
        processed = set()
        
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, 'r') as f:
                    checkpoint_data = json.load(f)
                    processed = set(checkpoint_data.get('processed_molecules', []))
                    
                logger.info(f"Loaded checkpoint with {len(processed)} processed molecules")
            except Exception as e:
                logger.warning(f"Failed to load checkpoint file: {str(e)}")
                
        return processed
        
    def _save_checkpoint(self):
        """Save checkpoint data to file."""
        try:
            os.makedirs(os.path.dirname(self.checkpoint_file), exist_ok=True)
            
            checkpoint_data = {
                'timestamp': datetime.now().isoformat(),
                'processed_molecules': list(self.processed_molecules),
                'progress': {
                    'total': self.results['total_molecules'],
                    'enhanced': self.results['molecules_enhanced'],
                    'skipped': self.results['molecules_skipped'],
                    'failed': self.results['molecules_failed']
                }
            }
            
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
                
            logger.debug(f"Saved checkpoint with {len(self.processed_molecules)} processed molecules")
        except Exception as e:
            logger.error(f"Failed to save checkpoint: {str(e)}")
    
    def find_molecules_to_enhance(self) -> List[Tuple]:
        """
        Find PubChem molecules with missing properties.
        
        Returns:
            List of molecule data tuples (id, name, pubchem_cid)
        """
        if self.dry_run:
            # Use sample data for dry run
            logger.info("Running in dry run mode with sample data")
            molecules_to_enhance = [
                ("CRYO0001", "Glycerol", "753"),
                ("CRYO0002", "DMSO", "679"),
                ("CRYO0003", "Ethylene glycol", "174"),
                ("CRYO0004", "Propylene glycol", "1030"),
                ("CRYO0005", "Trehalose", "7427")
            ]
        else:
            # Query database for molecules with missing properties in the normalized schema
            # Use the direct PostgreSQL connection
            query = """
                SELECT m.id, m.name, m.pubchem_cid
                FROM molecules m
                WHERE m.pubchem_cid IS NOT NULL
                AND (
                    NOT EXISTS (
                        SELECT 1 FROM molecular_properties mp
                        JOIN property_types pt ON mp.property_type_id = pt.id
                        WHERE mp.molecule_id = m.id AND pt.name = 'logP'
                    ) OR
                    NOT EXISTS (
                        SELECT 1 FROM molecular_properties mp
                        JOIN property_types pt ON mp.property_type_id = pt.id
                        WHERE mp.molecule_id = m.id AND pt.name = 'h_bond_donors'
                    ) OR
                    NOT EXISTS (
                        SELECT 1 FROM molecular_properties mp
                        JOIN property_types pt ON mp.property_type_id = pt.id
                        WHERE mp.molecule_id = m.id AND pt.name = 'h_bond_acceptors'
                    )
                )
            """
            
            # Filter out already processed molecules if checkpoint exists
            if self.processed_molecules:
                placeholders = ', '.join(['%s'] * len(self.processed_molecules))
                query += f" AND m.id NOT IN ({placeholders})"
                molecules_to_enhance = execute_query(query, list(self.processed_molecules))
            else:
                molecules_to_enhance = execute_query(query)
                
            # Log the data structure type for debugging
            if molecules_to_enhance and len(molecules_to_enhance) > 0:
                logger.debug(f"Molecule data structure type: {type(molecules_to_enhance[0])}")
                if isinstance(molecules_to_enhance[0], dict):
                    logger.debug(f"Molecule keys: {list(molecules_to_enhance[0].keys())}")
                logger.info(f"Found {len(molecules_to_enhance)} molecules to enhance")
            else:
                logger.info("No molecules found to enhance")
                molecules_to_enhance = []
                
        return molecules_to_enhance
    
    def _fetch_pubchem_properties(self, pubchem_cid: str) -> Dict:
        """
        Fetch properties for a compound from PubChem API with retry logic.
        
        Args:
            pubchem_cid: PubChem compound ID
            
        Returns:
            Dictionary of property data or empty dict if failed
        """
        # Implement retry logic manually instead of using decorator
        attempt = 0
        current_delay = 2  # Initial delay in seconds
        
        while True:
            try:
                # Wait for rate limiter
                self.rate_limiter.wait()
                
                # Get compound data from PubChem
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/property/MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount/JSON"
                response = requests.get(url, timeout=30)
                
                if not response.ok:
                    logger.warning(f"Failed to fetch data for PubChem CID: {pubchem_cid}, Status: {response.status_code}")
                    raise Exception(f"API request failed with status {response.status_code}")
                    
                try:
                    data = response.json()
                    compound_data = data.get("PropertyTable", {}).get("Properties", [{}])[0]
                    
                    if not compound_data:
                        logger.warning(f"No property data found for PubChem CID: {pubchem_cid}")
                        return {}
                        
                    return compound_data
                except Exception as e:
                    logger.warning(f"Failed to parse data for PubChem CID: {pubchem_cid}, Error: {str(e)}")
                    raise
                
            except Exception as e:
                attempt += 1
                if attempt >= MAX_RETRIES:
                    logger.error(f"All {MAX_RETRIES} attempts failed for PubChem CID {pubchem_cid}: {str(e)}")
                    return {}
                
                logger.warning(f"Attempt {attempt}/{MAX_RETRIES} failed for PubChem CID {pubchem_cid}: {str(e)}")
                logger.info(f"Retrying in {current_delay} seconds...")
                
                time.sleep(current_delay)
                current_delay *= 1.5  # Exponential backoff
        """
        Fetch properties for a compound from PubChem API with retry logic.
        
        Args:
            pubchem_cid: PubChem compound ID
            
        Returns:
            Dictionary of property data or empty dict if failed
        """
        # Wait for rate limiter
        self.rate_limiter.wait()
        
        # Get compound data from PubChem
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/property/MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount/JSON"
        response = requests.get(url, timeout=30)
        
        if not response.ok:
            logger.warning(f"Failed to fetch data for PubChem CID: {pubchem_cid}, Status: {response.status_code}")
            raise Exception(f"API request failed with status {response.status_code}")
            
        try:
            data = response.json()
            compound_data = data.get("PropertyTable", {}).get("Properties", [{}])[0]
            
            if not compound_data:
                logger.warning(f"No property data found for PubChem CID: {pubchem_cid}")
                return {}
                
            return compound_data
        except Exception as e:
            logger.warning(f"Failed to parse data for PubChem CID: {pubchem_cid}, Error: {str(e)}")
            raise
    
    def _enhance_molecule(self, molecule: Union[Dict, Tuple]) -> Tuple[bool, bool]:
        """
        Enhance a single molecule with missing properties.
        
        Args:
            molecule: Molecule data (dictionary with keys 'id', 'name', 'pubchem_cid' or tuple)
            
        Returns:
            Tuple of (success, is_enhanced)
        """
        # Handle both tuple and dictionary formats
        if isinstance(molecule, dict):
            molecule_id = molecule['id']
            molecule_name = molecule['name']
            pubchem_cid = molecule['pubchem_cid']
        elif isinstance(molecule, (list, tuple)) and len(molecule) >= 3:
            molecule_id = molecule[0]
            molecule_name = molecule[1]
            pubchem_cid = molecule[2]
        else:
            logger.error(f"Unsupported molecule data format: {type(molecule)}")
            return False, False
        
        logger.debug(f"Enhancing molecule {molecule_id} (CID {pubchem_cid})")
        
        # Skip if already processed (checkpoint)
        if molecule_id in self.processed_molecules:
            logger.debug(f"Skipping already processed molecule {molecule_id}")
            return True, False
        
        try:
            # Fetch properties from PubChem
            compound_data = self._fetch_pubchem_properties(pubchem_cid)
            
            if not compound_data:
                return False, False
            
            # Extract properties using the mapping
            properties = {}
            for pubchem_prop, our_prop in PROPERTY_MAPPING.items():
                if pubchem_prop in compound_data:
                    properties[our_prop] = compound_data[pubchem_prop]
            
            # If no properties extracted, skip
            if not properties:
                logger.warning(f"No properties found for PubChem CID: {pubchem_cid}")
                self.processed_molecules.add(molecule_id)
                return True, False
            
            if self.dry_run:
                # In dry run mode, just log the properties that would be updated
                logger.info(f"[DRY RUN] Would update molecule {molecule_id} with properties: {list(properties.keys())}")
                self.processed_molecules.add(molecule_id)
                return True, True
            else:
                # Use PropertyManager to update properties in the normalized schema
                success_count, total_count = self.property_manager.set_properties(
                    molecule_id, properties, created_by=None
                )
                
                # Update the molecule's updated_at timestamp
                execute_query("""
                    UPDATE molecules
                    SET updated_at = NOW()
                    WHERE id = %s
                """, (molecule_id,))
                
                logger.debug(f"Enhanced molecule {molecule_id} with {success_count}/{total_count} properties: {list(properties.keys())}")
                
                # Mark as processed for checkpoint
                self.processed_molecules.add(molecule_id)
                
                return True, True
                
        except Exception as e:
            logger.error(f"Error enhancing molecule {molecule_id}: {str(e)}")
            return False, False
    
    def process_batch(self, batch: List[Union[Dict, Tuple]]):
        """
        Process a batch of molecules in parallel.
        
        Args:
            batch: List of molecule data (dictionaries or tuples)
        """
        if not batch:
            logger.info("Empty batch received, skipping processing")
            return
            
        logger.info(f"Processing batch of {len(batch)} molecules")
        
        # Process batch in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.batch_size) as executor:
            future_to_molecule = {
                executor.submit(self._enhance_molecule, molecule): molecule
                for molecule in batch
            }
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(future_to_molecule):
                molecule = future_to_molecule[future]
                try:
                    success, is_enhanced = future.result()
                    
                    if success:
                        if is_enhanced:
                            self.results['molecules_enhanced'] += 1
                            # Handle both tuple and dictionary formats
                            if isinstance(molecule, dict):
                                molecule_id = molecule['id']
                                molecule_name = molecule['name']
                                pubchem_cid = molecule['pubchem_cid']
                            elif isinstance(molecule, (list, tuple)) and len(molecule) >= 3:
                                molecule_id = molecule[0]
                                molecule_name = molecule[1]
                                pubchem_cid = molecule[2]
                            else:
                                logger.error(f"Unsupported molecule data format: {type(molecule)}")
                                continue
                                
                            self.results['details'][molecule_id] = {
                                'status': 'enhanced',
                                'name': molecule_name,
                                'pubchem_cid': pubchem_cid
                            }
                        else:
                            self.results['molecules_skipped'] += 1
                            # Handle both tuple and dictionary formats
                            if isinstance(molecule, dict):
                                molecule_id = molecule['id']
                                molecule_name = molecule['name']
                                pubchem_cid = molecule['pubchem_cid']
                            else:
                                molecule_id = molecule[0]
                                molecule_name = molecule[1]
                                pubchem_cid = molecule[2]
                                
                            self.results['details'][molecule_id] = {
                                'status': 'skipped',
                                'name': molecule_name,
                                'pubchem_cid': pubchem_cid
                            }
                    else:
                        self.results['molecules_failed'] += 1
                        # Handle both tuple and dictionary formats
                        if isinstance(molecule, dict):
                            molecule_id = molecule['id']
                            molecule_name = molecule['name']
                            pubchem_cid = molecule['pubchem_cid']
                        else:
                            molecule_id = molecule[0]
                            molecule_name = molecule[1]
                            pubchem_cid = molecule[2]
                            
                        self.results['details'][molecule_id] = {
                            'status': 'failed',
                            'name': molecule_name,
                            'pubchem_cid': pubchem_cid
                        }
                        
                except Exception as e:
                    # Handle both tuple and dictionary formats
                    if isinstance(molecule, dict):
                        molecule_id = molecule['id']
                        molecule_name = molecule['name']
                        pubchem_cid = molecule['pubchem_cid']
                    elif isinstance(molecule, (list, tuple)) and len(molecule) >= 3:
                        molecule_id = molecule[0]
                        molecule_name = molecule[1]
                        pubchem_cid = molecule[2]
                    else:
                        logger.error(f"Unsupported molecule data format: {type(molecule)}")
                        self.results['molecules_failed'] += 1
                        continue
                        
                    logger.error(f"Error processing molecule {molecule_id}: {str(e)}")
                    self.results['molecules_failed'] += 1
                    self.results['details'][molecule_id] = {
                        'status': 'failed',
                        'name': molecule_name,
                        'pubchem_cid': pubchem_cid,
                        'error': str(e)
                    }
        
        # Save checkpoint after each batch
        self._save_checkpoint()
    
    def enhance_properties(self, output_report: Optional[str] = None) -> Dict:
        """
        Enhance PubChem molecules with missing properties.
        
        Args:
            output_report: Optional path to save the enhancement report
            
        Returns:
            Enhancement results dictionary
        """
        # Find molecules to enhance
        molecules_to_enhance = self.find_molecules_to_enhance()
        self.results['total_molecules'] = len(molecules_to_enhance)
        
        logger.info(f"Found {self.results['total_molecules']} PubChem molecules with missing properties")
        
        if not molecules_to_enhance:
            logger.info("No molecules to enhance, exiting")
            return self.results
        
        # Process molecules in batches
        process_in_batches(
            molecules_to_enhance,
            batch_size=self.batch_size,
            process_func=self.process_batch
        )
        
        # Generate report
        if output_report:
            os.makedirs(os.path.dirname(output_report), exist_ok=True)
            with open(output_report, 'w') as f:
                json.dump(self.results, f, indent=2)
            logger.info(f"Enhancement report saved to {output_report}")
        
        logger.info(f"Enhancement completed: {self.results['molecules_enhanced']} molecules enhanced, " +
                   f"{self.results['molecules_skipped']} molecules skipped, " +
                   f"{self.results['molecules_failed']} molecules failed")
        
        return self.results

def enhance_pubchem_properties(output_report: Optional[str] = None, batch_size: int = BATCH_SIZE, 
                              dry_run: bool = DRY_RUN, checkpoint_file: str = CHECKPOINT_FILE) -> Dict:
    """
    Enhance PubChem molecules with missing properties.
    
    Args:
        output_report: Optional path to save the enhancement report
        batch_size: Number of compounds to process in parallel
        dry_run: Whether to run without database updates
        checkpoint_file: Path to the checkpoint file for resumable operations
        
    Returns:
        Enhancement results dictionary
    """
    enhancer = PubChemEnhancer(batch_size, dry_run, checkpoint_file)
    return enhancer.enhance_properties(output_report)

def main():
    """CLI entry point for property enhancement."""
    parser = argparse.ArgumentParser(description='Enhance PubChem molecules with missing properties')
    parser.add_argument('--batch-size', type=int, default=BATCH_SIZE,
                      help='Number of molecules to process in parallel')
    parser.add_argument('--report',
                      default=f"reports/pubchem_enhancement_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                      help='Output file path for report')
    parser.add_argument('--dry-run', action='store_true',
                      help='Run in dry run mode without database updates')
    parser.add_argument('--checkpoint',
                      default=CHECKPOINT_FILE,
                      help='Checkpoint file for resumable operations')
    parser.add_argument('--reset-checkpoint', action='store_true',
                      help='Reset checkpoint and process all molecules')
    
    args = parser.parse_args()
    
    # Reset checkpoint if requested
    if args.reset_checkpoint and os.path.exists(args.checkpoint):
        os.remove(args.checkpoint)
        logger.info(f"Reset checkpoint file: {args.checkpoint}")
    
    try:
        enhance_pubchem_properties(args.report, args.batch_size, args.dry_run, args.checkpoint)
        return 0
    except Exception as e:
        logger.error(f"Enhancement failed: {str(e)}")
        return 1
    finally:
        # Ensure database connections are properly closed
        from database.connection import close_all_db_connections
        close_all_db_connections()

if __name__ == "__main__":
    sys.exit(main())