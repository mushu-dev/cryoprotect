# PubChem_CryoProtectants_Supabase_Enhanced.py
import os
import json
import logging
import time
import argparse
import sys
import concurrent.futures
import threading
from typing import Dict, List, Any, Optional, Set, Tuple
from datetime import datetime

# Import custom modules
from connection_pool_wrapper import ConnectionManager
from cryoprotectant_identifiers import CryoprotectantIdentifierManager
from pubchem.simple_client import PubChemClient
from pubchem.simple_rate_limiter import RateLimiter
from pubchem.cache import PubChemCache

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('pubchem_import.log')
    ]
)
logger = logging.getLogger(__name__)

# Constants
CHECKPOINT_DIR = "checkpoints"
CHECKPOINT_FILE = os.path.join(CHECKPOINT_DIR, "pubchem_import_enhanced.json")
BATCH_SIZE = 25  # Number of compounds to process in parallel
MAX_RETRIES = 3  # Maximum number of retries for failed operations

class PubChemImporter:
    """
    Enhanced PubChem importer with batch processing, robust error handling,
    and integration with connection pool and identifier manager.
    """
    
    def __init__(self, 
                checkpoint_file: str = CHECKPOINT_FILE,
                batch_size: int = BATCH_SIZE,
                max_retries: int = MAX_RETRIES):
        """
        Initialize the PubChem importer.
        
        Args:
            checkpoint_file: Path to checkpoint file
            batch_size: Number of compounds to process in parallel
            max_retries: Maximum number of retries for failed operations
        """
        self.checkpoint_file = checkpoint_file
        self.batch_size = batch_size
        self.max_retries = max_retries
        
        # Initialize components
        self.pubchem_client = PubChemClient()
        self.pubchem_cache = PubChemCache()
        self.rate_limiter = RateLimiter(requests_per_minute=5)
        self.id_manager = CryoprotectantIdentifierManager.get_instance()
        
        # Import state
        self.checkpoint = self._load_checkpoint()
        self.next_internal_id = self.checkpoint.get('next_internal_id', 1)
        
        # Statistics
        self.stats = {
            'compounds_processed': 0,
            'compounds_added': 0,
            'compounds_updated': 0,
            'compounds_skipped': 0,
            'compounds_failed': 0,
            'start_time': None,
            'end_time': None,
            'total_elapsed_time': 0
        }
    
    def _load_checkpoint(self) -> Dict[str, Any]:
        """
        Load checkpoint from file.
        
        Returns:
            Dictionary containing checkpoint data
        """
        if not os.path.exists(self.checkpoint_file):
            logger.info("No checkpoint file found at %s. Starting fresh import.", 
                       self.checkpoint_file)
            return {
                'processed_cids': [],
                'failed_cids': {},
                'next_internal_id': 1,
                'last_completed_batch': -1,
                'total_compounds': 0
            }
        
        try:
            with open(self.checkpoint_file, 'r') as f:
                checkpoint = json.load(f)
            
            logger.info("Loaded checkpoint with %d processed compounds, %d failed compounds",
                       len(checkpoint.get('processed_cids', [])), 
                       len(checkpoint.get('failed_cids', {})))
            
            return checkpoint
        except Exception as e:
            logger.error("Failed to load checkpoint: %s", str(e))
            return {
                'processed_cids': [],
                'failed_cids': {},
                'next_internal_id': 1,
                'last_completed_batch': -1,
                'total_compounds': 0
            }
    
    def _save_checkpoint(self) -> None:
        """Save checkpoint to file."""
        try:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.checkpoint_file), exist_ok=True)
            
            with open(self.checkpoint_file, 'w') as f:
                json.dump(self.checkpoint, f, indent=2)
            
            logger.debug("Saved checkpoint with %d processed compounds, %d failed compounds",
                        len(self.checkpoint.get('processed_cids', [])), 
                        len(self.checkpoint.get('failed_cids', {})))
        except Exception as e:
            logger.error("Failed to save checkpoint: %s", str(e))
    
    def _generate_internal_id(self) -> str:
        """
        Generate a new internal ID for a molecule.
        
        Returns:
            Internal ID string (e.g., CRYO0123)
        """
        internal_id = f"CRYO{self.next_internal_id:04d}"
        self.next_internal_id += 1
        self.checkpoint['next_internal_id'] = self.next_internal_id
        return internal_id
    
    def _get_or_create_internal_id(self, 
                                  pubchem_data: Dict[str, Any]) -> Tuple[str, bool]:
        """
        Get existing internal ID or create a new one.
        
        Args:
            pubchem_data: PubChem compound data
            
        Returns:
            Tuple of (internal_id, is_new)
        """
        pubchem_cid = str(pubchem_data.get('id', {}).get('id', {}).get('cid', ''))
        name = pubchem_data.get('name', '')
        inchi_key = pubchem_data.get('inchi_key', '')
        
        # Try to resolve with identifier manager
        internal_id, confidence = self.id_manager.resolve_identifier(
            pubchem_cid=pubchem_cid,
            name=name,
            inchi_key=inchi_key
        )
        
        if internal_id:
            return internal_id, False
        
        # Create new internal ID
        internal_id = self._generate_internal_id()
        
        # Add to identifier manager
        synonyms = []
        if name:
            synonyms.append(name)
        
        for synonym in pubchem_data.get('synonyms', [])[:5]:  # Limit to 5 synonyms
            if synonym not in synonyms:
                synonyms.append(synonym)
        
        self.id_manager.add_molecule(internal_id, {
            "internal_id": internal_id,
            "pubchem_cid": pubchem_cid,
            "names": synonyms,
            "inchi_key": inchi_key,
            "smiles": pubchem_data.get('smiles', ''),
            "formula": pubchem_data.get('molecular_formula', ''),
            "category": "imported"
        })
        
        # Save identifier list
        self.id_manager.save_identifiers()
        
        return internal_id, True
    
    def _fetch_compound_data(self, cid: str) -> Optional[Dict[str, Any]]:
        """
        Fetch compound data from PubChem.
        
        Args:
            cid: PubChem CID
            
        Returns:
            Dictionary containing compound data or None if failed
        """
        retry_count = 0
        while retry_count < self.max_retries:
            try:
                self.rate_limiter.wait()
                
                # Try to get from cache first
                cached_data = self.pubchem_cache.get_compound(cid)
                if cached_data:
                    logger.debug("Using cached data for CID %s", cid)
                    return cached_data
                
                # Fetch from PubChem
                logger.debug("Fetching data for CID %s from PubChem", cid)
                compound_data = self.pubchem_client.get_compound(cid)
                
                # Cache the result
                self.pubchem_cache.cache_compound(cid, compound_data)
                
                return compound_data
            
            except Exception as e:
                retry_count += 1
                logger.warning("Failed to fetch CID %s (attempt %d/%d): %s", 
                              cid, retry_count, self.max_retries, str(e))
                time.sleep(2 ** retry_count)  # Exponential backoff
        
        logger.error("Failed to fetch CID %s after %d attempts", cid, self.max_retries)
        return None
    
    def _process_compound(self, cid: str) -> bool:
        """
        Process a single compound.
        
        Args:
            cid: PubChem CID
            
        Returns:
            True if successful, False otherwise
        """
        if cid in self.checkpoint['processed_cids']:
            logger.debug("Skipping CID %s (already processed)", cid)
            self.stats['compounds_skipped'] += 1
            return True
        
        try:
            # Fetch compound data
            compound_data = self._fetch_compound_data(cid)
            if not compound_data:
                self.checkpoint['failed_cids'][cid] = "Failed to fetch data"
                self.stats['compounds_failed'] += 1
                return False
            
            # Get or create internal ID
            internal_id, is_new = self._get_or_create_internal_id(compound_data)
            
            # Insert or update in database
            self._insert_or_update_compound(internal_id, cid, compound_data, is_new)
            
            # Mark as processed
            if cid not in self.checkpoint['processed_cids']:
                self.checkpoint['processed_cids'].append(cid)
            
            if is_new:
                self.stats['compounds_added'] += 1
            else:
                self.stats['compounds_updated'] += 1
            
            self.stats['compounds_processed'] += 1
            
            return True
            
        except Exception as e:
            logger.error("Failed to process CID %s: %s", cid, str(e))
            self.checkpoint['failed_cids'][cid] = str(e)
            self.stats['compounds_failed'] += 1
            return False
    
    def _insert_or_update_compound(self, 
                                  internal_id: str, 
                                  cid: str, 
                                  compound_data: Dict[str, Any],
                                  is_new: bool) -> None:
        """
        Insert or update compound in database.
        
        Args:
            internal_id: Internal molecule ID
            cid: PubChem CID
            compound_data: PubChem compound data
            is_new: Whether this is a new compound
        """
        retry_count = 0
        
        while retry_count < self.max_retries:
            try:
                # Extract necessary data
                name = compound_data.get('name', '')
                synonyms = compound_data.get('synonyms', [])
                inchi = compound_data.get('inchi', '')
                inchi_key = compound_data.get('inchi_key', '')
                smiles = compound_data.get('smiles', '')
                molecular_formula = compound_data.get('molecular_formula', '')
                molecular_weight = compound_data.get('molecular_weight', 0.0)
                
                # Convert properties to JSON
                properties = {
                    'pubchem': {
                        'logs': {
                            'imported_at': datetime.now().isoformat(),
                            'source': 'PubChem',
                            'version': '2.0'
                        },
                        'basic': {
                            'name': name,
                            'synonyms': synonyms[:10],  # Limit to 10 synonyms
                            'formula': molecular_formula,
                            'molecular_weight': molecular_weight
                        },
                        'identifiers': {
                            'cid': cid,
                            'inchi': inchi,
                            'inchi_key': inchi_key,
                            'smiles': smiles
                        },
                        'properties': self._extract_properties(compound_data)
                    }
                }
                
                with ConnectionManager() as conn:
                    with conn.cursor() as cursor:
                        if is_new:
                            # Insert new molecule
                            insert_query = """
                            INSERT INTO molecules 
                            (id, name, properties, created_at, updated_at) 
                            VALUES (%s, %s, %s, NOW(), NOW())
                            """
                            cursor.execute(insert_query, 
                                          (internal_id, name, json.dumps(properties)))
                            
                            logger.debug("Inserted new molecule %s (CID %s)", internal_id, cid)
                        else:
                            # Update existing molecule
                            update_query = """
                            UPDATE molecules 
                            SET name = %s, properties = properties || %s, updated_at = NOW() 
                            WHERE id = %s
                            """
                            cursor.execute(update_query, 
                                          (name, json.dumps(properties), internal_id))
                            
                            logger.debug("Updated molecule %s (CID %s)", internal_id, cid)
                        
                        # Commit transaction
                        conn.commit()
                
                # Successfully processed
                return
                
            except Exception as e:
                retry_count += 1
                logger.warning("Database operation failed for CID %s (attempt %d/%d): %s", 
                              cid, retry_count, self.max_retries, str(e))
                time.sleep(1)
        
        # Failed after max retries
        logger.error("Failed to insert/update CID %s after %d attempts", cid, self.max_retries)
        raise Exception(f"Database operation failed for CID {cid}")
    
    def _extract_properties(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract useful properties from compound data.
        
        Args:
            compound_data: PubChem compound data
            
        Returns:
            Dictionary of extracted properties
        """
        properties = {}
        
        # Extract physical properties
        if 'props' in compound_data:
            for prop in compound_data.get('props', []):
                if prop.get('urn', {}).get('label') == 'LogP':
                    properties['logP'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Water Solubility':
                    properties['water_solubility'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Melting Point':
                    properties['melting_point'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Boiling Point':
                    properties['boiling_point'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Complexity':
                    properties['complexity'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'H-Bond Donor':
                    properties['h_bond_donor_count'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'H-Bond Acceptor':
                    properties['h_bond_acceptor_count'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'Rotatable Bond':
                    properties['rotatable_bond_count'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'Heavy Atom':
                    properties['heavy_atom_count'] = prop.get('value', {}).get('ival')
        
        return properties
    
    def _process_batch(self, cids: List[str]) -> None:
        """
        Process a batch of compounds in parallel.
        
        Args:
            cids: List of PubChem CIDs to process
        """
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.batch_size) as executor:
            # Submit all tasks
            future_to_cid = {executor.submit(self._process_compound, cid): cid for cid in cids}
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(future_to_cid):
                cid = future_to_cid[future]
                try:
                    success = future.result()
                    if success:
                        logger.debug("Successfully processed CID %s", cid)
                    else:
                        logger.warning("Failed to process CID %s", cid)
                except Exception as e:
                    logger.error("Exception processing CID %s: %s", cid, str(e))
                    self.checkpoint['failed_cids'][cid] = str(e)
                    self.stats['compounds_failed'] += 1
    
    def import_compounds(self, cid_list: List[str]) -> None:
        """
        Import compounds from a list of CIDs.
        
        Args:
            cid_list: List of PubChem CIDs to import
        """
        try:
            if not cid_list:
                logger.warning("Empty CID list provided, nothing to import")
                return
            
            # Update checkpoint with total compounds
            self.checkpoint['total_compounds'] = len(cid_list)
            
            # Start timing
            self.stats['start_time'] = time.time()
            
            # Process in batches
            total_batches = (len(cid_list) + self.batch_size - 1) // self.batch_size
            
            for batch_idx in range(total_batches):
                if batch_idx <= self.checkpoint.get('last_completed_batch', -1):
                    logger.info("Skipping batch %d/%d (already completed)", 
                                batch_idx + 1, total_batches)
                    continue
                
                # Get batch of CIDs
                start_idx = batch_idx * self.batch_size
                end_idx = min((batch_idx + 1) * self.batch_size, len(cid_list))
                batch_cids = cid_list[start_idx:end_idx]
                
                logger.info("Processing batch %d/%d (%d compounds)", 
                           batch_idx + 1, total_batches, len(batch_cids))
                
                # Process batch
                self._process_batch(batch_cids)
                
                # Update checkpoint
                self.checkpoint['last_completed_batch'] = batch_idx
                self._save_checkpoint()
                
                # Log progress
                progress = min(100.0, 100.0 * (batch_idx + 1) / total_batches)
                logger.info("Progress: %.1f%% (%d/%d batches, %d/%d compounds)", 
                           progress, batch_idx + 1, total_batches, 
                           len(self.checkpoint['processed_cids']), 
                           self.checkpoint['total_compounds'])
                
                # Sleep to avoid overtaxing PubChem API
                time.sleep(1)
            
            # Finish timing
            self.stats['end_time'] = time.time()
            self.stats['total_elapsed_time'] = self.stats['end_time'] - self.stats['start_time']
            
            # Final report
            logger.info("Import completed!")
            logger.info("Total compounds processed: %d", self.stats['compounds_processed'])
            logger.info("  - Added: %d", self.stats['compounds_added'])
            logger.info("  - Updated: %d", self.stats['compounds_updated'])
            logger.info("  - Skipped: %d", self.stats['compounds_skipped'])
            logger.info("  - Failed: %d", self.stats['compounds_failed'])
            logger.info("Total elapsed time: %.2f seconds", self.stats['total_elapsed_time'])
            
            # Save final checkpoint
            self._save_checkpoint()
            
        except Exception as e:
            logger.error("Import failed: %s", str(e))
            # Save checkpoint
            self._save_checkpoint()
    
    def generate_report(self, output_file: str) -> None:
        """
        Generate import report.
        
        Args:
            output_file: Path to output report file
        """
        try:
            # Create report data
            report = {
                'timestamp': datetime.now().isoformat(),
                'statistics': self.stats,
                'processed_compounds': len(self.checkpoint['processed_cids']),
                'failed_compounds': len(self.checkpoint['failed_cids']),
                'total_compounds': self.checkpoint['total_compounds'],
                'completion_percentage': (len(self.checkpoint['processed_cids']) / 
                                         max(1, self.checkpoint['total_compounds']) * 100),
                'failed_compound_details': self.checkpoint['failed_cids']
            }
            
            # Calculate rates
            if self.stats['total_elapsed_time'] > 0:
                report['processing_rate'] = (self.stats['compounds_processed'] / 
                                          self.stats['total_elapsed_time'])
            else:
                report['processing_rate'] = 0
                
            # Save report
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)
            
            logger.info("Report saved to %s", output_file)
            
        except Exception as e:
            logger.error("Failed to generate report: %s", str(e))

def main():
    """Main function for PubChem import."""
    parser = argparse.ArgumentParser(description='Import PubChem compounds to Supabase')
    parser.add_argument('--cid-file', type=str, help='File containing PubChem CIDs (one per line)')
    parser.add_argument('--cid-list', type=str, help='Comma-separated list of PubChem CIDs')
    parser.add_argument('--batch-size', type=int, default=BATCH_SIZE, 
                       help='Number of compounds to process in parallel')
    parser.add_argument('--checkpoint-file', type=str, default=CHECKPOINT_FILE,
                       help='Path to checkpoint file')
    parser.add_argument('--report-file', type=str, 
                       default=f"reports/pubchem_import_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                       help='Path to output report file')
    args = parser.parse_args()
    
    # Get CID list
    cid_list = []
    if args.cid_file:
        try:
            with open(args.cid_file, 'r') as f:
                cid_list = [line.strip() for line in f if line.strip()]
            logger.info("Loaded %d CIDs from %s", len(cid_list), args.cid_file)
        except Exception as e:
            logger.error("Failed to load CID file: %s", str(e))
            return 1
    elif args.cid_list:
        cid_list = [cid.strip() for cid in args.cid_list.split(',') if cid.strip()]
        logger.info("Using %d CIDs from command line", len(cid_list))
    else:
        # Use CIDs from identifier manager
        id_manager = CryoprotectantIdentifierManager.get_instance()
        cid_list = id_manager.get_all_pubchem_cids()
        logger.info("Using %d CIDs from identifier manager", len(cid_list))
    
    if not cid_list:
        logger.error("No CIDs provided. Use --cid-file or --cid-list option.")
        return 1
    
    # Create importer
    importer = PubChemImporter(
        checkpoint_file=args.checkpoint_file,
        batch_size=args.batch_size
    )
    
    # Import compounds
    importer.import_compounds(cid_list)
    
    # Generate report
    importer.generate_report(args.report_file)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
