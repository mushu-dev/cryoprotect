#!/usr/bin/env python3
"""
ChEMBL to PubChem cross-reference resolver.

This module resolves ChEMBL IDs to PubChem CIDs using multiple methods:
1. Direct API lookups using PubChem PUG REST API
2. InChIKey-based matching
3. SMILES-based similarity search
"""

import os
import sys
import uuid
import json
import logging
import time
import requests
from typing import Dict, List, Any, Optional, Tuple, Set
import psycopg2
from psycopg2.extras import RealDictCursor
from urllib.parse import quote

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Try importing RDKit for structure comparison
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    logger.warning("RDKit not available. Using fallback comparison methods.")
    RDKIT_AVAILABLE = False

# Database connection parameters (from environment or .env file)
DB_PARAMS = {
    'host': os.getenv('SUPABASE_DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
    'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
    'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
    'user': os.getenv('SUPABASE_DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
    'password': os.getenv('SUPABASE_DB_PASSWORD', 'LDHt$rkaM&Gmf3X@LQ37'),
    'sslmode': 'require'
}

# PubChem API Configuration
PUBCHEM_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW_BASE = "https://pubchem.ncbi.nlm.nih.gov/compound/"
REQUEST_DELAY = 0.5  # Time between API requests (seconds)
MAX_RETRIES = 3      # Maximum number of retries for API requests

class ChEMBLPubChemResolver:
    """
    Resolver for ChEMBL to PubChem cross-references.
    
    This class uses multiple methods to find PubChem CIDs for ChEMBL compounds.
    """
    
    def __init__(self, db_params: Dict[str, Any]):
        """
        Initialize the resolver.
        
        Args:
            db_params: Database connection parameters
        """
        self.db_params = db_params
        self.conn = None
        self.session = requests.Session()
        self.last_request_time = 0
        
        # Connect to database
        self._connect()
    
    def _connect(self):
        """Connect to the database."""
        try:
            self.conn = psycopg2.connect(**self.db_params)
            logger.info("Connected to database")
        except Exception as e:
            logger.error(f"Database connection error: {str(e)}")
            sys.exit(1)
    
    def _api_request(self, url: str, params: Dict[str, str] = None) -> Optional[Dict[str, Any]]:
        """
        Make a request to the PubChem API with rate limiting and retries.
        
        Args:
            url: URL to request
            params: Dictionary of URL parameters
            
        Returns:
            JSON response as dictionary or None if request failed
        """
        # Rate limiting
        now = time.time()
        time_since_last_request = now - self.last_request_time
        if time_since_last_request < REQUEST_DELAY:
            time.sleep(REQUEST_DELAY - time_since_last_request)
        
        self.last_request_time = time.time()
        
        # Make request with retries
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                response = self.session.get(url, params=params, timeout=30)
                response.raise_for_status()
                
                return response.json()
            except (requests.RequestException, json.JSONDecodeError) as e:
                logger.warning(f"API request failed (attempt {attempt}/{MAX_RETRIES}): {str(e)}")
                
                if attempt < MAX_RETRIES:
                    # Exponential backoff
                    sleep_time = REQUEST_DELAY * (2 ** (attempt - 1))
                    time.sleep(sleep_time)
                else:
                    logger.error(f"API request failed after {MAX_RETRIES} attempts: {str(e)}")
                    return None
    
    def find_pubchem_by_chembl_id(self, chembl_id: str) -> Optional[int]:
        """
        Find PubChem CID by ChEMBL ID using direct PubChem lookup.
        
        Args:
            chembl_id: ChEMBL ID
            
        Returns:
            PubChem CID as integer or None if not found
        """
        url = f"{PUBCHEM_API_BASE}/compound/xref/ChEMBL/{quote(chembl_id)}/cids/JSON"
        
        try:
            response = self._api_request(url)
            
            if response and 'IdentifierList' in response and 'CID' in response['IdentifierList']:
                cids = response['IdentifierList']['CID']
                if cids:
                    return cids[0]  # Return the first CID
            
            return None
        except Exception as e:
            logger.error(f"Error finding PubChem CID for ChEMBL ID {chembl_id}: {str(e)}")
            return None
    
    def find_pubchem_by_inchikey(self, inchikey: str) -> Optional[int]:
        """
        Find PubChem CID by InChIKey.
        
        Args:
            inchikey: InChIKey
            
        Returns:
            PubChem CID as integer or None if not found
        """
        url = f"{PUBCHEM_API_BASE}/compound/inchikey/{quote(inchikey)}/cids/JSON"
        
        try:
            response = self._api_request(url)
            
            if response and 'IdentifierList' in response and 'CID' in response['IdentifierList']:
                cids = response['IdentifierList']['CID']
                if cids:
                    return cids[0]  # Return the first CID
            
            return None
        except Exception as e:
            logger.error(f"Error finding PubChem CID for InChIKey {inchikey}: {str(e)}")
            return None
    
    def find_pubchem_by_smiles(self, smiles: str) -> Optional[int]:
        """
        Find PubChem CID by SMILES string.
        
        Args:
            smiles: SMILES string
            
        Returns:
            PubChem CID as integer or None if not found
        """
        url = f"{PUBCHEM_API_BASE}/compound/smiles/cids/JSON"
        
        try:
            response = self._api_request(url, params={'smiles': smiles})
            
            if response and 'IdentifierList' in response and 'CID' in response['IdentifierList']:
                cids = response['IdentifierList']['CID']
                if cids:
                    return cids[0]  # Return the first CID
            
            return None
        except Exception as e:
            logger.error(f"Error finding PubChem CID for SMILES {smiles}: {str(e)}")
            return None
    
    def fetch_molecule_details(self, pubchem_cid: int) -> Optional[Dict[str, Any]]:
        """
        Fetch molecule details from PubChem.
        
        Args:
            pubchem_cid: PubChem CID
            
        Returns:
            Dictionary of molecule details or None if request failed
        """
        url = f"{PUBCHEM_API_BASE}/compound/cid/{pubchem_cid}/JSON"
        
        try:
            response = self._api_request(url)
            
            if response and 'PC_Compounds' in response and len(response['PC_Compounds']) > 0:
                return response['PC_Compounds'][0]
            
            return None
        except Exception as e:
            logger.error(f"Error fetching details for PubChem CID {pubchem_cid}: {str(e)}")
            return None
    
    def get_molecules_without_pubchem(self, batch_size: int = 50) -> List[Dict[str, Any]]:
        """
        Get ChEMBL molecules without PubChem CIDs.
        
        Args:
            batch_size: Number of molecules to fetch at once
            
        Returns:
            List of molecules without PubChem CIDs
        """
        query = """
        SELECT id, name, chembl_id, inchikey, smiles
        FROM molecules
        WHERE chembl_id IS NOT NULL
        AND pubchem_cid IS NULL
        ORDER BY updated_at DESC
        LIMIT %s;
        """
        
        try:
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute(query, (batch_size,))
                return cursor.fetchall()
        except Exception as e:
            logger.error(f"Error fetching molecules without PubChem CIDs: {str(e)}")
            return []
    
    def update_pubchem_cid(self, molecule_id: str, pubchem_cid: int, method: str) -> bool:
        """
        Update the PubChem CID for a molecule.
        
        Args:
            molecule_id: Molecule ID
            pubchem_cid: PubChem CID
            method: Method used to find the CID
            
        Returns:
            True if update was successful, False otherwise
        """
        try:
            pubchem_link = f"{PUBCHEM_VIEW_BASE}{pubchem_cid}"
            
            # Start transaction
            with self.conn:
                with self.conn.cursor() as cursor:
                    # Update molecule
                    query = """
                    UPDATE molecules
                    SET pubchem_cid = %s, 
                        cid = %s,
                        pubchem_link = %s,
                        updated_at = NOW(),
                        modification_history = COALESCE(modification_history, '[]'::jsonb) || %s::jsonb
                    WHERE id = %s;
                    """
                    
                    # Create modification history entry
                    history_entry = json.dumps([{
                        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                        "operation": "update",
                        "field": "pubchem_cid",
                        "value": pubchem_cid,
                        "method": method
                    }])
                    
                    cursor.execute(
                        query, 
                        (
                            pubchem_cid,
                            pubchem_cid,
                            pubchem_link,
                            history_entry,
                            molecule_id
                        )
                    )
                    
                    # Commit transaction
                    self.conn.commit()
                    
                    logger.info(f"Updated PubChem CID for molecule {molecule_id} to {pubchem_cid} (method: {method})")
                    return True
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Error updating PubChem CID for molecule {molecule_id}: {str(e)}")
            return False
    
    def resolve_all_missing_pubchem_cids(self, batch_size: int = 50, limit: int = None) -> Dict[str, int]:
        """
        Resolve PubChem CIDs for all ChEMBL molecules without them.
        
        Args:
            batch_size: Number of molecules to process at once
            limit: Maximum number of molecules to process (None for all)
            
        Returns:
            Statistics about the resolution process
        """
        stats = {
            'total_molecules': 0,
            'resolved_by_chembl': 0,
            'resolved_by_inchikey': 0,
            'resolved_by_smiles': 0,
            'failed_to_resolve': 0
        }
        
        processed_count = 0
        
        while True:
            # Check if we've reached the limit
            if limit is not None and processed_count >= limit:
                logger.info(f"Reached limit of {limit} molecules")
                break
            
            # Get batch of molecules without PubChem CIDs
            molecules = self.get_molecules_without_pubchem(batch_size)
            
            if not molecules:
                logger.info("No more molecules without PubChem CIDs")
                break
            
            logger.info(f"Processing batch of {len(molecules)} molecules")
            stats['total_molecules'] += len(molecules)
            
            # Process each molecule
            for molecule in molecules:
                molecule_id = molecule['id']
                name = molecule['name']
                chembl_id = molecule['chembl_id']
                inchikey = molecule['inchikey']
                smiles = molecule['smiles']
                
                logger.info(f"Processing molecule: {name} (ChEMBL ID: {chembl_id})")
                
                # Try to find PubChem CID using different methods
                pubchem_cid = None
                method = ""
                
                # Method 1: Direct ChEMBL ID lookup
                if chembl_id:
                    pubchem_cid = self.find_pubchem_by_chembl_id(chembl_id)
                    if pubchem_cid:
                        method = "chembl_id"
                        stats['resolved_by_chembl'] += 1
                
                # Method 2: InChIKey lookup
                if not pubchem_cid and inchikey:
                    pubchem_cid = self.find_pubchem_by_inchikey(inchikey)
                    if pubchem_cid:
                        method = "inchikey"
                        stats['resolved_by_inchikey'] += 1
                
                # Method 3: SMILES lookup
                if not pubchem_cid and smiles:
                    pubchem_cid = self.find_pubchem_by_smiles(smiles)
                    if pubchem_cid:
                        method = "smiles"
                        stats['resolved_by_smiles'] += 1
                
                # Update PubChem CID if found
                if pubchem_cid:
                    success = self.update_pubchem_cid(molecule_id, pubchem_cid, method)
                    if not success:
                        logger.error(f"Failed to update PubChem CID for molecule {name} (ChEMBL ID: {chembl_id})")
                else:
                    logger.warning(f"Could not find PubChem CID for molecule {name} (ChEMBL ID: {chembl_id})")
                    stats['failed_to_resolve'] += 1
                
                processed_count += 1
                
                # Check if we've reached the limit
                if limit is not None and processed_count >= limit:
                    break
                
                # Small delay to avoid overloading the API
                time.sleep(0.1)
        
        return stats

def main():
    """Main function."""
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Resolve ChEMBL to PubChem cross-references')
    parser.add_argument('--batch-size', type=int, default=50, help='Batch size')
    parser.add_argument('--limit', type=int, default=None, help='Maximum number of molecules to process')
    parser.add_argument('--molecule-id', type=str, default=None, help='Process a specific molecule ID')
    parser.add_argument('--chembl-id', type=str, default=None, help='Look up a specific ChEMBL ID')
    args = parser.parse_args()
    
    # Initialize resolver
    resolver = ChEMBLPubChemResolver(DB_PARAMS)
    
    if args.chembl_id:
        # Look up a specific ChEMBL ID
        logger.info(f"Looking up ChEMBL ID: {args.chembl_id}")
        
        pubchem_cid = resolver.find_pubchem_by_chembl_id(args.chembl_id)
        
        if pubchem_cid:
            logger.info(f"Found PubChem CID: {pubchem_cid}")
            
            # Fetch molecule details
            details = resolver.fetch_molecule_details(pubchem_cid)
            
            if details:
                logger.info(f"Molecule details: {json.dumps(details, indent=2)}")
            else:
                logger.error(f"Could not fetch molecule details for PubChem CID {pubchem_cid}")
        else:
            logger.error(f"Could not find PubChem CID for ChEMBL ID {args.chembl_id}")
    elif args.molecule_id:
        # Process specific molecule
        logger.info(f"Processing molecule ID: {args.molecule_id}")
        
        # Get molecule info
        query = "SELECT name, chembl_id, inchikey, smiles FROM molecules WHERE id = %s;"
        with resolver.conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute(query, (args.molecule_id,))
            molecule = cursor.fetchone()
            
            if not molecule:
                logger.error(f"Molecule ID {args.molecule_id} not found")
                return
            
            name = molecule['name']
            chembl_id = molecule['chembl_id']
            inchikey = molecule['inchikey']
            smiles = molecule['smiles']
            
            logger.info(f"Molecule: {name}")
            logger.info(f"ChEMBL ID: {chembl_id}")
            
            # Try to find PubChem CID using different methods
            pubchem_cid = None
            method = ""
            
            # Method 1: Direct ChEMBL ID lookup
            if chembl_id:
                pubchem_cid = resolver.find_pubchem_by_chembl_id(chembl_id)
                if pubchem_cid:
                    method = "chembl_id"
            
            # Method 2: InChIKey lookup
            if not pubchem_cid and inchikey:
                pubchem_cid = resolver.find_pubchem_by_inchikey(inchikey)
                if pubchem_cid:
                    method = "inchikey"
            
            # Method 3: SMILES lookup
            if not pubchem_cid and smiles:
                pubchem_cid = resolver.find_pubchem_by_smiles(smiles)
                if pubchem_cid:
                    method = "smiles"
            
            # Update PubChem CID if found
            if pubchem_cid:
                logger.info(f"Found PubChem CID: {pubchem_cid} (method: {method})")
                
                success = resolver.update_pubchem_cid(args.molecule_id, pubchem_cid, method)
                
                if success:
                    logger.info(f"Successfully updated PubChem CID for molecule {name}")
                else:
                    logger.error(f"Failed to update PubChem CID for molecule {name}")
            else:
                logger.error(f"Could not find PubChem CID for molecule {name}")
    else:
        # Process all molecules without PubChem CIDs
        logger.info("Processing all ChEMBL molecules without PubChem CIDs")
        
        stats = resolver.resolve_all_missing_pubchem_cids(
            batch_size=args.batch_size,
            limit=args.limit
        )
        
        # Print statistics
        logger.info("=== PubChem CID Resolution Statistics ===")
        logger.info(f"Total molecules processed: {stats['total_molecules']}")
        logger.info(f"Resolved by ChEMBL ID: {stats['resolved_by_chembl']}")
        logger.info(f"Resolved by InChIKey: {stats['resolved_by_inchikey']}")
        logger.info(f"Resolved by SMILES: {stats['resolved_by_smiles']}")
        logger.info(f"Failed to resolve: {stats['failed_to_resolve']}")
        logger.info("========================================")

if __name__ == "__main__":
    main()