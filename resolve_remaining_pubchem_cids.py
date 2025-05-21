#!/usr/bin/env python3
"""
Resolve remaining molecules without PubChem CIDs.

This script attempts to find PubChem CIDs for the 7 remaining molecules
that couldn't be resolved through the standard process, using additional
lookup methods and APIs.
"""

import os
import sys
import json
import time
import uuid
import logging
import argparse
import psycopg2
import requests
from psycopg2.extras import RealDictCursor
from typing import Dict, List, Any, Set, Optional, Tuple

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/pubchem_resolution.log')
    ]
)

# Create logs directory if it doesn't exist
os.makedirs('logs', exist_ok=True)
os.makedirs('reports', exist_ok=True)

logger = logging.getLogger(__name__)

# Database connection parameters (from environment variables)
DB_PARAMS = {
    'host': os.getenv('SUPABASE_DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
    'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
    'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
    'user': os.getenv('SUPABASE_DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
    'password': os.getenv('SUPABASE_DB_PASSWORD', 'LDHt$rkaM&Gmf3X@LQ37'),
    'sslmode': 'require'
}

# API configuration
PUBCHEM_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
REQUEST_DELAY = 0.5  # Delay between API requests (seconds)
MAX_RETRIES = 3     # Maximum number of retries for API requests
USE_ALTERNATE_APIS = True  # Set to True to use additional APIs beyond PubChem

def connect_to_database():
    """
    Connect to the database.
    
    Returns:
        Database connection
    """
    logger.info("Connecting to database...")
    try:
        conn = psycopg2.connect(
            **DB_PARAMS,
            cursor_factory=RealDictCursor
        )
        logger.info("Connected to database")
        return conn
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        sys.exit(1)

def get_unresolved_molecules(conn) -> List[Dict[str, Any]]:
    """
    Get molecules without PubChem CIDs.
    
    Args:
        conn: Database connection
        
    Returns:
        List of molecules without PubChem CIDs
    """
    logger.info("Getting molecules without PubChem CIDs...")
    
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
            SELECT 
                id, name, smiles, inchikey, chembl_id
            FROM 
                molecules
            WHERE 
                chembl_id IS NOT NULL
                AND pubchem_cid IS NULL
            """)
            
            molecules = [dict(row) for row in cursor.fetchall()]
            logger.info(f"Found {len(molecules)} molecules without PubChem CIDs")
            return molecules
    except Exception as e:
        logger.error(f"Error getting unresolved molecules: {str(e)}")
        raise

def resolve_pubchem_cid_by_inchikey(inchikey: str) -> Optional[int]:
    """
    Resolve PubChem CID using InChIKey.
    
    Args:
        inchikey: InChIKey string
        
    Returns:
        PubChem CID or None if not found
    """
    if not inchikey:
        return None
        
    url = f"{PUBCHEM_API_BASE}/compound/inchikey/{inchikey}/cids/JSON"
    
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url)
            time.sleep(REQUEST_DELAY)  # Be nice to the PubChem API
            
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                    cids = data['IdentifierList']['CID']
                    if cids and len(cids) > 0:
                        return cids[0]
            elif response.status_code == 404:
                # Not found
                return None
                
        except Exception as e:
            logger.warning(f"Error resolving PubChem CID for InChIKey {inchikey} (attempt {attempt+1}/{MAX_RETRIES}): {str(e)}")
        
        # Wait longer between retries
        time.sleep(REQUEST_DELAY * (attempt + 1))
    
    return None

def resolve_pubchem_cid_by_chembl_id(chembl_id: str) -> Optional[int]:
    """
    Resolve PubChem CID using ChEMBL ID via ChEMBL API.
    
    Args:
        chembl_id: ChEMBL ID string
        
    Returns:
        PubChem CID or None if not found
    """
    if not chembl_id:
        return None
        
    url = f"{CHEMBL_API_BASE}/molecule/{chembl_id}.json"
    
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url)
            time.sleep(REQUEST_DELAY)  # Be nice to the ChEMBL API
            
            if response.status_code == 200:
                data = response.json()
                
                # Look for PubChem CID in cross-references
                if 'cross_references' in data and data['cross_references']:
                    for xref in data['cross_references']:
                        if xref['xref_src'].lower() == 'pubchem':
                            xref_id = xref['xref_id']
                            try:
                                return int(xref_id)
                            except ValueError:
                                logger.warning(f"Invalid PubChem CID format in ChEMBL API response: {xref_id}")
                
                # If no cross-references, try to use InChIKey from the response
                if 'molecule_structures' in data and data['molecule_structures']:
                    if 'standard_inchi_key' in data['molecule_structures']:
                        inchikey = data['molecule_structures']['standard_inchi_key']
                        if inchikey:
                            logger.info(f"Using InChIKey from ChEMBL API: {inchikey}")
                            return resolve_pubchem_cid_by_inchikey(inchikey)
                
                return None
            elif response.status_code == 404:
                # Not found
                return None
                
        except Exception as e:
            logger.warning(f"Error resolving PubChem CID for ChEMBL ID {chembl_id} (attempt {attempt+1}/{MAX_RETRIES}): {str(e)}")
        
        # Wait longer between retries
        time.sleep(REQUEST_DELAY * (attempt + 1))
    
    return None

def resolve_pubchem_cid_by_smiles(smiles: str) -> Optional[int]:
    """
    Resolve PubChem CID using SMILES via PubChem API.
    
    Args:
        smiles: SMILES string
        
    Returns:
        PubChem CID or None if not found
    """
    if not smiles:
        return None
        
    url = f"{PUBCHEM_API_BASE}/compound/smiles/cids/JSON"
    
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.post(url, data={'smiles': smiles})
            time.sleep(REQUEST_DELAY)  # Be nice to the PubChem API
            
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                    cids = data['IdentifierList']['CID']
                    if cids and len(cids) > 0:
                        return cids[0]
            elif response.status_code == 404:
                # Not found
                return None
                
        except Exception as e:
            logger.warning(f"Error resolving PubChem CID for SMILES (attempt {attempt+1}/{MAX_RETRIES}): {str(e)}")
        
        # Wait longer between retries
        time.sleep(REQUEST_DELAY * (attempt + 1))
    
    return None

def update_molecule_pubchem_cid(conn, molecule_id: str, pubchem_cid: int, method: str) -> bool:
    """
    Update molecule PubChem CID.
    
    Args:
        conn: Database connection
        molecule_id: Molecule ID
        pubchem_cid: PubChem CID
        method: Method used to resolve the CID
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Check if another molecule already has this PubChem CID
        with conn.cursor() as cursor:
            cursor.execute("""
            SELECT id, name FROM molecules 
            WHERE pubchem_cid = %s AND id != %s
            """, (pubchem_cid, molecule_id))
            
            existing = cursor.fetchone()
            
            if existing:
                logger.warning(f"PubChem CID {pubchem_cid} already assigned to molecule {existing['name']} (ID: {existing['id']})")
                
                # Add info to the modification_history JSONB field about the duplicate PubChem CID
                modification_note = {
                    "timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
                    "action": "duplicate_pubchem_cid_detected",
                    "details": {
                        "pubchem_cid": pubchem_cid,
                        "resolution_method": method,
                        "duplicate_with": {
                            "id": existing['id'],
                            "name": existing['name']
                        }
                    }
                }
                
                cursor.execute("""
                UPDATE molecules 
                SET modification_history = COALESCE(modification_history, '[]'::jsonb) || %s::jsonb,
                    updated_at = CURRENT_TIMESTAMP
                WHERE id = %s
                """, (
                    json.dumps([modification_note]),
                    molecule_id
                ))
                
                logger.info(f"Added duplicate CID note for molecule ID: {molecule_id}")
                conn.commit()
                return True
            else:
                # Update PubChem CID
                modification_note = {
                    "timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
                    "action": "pubchem_cid_updated",
                    "details": {
                        "pubchem_cid": pubchem_cid,
                        "resolution_method": method
                    }
                }
                
                cursor.execute("""
                UPDATE molecules 
                SET pubchem_cid = %s,
                    modification_history = COALESCE(modification_history, '[]'::jsonb) || %s::jsonb,
                    updated_at = CURRENT_TIMESTAMP
                WHERE id = %s
                """, (
                    pubchem_cid,
                    json.dumps([modification_note]),
                    molecule_id
                ))
                
                logger.info(f"Updated PubChem CID for molecule ID: {molecule_id} to {pubchem_cid} using {method}")
                conn.commit()
                return True
    except Exception as e:
        logger.error(f"Error updating molecule PubChem CID: {str(e)}")
        conn.rollback()
        return False

def resolve_molecules(conn, molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Attempt to resolve PubChem CIDs for molecules.
    
    Args:
        conn: Database connection
        molecules: List of molecules to resolve
        
    Returns:
        Dict with resolution results
    """
    logger.info(f"Attempting to resolve PubChem CIDs for {len(molecules)} molecules...")
    
    resolved_count = 0
    failed_count = 0
    resolution_methods = {}
    errors = []
    
    for idx, molecule in enumerate(molecules):
        try:
            molecule_id = molecule['id']
            name = molecule.get('name', 'Unknown')
            chembl_id = molecule.get('chembl_id')
            inchikey = molecule.get('inchikey')
            smiles = molecule.get('smiles')
            
            logger.info(f"[{idx+1}/{len(molecules)}] Processing molecule {name} (ID: {molecule_id})")
            
            # Try resolving by InChIKey (if available)
            if inchikey:
                pubchem_cid = resolve_pubchem_cid_by_inchikey(inchikey)
                if pubchem_cid:
                    if update_molecule_pubchem_cid(conn, molecule_id, pubchem_cid, "inchikey_lookup"):
                        resolved_count += 1
                        resolution_methods["inchikey_lookup"] = resolution_methods.get("inchikey_lookup", 0) + 1
                        continue
            
            # Try resolving by ChEMBL ID (if alternate APIs enabled)
            if USE_ALTERNATE_APIS and chembl_id:
                pubchem_cid = resolve_pubchem_cid_by_chembl_id(chembl_id)
                if pubchem_cid:
                    if update_molecule_pubchem_cid(conn, molecule_id, pubchem_cid, "chembl_api_lookup"):
                        resolved_count += 1
                        resolution_methods["chembl_api_lookup"] = resolution_methods.get("chembl_api_lookup", 0) + 1
                        continue
            
            # Try resolving by SMILES (last resort)
            if smiles:
                pubchem_cid = resolve_pubchem_cid_by_smiles(smiles)
                if pubchem_cid:
                    if update_molecule_pubchem_cid(conn, molecule_id, pubchem_cid, "smiles_lookup"):
                        resolved_count += 1
                        resolution_methods["smiles_lookup"] = resolution_methods.get("smiles_lookup", 0) + 1
                        continue
            
            # If all methods failed, record the error
            logger.warning(f"All resolution methods failed for molecule {name} (ID: {molecule_id})")
            failure_details = {
                "id": molecule_id,
                "name": name,
                "chembl_id": chembl_id,
                "has_inchikey": inchikey is not None,
                "has_smiles": smiles is not None
            }
            errors.append(failure_details)
            failed_count += 1
            
        except Exception as e:
            logger.error(f"Error processing molecule {molecule.get('name', 'Unknown')}: {str(e)}")
            errors.append({
                "id": molecule.get('id'),
                "name": molecule.get('name', 'Unknown'),
                "error": str(e)
            })
            failed_count += 1
    
    return {
        "resolved": resolved_count,
        "failed": failed_count,
        "resolution_methods": resolution_methods,
        "errors": errors
    }

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Resolve remaining molecules without PubChem CIDs')
    parser.add_argument('--use-alternate-apis', action='store_true', help='Use additional APIs beyond PubChem (ChEMBL API)')
    
    args = parser.parse_args()
    
    global USE_ALTERNATE_APIS
    USE_ALTERNATE_APIS = args.use_alternate_apis
    
    # Connect to database
    conn = connect_to_database()
    
    try:
        # Get unresolved molecules
        molecules = get_unresolved_molecules(conn)
        
        if not molecules:
            logger.info("No unresolved molecules found. All molecules have PubChem CIDs.")
            return 0
        
        # Resolve molecules
        results = resolve_molecules(conn, molecules)
        
        # Save results
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        results_file = f"reports/pubchem_resolution_{timestamp}.json"
        
        with open(results_file, 'w') as f:
            json.dump({
                "timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
                "total_molecules": len(molecules),
                "results": results
            }, f, indent=2)
        
        logger.info(f"Results saved to {results_file}")
        
        # Print summary
        print("\n=== PubChem CID Resolution Summary ===")
        print(f"Total molecules: {len(molecules)}")
        print(f"Successfully resolved: {results['resolved']}/{len(molecules)} ({results['resolved'] / len(molecules) * 100:.2f}%)")
        print(f"Failed: {results['failed']}/{len(molecules)} ({results['failed'] / len(molecules) * 100:.2f}%)")
        
        if results['resolution_methods']:
            print("\nResolution Methods:")
            for method, count in results['resolution_methods'].items():
                print(f"  - {method}: {count}")
        
        if results['failed'] > 0:
            print("\nFailed Molecules:")
            for error in results['errors']:
                print(f"  - {error.get('name', 'Unknown')} (ID: {error.get('id', 'Unknown')})")
        
        print(f"\nDetailed report saved to: {results_file}")
        print("========================================\n")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error in main function: {str(e)}")
        return 1
    
    finally:
        conn.close()
        logger.info("Database connection closed")

if __name__ == "__main__":
    sys.exit(main())