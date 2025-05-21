#!/usr/bin/env python3
"""
Fix data quality issues in ChEMBL molecules.

This script addresses:
1. Molecules with missing properties
2. Molecules with missing JSONB properties
3. Molecules missing PubChem cross-references
4. Molecules missing InChIKeys
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
from decimal import Decimal
from psycopg2.extras import RealDictCursor
from typing import Dict, List, Any, Set, Optional, Tuple

# Try importing RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, MolSurf, AllChem, DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    # Try to use mock_rdkit
    try:
        import mock_rdkit
        mock_rdkit.create_mock_rdkit()
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, MolSurf, AllChem, DataStructs
        RDKIT_AVAILABLE = True
        print("Using mock RDKit for basic property calculations")
    except ImportError:
        RDKIT_AVAILABLE = False
        print("WARNING: RDKit not available. Property calculation will be disabled.")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/fix_chembl_issues.log')
    ]
)

# Create logs directory if it doesn't exist
os.makedirs('logs', exist_ok=True)
os.makedirs('reports', exist_ok=True)

logger = logging.getLogger(__name__)

# Custom JSON encoder to handle Decimal types
class CustomJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Decimal):
            return float(obj)
        return super().default(obj)

# Database connection parameters (from environment variables)
DB_PARAMS = {
    'host': os.getenv('SUPABASE_DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
    'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
    'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
    'user': os.getenv('SUPABASE_DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
    'password': os.getenv('SUPABASE_DB_PASSWORD', 'LDHt$rkaM&Gmf3X@LQ37'),
    'sslmode': 'require'
}

# Define property functions
PROPERTY_DEFINITIONS = {
    'LogP': {
        'description': 'Octanol-water partition coefficient',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Descriptors.MolLogP if RDKIT_AVAILABLE else None
    },
    'TPSA': {
        'description': 'Topological polar surface area',
        'data_type': 'numeric',
        'unit': 'Å²',
        'rdkit_func': MolSurf.TPSA if RDKIT_AVAILABLE else None
    },
    'Molecular Weight': {
        'description': 'The molecular weight of the compound',
        'data_type': 'numeric',
        'unit': 'g/mol',
        'rdkit_func': Descriptors.MolWt if RDKIT_AVAILABLE else None
    },
    'Heavy Atom Count': {
        'description': 'Number of non-hydrogen atoms',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Lipinski.HeavyAtomCount if RDKIT_AVAILABLE else None
    },
    'Hydrogen Bond Donor Count': {
        'description': 'Number of hydrogen bond donors',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Lipinski.NumHDonors if RDKIT_AVAILABLE else None
    },
    'Hydrogen Bond Acceptor Count': {
        'description': 'Number of hydrogen bond acceptors',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Lipinski.NumHAcceptors if RDKIT_AVAILABLE else None
    },
    'Rotatable Bond Count': {
        'description': 'Number of rotatable bonds',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Descriptors.NumRotatableBonds if RDKIT_AVAILABLE else None
    },
    'Ring Count': {
        'description': 'Number of rings',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Descriptors.RingCount if RDKIT_AVAILABLE else None
    },
    'Aromatic Ring Count': {
        'description': 'Number of aromatic rings',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': lambda mol: rdMolDescriptors.CalcNumAromaticRings(mol) if RDKIT_AVAILABLE else None
    }
}

# PubChem API configuration
PUBCHEM_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
REQUEST_DELAY = 0.5  # Delay between API requests (seconds)
MAX_RETRIES = 3     # Maximum number of retries for API requests

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

def load_issues_report(file_path: str) -> Dict[str, Any]:
    """
    Load the issues report JSON file.
    
    Args:
        file_path: Path to the issues report file
        
    Returns:
        Dict containing the issues data
    """
    logger.info(f"Loading issues report from {file_path}")
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        logger.info("Issues report loaded successfully")
        return data
    except Exception as e:
        logger.error(f"Error loading issues report: {str(e)}")
        sys.exit(1)

def get_property_types(conn) -> Dict[str, Dict[str, Any]]:
    """
    Get property types from the database.
    
    Args:
        conn: Database connection
        
    Returns:
        Dict mapping property names to property type records
    """
    logger.info("Getting property types from database...")
    
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
            SELECT id, name, description, data_type, units
            FROM property_types
            """)
            
            property_types = {}
            for row in cursor.fetchall():
                property_types[row['name']] = dict(row)
            
            logger.info(f"Found {len(property_types)} property types in database")
            
            # Create any missing property types
            for prop_name, prop_def in PROPERTY_DEFINITIONS.items():
                if prop_name not in property_types:
                    logger.info(f"Creating missing property type: {prop_name}")
                    
                    cursor.execute("""
                    INSERT INTO property_types (id, name, description, data_type, units, created_at, updated_at)
                    VALUES (%s, %s, %s, %s, %s, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)
                    RETURNING id, name, description, data_type, units
                    """, (
                        str(uuid.uuid4()),
                        prop_name,
                        prop_def['description'],
                        prop_def['data_type'],
                        prop_def['unit']
                    ))
                    
                    new_prop = cursor.fetchone()
                    property_types[prop_name] = dict(new_prop)
                    logger.info(f"Created property type: {prop_name} (ID: {new_prop['id']})")
            
            conn.commit()
            return property_types
    except Exception as e:
        logger.error(f"Error getting property types: {str(e)}")
        conn.rollback()
        raise

def calculate_molecular_properties(smiles: str) -> Dict[str, float]:
    """
    Calculate molecular properties using RDKit.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dict mapping property names to calculated values
    """
    if not RDKIT_AVAILABLE:
        logger.error("RDKit is not available for property calculation")
        return {}
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return {}
        
        properties = {}
        
        for prop_name, prop_def in PROPERTY_DEFINITIONS.items():
            if prop_def['rdkit_func'] is not None:
                value = prop_def['rdkit_func'](mol)
                properties[prop_name] = value
        
        return properties
    except Exception as e:
        logger.error(f"Error calculating properties for {smiles}: {str(e)}")
        return {}

def fix_missing_properties(conn, molecules: List[Dict[str, Any]], property_types: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """
    Fix molecules with missing properties.
    
    Args:
        conn: Database connection
        molecules: List of molecules with missing properties
        property_types: Dict mapping property names to property type records
        
    Returns:
        Dict with results
    """
    if not RDKIT_AVAILABLE:
        logger.error("RDKit is not available for fixing missing properties")
        return {"fixed": 0, "failed": len(molecules), "errors": ["RDKit not available"]}
    
    logger.info(f"Fixing {len(molecules)} molecules with missing properties...")
    
    fixed_count = 0
    failed_count = 0
    errors = []
    
    for idx, molecule in enumerate(molecules):
        try:
            molecule_id = molecule['id']
            smiles = molecule['smiles']
            name = molecule.get('name', 'Unknown')
            
            logger.info(f"[{idx+1}/{len(molecules)}] Processing molecule {name} (ID: {molecule_id})")
            
            # Calculate properties
            properties = calculate_molecular_properties(smiles)
            
            if not properties:
                logger.warning(f"Could not calculate properties for {name} (ID: {molecule_id})")
                failed_count += 1
                errors.append(f"Failed to calculate properties for {name} (ID: {molecule_id})")
                continue
            
            # Update both the molecular_properties table and JSONB properties
            inserted_props = 0
            
            # First, insert properties into molecular_properties table
            for prop_name, value in properties.items():
                prop_type = property_types.get(prop_name)
                if not prop_type:
                    logger.warning(f"Property type {prop_name} not found in database")
                    continue
                
                with conn.cursor() as cursor:
                    # Check if property already exists
                    cursor.execute("""
                    SELECT id FROM molecular_properties 
                    WHERE molecule_id = %s AND property_type_id = %s
                    """, (molecule_id, prop_type['id']))
                    
                    existing = cursor.fetchone()
                    
                    if existing:
                        # Update existing property
                        cursor.execute("""
                        UPDATE molecular_properties 
                        SET numeric_value = %s, updated_at = CURRENT_TIMESTAMP
                        WHERE id = %s
                        """, (value, existing['id']))
                    else:
                        # Insert new property
                        cursor.execute("""
                        INSERT INTO molecular_properties (
                            id, molecule_id, property_type_id, numeric_value, 
                            text_value, property_value, created_at, updated_at
                        )
                        VALUES (%s, %s, %s, %s, NULL, %s, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)
                        """, (
                            str(uuid.uuid4()),
                            molecule_id,
                            prop_type['id'],
                            value,
                            str(value)
                        ))
                    
                    inserted_props += 1
            
            # Then, update the JSONB properties field
            with conn.cursor() as cursor:
                cursor.execute("""
                UPDATE molecules 
                SET properties = COALESCE(properties, '{}'::jsonb) || %s::jsonb,
                    updated_at = CURRENT_TIMESTAMP
                WHERE id = %s
                """, (json.dumps(properties, cls=CustomJSONEncoder), molecule_id))
            
            logger.info(f"Updated {inserted_props} properties for molecule {name} (ID: {molecule_id})")
            fixed_count += 1
            conn.commit()
            
        except Exception as e:
            logger.error(f"Error fixing properties for molecule {molecule.get('name', 'Unknown')}: {str(e)}")
            errors.append(f"Error for {molecule.get('name', 'Unknown')}: {str(e)}")
            failed_count += 1
            conn.rollback()
    
    return {
        "fixed": fixed_count,
        "failed": failed_count,
        "errors": errors
    }

def fix_missing_jsonb_properties(conn, molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Fix molecules with missing JSONB properties.
    
    Args:
        conn: Database connection
        molecules: List of molecules with missing JSONB properties
        
    Returns:
        Dict with results
    """
    logger.info(f"Fixing {len(molecules)} molecules with missing JSONB properties...")
    
    fixed_count = 0
    failed_count = 0
    errors = []
    
    for idx, molecule in enumerate(molecules):
        try:
            molecule_id = molecule['id']
            name = molecule.get('name', 'Unknown')
            
            logger.info(f"[{idx+1}/{len(molecules)}] Processing molecule {name} (ID: {molecule_id})")
            
            # Get properties from molecular_properties table
            with conn.cursor() as cursor:
                cursor.execute("""
                SELECT 
                    pt.name AS property_name,
                    mp.numeric_value,
                    mp.text_value
                FROM 
                    molecular_properties mp
                    JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE 
                    mp.molecule_id = %s
                """, (molecule_id,))
                
                prop_rows = cursor.fetchall()
            
            # Convert to JSONB format
            properties = {}
            for prop in prop_rows:
                prop_name = prop['property_name']
                value = prop['numeric_value'] if prop['numeric_value'] is not None else prop['text_value']
                properties[prop_name] = value
            
            # Update JSONB properties field
            if properties:
                with conn.cursor() as cursor:
                    cursor.execute("""
                    UPDATE molecules 
                    SET properties = %s::jsonb,
                        updated_at = CURRENT_TIMESTAMP
                    WHERE id = %s
                    """, (json.dumps(properties, cls=CustomJSONEncoder), molecule_id))
                
                logger.info(f"Updated JSONB properties for molecule {name} (ID: {molecule_id}) with {len(properties)} properties")
                fixed_count += 1
                conn.commit()
            else:
                logger.warning(f"No properties found for molecule {name} (ID: {molecule_id})")
                failed_count += 1
                errors.append(f"No properties found for {name} (ID: {molecule_id})")
        
        except Exception as e:
            logger.error(f"Error fixing JSONB properties for molecule {molecule.get('name', 'Unknown')}: {str(e)}")
            errors.append(f"Error for {molecule.get('name', 'Unknown')}: {str(e)}")
            failed_count += 1
            conn.rollback()
    
    return {
        "fixed": fixed_count,
        "failed": failed_count,
        "errors": errors
    }

def resolve_pubchem_cid_by_inchikey(inchikey: str) -> Optional[int]:
    """
    Resolve PubChem CID using InChIKey.
    
    Args:
        inchikey: InChIKey string
        
    Returns:
        PubChem CID or None if not found
    """
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

def fix_missing_pubchem_cids(conn, molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Fix molecules missing PubChem CIDs.
    
    Args:
        conn: Database connection
        molecules: List of molecules missing PubChem CIDs
        
    Returns:
        Dict with results
    """
    logger.info(f"Fixing {len(molecules)} molecules missing PubChem CIDs...")
    
    fixed_count = 0
    failed_count = 0
    errors = []
    
    for idx, molecule in enumerate(molecules):
        try:
            molecule_id = molecule['id']
            name = molecule.get('name', 'Unknown')
            chembl_id = molecule.get('chembl_id')
            inchikey = molecule.get('inchikey')
            
            logger.info(f"[{idx+1}/{len(molecules)}] Processing molecule {name} (ID: {molecule_id})")
            
            if not inchikey:
                logger.warning(f"Cannot resolve PubChem CID for {name} (ID: {molecule_id}): Missing InChIKey")
                failed_count += 1
                errors.append(f"Missing InChIKey for {name} (ID: {molecule_id})")
                continue
            
            # Resolve PubChem CID using InChIKey
            pubchem_cid = resolve_pubchem_cid_by_inchikey(inchikey)
            
            if pubchem_cid:
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
                            json.dumps([modification_note], cls=CustomJSONEncoder),
                            molecule_id
                        ))
                        
                        logger.info(f"Added duplicate CID note for molecule {name} (ID: {molecule_id}) referencing {existing['name']} (ID: {existing['id']})")
                        fixed_count += 1
                    else:
                        # Update PubChem CID
                        cursor.execute("""
                        UPDATE molecules 
                        SET pubchem_cid = %s,
                            updated_at = CURRENT_TIMESTAMP
                        WHERE id = %s
                        """, (pubchem_cid, molecule_id))
                        
                        logger.info(f"Updated PubChem CID for molecule {name} (ID: {molecule_id}) to {pubchem_cid}")
                        fixed_count += 1
                
                conn.commit()
            else:
                logger.warning(f"Could not resolve PubChem CID for {name} (ID: {molecule_id}, InChIKey: {inchikey})")
                failed_count += 1
                errors.append(f"Could not resolve PubChem CID for {name} (ID: {molecule_id})")
        
        except Exception as e:
            logger.error(f"Error fixing PubChem CID for molecule {molecule.get('name', 'Unknown')}: {str(e)}")
            errors.append(f"Error for {molecule.get('name', 'Unknown')}: {str(e)}")
            failed_count += 1
            conn.rollback()
    
    return {
        "fixed": fixed_count,
        "failed": failed_count,
        "errors": errors
    }

def generate_inchikey_from_smiles(smiles: str) -> Optional[str]:
    """
    Generate InChIKey from SMILES using RDKit.
    
    Args:
        smiles: SMILES string
        
    Returns:
        InChIKey or None if unable to generate
    """
    if not RDKIT_AVAILABLE:
        logger.error("RDKit is not available for InChIKey generation")
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None
        
        inchi = Chem.MolToInchi(mol)
        if inchi:
            inchikey = Chem.InchiToInchiKey(inchi)
            return inchikey
        
        return None
    except Exception as e:
        logger.error(f"Error generating InChIKey for {smiles}: {str(e)}")
        return None

def fix_missing_inchikeys(conn, molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Fix molecules missing InChIKeys.
    
    Args:
        conn: Database connection
        molecules: List of molecules missing InChIKeys
        
    Returns:
        Dict with results
    """
    if not RDKIT_AVAILABLE:
        logger.error("RDKit is not available for fixing missing InChIKeys")
        return {"fixed": 0, "failed": len(molecules), "errors": ["RDKit not available"]}
    
    logger.info(f"Fixing {len(molecules)} molecules missing InChIKeys...")
    
    fixed_count = 0
    failed_count = 0
    errors = []
    
    for idx, molecule in enumerate(molecules):
        try:
            molecule_id = molecule['id']
            name = molecule.get('name', 'Unknown')
            smiles = molecule.get('smiles')
            
            logger.info(f"[{idx+1}/{len(molecules)}] Processing molecule {name} (ID: {molecule_id})")
            
            if not smiles:
                logger.warning(f"Cannot generate InChIKey for {name} (ID: {molecule_id}): Missing SMILES")
                failed_count += 1
                errors.append(f"Missing SMILES for {name} (ID: {molecule_id})")
                continue
            
            # Generate InChIKey from SMILES
            inchikey = generate_inchikey_from_smiles(smiles)
            
            if inchikey:
                # Update InChIKey
                with conn.cursor() as cursor:
                    cursor.execute("""
                    UPDATE molecules 
                    SET inchikey = %s,
                        updated_at = CURRENT_TIMESTAMP
                    WHERE id = %s
                    """, (inchikey, molecule_id))
                
                logger.info(f"Updated InChIKey for molecule {name} (ID: {molecule_id}) to {inchikey}")
                fixed_count += 1
                conn.commit()
            else:
                logger.warning(f"Could not generate InChIKey for {name} (ID: {molecule_id}, SMILES: {smiles})")
                failed_count += 1
                errors.append(f"Could not generate InChIKey for {name} (ID: {molecule_id})")
        
        except Exception as e:
            logger.error(f"Error fixing InChIKey for molecule {molecule.get('name', 'Unknown')}: {str(e)}")
            errors.append(f"Error for {molecule.get('name', 'Unknown')}: {str(e)}")
            failed_count += 1
            conn.rollback()
    
    return {
        "fixed": fixed_count,
        "failed": failed_count,
        "errors": errors
    }

def process_in_batches(molecules, batch_size, process_func, conn, *args):
    """
    Process molecules in batches to avoid long-running transactions.
    
    Args:
        molecules: List of molecules to process
        batch_size: Size of each batch
        process_func: Function to apply to each batch
        conn: Database connection
        *args: Additional arguments to pass to process_func
        
    Returns:
        Dict with combined results from all batches
    """
    total_molecules = len(molecules)
    total_batches = (total_molecules + batch_size - 1) // batch_size
    
    logger.info(f"Processing {total_molecules} molecules in {total_batches} batches of {batch_size}")
    
    combined_results = {
        "fixed": 0,
        "failed": 0,
        "errors": []
    }
    
    for batch_idx in range(total_batches):
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, total_molecules)
        
        batch_molecules = molecules[start_idx:end_idx]
        
        logger.info(f"Processing batch {batch_idx + 1}/{total_batches} ({len(batch_molecules)} molecules)")
        
        batch_results = process_func(conn, batch_molecules, *args)
        
        # Combine results
        combined_results["fixed"] += batch_results.get("fixed", 0)
        combined_results["failed"] += batch_results.get("failed", 0)
        combined_results["errors"].extend(batch_results.get("errors", []))
        
        # Save intermediate results
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        intermediate_file = f"reports/fix_batch_{batch_idx+1}_of_{total_batches}_{timestamp}.json"
        
        with open(intermediate_file, 'w') as f:
            json.dump({
                "batch": batch_idx + 1,
                "total_batches": total_batches,
                "results": batch_results
            }, f, indent=2)
        
        logger.info(f"Batch {batch_idx + 1} completed. Intermediate results saved to {intermediate_file}")
        
        # Give the database a short break between batches
        if batch_idx < total_batches - 1:
            time.sleep(1)
    
    return combined_results

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Fix data quality issues in ChEMBL molecules')
    parser.add_argument('--issues-file', default='reports/molecule_issues.json', help='Path to the issues report file')
    parser.add_argument('--fix-properties', action='store_true', help='Fix missing molecular properties')
    parser.add_argument('--fix-jsonb', action='store_true', help='Fix missing JSONB properties')
    parser.add_argument('--fix-pubchem', action='store_true', help='Fix missing PubChem CIDs')
    parser.add_argument('--fix-inchikey', action='store_true', help='Fix missing InChIKeys')
    parser.add_argument('--fix-all', action='store_true', help='Fix all issues')
    parser.add_argument('--limit', type=int, help='Limit the number of molecules to process')
    parser.add_argument('--batch-size', type=int, default=50, help='Number of molecules to process in each batch')
    
    args = parser.parse_args()
    
    # Load issues report
    issues = load_issues_report(args.issues_file)
    
    # Connect to database
    conn = connect_to_database()
    
    try:
        results = {}
        
        # Get property types
        property_types = get_property_types(conn)
        
        # Fix missing properties
        if args.fix_properties or args.fix_all:
            molecules = issues.get('molecules_with_incomplete_properties', [])
            if args.limit:
                molecules = molecules[:args.limit]
            results['fix_properties'] = process_in_batches(molecules, args.batch_size, fix_missing_properties, conn, property_types)
        
        # Fix missing JSONB properties
        if args.fix_jsonb or args.fix_all:
            molecules = issues.get('molecules_missing_jsonb_properties', [])
            if args.limit:
                molecules = molecules[:args.limit]
            results['fix_jsonb'] = process_in_batches(molecules, args.batch_size, fix_missing_jsonb_properties, conn)
        
        # Fix missing PubChem CIDs
        if args.fix_pubchem or args.fix_all:
            molecules = issues.get('molecules_missing_pubchem_cid', [])
            if args.limit:
                molecules = molecules[:args.limit]
            results['fix_pubchem'] = process_in_batches(molecules, args.batch_size, fix_missing_pubchem_cids, conn)
        
        # Fix missing InChIKeys
        if args.fix_inchikey or args.fix_all:
            molecules = issues.get('molecules_missing_inchikey', [])
            if args.limit:
                molecules = molecules[:args.limit]
            results['fix_inchikey'] = process_in_batches(molecules, args.batch_size, fix_missing_inchikeys, conn)
        
        # Save results
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        results_file = f"reports/fix_chembl_issues_{timestamp}.json"
        
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"Results saved to {results_file}")
        
        # Print summary
        print("\n=== ChEMBL Molecule Fix Summary ===")
        for fix_type, fix_results in results.items():
            print(f"{fix_type}: Fixed {fix_results.get('fixed', 0)}, Failed {fix_results.get('failed', 0)}")
        print(f"Detailed report saved to: {results_file}")
        print("====================================\n")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error in main function: {str(e)}")
        return 1
    
    finally:
        conn.close()
        logger.info("Database connection closed")

if __name__ == "__main__":
    sys.exit(main())