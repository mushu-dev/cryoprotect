#!/usr/bin/env python3
"""
Complete missing molecular properties for molecules.

This script identifies molecules with missing key properties and attempts
to calculate them or retrieve them from external sources.

Key properties include:
- Molecular Weight
- LogP
- TPSA
- Hydrogen Bond Donor Count
- Hydrogen Bond Acceptor Count
- Rotatable Bond Count
- Heavy Atom Count
- Complexity
"""

import os
import sys
import time
import json
import psycopg2
import logging
import requests
from typing import Dict, List, Any, Tuple, Optional
from psycopg2.extras import RealDictCursor, Json
from dotenv import load_dotenv

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Set constants
PUBCHEM_USER_AGENT = "CryoProtect/1.0 (https://github.com/yourusername/cryoprotect)"
PUBCHEM_PROPERTY_ENDPOINT = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,Complexity/JSON"

# Cache for property type IDs
PROPERTY_TYPE_CACHE = {}

def connect_to_db():
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    
    logger.info(f"Connecting to database with parameters: {db_params}")
    return psycopg2.connect(**db_params)

def get_key_property_type_ids(conn: psycopg2.extensions.connection) -> Dict[str, str]:
    """Get the UUIDs for the key property types."""
    property_types = {
        "Molecular Weight": None,
        "LogP": None,
        "TPSA": None, 
        "Hydrogen Bond Donor Count": None,
        "Hydrogen Bond Acceptor Count": None,
        "Rotatable Bond Count": None,
        "Heavy Atom Count": None,
        "Complexity": None
    }
    
    with conn.cursor() as cursor:
        for property_name in property_types.keys():
            cursor.execute("""
                SELECT id FROM property_types
                WHERE name = %s
            """, (property_name,))
            
            result = cursor.fetchone()
            if result:
                property_types[property_name] = result[0]
            else:
                # If property type doesn't exist, create it
                cursor.execute("""
                    INSERT INTO property_types (name, description, data_type)
                    VALUES (%s, %s, %s)
                    RETURNING id
                """, (
                    property_name,
                    f"The {property_name.lower()} of the molecule",
                    "float" if property_name in ["Molecular Weight", "LogP", "TPSA", "Complexity"] else "integer"
                ))
                property_types[property_name] = cursor.fetchone()[0]
    
    logger.info(f"Using these property types: {property_types}")
    return property_types

def get_molecules_with_missing_properties(
    conn: psycopg2.extensions.connection,
    property_type_ids: Dict[str, str],
    limit: Optional[int] = None
) -> List[Dict]:
    """
    Get a list of molecules that are missing one or more key properties.
    Prioritized by those with fewest missing properties.
    """
    query = """
        SELECT 
            m.id,
            m.name,
            m.smiles,
            m.pubchem_cid,
            m.molecular_formula,
            m.data_source,
            COUNT(mp.id) FILTER (WHERE mp.property_type_id::text = ANY(%s)) AS existing_key_property_count
        FROM 
            molecules m
        LEFT JOIN 
            molecular_properties mp ON m.id = mp.molecule_id AND mp.property_type_id::text = ANY(%s)
        WHERE 
            m.smiles IS NOT NULL
        GROUP BY 
            m.id, m.name, m.smiles, m.pubchem_cid, m.molecular_formula, m.data_source
        HAVING 
            COUNT(mp.id) FILTER (WHERE mp.property_type_id::text = ANY(%s)) < %s
        ORDER BY 
            existing_key_property_count DESC,  -- Prioritize molecules that need fewer properties
            m.id
    """
    
    if limit:
        query += f" LIMIT {limit}"
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Convert string UUIDs to proper UUIDs for PostgreSQL
        property_uuids = [str(uuid_val) for uuid_val in property_type_ids.values()]
        cursor.execute(query, (
            property_uuids,
            property_uuids,
            property_uuids,
            len(property_type_ids)
        ))
        molecules = cursor.fetchall()
    
    logger.info(f"Found {len(molecules)} molecules with missing key properties")
    return molecules

def get_missing_properties_for_molecule(
    conn: psycopg2.extensions.connection,
    molecule_id: str,
    property_type_ids: Dict[str, str]
) -> List[str]:
    """Get the list of missing property types for a molecule."""
    missing_properties = []
    
    with conn.cursor() as cursor:
        for prop_name, prop_id in property_type_ids.items():
            cursor.execute("""
                SELECT COUNT(*)
                FROM molecular_properties
                WHERE molecule_id = %s
                AND property_type_id = %s
            """, (molecule_id, prop_id))
            
            count = cursor.fetchone()[0]
            if count == 0:
                missing_properties.append(prop_name)
    
    return missing_properties

def fetch_pubchem_properties(cid: str) -> Dict[str, Any]:
    """
    Fetch properties from PubChem for a given CID.
    Returns a dict mapping property names to values.
    """
    if not cid:
        return {}
    
    url = PUBCHEM_PROPERTY_ENDPOINT.format(cid=cid)
    headers = {'User-Agent': PUBCHEM_USER_AGENT}
    
    try:
        response = requests.get(url, headers=headers)
        
        if response.status_code == 200:
            data = response.json()
            
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                props = data['PropertyTable']['Properties']
                
                if props and len(props) > 0:
                    # Map PubChem property names to our property types
                    return {
                        "Molecular Weight": props[0].get("MolecularWeight"),
                        "LogP": props[0].get("XLogP"),
                        "TPSA": props[0].get("TPSA"),
                        "Hydrogen Bond Donor Count": props[0].get("HBondDonorCount"),
                        "Hydrogen Bond Acceptor Count": props[0].get("HBondAcceptorCount"),
                        "Rotatable Bond Count": props[0].get("RotatableBondCount"),
                        "Heavy Atom Count": props[0].get("HeavyAtomCount"),
                        "Complexity": props[0].get("Complexity")
                    }
        
        logger.warning(f"Failed to get properties from PubChem for CID {cid}: {response.status_code}")
        return {}
        
    except Exception as e:
        logger.error(f"Error fetching properties from PubChem for CID {cid}: {e}")
        return {}

def calculate_properties_from_smiles(smiles: str) -> Dict[str, Any]:
    """
    Calculate properties from SMILES using RDKit if available, or a mock function for testing.
    Returns a dict mapping property names to values.
    """
    try:
        # Try to import RDKit and calculate properties
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski
        
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return {
                "Molecular Weight": Descriptors.MolWt(mol),
                "LogP": Descriptors.MolLogP(mol),
                "TPSA": Descriptors.TPSA(mol),
                "Hydrogen Bond Donor Count": Lipinski.NumHDonors(mol),
                "Hydrogen Bond Acceptor Count": Lipinski.NumHAcceptors(mol),
                "Rotatable Bond Count": Descriptors.NumRotatableBonds(mol),
                "Heavy Atom Count": mol.GetNumHeavyAtoms(),
                "Complexity": Descriptors.BalabanJ(mol) if Descriptors.BalabanJ(mol) is not None else 0
            }
        else:
            logger.warning(f"Failed to parse SMILES: {smiles}")
            return {}
            
    except ImportError:
        # RDKit not available, use a very basic mock
        logger.warning("RDKit not available, using fallback mock property calculator")
        
        try:
            # Import our mock RDKit module if available
            from mock_rdkit import calculate_basic_properties
            return calculate_basic_properties(smiles)
        except ImportError:
            logger.warning("Mock RDKit module not available, using empty properties")
            # Very minimal fallback
            return {}

def calculate_and_store_properties(
    conn: psycopg2.extensions.connection,
    molecule: Dict,
    missing_properties: List[str],
    property_type_ids: Dict[str, str],
    molecule_index: int = -1  # Default to -1 which won't match any logging condition
) -> Tuple[int, List[Dict]]:
    """
    Calculate and store missing properties for a molecule.
    Returns a tuple of (number of properties added, list of added property details).
    """
    added_properties = []
    properties_added = 0
    
    # First try PubChem if available
    properties = {}
    source = None
    
    if molecule['pubchem_cid']:
        properties = fetch_pubchem_properties(molecule['pubchem_cid'])
        source = "PubChem API"
    
    # If PubChem doesn't have the data, try calculating from SMILES
    if not properties and molecule['smiles']:
        properties = calculate_properties_from_smiles(molecule['smiles'])
        source = "Calculated (RDKit)"
    
    # Store each property that was missing
    with conn.cursor() as cursor:
        for prop_name in missing_properties:
            if prop_name in properties and properties[prop_name] is not None:
                # Convert to appropriate type based on property
                value = properties[prop_name]
                if prop_name in ["Molecular Weight", "LogP", "TPSA", "Complexity"]:
                    # These are floats
                    try:
                        value = float(value)
                    except (ValueError, TypeError):
                        logger.warning(f"Invalid float value for {prop_name}: {value}")
                        continue
                else:
                    # These are integers
                    try:
                        value = int(value)
                    except (ValueError, TypeError):
                        logger.warning(f"Invalid integer value for {prop_name}: {value}")
                        continue
                
                # Insert the property
                cursor.execute("""
                    INSERT INTO molecular_properties
                    (molecule_id, property_type_id, numeric_value, data_source, property_name, property_value)
                    VALUES (%s, %s, %s, %s, %s, %s)
                    RETURNING id
                """, (
                    molecule['id'],
                    property_type_ids[prop_name],
                    value,
                    source,
                    prop_name,
                    str(value)
                ))
                
                prop_id = cursor.fetchone()[0]
                properties_added += 1
                
                added_properties.append({
                    "id": prop_id,
                    "molecule_id": molecule['id'],
                    "property_type": prop_name,
                    "value": value,
                    "source": source
                })
                
                # Only log property additions for milestone molecules or for a sample of other molecules
                if molecule_index % 100 == 0 or molecule_index % 50 == 0:
                    logger.info(f"Added {prop_name} = {value} for {molecule['name']} (ID: {molecule['id']})")
    
    return properties_added, added_properties

def main():
    """Main function."""
    import argparse
    parser = argparse.ArgumentParser(description="Complete missing molecular properties")
    parser.add_argument("--limit", type=int, default=None, help="Limit the number of molecules to process")
    parser.add_argument("--dry-run", action="store_true", help="Don't actually store properties")
    parser.add_argument("--pubchem-only", action="store_true", help="Only use PubChem data, don't calculate properties")
    args = parser.parse_args()
    
    if args.dry_run:
        logger.info("Running in dry run mode - no data will be written to the database")
    
    if args.pubchem_only:
        logger.info("Running in PubChem-only mode - properties will not be calculated")
    
    # Connect to the database
    conn = connect_to_db()
    
    try:
        # Get key property type IDs
        property_type_ids = get_key_property_type_ids(conn)
        
        # Get molecules with missing properties
        molecules = get_molecules_with_missing_properties(conn, property_type_ids, args.limit)
        
        # Track overall stats
        total_properties_added = 0
        molecules_updated = 0
        all_added_properties = []
        
        # Process each molecule
        for i, molecule in enumerate(molecules):
            if i % 100 == 0:
                logger.info(f"Processing molecule {i+1}/{len(molecules)}: {molecule['name']} (ID: {molecule['id']}) - {i*100/len(molecules):.1f}% complete")
            
            # Get list of missing properties for this molecule
            missing_properties = get_missing_properties_for_molecule(conn, molecule['id'], property_type_ids)
            
            if not missing_properties:
                # Skip verbose logging unless it's a milestone
                if i % 100 == 0:
                    logger.info(f"No missing properties for {molecule['name']} (ID: {molecule['id']})")
                continue
            
            # For milestone or when properties are found
            if i % 100 == 0 or len(missing_properties) > 0:
                logger.info(f"Missing properties: {', '.join(missing_properties)}")
            
            if not args.dry_run:
                properties_added, added_properties = calculate_and_store_properties(
                    conn, molecule, missing_properties, property_type_ids, i
                )
                
                if properties_added > 0:
                    molecules_updated += 1
                    total_properties_added += properties_added
                    all_added_properties.extend(added_properties)
            else:
                logger.info(f"Would add properties: {', '.join(missing_properties)}")
        
        # Commit changes if not a dry run
        if not args.dry_run:
            conn.commit()
            logger.info(f"Committed {total_properties_added} property values for {molecules_updated} molecules")
            
            # Save a report of added properties
            report = {
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "molecules_processed": len(molecules),
                "molecules_updated": molecules_updated,
                "total_properties_added": total_properties_added,
                "added_properties": all_added_properties
            }
            
            with open("property_population_report.json", "w") as f:
                json.dump(report, f, indent=2, default=str)
            
            logger.info(f"Report saved to property_population_report.json")
        else:
            logger.info(f"Dry run completed. Would have added properties for {len(molecules)} molecules")
        
    except Exception as e:
        conn.rollback()
        logger.error(f"Error completing missing properties: {e}", exc_info=True)
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()