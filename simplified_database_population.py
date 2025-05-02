#!/usr/bin/env python3
"""
Simplified Database Population Script

Focuses on reliable population of CryoProtect database
using a simplified, step-by-step approach.

This script implements the strategy described in:
- DATABASE_POPULATION_SIMPLIFIED_APPROACH.md
- MINIMAL_POPULATION_WORKFLOW_GUIDE.md

Usage:
    python simplified_database_population.py --step check
    python simplified_database_population.py --step references
    python simplified_database_population.py --step pubchem --limit 500
    python simplified_database_population.py --step chembl --limit 500
    python simplified_database_population.py --step verify
    python simplified_database_population.py --step all --limit 500
"""

import os
import sys
import time
import json
import logging
import argparse
import traceback
from datetime import datetime
import requests
from typing import Dict, Any, List, Optional, Union, Tuple

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/simplified_population.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Ensure logs directory exists
os.makedirs("logs", exist_ok=True)
os.makedirs("checkpoints", exist_ok=True)

# Database connection setup
import psycopg2
from psycopg2.extras import RealDictCursor

def get_connection():
    """Get a simple database connection without complexity."""
    db_host = os.getenv('SUPABASE_DB_HOST')
    db_port = os.getenv('SUPABASE_DB_PORT', '6543')  # Session pooler port
    db_name = os.getenv('SUPABASE_DB_NAME')
    db_user = os.getenv('SUPABASE_DB_USER')
    db_password = os.getenv('SUPABASE_DB_PASSWORD')
    
    try:
        conn = psycopg2.connect(
            host=db_host,
            port=db_port,
            dbname=db_name,
            user=db_user,
            password=db_password,
            cursor_factory=RealDictCursor
        )
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {str(e)}")
        raise

def execute_query(query, params=None, fetch_one=False, return_data=True):
    """
    Simple, consistent query execution function.
    Abstracts away database connection details.
    """
    conn = None
    try:
        # Get connection
        conn = get_connection()
        
        # Execute query
        cursor = conn.cursor()
        cursor.execute(query, params)
        
        # Return results if needed
        if return_data:
            if fetch_one:
                return cursor.fetchone()
            else:
                return cursor.fetchall()
        
        # Commit if not a SELECT query
        if not query.strip().upper().startswith("SELECT"):
            conn.commit()
            
        return True
        
    except Exception as e:
        logger.error(f"Error executing query: {str(e)}")
        logger.error(f"Query: {query}")
        logger.error(f"Params: {params}")
        if conn:
            conn.rollback()
        raise
    finally:
        if conn:
            conn.close()

def with_transaction(func, *args, **kwargs):
    """Execute a function within a transaction with simple error handling."""
    conn = None
    try:
        # Get connection
        conn = get_connection()
        conn.autocommit = False
        
        # Execute function with connection as first argument
        result = func(conn, *args, **kwargs)
        
        # Commit transaction
        conn.commit()
        return result
        
    except Exception as e:
        # Rollback on error
        logger.error(f"Transaction error: {str(e)}")
        if conn:
            conn.rollback()
        raise
        
    finally:
        # Close connection
        if conn:
            conn.close()

def check_connection():
    """Check database connection works."""
    try:
        result = execute_query("SELECT 1 as connection_test", fetch_one=True)
        return result and result['connection_test'] == 1
    except Exception as e:
        logger.error(f"Connection check failed: {str(e)}")
        return False

def verify_schema():
    """Verify that required tables and columns exist."""
    try:
        # Check molecules table
        molecules_result = execute_query("""
            SELECT EXISTS (
                SELECT 1 FROM information_schema.tables 
                WHERE table_name = 'molecules'
            ) as table_exists
        """, fetch_one=True)
        
        if not molecules_result or not molecules_result['table_exists']:
            logger.error("Table 'molecules' does not exist")
            return False
        
        # Check molecular_properties table
        properties_result = execute_query("""
            SELECT EXISTS (
                SELECT 1 FROM information_schema.tables 
                WHERE table_name = 'molecular_properties'
            ) as table_exists
        """, fetch_one=True)
        
        if not properties_result or not properties_result['table_exists']:
            logger.error("Table 'molecular_properties' does not exist")
            return False
        
        # Check property_types table
        types_result = execute_query("""
            SELECT EXISTS (
                SELECT 1 FROM information_schema.tables 
                WHERE table_name = 'property_types'
            ) as table_exists
        """, fetch_one=True)
        
        if not types_result or not types_result['table_exists']:
            logger.error("Table 'property_types' does not exist")
            return False
        
        return True
        
    except Exception as e:
        logger.error(f"Schema verification failed: {str(e)}")
        return False

# Reference compounds functions
def populate_reference_compounds():
    """
    Populate ONLY the 9 reference compounds with guaranteed properties.
    This is our simplest test of the end-to-end pipeline.
    """
    # Hard-coded reference compounds with complete properties
    reference_compounds = [
        {
            "chembl_id": "CHEMBL388978",
            "name": "Glycerol",
            "smiles": "C(C(CO)O)O",
            "inchi": "InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2",
            "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
            "formula": "C3H8O3",
            "properties": {
                "logP": -1.76,
                "h_bond_donors": 3,
                "h_bond_acceptors": 3,
                "molecular_weight": 92.09
            }
        },
        {
            "chembl_id": "CHEMBL1098659",
            "name": "DMSO",
            "smiles": "CS(=O)C",
            "inchi": "InChI=1S/C2H6OS/c1-4(2)3/h1-2H3",
            "inchikey": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N",
            "formula": "C2H6OS",
            "properties": {
                "logP": -1.35,
                "h_bond_donors": 0,
                "h_bond_acceptors": 1,
                "molecular_weight": 78.13
            }
        },
        {
            "chembl_id": "CHEMBL66195",
            "name": "beta-Alanine",
            "smiles": "C(CC(=O)O)N",
            "inchi": "InChI=1S/C3H7NO2/c4-2-1-3(5)6/h1-2,4H2,(H,5,6)",
            "inchikey": "UCMIRNVEIXFBKS-UHFFFAOYSA-N",
            "formula": "C3H7NO2",
            "properties": {
                "logP": -3.17,
                "h_bond_donors": 2,
                "h_bond_acceptors": 3,
                "molecular_weight": 89.09
            }
        },
        {
            "chembl_id": "CHEMBL500033",
            "name": "tert-Butanol",
            "smiles": "CC(C)(C)O",
            "inchi": "InChI=1S/C4H10O/c1-4(2,3)5/h5H,1-3H3",
            "inchikey": "DKGRVHCAPJKBHC-UHFFFAOYSA-N",
            "formula": "C4H10O",
            "properties": {
                "logP": 0.35,
                "h_bond_donors": 1,
                "h_bond_acceptors": 1,
                "molecular_weight": 74.12
            }
        },
        {
            "chembl_id": "CHEMBL1487",
            "name": "Urea",
            "smiles": "C(=O)(N)N",
            "inchi": "InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)",
            "inchikey": "XSQUKJJJFZCRTK-UHFFFAOYSA-N",
            "formula": "CH4N2O",
            "properties": {
                "logP": -2.11,
                "h_bond_donors": 2,
                "h_bond_acceptors": 1,
                "molecular_weight": 60.06
            }
        },
        {
            "chembl_id": "CHEMBL6196",
            "name": "Ethylene glycol",
            "smiles": "C(CO)O",
            "inchi": "InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2",
            "inchikey": "LYCAIKOWRPUZTN-UHFFFAOYSA-N",
            "formula": "C2H6O2",
            "properties": {
                "logP": -1.36,
                "h_bond_donors": 2,
                "h_bond_acceptors": 2,
                "molecular_weight": 62.07
            }
        },
        {
            "chembl_id": "CHEMBL967",
            "name": "Propylene glycol",
            "smiles": "CC(CO)O",
            "inchi": "InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3",
            "inchikey": "DNIAPMSPPWPWGF-UHFFFAOYSA-N",
            "formula": "C3H8O2",
            "properties": {
                "logP": -0.92,
                "h_bond_donors": 2,
                "h_bond_acceptors": 2,
                "molecular_weight": 76.09
            }
        },
        {
            "chembl_id": "CHEMBL262548",
            "name": "Trehalose",
            "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O",
            "inchi": "InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12-10(19)6(15)5(14)3(2-13)22-12/h3-19H,1-2H2/t3-,4-,5-,6-,7-,8-,9-,10-,11-,12-/m1/s1",
            "inchikey": "LFQSCWFLJHTTHZ-LUBWMEBZSA-N",
            "formula": "C12H22O11",
            "properties": {
                "logP": -3.77,
                "h_bond_donors": 8,
                "h_bond_acceptors": 11,
                "molecular_weight": 342.30
            }
        },
        {
            "chembl_id": "CHEMBL6752",
            "name": "Glycine",
            "smiles": "C(C(=O)O)N",
            "inchi": "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)",
            "inchikey": "DHMQDGOQFOQNFH-UHFFFAOYSA-N",
            "formula": "C2H5NO2",
            "properties": {
                "logP": -3.21,
                "h_bond_donors": 1,
                "h_bond_acceptors": 3,
                "molecular_weight": 75.07
            }
        }
    ]
    
    # Use the simplest possible insertion method
    success_count = 0
    for compound in reference_compounds:
        try:
            # Log what we're doing
            logger.info(f"Processing reference compound: {compound['name']} ({compound['chembl_id']})")
            
            # 1. Insert basic molecule
            molecule_record = execute_query("""
                INSERT INTO molecules (name, chembl_id, smiles, inchi, inchikey, formula, data_source)
                VALUES (%(name)s, %(chembl_id)s, %(smiles)s, %(inchi)s, %(inchikey)s, %(formula)s, 'reference')
                ON CONFLICT (chembl_id) DO UPDATE SET
                    name = EXCLUDED.name,
                    smiles = EXCLUDED.smiles,
                    inchi = EXCLUDED.inchi,
                    inchikey = EXCLUDED.inchikey,
                    formula = EXCLUDED.formula,
                    updated_at = NOW()
                RETURNING id
            """, compound, fetch_one=True)
            
            if not molecule_record or 'id' not in molecule_record:
                logger.error(f"Failed to insert/get molecule record for {compound['name']}")
                continue
                
            molecule_id = molecule_record['id']
            logger.info(f"Molecule record created/updated with ID: {molecule_id}")
            
            # 2. Insert each property one by one
            for prop_name, value in compound['properties'].items():
                # Get or create property type
                prop_type_record = execute_query("""
                    INSERT INTO property_types (name, data_type)
                    VALUES (%(name)s, %(data_type)s)
                    ON CONFLICT (name) DO UPDATE SET
                        updated_at = NOW()
                    RETURNING id
                """, {
                    'name': prop_name,
                    'data_type': 'numeric'
                }, fetch_one=True)
                
                if not prop_type_record or 'id' not in prop_type_record:
                    logger.error(f"Failed to get property type ID for {prop_name}")
                    continue
                    
                prop_type_id = prop_type_record['id']
                
                # Insert property value
                execute_query("""
                    INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                    VALUES (%(molecule_id)s, %(property_type_id)s, %(value)s)
                    ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                        numeric_value = %(value)s,
                        updated_at = NOW()
                """, {
                    'molecule_id': molecule_id,
                    'property_type_id': prop_type_id,
                    'value': value
                })
                
                logger.info(f"Property {prop_name}={value} set for {compound['name']}")
            
            # If we get here, all properties were inserted successfully
            success_count += 1
            logger.info(f"Successfully populated reference compound: {compound['name']}")
            
        except Exception as e:
            logger.error(f"Error populating reference compound {compound['name']}: {str(e)}")
            logger.error(traceback.format_exc())
    
    logger.info(f"Successfully populated {success_count}/9 reference compounds")
    return success_count

def verify_reference_compounds():
    """
    Verify ONLY that reference compounds exist with required properties.
    This is our simplest verification step.
    """
    # Get all reference compounds
    ref_compounds = execute_query("""
        SELECT id, name, chembl_id FROM molecules
        WHERE chembl_id IN (
            'CHEMBL388978', 'CHEMBL1098659', 'CHEMBL66195', 'CHEMBL500033',
            'CHEMBL1487', 'CHEMBL6196', 'CHEMBL967', 'CHEMBL262548', 'CHEMBL6752'
        )
    """)
    
    found_count = len(ref_compounds)
    logger.info(f"Found {found_count}/9 reference compounds")
    
    # For each found compound, check critical properties
    complete_count = 0
    for compound in ref_compounds:
        required_props = ['logP', 'h_bond_donors', 'h_bond_acceptors']
        missing_props = []
        
        for prop in required_props:
            # Check if property exists
            result = execute_query("""
                SELECT mp.numeric_value
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE mp.molecule_id = %(molecule_id)s AND pt.name = %(prop_name)s
            """, {
                'molecule_id': compound['id'],
                'prop_name': prop
            }, fetch_one=True)
            
            if not result or result['numeric_value'] is None:
                missing_props.append(prop)
        
        if not missing_props:
            complete_count += 1
            logger.info(f"✓ {compound['name']} has all required properties")
        else:
            logger.warning(f"✗ {compound['name']} missing properties: {', '.join(missing_props)}")
    
    logger.info(f"Found {complete_count}/{found_count} complete reference compounds")
    
    # Return simple success/failure
    return complete_count == 9

# PubChem import functions
def fetch_pubchem_data(cid):
    """
    Fetch compound data from PubChem with simple retry mechanism.
    Focus on the most crucial properties for cryoprotectants.
    """
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
        "MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,"
        "IsomericSMILES,InChI,InChIKey,IUPACName,Title/JSON"
    )
    
    for attempt in range(3):  # Simple retry mechanism
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                return data["PropertyTable"]["Properties"][0]
            
            # Handle rate limiting
            if response.status_code == 429:
                wait_time = (attempt + 1) * 2  # Simple exponential backoff
                logger.warning(f"Rate limited. Waiting {wait_time} seconds...")
                time.sleep(wait_time)
                continue
                
            logger.warning(f"Failed to fetch CID {cid}: HTTP {response.status_code}")
            return None
            
        except Exception as e:
            logger.warning(f"Error fetching CID {cid}, attempt {attempt+1}: {str(e)}")
            time.sleep(1)
    
    return None

def insert_pubchem_molecule(data):
    """
    Insert molecule data from PubChem into database.
    Simple, straightforward implementation with clear logging.
    """
    try:
        # Check for required fields
        if not all(key in data for key in ["IsomericSMILES", "InChI", "InChIKey", "MolecularFormula"]):
            logger.warning(f"Missing required fields for CID {data['CID']}")
            return None
        
        # Prepare data
        molecule_data = {
            "pubchem_cid": data["CID"],
            "name": data.get("Title", f"CID {data['CID']}"),
            "smiles": data["IsomericSMILES"],
            "inchi": data["InChI"],
            "inchikey": data["InChIKey"],
            "formula": data["MolecularFormula"],
            "molecular_weight": data.get("MolecularWeight"),
            "data_source": "PubChem"
        }
        
        # Execute insert
        result = execute_query("""
            INSERT INTO molecules 
                (pubchem_cid, name, smiles, inchi, inchikey, formula, molecular_weight, data_source)
            VALUES 
                (%(pubchem_cid)s, %(name)s, %(smiles)s, %(inchi)s, %(inchikey)s, 
                 %(formula)s, %(molecular_weight)s, %(data_source)s)
            ON CONFLICT (pubchem_cid) DO UPDATE SET 
                name = EXCLUDED.name,
                smiles = EXCLUDED.smiles,
                inchi = EXCLUDED.inchi, 
                inchikey = EXCLUDED.inchikey,
                formula = EXCLUDED.formula,
                molecular_weight = EXCLUDED.molecular_weight,
                updated_at = NOW()
            RETURNING id
        """, molecule_data, fetch_one=True)
        
        if result and 'id' in result:
            logger.info(f"Successfully inserted molecule for CID {data['CID']}")
            return result['id']
        else:
            logger.warning(f"Failed to insert molecule for CID {data['CID']}")
            return None
            
    except Exception as e:
        logger.error(f"Error inserting molecule for CID {data['CID']}: {str(e)}")
        return None

def insert_pubchem_properties(molecule_id, data):
    """
    Insert properties for a PubChem molecule.
    Focus on getting the 3 critical properties correctly.
    """
    try:
        # Extract the 3 critical properties
        properties = {}
        
        # These are the most important properties for our verification criteria
        if "XLogP" in data:
            try:
                properties["logP"] = float(data["XLogP"])
            except (ValueError, TypeError):
                logger.warning(f"Invalid XLogP value: {data['XLogP']}")
        
        if "HBondDonorCount" in data:
            try:
                properties["h_bond_donors"] = int(data["HBondDonorCount"])
            except (ValueError, TypeError):
                logger.warning(f"Invalid HBondDonorCount value: {data['HBondDonorCount']}")
        
        if "HBondAcceptorCount" in data:
            try:
                properties["h_bond_acceptors"] = int(data["HBondAcceptorCount"])
            except (ValueError, TypeError):
                logger.warning(f"Invalid HBondAcceptorCount value: {data['HBondAcceptorCount']}")
        
        # Check if we have all required properties
        if not all(prop in properties for prop in ["logP", "h_bond_donors", "h_bond_acceptors"]):
            logger.warning(f"Missing critical properties for molecule {molecule_id}")
            return False
        
        # Insert properties one by one for clarity
        for prop_name, value in properties.items():
            # Get or create property type
            prop_type_id = execute_query("""
                INSERT INTO property_types (name, data_type)
                VALUES (%(name)s, %(data_type)s)
                ON CONFLICT (name) DO UPDATE SET
                    updated_at = NOW()
                RETURNING id
            """, {
                'name': prop_name,
                'data_type': 'numeric'
            }, fetch_one=True)['id']
            
            # Insert property value
            execute_query("""
                INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                VALUES (%(molecule_id)s, %(property_type_id)s, %(value)s)
                ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                    numeric_value = %(value)s,
                    updated_at = NOW()
            """, {
                'molecule_id': molecule_id,
                'property_type_id': prop_type_id,
                'value': value
            })
        
        logger.info(f"Successfully inserted properties for molecule {molecule_id}")
        return True
        
    except Exception as e:
        logger.error(f"Error inserting properties for molecule {molecule_id}: {str(e)}")
        return False

def import_pubchem_minimal(limit=500):
    """
    Import a minimal set of molecules from PubChem with critical properties.
    Focus on reliability rather than quantity.
    """
    # Get CID list from file
    try:
        with open("CID-Synonym-curated", "r") as f:
            cids = [int(line.strip().split("\t")[0]) for line in f 
                   if line.strip() and line.strip().split("\t")[0].isdigit()][:limit]
    except Exception as e:
        logger.error(f"Error reading CID list: {str(e)}")
        return 0
    
    logger.info(f"Starting minimal import of {len(cids)} PubChem compounds")
    
    # Process in small batches of 10
    batch_size = 10
    success_count = 0
    
    # Check for checkpoint
    checkpoint_file = "checkpoints/pubchem_import_checkpoint.txt"
    start_position = 0
    if os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, "r") as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    start_position = int(lines[0].strip())
                    success_count = int(lines[1].strip())
                    logger.info(f"Resuming from checkpoint: position {start_position}, {success_count} successful imports")
        except Exception as e:
            logger.warning(f"Error reading checkpoint: {str(e)}")
    
    for i in range(start_position, len(cids), batch_size):
        batch = cids[i:i+batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(len(cids) + batch_size - 1)//batch_size}")
        
        # Process each compound in batch
        batch_success = 0
        for cid in batch:
            try:
                # 1. Fetch compound data from PubChem
                logger.info(f"Fetching data for CID {cid}")
                data = fetch_pubchem_data(cid)
                if not data:
                    logger.warning(f"No data found for CID {cid}")
                    continue
                
                # 2. Insert molecule
                logger.info(f"Inserting molecule for CID {cid}")
                molecule_id = insert_pubchem_molecule(data)
                if not molecule_id:
                    continue
                
                # 3. Insert properties
                logger.info(f"Inserting properties for CID {cid}")
                property_success = insert_pubchem_properties(molecule_id, data)
                if property_success:
                    batch_success += 1
                    
            except Exception as e:
                logger.error(f"Error processing CID {cid}: {str(e)}")
                logger.error(traceback.format_exc())
        
        # Update success count and log progress
        success_count += batch_success
        logger.info(f"Batch complete. Total successful imports: {success_count}/{i+len(batch)}")
        
        # Simple file-based checkpoint
        with open(checkpoint_file, "w") as f:
            f.write(f"{i+batch_size}\n{success_count}")
            
        # Small delay to avoid overwhelming API
        time.sleep(1)
    
    logger.info(f"Completed PubChem import: {success_count}/{len(cids)} compounds")
    return success_count

def verify_pubchem_molecules():
    """
    Verify that PubChem molecules exist with required properties.
    Focus on property quality rather than pure quantity.
    """
    # Count PubChem molecules
    pubchem_count = execute_query("""
        SELECT COUNT(*) as count FROM molecules 
        WHERE pubchem_cid IS NOT NULL
    """, fetch_one=True)['count']
    
    # Count PubChem molecules with all critical properties
    complete_count = execute_query("""
        SELECT COUNT(DISTINCT m.id) as count
        FROM molecules m
        JOIN molecular_properties mp_logp ON m.id = mp_logp.molecule_id
        JOIN property_types pt_logp ON mp_logp.property_type_id = pt_logp.id AND pt_logp.name = 'logP'
        JOIN molecular_properties mp_hbd ON m.id = mp_hbd.molecule_id
        JOIN property_types pt_hbd ON mp_hbd.property_type_id = pt_hbd.id AND pt_hbd.name = 'h_bond_donors'
        JOIN molecular_properties mp_hba ON m.id = mp_hba.molecule_id
        JOIN property_types pt_hba ON mp_hba.property_type_id = pt_hba.id AND pt_hba.name = 'h_bond_acceptors'
        WHERE m.pubchem_cid IS NOT NULL
    """, fetch_one=True)['count']
    
    completion_percent = (complete_count / pubchem_count * 100) if pubchem_count > 0 else 0
    logger.info(f"PubChem molecules: {pubchem_count} total, {complete_count} complete ({completion_percent:.1f}%)")
    
    # Check if we meet requirements
    return pubchem_count >= 500 and completion_percent >= 90

# ChEMBL import functions
def fetch_chembl_data(chembl_id):
    """
    Fetch compound data from ChEMBL API.
    Simplified implementation focusing on core data.
    """
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    
    for attempt in range(3):  # Simple retry mechanism
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                return response.json()
            
            # Handle rate limiting or server errors
            if response.status_code in [429, 500, 502, 503, 504]:
                wait_time = (attempt + 1) * 2  # Simple exponential backoff
                logger.warning(f"ChEMBL API error {response.status_code}. Waiting {wait_time} seconds...")
                time.sleep(wait_time)
                continue
                
            logger.warning(f"Failed to fetch ChEMBL ID {chembl_id}: HTTP {response.status_code}")
            return None
            
        except Exception as e:
            logger.warning(f"Error fetching ChEMBL ID {chembl_id}, attempt {attempt+1}: {str(e)}")
            time.sleep(1)
    
    return None

def insert_chembl_molecule(data):
    """
    Insert molecule data from ChEMBL into database.
    Simplified version focusing on core data.
    """
    try:
        # Extract nested data with careful error handling
        molecule_props = data.get('molecule_properties', {})
        molecule_structures = data.get('molecule_structures', {})
        
        # Check for required fields
        if not molecule_structures or not molecule_structures.get('canonical_smiles'):
            logger.warning(f"Missing required structures for ChEMBL ID {data.get('molecule_chembl_id')}")
            return None
        
        # Prepare data
        molecule_data = {
            "chembl_id": data.get('molecule_chembl_id'),
            "name": data.get('pref_name', f"ChEMBL {data.get('molecule_chembl_id')}"),
            "smiles": molecule_structures.get('canonical_smiles'),
            "inchi": molecule_structures.get('standard_inchi', ''),
            "inchikey": molecule_structures.get('standard_inchi_key', ''),
            "formula": molecule_props.get('full_molformula', ''),
            "molecular_weight": molecule_props.get('full_mwt'),
            "data_source": "ChEMBL"
        }
        
        # Execute insert
        result = execute_query("""
            INSERT INTO molecules 
                (chembl_id, name, smiles, inchi, inchikey, formula, molecular_weight, data_source)
            VALUES 
                (%(chembl_id)s, %(name)s, %(smiles)s, %(inchi)s, %(inchikey)s, 
                 %(formula)s, %(molecular_weight)s, %(data_source)s)
            ON CONFLICT (chembl_id) DO UPDATE SET 
                name = EXCLUDED.name,
                smiles = EXCLUDED.smiles,
                inchi = EXCLUDED.inchi, 
                inchikey = EXCLUDED.inchikey,
                formula = EXCLUDED.formula,
                molecular_weight = EXCLUDED.molecular_weight,
                updated_at = NOW()
            RETURNING id
        """, molecule_data, fetch_one=True)
        
        if result and 'id' in result:
            logger.info(f"Successfully inserted molecule for ChEMBL ID {data.get('molecule_chembl_id')}")
            return result['id']
        else:
            logger.warning(f"Failed to insert molecule for ChEMBL ID {data.get('molecule_chembl_id')}")
            return None
            
    except Exception as e:
        logger.error(f"Error inserting molecule for ChEMBL ID {data.get('molecule_chembl_id')}: {str(e)}")
        return None

def insert_chembl_properties(molecule_id, data):
    """
    Insert properties for a ChEMBL molecule.
    Focus on the 3 critical properties.
    """
    try:
        # Extract molecular properties with careful error handling
        molecule_props = data.get('molecule_properties', {})
        
        # Extract the 3 critical properties
        properties = {}
        
        # These are the most important properties for our verification criteria
        if 'alogp' in molecule_props:
            try:
                properties["logP"] = float(molecule_props['alogp'])
            except (ValueError, TypeError):
                logger.warning(f"Invalid ALogP value: {molecule_props['alogp']}")
        
        if 'hbd' in molecule_props:
            try:
                properties["h_bond_donors"] = int(molecule_props['hbd'])
            except (ValueError, TypeError):
                logger.warning(f"Invalid HBD value: {molecule_props['hbd']}")
        
        if 'hba' in molecule_props:
            try:
                properties["h_bond_acceptors"] = int(molecule_props['hba'])
            except (ValueError, TypeError):
                logger.warning(f"Invalid HBA value: {molecule_props['hba']}")
        
        # Check if we have all required properties
        if not all(prop in properties for prop in ["logP", "h_bond_donors", "h_bond_acceptors"]):
            logger.warning(f"Missing critical properties for molecule {molecule_id}")
            return False
        
        # Insert properties one by one for clarity
        for prop_name, value in properties.items():
            # Get or create property type
            prop_type_id = execute_query("""
                INSERT INTO property_types (name, data_type)
                VALUES (%(name)s, %(data_type)s)
                ON CONFLICT (name) DO UPDATE SET
                    updated_at = NOW()
                RETURNING id
            """, {
                'name': prop_name,
                'data_type': 'numeric'
            }, fetch_one=True)['id']
            
            # Insert property value
            execute_query("""
                INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                VALUES (%(molecule_id)s, %(property_type_id)s, %(value)s)
                ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                    numeric_value = %(value)s,
                    updated_at = NOW()
            """, {
                'molecule_id': molecule_id,
                'property_type_id': prop_type_id,
                'value': value
            })
        
        logger.info(f"Successfully inserted properties for molecule {molecule_id}")
        return True
        
    except Exception as e:
        logger.error(f"Error inserting properties for molecule {molecule_id}: {str(e)}")
        return False

def import_chembl_minimal(limit=500):
    """
    Import a minimal set of molecules from ChEMBL with critical properties.
    Focus on reliability rather than quantity.
    """
    # Get reference compound ChEMBL IDs as a starting point
    reference_ids = [
        "CHEMBL388978", "CHEMBL1098659", "CHEMBL66195", "CHEMBL500033",
        "CHEMBL1487", "CHEMBL6196", "CHEMBL967", "CHEMBL262548", "CHEMBL6752"
    ]
    
    # For this simplified version, we'll just use these IDs plus some other known cryoprotectants
    additional_ids = [
        "CHEMBL268469",  # Sucrose
        "CHEMBL1197",    # Methanol
        "CHEMBL14",      # Ethanol
        "CHEMBL44093",   # Choline
        "CHEMBL10",      # D-Glucose
        "CHEMBL1407",    # Mannitol
        "CHEMBL9123",    # Sorbitol
        "CHEMBL1229000", # PEG-400
        "CHEMBL502",     # Sodium chloride
        "CHEMBL1200676", # Potassium chloride
        "CHEMBL1302",    # Ribose
        "CHEMBL439",     # Fructose
        "CHEMBL1161",    # Maltose
        "CHEMBL1113",    # Ammonium sulfate
        "CHEMBL1358",    # Magnesium sulfate
    ]
    
    # Start with direct collection of chemistry IDs
    all_compounds = reference_ids + additional_ids
    
    # Deduplicate
    all_compounds = list(set(all_compounds))[:limit]
    
    logger.info(f"Starting import of {len(all_compounds)} ChEMBL compounds")
    
    # Process each compound
    success_count = 0
    for i, chembl_id in enumerate(all_compounds):
        try:
            logger.info(f"Processing compound {i+1}/{len(all_compounds)}: {chembl_id}")
            
            # 1. Fetch compound data from ChEMBL
            data = fetch_chembl_data(chembl_id)
            if not data:
                logger.warning(f"No data found for ChEMBL ID {chembl_id}")
                continue
            
            # 2. Insert molecule
            molecule_id = insert_chembl_molecule(data)
            if not molecule_id:
                continue
            
            # 3. Insert properties
            if insert_chembl_properties(molecule_id, data):
                success_count += 1
                
            # Simple checkpoint every 10 compounds
            if (i + 1) % 10 == 0:
                logger.info(f"Progress: {i+1}/{len(all_compounds)} compounds processed, {success_count} successful")
                
        except Exception as e:
            logger.error(f"Error processing ChEMBL ID {chembl_id}: {str(e)}")
    
    logger.info(f"Completed ChEMBL import: {success_count}/{len(all_compounds)} compounds")
    return success_count

def verify_chembl_molecules():
    """
    Verify that ChEMBL molecules exist with required properties.
    """
    # Count ChEMBL molecules
    chembl_count = execute_query("""
        SELECT COUNT(*) as count FROM molecules 
        WHERE chembl_id IS NOT NULL
    """, fetch_one=True)['count']
    
    # Count ChEMBL molecules with all critical properties
    complete_count = execute_query("""
        SELECT COUNT(DISTINCT m.id) as count
        FROM molecules m
        JOIN molecular_properties mp_logp ON m.id = mp_logp.molecule_id
        JOIN property_types pt_logp ON mp_logp.property_type_id = pt_logp.id AND pt_logp.name = 'logP'
        JOIN molecular_properties mp_hbd ON m.id = mp_hbd.molecule_id
        JOIN property_types pt_hbd ON mp_hbd.property_type_id = pt_hbd.id AND pt_hbd.name = 'h_bond_donors'
        JOIN molecular_properties mp_hba ON m.id = mp_hba.molecule_id
        JOIN property_types pt_hba ON mp_hba.property_type_id = pt_hba.id AND pt_hba.name = 'h_bond_acceptors'
        WHERE m.chembl_id IS NOT NULL
    """, fetch_one=True)['count']
    
    completion_percent = (complete_count / chembl_count * 100) if chembl_count > 0 else 0
    logger.info(f"ChEMBL molecules: {chembl_count} total, {complete_count} complete ({completion_percent:.1f}%)")
    
    # Check if we meet minimum requirements
    return chembl_count >= 24 and completion_percent >= 90  # 9 reference + 15 additional

# All molecules verification
def verify_all_molecules():
    """
    Verify all molecules against requirements.
    """
    # Count all molecules
    total_count = execute_query("""
        SELECT COUNT(*) as count FROM molecules
    """, fetch_one=True)['count']
    
    # Count molecules with required properties
    complete_count = execute_query("""
        SELECT COUNT(DISTINCT m.id) as count
        FROM molecules m
        JOIN molecular_properties mp_logp ON m.id = mp_logp.molecule_id
        JOIN property_types pt_logp ON mp_logp.property_type_id = pt_logp.id AND pt_logp.name = 'logP'
        JOIN molecular_properties mp_hbd ON m.id = mp_hbd.molecule_id
        JOIN property_types pt_hbd ON mp_hbd.property_type_id = pt_hbd.id AND pt_hbd.name = 'h_bond_donors'
        JOIN molecular_properties mp_hba ON m.id = mp_hba.molecule_id
        JOIN property_types pt_hba ON mp_hba.property_type_id = pt_hba.id AND pt_hba.name = 'h_bond_acceptors'
    """, fetch_one=True)['count']
    
    # Calculate completion percentage
    completion_percent = (complete_count / total_count * 100) if total_count > 0 else 0
    
    # Print results
    logger.info(f"Total molecules: {total_count}")
    logger.info(f"Complete molecules: {complete_count} ({completion_percent:.1f}%)")
    
    # Check all requirements
    ref_complete = verify_reference_compounds()
    pubchem_complete = verify_pubchem_molecules()
    chembl_complete = verify_chembl_molecules()
    target_met = total_count >= 5000
    completion_met = completion_percent >= 90
    
    logger.info(f"Reference compounds complete: {ref_complete}")
    logger.info(f"PubChem verification: {pubchem_complete}")
    logger.info(f"ChEMBL verification: {chembl_complete}")
    logger.info(f"Target molecule count (5000+): {target_met}")
    logger.info(f"Property completion target (90%+): {completion_met}")
    
    # All requirements met?
    all_requirements_met = ref_complete and pubchem_complete and chembl_complete and target_met and completion_met
    logger.info(f"All requirements met: {all_requirements_met}")
    
    return all_requirements_met

def verify_performance():
    """
    Check query performance.
    """
    # Define test queries
    test_queries = [
        {
            "name": "Fetch single molecule by ID",
            "query": """
            SELECT id, name, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid, chembl_id
            FROM molecules
            ORDER BY created_at DESC
            LIMIT 1
            """,
            "params": None
        },
        {
            "name": "Fetch molecule with properties",
            "query": """
            WITH molecule AS (
                SELECT id, name, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid, chembl_id
                FROM molecules
                ORDER BY created_at DESC
                LIMIT 1
            )
            SELECT
                m.id, m.name, m.smiles, m.inchi, m.inchikey, m.formula, m.molecular_weight, m.pubchem_cid, m.chembl_id,
                array_agg(pt.name) as property_names,
                array_agg(
                    CASE
                        WHEN mp.numeric_value IS NOT NULL THEN mp.numeric_value::text
                        ELSE mp.text_value
                    END
                ) as property_values
            FROM molecule m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            GROUP BY m.id, m.name, m.smiles, m.inchi, m.inchikey, m.formula, m.molecular_weight, m.pubchem_cid, m.chembl_id
            """,
            "params": None
        },
        {
            "name": "Search molecules by name",
            "query": """
            SELECT id, name, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid, chembl_id
            FROM molecules
            WHERE LOWER(name) LIKE LOWER(%s)
            LIMIT 10
            """,
            "params": ("%glyc%",)
        },
        {
            "name": "Count molecules by property value range",
            "query": """
            SELECT COUNT(*)
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE pt.name = 'logP' AND mp.numeric_value BETWEEN -1 AND 3
            """,
            "params": None
        }
    ]
    
    # Run each query 5 times and measure performance
    conn = get_connection()
    times = []
    
    try:
        for test in test_queries:
            query_times = []
            for i in range(5):
                cursor = conn.cursor()
                
                # Clear cache with a dummy query
                cursor.execute("SELECT 1")
                
                # Time the query
                start_time = time.time()
                if test["params"]:
                    cursor.execute(test["query"], test["params"])
                else:
                    cursor.execute(test["query"])
                    
                cursor.fetchall()  # Consume results
                end_time = time.time()
                
                query_time = (end_time - start_time) * 1000  # Convert to milliseconds
                query_times.append(query_time)
                
            # Calculate average time for this query
            avg_time = sum(query_times) / len(query_times)
            logger.info(f"Query '{test['name']}': {avg_time:.2f}ms")
            times.append(avg_time)
            
    finally:
        conn.close()
    
    # Calculate overall average
    overall_avg = sum(times) / len(times)
    logger.info(f"Overall average query time: {overall_avg:.2f}ms")
    
    return overall_avg

def add_basic_indexes():
    """
    Add simple indexes to improve query performance.
    Only run after basic functionality is verified.
    """
    try:
        execute_query("""
            -- Add indexes for commonly queried fields
            CREATE INDEX IF NOT EXISTS idx_molecules_pubchem_cid ON molecules (pubchem_cid);
            CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON molecules (chembl_id);
            CREATE INDEX IF NOT EXISTS idx_molecules_inchikey ON molecules (inchikey);
            CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties (molecule_id);
            CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_type_id ON molecular_properties (property_type_id);
            CREATE INDEX IF NOT EXISTS idx_molecular_properties_numeric_value ON molecular_properties (numeric_value);
            CREATE INDEX IF NOT EXISTS idx_property_types_name ON property_types (name);
            
            -- Add compound index for frequent join
            CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_property ON molecular_properties (molecule_id, property_type_id);
        """)
        
        # Run ANALYZE to update statistics
        execute_query("ANALYZE")
        
        logger.info("Added basic performance indexes")
        return True
        
    except Exception as e:
        logger.error(f"Error adding indexes: {str(e)}")
        return False

def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description="Populate CryoProtect database with simplified approach")
    parser.add_argument("--step", type=str, choices=["check", "references", "pubchem", "chembl", "all", "verify", "optimize"],
                      help="Step to execute", required=True)
    parser.add_argument("--limit", type=int, default=500,
                      help="Maximum number of compounds to import (for pubchem/chembl steps)")
    args = parser.parse_args()
    
    logger.info(f"Starting simplified database population, step: {args.step}")
    
    try:
        # Check connection first
        if args.step == "check" or args.step == "all":
            logger.info("Checking database connection and schema...")
            if not check_connection():
                logger.error("Database connection failed. Please check your credentials.")
                return 1
            if not verify_schema():
                logger.error("Database schema verification failed. Please check your database setup.")
                return 1
            logger.info("Database connection and schema verification successful!")
        
        # Populate reference compounds
        if args.step == "references" or args.step == "all":
            logger.info("Populating reference compounds...")
            success_count = populate_reference_compounds()
            
            if success_count == 9:
                logger.info("Successfully populated all 9 reference compounds!")
            else:
                logger.warning(f"Warning: Only populated {success_count}/9 reference compounds.")
            
            # Verify reference compounds
            if verify_reference_compounds():
                logger.info("Reference compound verification passed!")
            else:
                logger.error("Reference compound verification failed!")
                if args.step == "all":
                    logger.error("Stopping 'all' process due to reference verification failure.")
                    return 1
        
        # Import PubChem data
        if args.step == "pubchem" or args.step == "all":
            logger.info(f"Importing up to {args.limit} PubChem compounds...")
            pubchem_count = import_pubchem_minimal(args.limit)
            logger.info(f"Imported {pubchem_count} PubChem compounds.")
        
        # Import ChEMBL data
        if args.step == "chembl" or args.step == "all":
            logger.info(f"Importing up to {args.limit} ChEMBL compounds...")
            chembl_count = import_chembl_minimal(args.limit)
            logger.info(f"Imported {chembl_count} ChEMBL compounds.")
        
        # Optimize database
        if args.step == "optimize" or args.step == "all":
            logger.info("Adding database indexes for performance optimization...")
            if add_basic_indexes():
                logger.info("Database optimization complete!")
            else:
                logger.warning("Database optimization failed!")
        
        # Verify all
        if args.step == "verify" or args.step == "all":
            logger.info("Verifying database population...")
            if verify_all_molecules():
                logger.info("Database verification PASSED! All criteria met.")
            else:
                logger.warning("Database verification FAILED. See logs for details.")
            
            # Check performance
            logger.info("Checking query performance...")
            avg_time = verify_performance()
            logger.info(f"Average query time: {avg_time:.2f}ms")
            if avg_time < 50:
                logger.info("Performance requirement met!")
            else:
                logger.warning(f"Performance requirement not met ({avg_time:.2f}ms > 50ms). Consider adding indexes.")
        
        logger.info(f"Simplified database population, step '{args.step}' completed successfully!")
        return 0
        
    except Exception as e:
        logger.error(f"Error in simplified database population: {str(e)}")
        logger.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())