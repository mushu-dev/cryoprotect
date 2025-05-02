#!/usr/bin/env python3
"""
CryoProtect v2 - Molecule Population Script

This script populates the molecules and molecular_properties tables with scientifically accurate
cryoprotectant data. It ensures proper chemical identifiers (SMILES, InChI, InChIKey) and
physical properties are included.

Usage:
    python populate_molecules.py [--dry-run]

Environment variables required (from .env):
    SUPABASE_URL, SUPABASE_KEY, SUPABASE_USER, SUPABASE_PASSWORD
"""

import os
import json
import uuid
import argparse
import logging
from datetime import datetime
import time
from dotenv import load_dotenv
from supabase import create_client, Client
from service_role_helper import get_supabase_client, get_user_id, ensure_user_profile

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("molecule_population.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Supabase connection
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
SUPABASE_USER = os.getenv("SUPABASE_USER")
SUPABASE_PASSWORD = os.getenv("SUPABASE_PASSWORD")

if not SUPABASE_URL or not SUPABASE_KEY:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

# Common cryoprotectants with accurate scientific data
CRYOPROTECTANTS = [
    {
        "name": "Dimethyl sulfoxide",
        "smiles": "CS(=O)C",
        "inchi": "InChI=1S/C2H6OS/c1-4(2)3/h1-2H3",
        "inchikey": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N",
        "formula": "C2H6OS",
        "molecular_weight": 78.13,
        "properties": [
            {"property_type": "LogP", "value": -1.35, "unit": ""},
            {"property_type": "Melting Point", "value": 18.5, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 189, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -137, "unit": "°C"}
        ]
    },
    {
        "name": "Glycerol",
        "smiles": "C(C(CO)O)O",
        "inchi": "InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2",
        "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
        "formula": "C3H8O3",
        "molecular_weight": 92.09,
        "properties": [
            {"property_type": "LogP", "value": -1.76, "unit": ""},
            {"property_type": "Melting Point", "value": 17.8, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 290, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -93, "unit": "°C"}
        ]
    },
    {
        "name": "Ethylene glycol",
        "smiles": "C(CO)O",
        "inchi": "InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2",
        "inchikey": "LYCAIKOWRPUZTN-UHFFFAOYSA-N",
        "formula": "C2H6O2",
        "molecular_weight": 62.07,
        "properties": [
            {"property_type": "LogP", "value": -1.36, "unit": ""},
            {"property_type": "Melting Point", "value": -12.9, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 197.3, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -128, "unit": "°C"}
        ]
    },
    {
        "name": "Propylene glycol",
        "smiles": "CC(CO)O",
        "inchi": "InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3",
        "inchikey": "DNIAPMSPPWPWGF-UHFFFAOYSA-N",
        "formula": "C3H8O2",
        "molecular_weight": 76.09,
        "properties": [
            {"property_type": "LogP", "value": -0.92, "unit": ""},
            {"property_type": "Melting Point", "value": -59, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 188.2, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -108, "unit": "°C"}
        ]
    },
    {
        "name": "Trehalose",
        "smiles": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)CO)O",
        "inchi": "InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12-10(19)6(15)5(14)3(2-13)22-12/h3-19H,1-2H2/t3-,4-,5-,6-,7-,8+,9-,10-,11-,12+/m1/s1",
        "inchikey": "OHCBMWOFNXFUKL-JCCZQYLRSA-N",
        "formula": "C12H22O11",
        "molecular_weight": 342.30,
        "properties": [
            {"property_type": "LogP", "value": -4.23, "unit": ""},
            {"property_type": "Melting Point", "value": 203, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": 115, "unit": "°C"}
        ]
    },
    {
        "name": "Sucrose",
        "smiles": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O",
        "inchi": "InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12(3-15)10(19)6(2-14)22-5(12)20/h4-11,13-20H,1-3H2/t4-,5+,6-,7-,8+,9-,10+,11-,12+/m1/s1",
        "inchikey": "CZMRCDWAGMRECN-UGDNZRGBSA-N",
        "formula": "C12H22O11",
        "molecular_weight": 342.30,
        "properties": [
            {"property_type": "LogP", "value": -3.76, "unit": ""},
            {"property_type": "Melting Point", "value": 186, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": 65, "unit": "°C"}
        ]
    },
    {
        "name": "Methanol",
        "smiles": "CO",
        "inchi": "InChI=1S/CH4O/c1-2/h2H,1H3",
        "inchikey": "OKKJLVBELUTLKV-UHFFFAOYSA-N",
        "formula": "CH4O",
        "molecular_weight": 32.04,
        "properties": [
            {"property_type": "LogP", "value": -0.77, "unit": ""},
            {"property_type": "Melting Point", "value": -97.6, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 64.7, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -175, "unit": "°C"}
        ]
    },
    {
        "name": "Formamide",
        "smiles": "C(=O)N",
        "inchi": "InChI=1S/CH3NO/c2-1-3/h1H,(H2,2,3)",
        "inchikey": "ZHNUHDYFZUAESO-UHFFFAOYSA-N",
        "formula": "CH3NO",
        "molecular_weight": 45.04,
        "properties": [
            {"property_type": "LogP", "value": -1.51, "unit": ""},
            {"property_type": "Melting Point", "value": 2.55, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 210, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -113, "unit": "°C"}
        ]
    },
    {
        "name": "Acetamide",
        "smiles": "CC(=O)N",
        "inchi": "InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)",
        "inchikey": "DLFVBJFMPXGRIB-UHFFFAOYSA-N",
        "formula": "C2H5NO",
        "molecular_weight": 59.07,
        "properties": [
            {"property_type": "LogP", "value": -1.26, "unit": ""},
            {"property_type": "Melting Point", "value": 79.5, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 221.2, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -73, "unit": "°C"}
        ]
    },
    {
        "name": "1,2-Propanediol",
        "smiles": "CC(CO)O",
        "inchi": "InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3",
        "inchikey": "DNIAPMSPPWPWGF-UHFFFAOYSA-N",
        "formula": "C3H8O2",
        "molecular_weight": 76.09,
        "properties": [
            {"property_type": "LogP", "value": -0.92, "unit": ""},
            {"property_type": "Melting Point", "value": -59, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 188.2, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -108, "unit": "°C"}
        ]
    }
]

def connect_to_supabase():
    """Connect to Supabase using service role helper."""
    try:
        # Use the service role helper to get a client
        supabase = get_supabase_client()
        logger.info("Connected to Supabase using service role")
        return supabase
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        
        # Fall back to original method if service role helper fails
        supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
        
        # Authenticate if credentials provided
        if SUPABASE_USER and SUPABASE_PASSWORD:
            try:
                response = supabase.auth.sign_in_with_password({
                    "email": SUPABASE_USER,
                    "password": SUPABASE_PASSWORD
                })
                if hasattr(response, 'error') and response.error:
                    logger.warning(f"Authentication error: {response.error}")
                    logger.warning("Continuing without authentication. Some operations may fail.")
                else:
                    logger.info(f"Authenticated as {SUPABASE_USER}")
            except Exception as e:
                logger.warning(f"Authentication error: {str(e)}")
                logger.warning("Continuing without authentication. Some operations may fail.")
        else:
            logger.warning("No authentication credentials provided. Continuing without authentication.")
        
        return supabase

def get_authenticated_user_id(supabase):
    """Get the user ID of the authenticated user."""
    try:
        # Try to get user ID from service role helper first
        user_id = get_user_id()
        if user_id:
            logger.info(f"Using service role user ID: {user_id}")
            return user_id
            
        # Fall back to original method
        user = supabase.auth.get_user()
        if hasattr(user, 'user') and user.user:
            return user.user.id
        return None
    except Exception as e:
        logger.error(f"Error getting user ID: {str(e)}")
        return None

def create_user_profile(supabase, auth_user_id, dry_run=False):
    """Create a user profile if it doesn't exist."""
    try:
        # Try to use service role helper first
        profile_id = ensure_user_profile(supabase)
        if profile_id:
            logger.info(f"User profile ensured with ID: {profile_id}")
            return profile_id
    except Exception as e:
        logger.warning(f"Error using service role helper for profile: {str(e)}")
        # Fall back to original method
    
    if not auth_user_id:
        logger.warning("No auth user ID provided. Skipping user profile creation.")
        return None
    
    # Check if profile already exists
    response = supabase.table("user_profile").select("*").eq("auth_user_id", auth_user_id).execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"User profile already exists for {auth_user_id}")
        return response.data[0]["id"]
    
    # Create new profile
    profile_id = str(uuid.uuid4())
    profile_data = {
        "id": profile_id,
        "auth_user_id": auth_user_id,
        "display_name": SUPABASE_USER.split('@')[0] if SUPABASE_USER else "CryoProtect User",
        "email": SUPABASE_USER,
        "affiliation": "CryoProtect Project",
        "created_at": datetime.now().isoformat(),
        "updated_at": datetime.now().isoformat()
    }
    
    if not dry_run:
        response = supabase.table("user_profile").insert(profile_data).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error creating user profile: {response.error}")
            return None
        logger.info(f"Created user profile with ID: {profile_id}")
    else:
        logger.info(f"DRY RUN: Would create user profile with ID: {profile_id}")
    
    return profile_id

def populate_molecules(supabase, user_profile_id, dry_run=False):
    """Populate the molecules table with scientifically accurate data."""
    if not user_profile_id:
        logger.warning("No user profile ID provided. Skipping molecules population.")
        return {}

    # Check if molecules already exist
    response = supabase.table("molecules").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Molecules already exist ({len(response.data)} found)")
        return {mol["name"].lower(): mol["id"] for mol in response.data if "name" in mol}

    # Prepare molecule data
    molecule_map = {}
    batch_molecules = []
    for cryo in CRYOPROTECTANTS:
        molecule_id = str(uuid.uuid4())
        molecule_data = {
            "id": molecule_id,
            "name": cryo["name"],
            "smiles": cryo["smiles"],
            "inchi": cryo["inchi"],
            "inchikey": cryo["inchikey"],
            "formula": cryo["formula"],
            "molecular_weight": cryo["molecular_weight"],
            "created_by": user_profile_id,
            "data_source": "CryoProtect v2 Database Population Script",
            "version": 1,
            "modification_history": json.dumps([{
                "timestamp": datetime.now().isoformat(),
                "action": "created",
                "user_id": user_profile_id
            }]),
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat()
        }
        batch_molecules.append(molecule_data)
        molecule_map[cryo["name"].lower()] = molecule_id

    # Performance Optimization: Batch insert all molecules at once to reduce API calls and improve speed.
    if not dry_run and batch_molecules:
        response = supabase.table("molecules").insert(batch_molecules).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error batch inserting molecules: {response.error}")
        else:
            logger.info(f"Batch inserted {len(batch_molecules)} molecules")
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_molecules)} molecules")

    logger.info(f"Populated {len(molecule_map)} molecules")
    return molecule_map

def populate_molecular_properties(supabase, molecule_map, user_profile_id, dry_run=False):
    """Populate the molecular_properties table with properties for each molecule."""
    if not molecule_map or not user_profile_id:
        logger.warning("Missing molecule_map or user_profile_id. Skipping property population.")
        return
    
    # Check if properties already exist
    response = supabase.table("molecular_properties").select("*").execute()
    if hasattr(response, 'data') and response.data and len(response.data) > 10:
        logger.info(f"Molecular properties already exist ({len(response.data)} found)")
        return
    
    # Populate properties for each molecule
    property_count = 0
    batch_properties = []
    for cryo in CRYOPROTECTANTS:
        molecule_id = molecule_map.get(cryo["name"].lower())
        if not molecule_id:
            logger.warning(f"Molecule ID not found for {cryo['name']}. Skipping properties.")
            continue
        
        for prop in cryo["properties"]:
            property_id = str(uuid.uuid4())
            property_data = {
                "id": property_id,
                "molecule_id": molecule_id,
                "property_type": prop["property_type"],
                "value": prop["value"],
                "unit": prop["unit"],
                "created_by": user_profile_id,
                "data_source": "CryoProtect v2 Database Population Script",
                "version": 1,
                "modification_history": json.dumps([{
                    "timestamp": datetime.now().isoformat(),
                    "action": "created",
                    "user_id": user_profile_id
                }]),
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat()
            }
            batch_properties.append(property_data)
            property_count += 1

    # Performance Optimization: Batch insert all properties at once to reduce API calls and improve speed.
    if not dry_run and batch_properties:
        response = supabase.table("molecular_properties").insert(batch_properties).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error batch inserting properties: {response.error}")
        else:
            logger.info(f"Batch inserted {len(batch_properties)} molecular properties")
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_properties)} molecular properties")

    logger.info(f"Populated {property_count} molecular properties")

def main():
    parser = argparse.ArgumentParser(description="Populate molecules and molecular_properties tables with scientifically accurate data.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 Molecule Population")
    
    # Connect to Supabase
    supabase = connect_to_supabase()
    
    # Get authenticated user ID
    auth_user_id = get_authenticated_user_id(supabase)
    
    # Create user profile if needed
    user_profile_id = create_user_profile(supabase, auth_user_id, args.dry_run)
    
    # Populate molecules
    molecule_map = populate_molecules(supabase, user_profile_id, args.dry_run)
    
    # Populate molecular properties
    populate_molecular_properties(supabase, molecule_map, user_profile_id, args.dry_run)
    
    logger.info("Molecule population complete")

if __name__ == "__main__":
    main()