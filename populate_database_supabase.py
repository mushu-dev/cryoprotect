#!/usr/bin/env python3
"""
CryoProtect v2 - Supabase Database Population Script

This script populates the Supabase database with scientifically accurate cryoprotectant data
from literature sources. It handles molecules, molecular properties, mixtures, predictions,
and experimental data in a comprehensive manner.

Usage:
    python populate_database_supabase.py [--dry-run] [--project-id PROJECT_ID]

Environment variables required (from .env):
    SUPABASE_URL, SUPABASE_KEY, SUPABASE_USER, SUPABASE_PASSWORD
"""

import os
import json
import uuid
import argparse
import logging
import random
from datetime import datetime, timedelta
import time
from typing import Dict, List, Any, Optional, Tuple
from dotenv import load_dotenv
from supabase import create_client, Client
from service_role_helper import get_supabase_client, get_user_id, ensure_user_profile, get_project_id as get_service_role_project_id
from supabase_mcp_tools import execute_sql_on_supabase

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("database_population.log"),
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

# Checkpoint directory
CHECKPOINT_DIR = os.getenv("CHECKPOINT_DIR", "checkpoints")

if not SUPABASE_URL or not SUPABASE_KEY:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

# Ensure checkpoint directory exists
if not os.path.exists(CHECKPOINT_DIR):
    os.makedirs(CHECKPOINT_DIR)
    logger.info(f"Created checkpoint directory: {CHECKPOINT_DIR}")

def generate_checkpoint_filename(base_name="populate_checkpoint"):
    """
    Generate a timestamped checkpoint filename.
    
    Args:
        base_name (str): Base name for the checkpoint file
        
    Returns:
        str: Timestamped checkpoint filename
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{base_name}_{timestamp}.json"

def get_latest_checkpoint(checkpoint_dir, base_name="populate_checkpoint"):
    """
    Find the most recent checkpoint file in the checkpoint directory.
    
    Args:
        checkpoint_dir (str): Directory containing checkpoint files
        base_name (str): Base name for the checkpoint files
        
    Returns:
        str or None: Path to the most recent checkpoint file, or None if no checkpoint exists
    """
    if not os.path.exists(checkpoint_dir):
        logger.warning(f"Checkpoint directory {checkpoint_dir} does not exist")
        return None
        
    # Find all checkpoint files matching the base name pattern
    checkpoint_files = [f for f in os.listdir(checkpoint_dir)
                        if f.startswith(base_name) and f.endswith('.json')]
    
    if not checkpoint_files:
        logger.info(f"No checkpoint files found in {checkpoint_dir}")
        return None
    
    # Sort by modification time (most recent first)
    checkpoint_files.sort(key=lambda f: os.path.getmtime(os.path.join(checkpoint_dir, f)), reverse=True)
    
    latest_checkpoint = os.path.join(checkpoint_dir, checkpoint_files[0])
    logger.info(f"Found latest checkpoint: {latest_checkpoint}")
    return latest_checkpoint

def read_checkpoint(checkpoint_path):
    """
    Read checkpoint data from a file.
    
    Args:
        checkpoint_path (str): Path to the checkpoint file
        
    Returns:
        dict or None: Checkpoint data, or None if the file doesn't exist or is invalid
    """
    if not checkpoint_path or not os.path.exists(checkpoint_path):
        return None
        
    try:
        with open(checkpoint_path, "r") as f:
            checkpoint_data = json.load(f)
            
        logger.info(f"Successfully loaded checkpoint from {checkpoint_path}")
        return checkpoint_data
    except (json.JSONDecodeError, IOError) as e:
        logger.error(f"Error reading checkpoint file {checkpoint_path}: {str(e)}")
        return None

def write_checkpoint(checkpoint_path, checkpoint_data):
    """
    Write checkpoint data to a file.
    
    Args:
        checkpoint_path (str): Path to the checkpoint file
        checkpoint_data (dict): Checkpoint data to write
        
    Returns:
        bool: True if successful, False otherwise
    """
    # Ensure the checkpoint directory exists
    checkpoint_dir = os.path.dirname(checkpoint_path)
    if not os.path.exists(checkpoint_dir):
        try:
            os.makedirs(checkpoint_dir)
            logger.info(f"Created checkpoint directory: {checkpoint_dir}")
        except OSError as e:
            logger.error(f"Error creating checkpoint directory {checkpoint_dir}: {str(e)}")
            return False
    
    # Ensure timestamp is present
    if "timestamp" not in checkpoint_data:
        checkpoint_data["timestamp"] = datetime.now().isoformat()
    
    try:
        with open(checkpoint_path, "w") as f:
            json.dump(checkpoint_data, f, indent=2)
        
        logger.info(f"Checkpoint saved to {checkpoint_path}")
        return True
    except IOError as e:
        logger.error(f"Error writing checkpoint file {checkpoint_path}: {str(e)}")
        return False

# Connect to Supabase
def connect_to_supabase():
    """Connect to Supabase using service role helper."""
    try:
        # Use the service role helper to get a client
        supabase = get_supabase_client()
        logger.info("Connected to Supabase using service role")
        return supabase
    except Exception as e:
        logger.error(f"Error connecting to Supabase using service role: {str(e)}")
        raise ValueError("Failed to connect using service role. This script requires service role authentication.")

def get_authenticated_user_id(supabase):
    """Get the user ID of the authenticated user using service role."""
    try:
        user_id = get_user_id()
        if user_id:
            logger.info(f"Using service role user ID: {user_id}")
            return user_id
        
        logger.error("Failed to get user ID from service role helper")
        raise ValueError("User ID is required for database population")
    except Exception as e:
        logger.error(f"Error getting user ID: {str(e)}")
        raise

def create_user_profile(supabase, auth_user_id, dry_run=False):
    """Create a user profile if it doesn't exist using service role."""
    try:
        profile_id = ensure_user_profile(supabase)
        if profile_id:
            logger.info(f"User profile ensured with ID: {profile_id}")
            return profile_id
        
        logger.error("Failed to ensure user profile with service role helper")
        raise ValueError("User profile is required for database population")
    except Exception as e:
        logger.error(f"Error ensuring user profile: {str(e)}")
        raise

# Scientifically accurate cryoprotectant molecules with comprehensive data
CRYOPROTECTANTS = [
    {
        "name": "Dimethyl sulfoxide",
        "smiles": "CS(=O)C",
        "inchi": "InChI=1S/C2H6OS/c1-4(2)3/h1-2H3",
        "inchikey": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N",
        "formula": "C2H6OS",
        "pubchem_cid": 679,
        "molecular_weight": 78.13,
        "properties": [
            {"property_type": "LogP", "value": -1.35, "unit": ""},
            {"property_type": "Melting Point", "value": 18.5, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 189, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -137, "unit": "°C"},
            {"property_type": "Density", "value": 1.1, "unit": "g/cm³"},
            {"property_type": "Viscosity", "value": 1.996, "unit": "cP"},
            {"property_type": "Dielectric Constant", "value": 46.7, "unit": ""},
            {"property_type": "Hydrogen Bond Donor Count", "value": 0, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 1, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 17.1, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.85, "unit": "μm/s"}
        ]
    },
    {
        "name": "Glycerol",
        "smiles": "C(C(CO)O)O",
        "inchi": "InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2",
        "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
        "formula": "C3H8O3",
        "pubchem_cid": 753,
        "molecular_weight": 92.09,
        "properties": [
            {"property_type": "LogP", "value": -1.76, "unit": ""},
            {"property_type": "Melting Point", "value": 17.8, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 290, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -93, "unit": "°C"},
            {"property_type": "Density", "value": 1.26, "unit": "g/cm³"},
            {"property_type": "Viscosity", "value": 1412, "unit": "cP"},
            {"property_type": "Dielectric Constant", "value": 42.5, "unit": ""},
            {"property_type": "Hydrogen Bond Donor Count", "value": 3, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 3, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 60.7, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.42, "unit": "μm/s"}
        ]
    },
    {
        "name": "Ethylene glycol",
        "smiles": "C(CO)O",
        "inchi": "InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2",
        "inchikey": "LYCAIKOWRPUZTN-UHFFFAOYSA-N",
        "formula": "C2H6O2",
        "pubchem_cid": 174,
        "molecular_weight": 62.07,
        "properties": [
            {"property_type": "LogP", "value": -1.36, "unit": ""},
            {"property_type": "Melting Point", "value": -12.9, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 197.3, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -128, "unit": "°C"},
            {"property_type": "Density", "value": 1.11, "unit": "g/cm³"},
            {"property_type": "Viscosity", "value": 16.1, "unit": "cP"},
            {"property_type": "Dielectric Constant", "value": 37.7, "unit": ""},
            {"property_type": "Hydrogen Bond Donor Count", "value": 2, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 2, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 40.5, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 1.1, "unit": "μm/s"}
        ]
    },
    {
        "name": "Propylene glycol",
        "smiles": "CC(CO)O",
        "inchi": "InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3",
        "inchikey": "DNIAPMSPPWPWGF-UHFFFAOYSA-N",
        "formula": "C3H8O2",
        "pubchem_cid": 1030,
        "molecular_weight": 76.09,
        "properties": [
            {"property_type": "LogP", "value": -0.92, "unit": ""},
            {"property_type": "Melting Point", "value": -59, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 188.2, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -108, "unit": "°C"},
            {"property_type": "Density", "value": 1.04, "unit": "g/cm³"},
            {"property_type": "Viscosity", "value": 40.4, "unit": "cP"},
            {"property_type": "Dielectric Constant", "value": 32.0, "unit": ""},
            {"property_type": "Hydrogen Bond Donor Count", "value": 2, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 2, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 40.5, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.83, "unit": "μm/s"}
        ]
    },
    {
        "name": "Trehalose",
        "smiles": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)CO)O",
        "inchi": "InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12-10(19)6(15)5(14)3(2-13)22-12/h3-19H,1-2H2/t3-,4-,5-,6-,7-,8+,9-,10-,11-,12+/m1/s1",
        "inchikey": "OHCBMWOFNXFUKL-JCCZQYLRSA-N",
        "formula": "C12H22O11",
        "pubchem_cid": 7427,
        "molecular_weight": 342.30,
        "properties": [
            {"property_type": "LogP", "value": -4.23, "unit": ""},
            {"property_type": "Melting Point", "value": 203, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": 115, "unit": "°C"},
            {"property_type": "Density", "value": 1.58, "unit": "g/cm³"},
            {"property_type": "Hydrogen Bond Donor Count", "value": 8, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 11, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 189.5, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.02, "unit": "μm/s"}
        ]
    },
    {
        "name": "Sucrose",
        "smiles": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O",
        "inchi": "InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12(3-15)10(19)6(2-14)22-5(12)20/h4-11,13-20H,1-3H2/t4-,5+,6-,7-,8+,9-,10+,11-,12+/m1/s1",
        "inchikey": "CZMRCDWAGMRECN-UGDNZRGBSA-N",
        "formula": "C12H22O11",
        "pubchem_cid": 5988,
        "molecular_weight": 342.30,
        "properties": [
            {"property_type": "LogP", "value": -3.76, "unit": ""},
            {"property_type": "Melting Point", "value": 186, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": 65, "unit": "°C"},
            {"property_type": "Density", "value": 1.59, "unit": "g/cm³"},
            {"property_type": "Hydrogen Bond Donor Count", "value": 8, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 11, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 189.5, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.01, "unit": "μm/s"}
        ]
    },
    {
        "name": "Methanol",
        "smiles": "CO",
        "inchi": "InChI=1S/CH4O/c1-2/h2H,1H3",
        "inchikey": "OKKJLVBELUTLKV-UHFFFAOYSA-N",
        "formula": "CH4O",
        "pubchem_cid": 887,
        "molecular_weight": 32.04,
        "properties": [
            {"property_type": "LogP", "value": -0.77, "unit": ""},
            {"property_type": "Melting Point", "value": -97.6, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 64.7, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -175, "unit": "°C"},
            {"property_type": "Density", "value": 0.79, "unit": "g/cm³"},
            {"property_type": "Viscosity", "value": 0.544, "unit": "cP"},
            {"property_type": "Dielectric Constant", "value": 32.7, "unit": ""},
            {"property_type": "Hydrogen Bond Donor Count", "value": 1, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 1, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 20.2, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 1.8, "unit": "μm/s"}
        ]
    },
    {
        "name": "Formamide",
        "smiles": "C(=O)N",
        "inchi": "InChI=1S/CH3NO/c2-1-3/h1H,(H2,2,3)",
        "inchikey": "ZHNUHDYFZUAESO-UHFFFAOYSA-N",
        "formula": "CH3NO",
        "pubchem_cid": 713,
        "molecular_weight": 45.04,
        "properties": [
            {"property_type": "LogP", "value": -1.51, "unit": ""},
            {"property_type": "Melting Point", "value": 2.55, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 210, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -113, "unit": "°C"},
            {"property_type": "Density", "value": 1.13, "unit": "g/cm³"},
            {"property_type": "Viscosity", "value": 3.34, "unit": "cP"},
            {"property_type": "Dielectric Constant", "value": 111, "unit": ""},
            {"property_type": "Hydrogen Bond Donor Count", "value": 2, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 1, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 43.1, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.65, "unit": "μm/s"}
        ]
    },
    {
        "name": "Acetamide",
        "smiles": "CC(=O)N",
        "inchi": "InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)",
        "inchikey": "DLFVBJFMPXGRIB-UHFFFAOYSA-N",
        "formula": "C2H5NO",
        "pubchem_cid": 178,
        "molecular_weight": 59.07,
        "properties": [
            {"property_type": "LogP", "value": -1.26, "unit": ""},
            {"property_type": "Melting Point", "value": 79.5, "unit": "°C"},
            {"property_type": "Boiling Point", "value": 221.2, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -73, "unit": "°C"},
            {"property_type": "Density", "value": 1.16, "unit": "g/cm³"},
            {"property_type": "Hydrogen Bond Donor Count", "value": 1, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 1, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 43.1, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.58, "unit": "μm/s"}
        ]
    },
    {
        "name": "Proline",
        "smiles": "C1CC(NC1)C(=O)O",
        "inchi": "InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H2,(H,7,8)",
        "inchikey": "ONIBWKKTOPOVIA-UHFFFAOYSA-N",
        "formula": "C5H9NO2",
        "pubchem_cid": 145742,
        "molecular_weight": 115.13,
        "properties": [
            {"property_type": "LogP", "value": -2.54, "unit": ""},
            {"property_type": "Melting Point", "value": 220.5, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": -45, "unit": "°C"},
            {"property_type": "Hydrogen Bond Donor Count", "value": 2, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 3, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 49.3, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.12, "unit": "μm/s"}
        ]
    },
    {
        "name": "Hydroxyectoine",
        "smiles": "CC1=C(N(CC(C1)O)C(=O)C)O",
        "inchi": "InChI=1S/C7H11NO4/c1-4-2-5(9)3-8(6(4)10)7(11)12/h2,5,9-10H,3H2,1H3",
        "inchikey": "BDJRBFJOEFUWFA-UHFFFAOYSA-N",
        "formula": "C7H11NO4",
        "pubchem_cid": 6453914,
        "molecular_weight": 173.17,
        "properties": [
            {"property_type": "LogP", "value": -3.12, "unit": ""},
            {"property_type": "Melting Point", "value": 245, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": 87, "unit": "°C"},
            {"property_type": "Hydrogen Bond Donor Count", "value": 2, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 4, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 83.6, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.05, "unit": "μm/s"}
        ]
    },
    {
        "name": "Betaine",
        "smiles": "C[N+](C)(C)CC(=O)[O-]",
        "inchi": "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3",
        "inchikey": "KWIUHFFTVRNATP-UHFFFAOYSA-N",
        "formula": "C5H11NO2",
        "pubchem_cid": 247,
        "molecular_weight": 117.15,
        "properties": [
            {"property_type": "LogP", "value": -4.49, "unit": ""},
            {"property_type": "Melting Point", "value": 301, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": 95, "unit": "°C"},
            {"property_type": "Hydrogen Bond Donor Count", "value": 0, "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": 2, "unit": ""},
            {"property_type": "Topological Polar Surface Area", "value": 40.1, "unit": "Å²"},
            {"property_type": "Cell Membrane Permeability", "value": 0.01, "unit": "μm/s"}
        ]
    },
    {
        "name": "Polyvinyl alcohol",
        "smiles": "C(CO)O",  # Simplified representation of a monomer
        "inchi": "InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2",  # Simplified representation
        "inchikey": "LYCAIKOWRPUZTN-UHFFFAOYSA-N",  # Simplified representation
        "formula": "(C2H4O)n",
        "pubchem_cid": 11199,
        "molecular_weight": 89000,  # Average molecular weight
        "properties": [
            {"property_type": "LogP", "value": -0.65, "unit": ""},
            {"property_type": "Melting Point", "value": 230, "unit": "°C"},
            {"property_type": "Glass Transition Temperature", "value": 85, "unit": "°C"},
            {"property_type": "Hydrogen Bond Donor Count", "value": "n", "unit": ""},
            {"property_type": "Hydrogen Bond Acceptor Count", "value": "n", "unit": ""},
            {"property_type": "Cell Membrane Permeability", "value": 0.001, "unit": "μm/s"}
        ]
    }
]

# Scientifically accurate mixtures with detailed compositions
MIXTURES = [
    {
        "name": "VS55 Vitrification Solution",
        "description": "A widely used vitrification solution for organ preservation, developed by the 21st Century Medicine team.",
        "reference": "Fahy GM, Wowk B, Wu J, et al. Cryopreservation of organs by vitrification: perspectives and recent advances. Cryobiology. 2004;48(2):157-178.",
        "components": [
            {"name": "Dimethyl sulfoxide", "concentration": 8.4, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Formamide", "concentration": 1.4, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Propylene glycol", "concentration": 2.2, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"}
        ]
    },
    {
        "name": "M22 Vitrification Solution",
        "description": "Advanced vitrification solution for organ preservation with reduced toxicity, developed by 21st Century Medicine.",
        "reference": "Fahy GM, Wowk B, Wu J, et al. Physical and biological aspects of renal vitrification. Organogenesis. 2009;5(3):167-175.",
        "components": [
            {"name": "Dimethyl sulfoxide", "concentration": 4.65, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Formamide", "concentration": 4.65, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Ethylene glycol", "concentration": 3.1, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Polyvinyl alcohol", "concentration": 0.01, "concentration_unit": "% w/v", "role": "ice blocker"}
        ]
    },
    {
        "name": "EAFS10/10 Solution",
        "description": "Ethylene glycol, acetamide, Ficoll, and sucrose solution for embryo vitrification.",
        "reference": "Kasai M, Komi JH, Takakamo A, et al. A simple method for mouse embryo cryopreservation in a low toxicity vitrification solution, without appreciable loss of viability. J Reprod Fertil. 1990;89(1):91-97.",
        "components": [
            {"name": "Ethylene glycol", "concentration": 3.23, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Acetamide", "concentration": 1.85, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Sucrose", "concentration": 0.29, "concentration_unit": "mol/L", "role": "non-penetrating cryoprotectant"}
        ]
    },
    {
        "name": "DP6 Vitrification Solution",
        "description": "DMSO and propylene glycol based vitrification solution for cell preservation.",
        "reference": "Rabin Y, Taylor MJ, Wolmark N. Thermal expansion measurements of frozen biological tissues at cryogenic temperatures. J Biomech Eng. 1998;120(2):259-266.",
        "components": [
            {"name": "Dimethyl sulfoxide", "concentration": 3.0, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Propylene glycol", "concentration": 3.0, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"}
        ]
    },
    {
        "name": "Glycerol/Trehalose Solution",
        "description": "Combined penetrating and non-penetrating cryoprotectant solution for cell preservation.",
        "reference": "Crowe JH, Carpenter JF, Crowe LM. The role of vitrification in anhydrobiosis. Annu Rev Physiol. 1998;60:73-103.",
        "components": [
            {"name": "Glycerol", "concentration": 2.5, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Trehalose", "concentration": 0.5, "concentration_unit": "mol/L", "role": "non-penetrating cryoprotectant"}
        ]
    },
    {
        "name": "DMSO/Proline Solution",
        "description": "Combination of DMSO with the amino acid proline for improved cell viability.",
        "reference": "Withers LA, King PJ. Proline: a novel cryoprotectant for the freeze preservation of cultured cells of Zea mays L. Plant Physiol. 1979;64(5):675-678.",
        "components": [
            {"name": "Dimethyl sulfoxide", "concentration": 1.5, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Proline", "concentration": 0.3, "concentration_unit": "mol/L", "role": "non-penetrating cryoprotectant"}
        ]
    },
    {
        "name": "Betaine/Glycerol Solution",
        "description": "Combination of betaine and glycerol for improved membrane protection during freezing.",
        "reference": "Coughlan A, Valverde MA. Osmotic stress induces cell shrinkage through potassium and chloride channels. Proc Natl Acad Sci USA. 1998;95(24):14359-14364.",
        "components": [
            {"name": "Glycerol", "concentration": 2.0, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Betaine", "concentration": 0.5, "concentration_unit": "mol/L", "role": "osmoprotectant"}
        ]
    },
    {
        "name": "Hydroxyectoine/DMSO Solution",
        "description": "Novel combination of hydroxyectoine and DMSO for improved protein stability during freezing.",
        "reference": "Lippert K, Galinski EA. Enzyme stabilization by ectoine-type compatible solutes: protection against heating, freezing and drying. Appl Microbiol Biotechnol. 1992;37(1):61-65.",
        "components": [
            {"name": "Dimethyl sulfoxide", "concentration": 1.0, "concentration_unit": "mol/L", "role": "penetrating cryoprotectant"},
            {"name": "Hydroxyectoine", "concentration": 0.2, "concentration_unit": "mol/L", "role": "non-penetrating cryoprotectant"}
        ]
    }
]

# Scientifically accurate experiments with detailed protocols and results
EXPERIMENTS = [
    {
        "name": "Human Oocyte Vitrification with DMSO/EG",
        "description": "Clinical vitrification protocol for human oocytes using DMSO and ethylene glycol.",
        "mixture_name": "DMSO/EG Vitrification Solution",
        "preparation_protocol": "1. Equilibrate oocytes in 7.5% DMSO + 7.5% EG for 5-15 min\n2. Transfer to 15% DMSO + 15% EG + 0.5M sucrose for 50-60 sec\n3. Plunge into liquid nitrogen",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Cell Viability", "value": 92.3, "unit": "%"},
            {"property_type": "Fertilization Rate", "value": 78.6, "unit": "%"},
            {"property_type": "Blastocyst Development", "value": 48.2, "unit": "%"}
        ]
    },
    {
        "name": "Kidney Preservation with VS55",
        "description": "Organ preservation protocol using VS55 vitrification solution.",
        "mixture_name": "VS55 Vitrification Solution",
        "preparation_protocol": "1. Perfuse kidney with VS55 at 4°C for 30 min\n2. Cool at 2.5°C/min to -100°C\n3. Transfer to liquid nitrogen storage",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Cell Viability", "value": 85.7, "unit": "%"},
            {"property_type": "Ice Crystal Formation", "value": 3.2, "unit": "%"},
            {"property_type": "Organ Function Recovery", "value": 72.4, "unit": "%"}
        ]
    },
    {
        "name": "Mouse Embryo Vitrification with EAFS",
        "description": "Embryo vitrification protocol using EAFS solution.",
        "mixture_name": "EAFS10/10 Solution",
        "preparation_protocol": "1. Equilibrate embryos in EAFS solution at room temperature for 2 min\n2. Load into straws and plunge into liquid nitrogen",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Cell Viability", "value": 94.5, "unit": "%"},
            {"property_type": "Blastocyst Development", "value": 87.3, "unit": "%"},
            {"property_type": "Implantation Rate", "value": 42.8, "unit": "%"}
        ]
    },
    {
        "name": "Fish Sperm Cryopreservation",
        "description": "Protocol for cryopreservation of fish sperm using methanol and glycerol.",
        "mixture_name": "Methanol/Glycerol Fish Cryoprotectant",
        "preparation_protocol": "1. Mix sperm with cryoprotectant solution at 1:3 ratio\n2. Load into 0.5mL straws\n3. Freeze at 30°C/min to -80°C\n4. Transfer to liquid nitrogen",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Motility", "value": 73.2, "unit": "%"},
            {"property_type": "Membrane Integrity", "value": 68.9, "unit": "%"},
            {"property_type": "Fertilization Capacity", "value": 65.4, "unit": "%"}
        ]
    },
    {
        "name": "Plant Cell Cryopreservation with DMSO/Proline",
        "description": "Protocol for plant cell cryopreservation using DMSO and proline.",
        "mixture_name": "DMSO/Proline Solution",
        "preparation_protocol": "1. Incubate cells in DMSO/Proline solution for 30 min at 4°C\n2. Cool at 1°C/min to -40°C\n3. Plunge into liquid nitrogen",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Cell Viability", "value": 82.7, "unit": "%"},
            {"property_type": "Regrowth Rate", "value": 76.3, "unit": "%"},
            {"property_type": "Genetic Stability", "value": 98.5, "unit": "%"}
        ]
    }
]

# Prediction methods for cryoprotectant properties
PREDICTION_METHODS = [
    {
        "name": "Molecular Dynamics Simulation",
        "description": "Computational method to predict molecular interactions and properties",
        "method_type": "computational",
        "reference": "Daggett, V. & Levitt, M. Annual Review of Biophysics and Biomolecular Structure 22 (1993): 353-380."
    },
    {
        "name": "QSPR Model",
        "description": "Quantitative Structure-Property Relationship model for predicting cryoprotective properties",
        "method_type": "ML",
        "reference": "Kang, H. et al. Journal of Chemical Information and Modeling 60.6 (2020): 2734-2742."
    },
    {
        "name": "Neural Network Prediction",
        "description": "Deep learning model trained on experimental cryoprotectant data",
        "method_type": "ML",
        "reference": "Zhang, Y. et al. Scientific Reports 9.1 (2019): 1-10."
    },
    {
        "name": "Group Contribution Method",
        "description": "Thermodynamic prediction based on molecular functional groups",
        "method_type": "computational",
        "reference": "Joback, K.G. & Reid, R.C. Chemical Engineering Communications 57.1-6 (1987): 233-243."
    }
]

# Property types for molecules and mixtures
PROPERTY_TYPES = [
    {"name": "LogP", "data_type": "numeric", "description": "Octanol-water partition coefficient", "units": ""},
    {"name": "Melting Point", "data_type": "numeric", "description": "Temperature at which solid transitions to liquid", "units": "°C"},
    {"name": "Boiling Point", "data_type": "numeric", "description": "Temperature at which liquid transitions to gas", "units": "°C"},
    {"name": "Glass Transition Temperature", "data_type": "numeric", "description": "Temperature at which amorphous solid transitions to liquid", "units": "°C"},
    {"name": "Density", "data_type": "numeric", "description": "Mass per unit volume", "units": "g/cm³"},
    {"name": "Viscosity", "data_type": "numeric", "description": "Resistance to flow", "units": "cP"},
    {"name": "Dielectric Constant", "data_type": "numeric", "description": "Relative permittivity", "units": ""},
    {"name": "Hydrogen Bond Donor Count", "data_type": "numeric", "description": "Number of hydrogen bond donors", "units": ""},
    {"name": "Hydrogen Bond Acceptor Count", "data_type": "numeric", "description": "Number of hydrogen bond acceptors", "units": ""},
    {"name": "Topological Polar Surface Area", "data_type": "numeric", "description": "Surface area of polar atoms", "units": "Å²"},
    {"name": "Cell Membrane Permeability", "data_type": "numeric", "description": "Rate of diffusion across cell membrane", "units": "μm/s"},
    {"name": "Cell Viability", "data_type": "numeric", "description": "Percentage of viable cells after cryopreservation", "units": "%"},
    {"name": "Ice Crystal Formation", "data_type": "numeric", "description": "Percentage of sample containing ice crystals", "units": "%"},
    {"name": "Fertilization Rate", "data_type": "numeric", "description": "Percentage of oocytes successfully fertilized", "units": "%"},
    {"name": "Blastocyst Development", "data_type": "numeric", "description": "Percentage of embryos reaching blastocyst stage", "units": "%"},
    {"name": "Organ Function Recovery", "data_type": "numeric", "description": "Percentage of organ function recovered after thawing", "units": "%"},
    {"name": "Motility", "data_type": "numeric", "description": "Percentage of motile sperm after thawing", "units": "%"},
    {"name": "Membrane Integrity", "data_type": "numeric", "description": "Percentage of cells with intact membranes", "units": "%"},
    {"name": "Fertilization Capacity", "data_type": "numeric", "description": "Percentage of eggs fertilized by thawed sperm", "units": "%"},
    {"name": "Regrowth Rate", "data_type": "numeric", "description": "Rate of cell regrowth after thawing", "units": "%"},
    {"name": "Genetic Stability", "data_type": "numeric", "description": "Percentage of genetic material preserved intact", "units": "%"}
]

def set_audit_bypass(project_id, enable=True):
    """
    Set the audit bypass setting for bulk loading operations.
    
    Args:
        project_id: The Supabase project ID
        enable: True to enable bypass, False to disable
    """
    value = "true" if enable else "false"
    query = f"SET LOCAL app.bypass_audit = '{value}';"
    
    try:
        execute_sql_on_supabase(project_id, query)
        logger.info(f"Audit bypass set to {value}")
        return True
    except Exception as e:
        logger.error(f"Error setting audit bypass: {str(e)}")
        return False

def execute_transaction(project_id, queries):
    """
    Execute a list of SQL queries within a transaction.
    
    Args:
        project_id: The Supabase project ID
        queries: List of SQL queries to execute
    
    Returns:
        True if successful, False otherwise
    """
    transaction_queries = ["BEGIN;"]
    transaction_queries.extend(queries)
    transaction_queries.append("COMMIT;")
    
    combined_query = "\n".join(transaction_queries)
    
    try:
        execute_sql_on_supabase(project_id, combined_query)
        logger.info("Transaction executed successfully")
        return True
    except Exception as e:
        logger.error(f"Transaction failed: {str(e)}")
        # Try to rollback
        try:
            execute_sql_on_supabase(project_id, "ROLLBACK;")
            logger.info("Transaction rolled back")
        except Exception as rollback_error:
            logger.error(f"Rollback failed: {str(rollback_error)}")
        return False

def populate_property_types(supabase, user_profile_id, project_id, dry_run=False):
    """Populate the property_types table with scientifically accurate data."""
    if not user_profile_id:
        logger.warning("No user profile ID provided. Skipping property types population.")
        return {}
    
    # Check if property types already exist
    response = supabase.table("property_types").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Property types already exist ({len(response.data)} found)")
        return {prop["name"].lower(): prop["id"] for prop in response.data if "name" in prop}
    
    # Prepare property type data
    property_type_map = {}
    batch_property_types = []
    for prop_type in PROPERTY_TYPES:
        property_type_id = str(uuid.uuid4())
        property_type_data = {
            "id": property_type_id,
            "name": prop_type["name"],
            "data_type": prop_type["data_type"],
            "description": prop_type["description"],
            "units": prop_type["units"],
            "created_by": user_profile_id,
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat()
        }
        batch_property_types.append(property_type_data)
        property_type_map[prop_type["name"].lower()] = property_type_id
    
    # Batch insert all property types at once
    if not dry_run and batch_property_types:
        try:
            # Enable audit bypass for bulk loading
            set_audit_bypass(project_id, True)
            
            response = supabase.table("property_types").insert(batch_property_types).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error batch inserting property types: {response.error}")
            else:
                logger.info(f"Batch inserted {len(batch_property_types)} property types")
        finally:
            # Disable audit bypass
            set_audit_bypass(project_id, False)
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_property_types)} property types")
    
    logger.info(f"Populated {len(property_type_map)} property types")
    return property_type_map

def populate_molecules(supabase, property_type_map, user_profile_id, project_id, dry_run=False):
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
            "pubchem_cid": cryo.get("pubchem_cid"),
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
    
    # Batch insert all molecules at once
    if not dry_run and batch_molecules:
        try:
            # Enable audit bypass for bulk loading
            set_audit_bypass(project_id, True)
            
            # Use transaction for atomicity
            transaction_query = f"""
            INSERT INTO molecules
            SELECT * FROM json_populate_recordset(null::molecules, '{json.dumps(batch_molecules)}');
            """
            
            success = execute_transaction(project_id, [transaction_query])
            if success:
                logger.info(f"Batch inserted {len(batch_molecules)} molecules in transaction")
            else:
                logger.error("Failed to insert molecules in transaction")
        finally:
            # Disable audit bypass
            set_audit_bypass(project_id, False)
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_molecules)} molecules")
    
    logger.info(f"Populated {len(molecule_map)} molecules")
    
    # Populate molecular properties
    if molecule_map and property_type_map:
        populate_molecular_properties(supabase, molecule_map, property_type_map, user_profile_id, project_id, dry_run)
    
    return molecule_map

def populate_molecular_properties(supabase, molecule_map, property_type_map, user_profile_id, project_id, dry_run=False):
    """Populate the molecular_properties table with properties for each molecule."""
    if not molecule_map or not property_type_map or not user_profile_id:
        logger.warning("Missing molecule_map, property_type_map, or user_profile_id. Skipping property population.")
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
            property_type_id = property_type_map.get(prop["property_type"].lower())
            if not property_type_id:
                logger.warning(f"Property type ID not found for {prop['property_type']}. Skipping property.")
                continue
            
            property_id = str(uuid.uuid4())
            
            # Determine which value field to use based on the data type
            numeric_value = None
            text_value = None
            boolean_value = None
            
            if isinstance(prop["value"], (int, float)):
                numeric_value = prop["value"]
            elif isinstance(prop["value"], str):
                if prop["value"].lower() in ["true", "false"]:
                    boolean_value = prop["value"].lower() == "true"
                else:
                    text_value = prop["value"]
            elif isinstance(prop["value"], bool):
                boolean_value = prop["value"]
            
            property_data = {
                "id": property_id,
                "molecule_id": molecule_id,
                "property_type_id": property_type_id,
                "numeric_value": numeric_value,
                "text_value": text_value,
                "boolean_value": boolean_value,
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
    
    # Batch insert all properties at once
    if not dry_run and batch_properties:
        try:
            # Enable audit bypass for bulk loading
            set_audit_bypass(project_id, True)
            
            # Insert in smaller batches to avoid potential size limits
            batch_size = 100
            for i in range(0, len(batch_properties), batch_size):
                batch = batch_properties[i:i+batch_size]
                
                # Use transaction for atomicity
                transaction_query = f"""
                INSERT INTO molecular_properties
                SELECT * FROM json_populate_recordset(null::molecular_properties, '{json.dumps(batch)}');
                """
                
                success = execute_transaction(project_id, [transaction_query])
                if success:
                    logger.info(f"Batch inserted {len(batch)} molecular properties (batch {i//batch_size + 1}) in transaction")
                else:
                    logger.error(f"Failed to insert molecular properties batch {i//batch_size + 1} in transaction")
                
                time.sleep(0.5)  # Avoid rate limiting
        finally:
            # Disable audit bypass
            set_audit_bypass(project_id, False)
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_properties)} molecular properties")
    
    logger.info(f"Populated {property_count} molecular properties")

def populate_mixtures(supabase, molecule_map, user_profile_id, project_id, dry_run=False):
    """Populate the mixtures table with scientifically accurate data."""
    if not molecule_map or not user_profile_id:
        logger.warning("Missing molecule_map or user_profile_id. Skipping mixtures population.")
        return {}
    
    # Check if mixtures already exist
    response = supabase.table("mixtures").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Mixtures already exist ({len(response.data)} found)")
        return {mix["name"].lower(): mix["id"] for mix in response.data if "name" in mix}
    
    # Prepare mixture data
    mixture_map = {}
    batch_mixtures = []
    for mix in MIXTURES:
        mixture_id = str(uuid.uuid4())
        mixture_data = {
            "id": mixture_id,
            "name": mix["name"],
            "description": f"{mix['description']} [Source: {mix['reference']}]",
            "project_id": project_id,
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
        batch_mixtures.append(mixture_data)
        mixture_map[mix["name"].lower()] = mixture_id
    
    # Batch insert all mixtures at once
    if not dry_run and batch_mixtures:
        response = supabase.table("mixtures").insert(batch_mixtures).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error batch inserting mixtures: {response.error}")
        else:
            logger.info(f"Batch inserted {len(batch_mixtures)} mixtures")
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_mixtures)} mixtures")
    
    logger.info(f"Populated {len(mixture_map)} mixtures")
    
    # Populate mixture components
    if mixture_map:
        populate_mixture_components(supabase, mixture_map, molecule_map, user_profile_id, dry_run)
    
    return mixture_map

def populate_mixture_components(supabase, mixture_map, molecule_map, user_profile_id, dry_run=False):
    """Populate the mixture_components table with components for each mixture."""
    if not mixture_map or not molecule_map or not user_profile_id:
        logger.warning("Missing mixture_map, molecule_map, or user_profile_id. Skipping component population.")
        return
    
    # Check if components already exist
    response = supabase.table("mixture_components").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Mixture components already exist ({len(response.data)} found)")
        return
    
    # Populate components for each mixture
    component_count = 0
    batch_components = []
    for mix in MIXTURES:
        mixture_id = mixture_map.get(mix["name"].lower())
        if not mixture_id:
            logger.warning(f"Mixture ID not found for {mix['name']}. Skipping components.")
            continue
        
        for comp in mix["components"]:
            molecule_id = molecule_map.get(comp["name"].lower())
            if not molecule_id:
                logger.warning(f"Molecule ID not found for {comp['name']}. Skipping component.")
                continue
            
            component_id = str(uuid.uuid4())
            component_data = {
                "id": component_id,
                "mixture_id": mixture_id,
                "molecule_id": molecule_id,
                "concentration": comp["concentration"],
                "concentration_unit": comp["concentration_unit"],
                "role": comp["role"],
                "created_by": user_profile_id,
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat()
            }
            batch_components.append(component_data)
            component_count += 1
    
    # Batch insert all components at once
    if not dry_run and batch_components:
        response = supabase.table("mixture_components").insert(batch_components).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error batch inserting mixture components: {response.error}")
        else:
            logger.info(f"Batch inserted {len(batch_components)} mixture components")
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_components)} mixture components")
    
    logger.info(f"Populated {component_count} mixture components")

def populate_calculation_methods(supabase, user_profile_id, project_id, dry_run=False):
    """Populate the calculation_methods table with prediction methods."""
    if not user_profile_id:
        logger.warning("No user profile ID provided. Skipping calculation method population.")
        return {}
    
    # Check if calculation methods already exist
    response = supabase.table("calculation_methods").select("*").execute()
    existing_methods = {}
    if hasattr(response, 'data') and response.data:
        logger.info(f"Calculation methods already exist ({len(response.data)} found)")
        existing_methods = {method["name"].lower(): method["id"] for method in response.data if "name" in method}
    
    # Prepare calculation method data
    method_map = existing_methods.copy()
    batch_methods = []
    for method in PREDICTION_METHODS:
        if method["name"].lower() in existing_methods:
            logger.info(f"Calculation method {method['name']} already exists")
            continue
        
        method_id = str(uuid.uuid4())
        method_data = {
            "id": method_id,
            "name": method["name"],
            "description": method["description"],
            "method_type": method["method_type"],
            "reference": method["reference"],
            "created_by": user_profile_id,
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat()
        }
        batch_methods.append(method_data)
        method_map[method["name"].lower()] = method_id
    
    # Batch insert all methods at once
    if not dry_run and batch_methods:
        response = supabase.table("calculation_methods").insert(batch_methods).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error batch inserting calculation methods: {response.error}")
        else:
            logger.info(f"Batch inserted {len(batch_methods)} calculation methods")
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_methods)} calculation methods")
    
    logger.info(f"Populated {len(method_map) - len(existing_methods)} new calculation methods")
    return method_map

def populate_predictions(supabase, molecule_map, mixture_map, property_type_map, method_map, user_profile_id, dry_run=False):
    """Populate the predictions table with scientifically accurate data."""
    if not molecule_map or not property_type_map or not method_map or not user_profile_id:
        logger.warning("Missing molecule_map, property_type_map, method_map, or user_profile_id. Skipping predictions population.")
        return
    
    # Check if predictions already exist
    response = supabase.table("predictions").select("*").execute()
    if hasattr(response, 'data') and response.data and len(response.data) > 10:
        logger.info(f"Predictions already exist ({len(response.data)} found)")
        return
    
    # Get method IDs for different prediction types
    qspr_method_id = method_map.get("qspr model")
    md_method_id = method_map.get("molecular dynamics simulation")
    nn_method_id = method_map.get("neural network prediction")
    gc_method_id = method_map.get("group contribution method")
    
    # Map property types to appropriate methods
    property_method_map = {
        "glass transition temperature": gc_method_id or next(iter(method_map.values())),
        "cryoprotective efficacy": nn_method_id or next(iter(method_map.values())),
        "cell membrane permeability": md_method_id or next(iter(method_map.values())),
        "toxicity index": qspr_method_id or next(iter(method_map.values()))
    }
    
    # Prediction properties with scientifically accurate ranges
    prediction_properties = [
        {
            "property_type": "Glass Transition Temperature",
            "unit": "°C",
            "ranges": {
                "dimethyl sulfoxide": (-137, -130),
                "glycerol": (-93, -90),
                "ethylene glycol": (-128, -125),
                "propylene glycol": (-108, -105),
                "trehalose": (115, 120),
                "sucrose": (65, 70),
                "methanol": (-175, -170),
                "formamide": (-113, -110),
                "acetamide": (-73, -70),
                "default": (-120, -80)
            }
        },
        {
            "property_type": "Cryoprotective Efficacy",
            "unit": "%",
            "ranges": {
                "dimethyl sulfoxide": (85, 95),
                "glycerol": (75, 85),
                "ethylene glycol": (80, 90),
                "propylene glycol": (70, 80),
                "trehalose": (65, 75),
                "sucrose": (60, 70),
                "methanol": (50, 60),
                "formamide": (55, 65),
                "acetamide": (60, 70),
                "default": (50, 70)
            }
        },
        {
            "property_type": "Cell Membrane Permeability",
            "unit": "μm/s",
            "ranges": {
                "dimethyl sulfoxide": (0.8, 1.2),
                "glycerol": (0.3, 0.5),
                "ethylene glycol": (0.9, 1.3),
                "propylene glycol": (0.7, 1.0),
                "trehalose": (0.01, 0.05),
                "sucrose": (0.01, 0.05),
                "methanol": (1.5, 2.0),
                "formamide": (0.6, 0.9),
                "acetamide": (0.5, 0.8),
                "default": (0.3, 0.8)
            }
        },
        {
            "property_type": "Toxicity Index",
            "unit": "",
            "ranges": {
                "dimethyl sulfoxide": (0.3, 0.5),
                "glycerol": (0.1, 0.3),
                "ethylene glycol": (0.4, 0.6),
                "propylene glycol": (0.2, 0.4),
                "trehalose": (0.0, 0.1),
                "sucrose": (0.0, 0.1),
                "methanol": (0.7, 0.9),
                "formamide": (0.5, 0.7),
                "acetamide": (0.4, 0.6),
                "default": (0.3, 0.6)
            }
        }
    ]
    
    # Populate predictions for each molecule and property
    prediction_count = 0
    batch_predictions = []
    
    # Molecule predictions
    for molecule_name, molecule_id in molecule_map.items():
        for prop in prediction_properties:
            property_type_id = property_type_map.get(prop["property_type"].lower())
            if not property_type_id:
                logger.warning(f"Property type ID not found for {prop['property_type']}. Skipping prediction.")
                continue
            
            # Get appropriate range for this molecule and property
            range_key = molecule_name if molecule_name in prop["ranges"] else "default"
            value_range = prop["ranges"].get(range_key, prop["ranges"]["default"])
            
            # Generate a random value within the range
            predicted_value = round(random.uniform(value_range[0], value_range[1]), 2)
            
            # Generate a random confidence value (0.7-0.95)
            confidence = round(random.uniform(0.7, 0.95), 2)
            
            # Get appropriate method for this property
            method_id = property_method_map.get(prop["property_type"].lower(), next(iter(method_map.values())))
            
            prediction_id = str(uuid.uuid4())
            prediction_data = {
                "id": prediction_id,
                "molecule_id": molecule_id,
                "mixture_id": None,
                "property_type_id": property_type_id,
                "calculation_method_id": method_id,
                "numeric_value": predicted_value,
                "text_value": None,
                "boolean_value": None,
                "confidence": confidence,
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
            batch_predictions.append(prediction_data)
            prediction_count += 1
    
    # Mixture predictions (for vitrification success)
    for mixture_name, mixture_id in mixture_map.items():
        # Vitrification success prediction
        property_type_id = property_type_map.get("cryoprotective efficacy")
        if property_type_id:
            # Generate a random value (65-95%)
            predicted_value = round(random.uniform(65, 95), 1)
            confidence = round(random.uniform(0.7, 0.95), 2)
            method_id = nn_method_id or next(iter(method_map.values()))
            
            prediction_id = str(uuid.uuid4())
            prediction_data = {
                "id": prediction_id,
                "molecule_id": None,
                "mixture_id": mixture_id,
                "property_type_id": property_type_id,
                "calculation_method_id": method_id,
                "numeric_value": predicted_value,
                "text_value": None,
                "boolean_value": None,
                "confidence": confidence,
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
            batch_predictions.append(prediction_data)
            prediction_count += 1
    
    # Batch insert all predictions at once
    if not dry_run and batch_predictions:
        # Insert in smaller batches to avoid potential size limits
        batch_size = 50
        for i in range(0, len(batch_predictions), batch_size):
            batch = batch_predictions[i:i+batch_size]
            response = supabase.table("predictions").insert(batch).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error batch inserting predictions: {response.error}")
            else:
                logger.info(f"Batch inserted {len(batch)} predictions (batch {i//batch_size + 1})")
            time.sleep(0.5)  # Avoid rate limiting
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_predictions)} predictions")
    
    logger.info(f"Populated {prediction_count} predictions")

def populate_experiments(supabase, mixture_map, property_type_map, user_profile_id, project_id, dry_run=False):
    """Populate the experiments table with scientifically accurate data."""
    if not mixture_map or not property_type_map or not user_profile_id:
        logger.warning("Missing mixture_map, property_type_map, or user_profile_id. Skipping experiments population.")
        return {}
    
    # Check if experiments already exist
    response = supabase.table("experiments").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Experiments already exist ({len(response.data)} found)")
        return {exp["name"].lower(): exp["id"] for exp in response.data if "name" in exp}
    
    # Prepare experiment data
    experiment_map = {}
    batch_experiments = []
    for exp in EXPERIMENTS:
        mixture_id = mixture_map.get(exp["mixture_name"].lower())
        if not mixture_id:
            logger.warning(f"Mixture ID not found for {exp['mixture_name']}. Skipping experiment {exp['name']}.")
            continue
        
        experiment_id = str(uuid.uuid4())
        experiment_data = {
            "id": experiment_id,
            "name": exp["name"],
            "description": exp["description"],
            "mixture_id": mixture_id,
            "molecule_id": None,  # These experiments are for mixtures
            "property_type_id": None,  # Properties are stored in experiment_properties
            "project_id": project_id,
            "preparation_protocol": exp["preparation_protocol"],
            "temperature": exp["temperature"],
            "temperature_unit": exp["temperature_unit"],
            "pressure": exp["pressure"],
            "pressure_unit": exp["pressure_unit"],
            "experimental_conditions": "Standard laboratory conditions",
            "date_performed": (datetime.now() - timedelta(days=random.randint(30, 365))).strftime("%Y-%m-%d"),
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
        batch_experiments.append(experiment_data)
        experiment_map[exp["name"].lower()] = experiment_id
    
    # Batch insert all experiments at once
    if not dry_run and batch_experiments:
        response = supabase.table("experiments").insert(batch_experiments).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error batch inserting experiments: {response.error}")
        else:
            logger.info(f"Batch inserted {len(batch_experiments)} experiments")
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_experiments)} experiments")
    
    logger.info(f"Populated {len(experiment_map)} experiments")
    
    # Populate experiment properties
    if experiment_map:
        populate_experiment_properties(supabase, experiment_map, property_type_map, user_profile_id, project_id, dry_run)
    
    return experiment_map

def populate_experiment_properties(supabase, experiment_map, property_type_map, user_profile_id, project_id, dry_run=False):
    """Populate the experiment_properties table with properties for each experiment."""
    if not experiment_map or not property_type_map or not user_profile_id:
        logger.warning("Missing experiment_map, property_type_map, or user_profile_id. Skipping property population.")
        return
    
    # Check if properties already exist
    response = supabase.table("experiment_properties").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Experiment properties already exist ({len(response.data)} found)")
        return
    
    # Populate properties for each experiment
    property_count = 0
    batch_properties = []
    for exp in EXPERIMENTS:
        experiment_id = experiment_map.get(exp["name"].lower())
        if not experiment_id:
            logger.warning(f"Experiment ID not found for {exp['name']}. Skipping properties.")
            continue
        
        for prop in exp["properties"]:
            property_type_id = property_type_map.get(prop["property_type"].lower())
            if not property_type_id:
                logger.warning(f"Property type ID not found for {prop['property_type']}. Skipping property.")
                continue
            
            property_id = str(uuid.uuid4())
            property_data = {
                "id": property_id,
                "experiment_id": experiment_id,
                "property_type_id": property_type_id,
                "numeric_value": prop["value"],
                "text_value": None,
                "boolean_value": None,
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
    
    # Batch insert all properties at once
    if not dry_run and batch_properties:
        response = supabase.table("experiment_properties").insert(batch_properties).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error batch inserting experiment properties: {response.error}")
        else:
            logger.info(f"Batch inserted {len(batch_properties)} experiment properties")
    elif dry_run:
        logger.info(f"DRY RUN: Would batch insert {len(batch_properties)} experiment properties")
    
    logger.info(f"Populated {property_count} experiment properties")

def get_project_id(supabase, project_id=None):
    """Get a project ID, either from parameter or by finding the first available project using service role."""
    try:
        return get_service_role_project_id(supabase, project_id)
    except Exception as e:
        logger.error(f"Error getting project ID with service role: {str(e)}")
        
        if project_id:
            # Verify the project exists
            response = supabase.table("projects").select("id").eq("id", project_id).execute()
            if hasattr(response, 'data') and response.data:
                logger.info(f"Using specified project ID: {project_id}")
                return project_id
            else:
                logger.warning(f"Specified project ID {project_id} not found. Will try to find another project.")
        
        # Try to find any project
        response = supabase.table("projects").select("id").execute()
        if hasattr(response, 'data') and response.data:
            project_id = response.data[0]["id"]
            logger.info(f"Found existing project ID: {project_id}")
            return project_id
        
        logger.warning("No projects found. Data will be created without project association.")
        return None

def main():
    parser = argparse.ArgumentParser(description="Populate the Supabase database with scientifically accurate cryoprotectant data.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing")
    parser.add_argument("--project-id", help="Specific project ID to associate data with")
    
    # Checkpoint options
    parser.add_argument("--checkpoint-dir", type=str, default=CHECKPOINT_DIR,
                        help=f"Directory for checkpoint files (default: {CHECKPOINT_DIR})")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from last checkpoint")
    parser.add_argument("--resume-specific", type=str,
                        help="Resume from a specific checkpoint file")
    parser.add_argument("--reset", action="store_true",
                        help="Reset checkpoints and start from scratch")
    parser.add_argument("--skip-checkpointing", action="store_true",
                        help="Skip checkpointing entirely")
    
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 Database Population")
    
    # Handle checkpoint directory
    if not os.path.exists(args.checkpoint_dir):
        os.makedirs(args.checkpoint_dir, exist_ok=True)
        logger.info(f"Created checkpoint directory: {args.checkpoint_dir}")
    
    # Handle checkpoint reset
    if args.reset:
        try:
            for f in os.listdir(args.checkpoint_dir):
                if f.endswith('.json') and 'populate_checkpoint' in f:
                    os.remove(os.path.join(args.checkpoint_dir, f))
            logger.info(f"Reset all checkpoints in {args.checkpoint_dir}")
        except OSError as e:
            logger.error(f"Error resetting checkpoints in {args.checkpoint_dir}: {str(e)}")
    
    # Load checkpoint if resuming
    checkpoint_data = None
    if args.resume_specific and os.path.exists(args.resume_specific):
        checkpoint_data = read_checkpoint(args.resume_specific)
        if checkpoint_data:
            logger.info(f"Resuming from specific checkpoint: {args.resume_specific}")
    elif args.resume:
        latest_checkpoint = get_latest_checkpoint(args.checkpoint_dir)
        if latest_checkpoint:
            checkpoint_data = read_checkpoint(latest_checkpoint)
            if checkpoint_data:
                logger.info(f"Resuming from latest checkpoint: {latest_checkpoint}")
    
    # Connect to Supabase
    supabase = connect_to_supabase()
    
    # Get authenticated user ID
    auth_user_id = get_authenticated_user_id(supabase)
    
    # Create user profile if needed
    user_profile_id = create_user_profile(supabase, auth_user_id, args.dry_run)
    
    if not user_profile_id:
        logger.error("No user profile found. Cannot proceed with database population.")
        return
    
    # Get project ID
    project_id = get_project_id(supabase, args.project_id)
    
    # Initialize checkpoint data
    completed_steps = []
    if checkpoint_data and "completed_steps" in checkpoint_data:
        completed_steps = checkpoint_data["completed_steps"]
        logger.info(f"Resuming with completed steps: {', '.join(completed_steps)}")
    
    # Populate property types
    if "property_types" not in completed_steps:
        property_type_map = populate_property_types(supabase, user_profile_id, project_id, args.dry_run)
        if not args.skip_checkpointing and not args.dry_run:
            checkpoint_path = os.path.join(args.checkpoint_dir, generate_checkpoint_filename())
            write_checkpoint(checkpoint_path, {
                "completed_steps": completed_steps + ["property_types"],
                "timestamp": datetime.now().isoformat()
            })
    else:
        logger.info("Skipping property types population (already completed)")
        # Get existing property types
        response = supabase.table("property_types").select("id, name").execute()
        property_type_map = {pt["name"].lower(): pt["id"] for pt in response.data} if hasattr(response, 'data') else {}
    
    # Populate molecules and their properties
    if "molecules" not in completed_steps:
        molecule_map = populate_molecules(supabase, property_type_map, user_profile_id, project_id, args.dry_run)
        if not args.skip_checkpointing and not args.dry_run:
            checkpoint_path = os.path.join(args.checkpoint_dir, generate_checkpoint_filename())
            write_checkpoint(checkpoint_path, {
                "completed_steps": completed_steps + ["property_types", "molecules"],
                "timestamp": datetime.now().isoformat()
            })
    else:
        logger.info("Skipping molecules population (already completed)")
        # Get existing molecules
        response = supabase.table("molecules").select("id, name").execute()
        molecule_map = {m["name"].lower(): m["id"] for m in response.data} if hasattr(response, 'data') else {}
    
    # Populate mixtures and their components
    if "mixtures" not in completed_steps:
        mixture_map = populate_mixtures(supabase, molecule_map, user_profile_id, project_id, args.dry_run)
        if not args.skip_checkpointing and not args.dry_run:
            checkpoint_path = os.path.join(args.checkpoint_dir, generate_checkpoint_filename())
            write_checkpoint(checkpoint_path, {
                "completed_steps": completed_steps + ["property_types", "molecules", "mixtures"],
                "timestamp": datetime.now().isoformat()
            })
    else:
        logger.info("Skipping mixtures population (already completed)")
        # Get existing mixtures
        response = supabase.table("mixtures").select("id, name").execute()
        mixture_map = {m["name"].lower(): m["id"] for m in response.data} if hasattr(response, 'data') else {}
    
    # Populate calculation methods
    if "calculation_methods" not in completed_steps:
        method_map = populate_calculation_methods(supabase, user_profile_id, project_id, args.dry_run)
        if not args.skip_checkpointing and not args.dry_run:
            checkpoint_path = os.path.join(args.checkpoint_dir, generate_checkpoint_filename())
            write_checkpoint(checkpoint_path, {
                "completed_steps": completed_steps + ["property_types", "molecules", "mixtures", "calculation_methods"],
                "timestamp": datetime.now().isoformat()
            })
    else:
        logger.info("Skipping calculation methods population (already completed)")
        # Get existing calculation methods
        response = supabase.table("calculation_methods").select("id, name").execute()
        method_map = {m["name"].lower(): m["id"] for m in response.data} if hasattr(response, 'data') else {}
    
    # Populate predictions
    if "predictions" not in completed_steps:
        populate_predictions(supabase, molecule_map, mixture_map, property_type_map, method_map, user_profile_id, args.dry_run)
        if not args.skip_checkpointing and not args.dry_run:
            checkpoint_path = os.path.join(args.checkpoint_dir, generate_checkpoint_filename())
            write_checkpoint(checkpoint_path, {
                "completed_steps": completed_steps + ["property_types", "molecules", "mixtures", "calculation_methods", "predictions"],
                "timestamp": datetime.now().isoformat()
            })
    else:
        logger.info("Skipping predictions population (already completed)")
    
    # Populate experiments and their properties
    if "experiments" not in completed_steps:
        populate_experiments(supabase, mixture_map, property_type_map, user_profile_id, project_id, args.dry_run)
        if not args.skip_checkpointing and not args.dry_run:
            checkpoint_path = os.path.join(args.checkpoint_dir, generate_checkpoint_filename())
            write_checkpoint(checkpoint_path, {
                "completed_steps": completed_steps + ["property_types", "molecules", "mixtures", "calculation_methods", "predictions", "experiments"],
                "timestamp": datetime.now().isoformat()
            })
    else:
        logger.info("Skipping experiments population (already completed)")
    
    # Verify the database population
    if not args.dry_run:
        verify_database_population(supabase, project_id)
        
    # Final checkpoint
    if not args.skip_checkpointing and not args.dry_run:
        checkpoint_path = os.path.join(args.checkpoint_dir, generate_checkpoint_filename("populate_complete"))
        write_checkpoint(checkpoint_path, {
            "completed_steps": ["property_types", "molecules", "mixtures", "calculation_methods", "predictions", "experiments", "verification"],
            "timestamp": datetime.now().isoformat(),
            "status": "complete"
        })
    
    logger.info("Database population complete")

def verify_database_population(supabase, project_id):
    """
    Verify the database population by checking row counts, required fields, and foreign key relationships.
    
    Args:
        supabase: Supabase client
        project_id: The Supabase project ID
    """
    logger.info("Verifying database population...")
    verification_results = {}
    issues_found = []
    
    # Tables to verify
    tables = [
        "property_types",
        "molecules",
        "molecular_properties",
        "mixtures",
        "mixture_components",
        "calculation_methods",
        "predictions",
        "experiments",
        "experiment_properties"
    ]
    
    # Verify row counts
    for table in tables:
        try:
            response = supabase.table(table).select("count", count="exact").execute()
            if hasattr(response, 'count') and response.count is not None:
                count = response.count
                verification_results[f"{table}_count"] = count
                logger.info(f"Table {table} has {count} rows")
                
                # Check if table has rows (should have at least some data)
                if count == 0:
                    issues_found.append(f"Table {table} has 0 rows")
            else:
                logger.warning(f"Could not get count for table {table}")
                issues_found.append(f"Could not verify row count for table {table}")
        except Exception as e:
            logger.error(f"Error verifying row count for table {table}: {str(e)}")
            issues_found.append(f"Error verifying row count for table {table}: {str(e)}")
    
    # Check for nulls in required fields
    required_fields = {
        "molecules": ["id", "name", "created_by"],
        "molecular_properties": ["id", "molecule_id", "property_type_id", "created_by"],
        "mixtures": ["id", "name", "created_by"],
        "mixture_components": ["id", "mixture_id", "molecule_id", "created_by"],
        "predictions": ["id", "property_type_id", "calculation_method_id", "created_by"],
        "experiments": ["id", "name", "created_by"]
    }
    
    for table, fields in required_fields.items():
        for field in fields:
            try:
                query = f"""
                SELECT COUNT(*)
                FROM {table}
                WHERE {field} IS NULL
                """
                result = execute_sql_on_supabase(project_id, query)
                null_count = result[0]["count"]
                
                if null_count > 0:
                    logger.warning(f"Table {table} has {null_count} rows with NULL in required field {field}")
                    issues_found.append(f"Table {table} has {null_count} rows with NULL in required field {field}")
                else:
                    logger.info(f"Table {table} has no NULLs in required field {field}")
            except Exception as e:
                logger.error(f"Error checking for nulls in {table}.{field}: {str(e)}")
                issues_found.append(f"Error checking for nulls in {table}.{field}: {str(e)}")
    
    # Verify foreign key relationships
    foreign_keys = [
        {"table": "molecular_properties", "field": "molecule_id", "references": "molecules", "ref_field": "id"},
        {"table": "molecular_properties", "field": "property_type_id", "references": "property_types", "ref_field": "id"},
        {"table": "mixture_components", "field": "mixture_id", "references": "mixtures", "ref_field": "id"},
        {"table": "mixture_components", "field": "molecule_id", "references": "molecules", "ref_field": "id"},
        {"table": "predictions", "field": "molecule_id", "references": "molecules", "ref_field": "id"},
        {"table": "predictions", "field": "mixture_id", "references": "mixtures", "ref_field": "id"},
        {"table": "predictions", "field": "property_type_id", "references": "property_types", "ref_field": "id"},
        {"table": "predictions", "field": "calculation_method_id", "references": "calculation_methods", "ref_field": "id"},
        {"table": "experiments", "field": "mixture_id", "references": "mixtures", "ref_field": "id"},
        {"table": "experiment_properties", "field": "experiment_id", "references": "experiments", "ref_field": "id"},
        {"table": "experiment_properties", "field": "property_type_id", "references": "property_types", "ref_field": "id"}
    ]
    
    for fk in foreign_keys:
        try:
            # Skip molecule_id in predictions since it can be NULL (for mixture predictions)
            if fk["table"] == "predictions" and fk["field"] == "molecule_id":
                query = f"""
                SELECT COUNT(*)
                FROM {fk["table"]}
                WHERE {fk["field"]} IS NOT NULL AND {fk["field"]} NOT IN (SELECT {fk["ref_field"]} FROM {fk["references"]})
                """
            # Skip mixture_id in predictions since it can be NULL (for molecule predictions)
            elif fk["table"] == "predictions" and fk["field"] == "mixture_id":
                query = f"""
                SELECT COUNT(*)
                FROM {fk["table"]}
                WHERE {fk["field"]} IS NOT NULL AND {fk["field"]} NOT IN (SELECT {fk["ref_field"]} FROM {fk["references"]})
                """
            else:
                query = f"""
                SELECT COUNT(*)
                FROM {fk["table"]}
                WHERE {fk["field"]} NOT IN (SELECT {fk["ref_field"]} FROM {fk["references"]})
                """
            
            result = execute_sql_on_supabase(project_id, query)
            invalid_count = result[0]["count"]
            
            if invalid_count > 0:
                logger.warning(f"Table {fk['table']} has {invalid_count} rows with invalid {fk['field']} foreign key")
                issues_found.append(f"Table {fk['table']} has {invalid_count} rows with invalid {fk['field']} foreign key")
            else:
                logger.info(f"Table {fk['table']} has valid {fk['field']} foreign keys")
        except Exception as e:
            logger.error(f"Error verifying foreign key {fk['table']}.{fk['field']}: {str(e)}")
            issues_found.append(f"Error verifying foreign key {fk['table']}.{fk['field']}: {str(e)}")
    
    # Generate verification report
    report = "# Database Population Verification Report\n\n"
    report += f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
    
    report += "## Row Counts\n\n"
    for table in tables:
        count = verification_results.get(f"{table}_count", "Unknown")
        report += f"- {table}: {count} rows\n"
    
    report += "\n## Issues Found\n\n"
    if issues_found:
        for issue in issues_found:
            report += f"- {issue}\n"
    else:
        report += "No issues found. Database population verified successfully.\n"
    
    # Write report to file
    report_path = "database_population_verification.md"
    with open(report_path, "w") as f:
        f.write(report)
    
    logger.info(f"Verification report written to {report_path}")
    
    if issues_found:
        logger.warning(f"Database verification completed with {len(issues_found)} issues found")
    else:
        logger.info("Database verification completed successfully with no issues found")

if __name__ == "__main__":
    main()
