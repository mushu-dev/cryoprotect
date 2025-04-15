#!/usr/bin/env python3
"""
CryoProtect Analyzer - PubChem Cryoprotectant Data Importer (Supabase Version)

This script retrieves cryoprotectant data from PubChem, filters molecules based on
predefined criteria, scores them, and stores the results in a Supabase database.

Prerequisites:
- Python 3.6+ installed
- Supabase project with the CryoProtect schema applied
- supabase-py package installed (pip install supabase)
- python-dotenv package installed (pip install python-dotenv)
"""

import os
import time
import requests
import json
from datetime import datetime
import logging
from dotenv import load_dotenv
from supabase import create_client, Client

# Load environment variables from .env file
load_dotenv()

# Supabase connection
supabase_url = os.getenv("SUPABASE_URL")
supabase_key = os.getenv("SUPABASE_KEY")

if not supabase_url or not supabase_key:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

supabase: Client = create_client(supabase_url, supabase_key)

# Scoring Weights (Total = 200)
WEIGHTS = {
    "hydrogen_bonding": 50,
    "solubility_polarity": 40,
    "membrane_permeability": 40,
    "toxicity_biocompatibility": 30,
    "protein_stabilization": 20,
    "stability_reactivity": 10,
    "environmental_safety": 10
}

# Stage 1: Core Filtering Criteria
CORE_CRITERIA = {
    "logP_range": (-1.5, 0),
    "mw_range": (50, 150),
    "TPSA_range": (40, 100),
    "functional_groups": ["OH", "CONH2", "S=O"]
}

# PubChem API delay
PUBCHEM_API_DELAY = float(os.getenv("PUBCHEM_API_DELAY", "0.2"))

# CID File
CID_FILE = "CID-Synonym-filtered"  # Updated CID list location

# Set up logging
LOG_FILE = "cryoprotectant_analysis.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def get_cid_list():
    """Retrieve all CIDs from the downloaded PubChem CID list."""
    if not os.path.exists(CID_FILE):
        logger.warning(f"⚠️ CID file '{CID_FILE}' not found. Download it first.")
        return []

    with open(CID_FILE, "r") as file:
        cids = [int(line.strip().split("\t")[0]) for line in file if line.strip().split("\t")[0].isdigit()]

    logger.info(f"✅ Loaded {len(cids)} CIDs from PubChem's CID list.")
    return cids


def get_molecule_properties(cid):
    """Fetch molecular properties from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,IsomericSMILES/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        properties = data["PropertyTable"]["Properties"][0]
        return {
            "CID": cid,
            "Molecular Formula": properties.get("MolecularFormula"),
            "Molecular Weight": properties.get("MolecularWeight"),
            "LogP": properties.get("XLogP"),
            "TPSA": properties.get("TPSA"),
            "H-Bond Donors": properties.get("HBondDonorCount"),
            "H-Bond Acceptors": properties.get("HBondAcceptorCount"),
            "SMILES": properties.get("IsomericSMILES"),
            "PubChem Link": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        }
    logger.warning(f"⚠️ Warning: No molecular properties found for CID {cid}.")
    return {"CID": cid, "Error": "No data found"}


def get_additional_properties(cid):
    """Fetch toxicity, stability, and environmental impact from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        props = data["PC_Compounds"][0]["props"]
        
        extra_data = {
            "Toxicity": None,
            "Stability": None,
            "Environmental Safety": None
        }

        for prop in props:
            label = prop.get("urn", {}).get("label", "")
            value = prop.get("value", {}).get("sval", "")

            if "Toxicity" in label:
                extra_data["Toxicity"] = value
            elif "Stability" in label:
                extra_data["Stability"] = value
            elif "Environmental" in label:
                extra_data["Environmental Safety"] = value

        return extra_data
    logger.warning(f"⚠️ Warning: No additional properties found for CID {cid}.")
    return {"Error": "No additional data found"}


def filter_molecule(molecule):
    """Initial filtering based on core cryoprotectant properties, ensuring numerical values."""
    if "Error" in molecule:
        return False

    try:
        mw = float(molecule["Molecular Weight"]) if molecule["Molecular Weight"] else None
        logp = float(molecule["LogP"]) if molecule["LogP"] else None
        tpsa = float(molecule["TPSA"]) if molecule["TPSA"] else None
    except (ValueError, TypeError):
        logger.warning(f"⚠️ Skipping CID {molecule['CID']} due to invalid numerical values.")
        return False

    smiles = molecule["SMILES"]

    if mw is None or not (CORE_CRITERIA["mw_range"][0] <= mw <= CORE_CRITERIA["mw_range"][1]):
        return False
    if logp is None or not (CORE_CRITERIA["logP_range"][0] <= logp <= CORE_CRITERIA["logP_range"][1]):
        return False
    if tpsa is None or not (CORE_CRITERIA["TPSA_range"][0] <= tpsa <= CORE_CRITERIA["TPSA_range"][1]):
        return False
    if not any(group in smiles for group in CORE_CRITERIA["functional_groups"]):
        return False

    return True


def score_molecule(molecule, extra_properties):
    """Compute final score out of 200 based on all properties."""
    score = 0

    # Hydrogen Bonding
    score += WEIGHTS["hydrogen_bonding"]

    # Solubility & Permeability
    score += WEIGHTS["solubility_polarity"]
    score += WEIGHTS["membrane_permeability"]

    # Toxicity & Stability
    if extra_properties.get("Toxicity"):
        score += WEIGHTS["toxicity_biocompatibility"]
    if extra_properties.get("Stability"):
        score += WEIGHTS["stability_reactivity"]

    # Environmental Safety
    if extra_properties.get("Environmental Safety"):
        score += WEIGHTS["environmental_safety"]

    return score


def import_molecule_to_supabase(molecule, extra_properties, score):
    """Import molecule and its properties to Supabase database."""
    try:
        # Get the current user ID (if authenticated)
        user_id = supabase.auth.current_user.id if hasattr(supabase.auth, 'current_user') and supabase.auth.current_user else None

        # 1. First, insert the molecule
        response = supabase.table("molecules").insert({
            "cid": molecule["CID"],
            "name": f"Compound {molecule['CID']}",  # Default name
            "molecular_formula": molecule.get("Molecular Formula"),
            "smiles": molecule.get("SMILES"),
            "created_by": user_id
        }).execute()

        if response.error:
            logger.error(f"Error inserting molecule: {response.error}")
            return None

        molecule_id = response.data[0]["id"]
        logger.info(f"Inserted molecule with ID: {molecule_id}")

        # 2. Get property types
        response = supabase.table("property_types").select("id, name, data_type").execute()
        if response.error:
            logger.error(f"Error fetching property types: {response.error}")
            return None

        property_types = response.data

        # 3. Prepare property inserts
        property_inserts = []
        
        # Basic properties
        properties_to_insert = {
            "Molecular Weight": molecule.get("Molecular Weight"),
            "LogP": molecule.get("LogP"),
            "TPSA": molecule.get("TPSA"),
            "H-Bond Donors": molecule.get("H-Bond Donors"),
            "H-Bond Acceptors": molecule.get("H-Bond Acceptors"),
            "Toxicity": extra_properties.get("Toxicity"),
            "Stability": extra_properties.get("Stability"),
            "Environmental Safety": extra_properties.get("Environmental Safety"),
            "Total Score": score
        }
        
        for property_name, value in properties_to_insert.items():
            if value is None:
                continue
                
            property_type = next((pt for pt in property_types if pt["name"] == property_name), None)
            if not property_type:
                logger.warning(f"Property type '{property_name}' not found, skipping")
                continue
            
            property_insert = {
                "molecule_id": molecule_id,
                "property_type_id": property_type["id"],
                "created_by": user_id
            }
            
            # Set the appropriate value field based on data type
            if property_type["data_type"] == "numeric":
                try:
                    property_insert["numeric_value"] = float(value) if value is not None else None
                except (ValueError, TypeError):
                    logger.warning(f"Could not convert '{value}' to float for property '{property_name}', skipping")
                    continue
            elif property_type["data_type"] == "text":
                property_insert["text_value"] = str(value) if value is not None else None
            elif property_type["data_type"] == "boolean":
                property_insert["boolean_value"] = bool(value) if value is not None else None
            
            property_inserts.append(property_insert)
        
        # 4. Insert properties if we have any
        if property_inserts:
            response = supabase.table("molecular_properties").insert(property_inserts).execute()
            
            if response.error:
                logger.error(f"Error inserting properties: {response.error}")
                return None
            
            logger.info(f"Inserted {len(property_inserts)} properties for molecule {molecule_id}")
        
        return molecule_id
    
    except Exception as e:
        logger.error(f"Error in import_molecule_to_supabase: {str(e)}")
        return None


def process_molecules():
    """Process molecules from PubChem and store in Supabase."""
    logger.info("🚀 Starting full dataset processing...")

    # Fetch all CIDs from the file
    cids = get_cid_list()
    
    # Track statistics
    total_processed = 0
    total_imported = 0
    
    for cid in cids:
        try:
            # Get molecule properties
            molecule = get_molecule_properties(cid)
            
            # Apply initial filtering
            if filter_molecule(molecule):
                # Get additional properties
                extra_properties = get_additional_properties(cid)
                
                # Calculate score
                score = score_molecule(molecule, extra_properties)
                
                # Import to Supabase
                molecule_id = import_molecule_to_supabase(molecule, extra_properties, score)
                
                if molecule_id:
                    logger.info(f"✅ Processed and imported CID {cid}: Score = {score}")
                    total_imported += 1
                else:
                    logger.warning(f"⚠️ Failed to import CID {cid} to database")
            
            total_processed += 1
            
            # Sleep to avoid overwhelming the PubChem API
            time.sleep(PUBCHEM_API_DELAY)
            
        except Exception as e:
            logger.error(f"❌ Error processing CID {cid}: {str(e)}")
    
    logger.info(f"✅ Processing complete! Processed {total_processed} molecules, imported {total_imported} to database.")


if __name__ == "__main__":
    try:
        # Check if we need to authenticate
        supabase_user = os.getenv("SUPABASE_USER")
        supabase_password = os.getenv("SUPABASE_PASSWORD")
        
        if supabase_user and supabase_password:
            try:
                response = supabase.auth.sign_in_with_password({
                    "email": supabase_user,
                    "password": supabase_password
                })
                
                if response.error:
                    logger.warning(f"Authentication error: {response.error}")
                    logger.warning("Continuing without authentication. Some database operations may fail.")
                else:
                    logger.info(f"Authenticated as {supabase_user}")
            except Exception as e:
                logger.warning(f"Authentication error: {str(e)}")
                logger.warning("Continuing without authentication. Some database operations may fail.")
        else:
            logger.warning("No authentication credentials provided. Continuing without authentication.")
            logger.warning("Some database operations may fail due to Row Level Security (RLS) policies.")
        
        # Process molecules
        process_molecules()
        
    except Exception as e:
        logger.error(f"❌ Fatal error: {str(e)}")