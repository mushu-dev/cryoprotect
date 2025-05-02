#!/usr/bin/env python
"""
CryoProtect v2 - Test Data Loader

This script loads the test datasets (core_cryoprotectants.json, mixtures.json, and edge_cases.json)
into the CryoProtect database. It can be used for integration testing and development.

Usage:
    python load_test_data.py [--dataset core|mixtures|edge|all] [--clear]

Options:
    --dataset: Specify which dataset to load (default: all)
    --clear: Clear existing data before loading (use with caution)
"""

import os
import sys
import json
import argparse
import logging
import uuid
from datetime import datetime
from typing import Dict, List, Any, Optional

# Add parent directory to path to import from api
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

try:
    from api.models import Molecule, MolecularProperty, Mixture, MixtureComponent
    from api.rdkit_utils import calculate_all_properties
except ImportError:
    print("Error: Could not import required modules from CryoProtect API.")
    print("Make sure you're running this script from the CryoProtect project directory.")
    sys.exit(1)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constants
DATASETS_DIR = os.path.dirname(os.path.abspath(__file__))
CORE_DATASET = os.path.join(DATASETS_DIR, 'core_cryoprotectants.json')
MIXTURES_DATASET = os.path.join(DATASETS_DIR, 'mixtures.json')
EDGE_CASES_DATASET = os.path.join(DATASETS_DIR, 'edge_cases.json')

# Default user ID for test data
DEFAULT_USER_ID = "00000000-0000-0000-0000-000000000001"

def load_json_file(file_path: str) -> Dict[str, Any]:
    """Load a JSON file and return its contents."""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error loading {file_path}: {str(e)}")
        return {}

def clear_existing_data() -> None:
    """Clear existing test data from the database."""
    logger.info("Clearing existing test data...")
    
    # Delete test mixtures first (due to foreign key constraints)
    try:
        # Get test mixtures (those created by this script)
        test_mixtures = Mixture.get_all_by_creator(DEFAULT_USER_ID)
        for mixture in test_mixtures:
            # Delete mixture components
            MixtureComponent.delete_by_mixture_id(mixture["id"])
            # Delete mixture
            Mixture.delete(mixture["id"])
        logger.info(f"Deleted {len(test_mixtures)} test mixtures and their components.")
    except Exception as e:
        logger.error(f"Error deleting test mixtures: {str(e)}")
    
    # Delete test molecules
    try:
        # Get test molecules (those created by this script)
        test_molecules = Molecule.get_all_by_creator(DEFAULT_USER_ID)
        for molecule in test_molecules:
            # Delete molecular properties
            MolecularProperty.delete_by_molecule_id(molecule["id"])
            # Delete molecule
            Molecule.delete(molecule["id"])
        logger.info(f"Deleted {len(test_molecules)} test molecules and their properties.")
    except Exception as e:
        logger.error(f"Error deleting test molecules: {str(e)}")

def load_core_cryoprotectants() -> Dict[str, str]:
    """
    Load the core cryoprotectants dataset.
    
    Returns:
        Dictionary mapping molecule names to their IDs
    """
    logger.info("Loading core cryoprotectants dataset...")
    data = load_json_file(CORE_DATASET)
    
    if not data or "molecules" not in data:
        logger.error("Invalid core cryoprotectants dataset format.")
        return {}
    
    molecule_id_map = {}
    
    for molecule_data in data["molecules"]:
        try:
            # Create molecule record
            molecule_id = str(uuid.uuid4())
            
            molecule = {
                "id": molecule_id,
                "name": molecule_data["name"],
                "smiles": molecule_data["smiles"],
                "molecular_formula": molecule_data["formula"],
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat(),
                "created_by": DEFAULT_USER_ID
            }
            
            # Insert molecule
            Molecule.insert(molecule)
            molecule_id_map[molecule_data["name"]] = molecule_id
            
            # Calculate additional properties using RDKit
            try:
                rdkit_props = calculate_all_properties(molecule_data["smiles"])
            except Exception as e:
                logger.warning(f"Error calculating RDKit properties for {molecule_data['name']}: {str(e)}")
                rdkit_props = {}
            
            # Insert molecular properties
            properties_to_insert = [
                # Basic properties from dataset
                {"name": "Molecular Weight", "value": molecule_data["properties"]["molecular_weight"], "data_type": "numeric"},
                {"name": "LogP", "value": molecule_data["properties"]["logp"], "data_type": "numeric"},
                {"name": "TPSA", "value": molecule_data["properties"].get("tpsa", 0), "data_type": "numeric"},
                {"name": "H-Bond Donors", "value": molecule_data["properties"]["h_bond_donors"], "data_type": "numeric"},
                {"name": "H-Bond Acceptors", "value": molecule_data["properties"]["h_bond_acceptors"], "data_type": "numeric"},
                {"name": "Melting Point", "value": molecule_data["properties"].get("melting_point"), "data_type": "numeric"},
                {"name": "Boiling Point", "value": molecule_data["properties"].get("boiling_point"), "data_type": "numeric"},
                {"name": "Glass Transition Temp", "value": molecule_data["properties"].get("glass_transition_temp"), "data_type": "numeric"},
                
                # Cryoprotection scores (calculated based on properties)
                {"name": "Cryoprotection Hydrogen Bonding Score", "value": score_hydrogen_bonding(molecule_data), "data_type": "numeric"},
                {"name": "Cryoprotection Logp Score", "value": score_logp(molecule_data), "data_type": "numeric"},
                {"name": "Cryoprotection Molecular Size Score", "value": score_molecular_size(molecule_data), "data_type": "numeric"},
                {"name": "Cryoprotection Tpsa Score", "value": score_tpsa(molecule_data), "data_type": "numeric"},
                {"name": "Cryoprotection Score", "value": calculate_total_score(molecule_data), "data_type": "numeric"}
            ]
            
            for prop in properties_to_insert:
                if prop["value"] is not None:
                    property_id = str(uuid.uuid4())
                    
                    # Get property type ID
                    property_type = MolecularProperty.get_property_type_by_name(prop["name"])
                    if not property_type:
                        logger.warning(f"Property type {prop['name']} not found. Skipping.")
                        continue
                    
                    property_data = {
                        "id": property_id,
                        "molecule_id": molecule_id,
                        "property_type_id": property_type["id"],
                        "numeric_value": prop["value"] if prop["data_type"] == "numeric" else None,
                        "text_value": prop["value"] if prop["data_type"] == "text" else None,
                        "boolean_value": prop["value"] if prop["data_type"] == "boolean" else None,
                        "created_at": datetime.now().isoformat(),
                        "updated_at": datetime.now().isoformat(),
                        "created_by": DEFAULT_USER_ID
                    }
                    
                    MolecularProperty.insert(property_data)
            
            logger.info(f"Added molecule: {molecule_data['name']}")
            
        except Exception as e:
            logger.error(f"Error adding molecule {molecule_data['name']}: {str(e)}")
    
    logger.info(f"Successfully loaded {len(molecule_id_map)} core cryoprotectants.")
    return molecule_id_map

def load_edge_cases() -> Dict[str, str]:
    """
    Load the edge cases dataset.
    
    Returns:
        Dictionary mapping molecule names to their IDs
    """
    logger.info("Loading edge cases dataset...")
    data = load_json_file(EDGE_CASES_DATASET)
    
    if not data or "molecules" not in data:
        logger.error("Invalid edge cases dataset format.")
        return {}
    
    molecule_id_map = {}
    
    for molecule_data in data["molecules"]:
        try:
            # Create molecule record
            molecule_id = str(uuid.uuid4())
            
            molecule = {
                "id": molecule_id,
                "name": molecule_data["name"],
                "smiles": molecule_data["smiles"],
                "molecular_formula": molecule_data["formula"],
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat(),
                "created_by": DEFAULT_USER_ID
            }
            
            # Insert molecule
            Molecule.insert(molecule)
            molecule_id_map[molecule_data["name"]] = molecule_id
            
            # Calculate additional properties using RDKit
            try:
                rdkit_props = calculate_all_properties(molecule_data["smiles"])
            except Exception as e:
                logger.warning(f"Error calculating RDKit properties for {molecule_data['name']}: {str(e)}")
                rdkit_props = {}
            
            # Insert molecular properties
            properties_to_insert = [
                # Basic properties from dataset
                {"name": "Molecular Weight", "value": molecule_data["properties"]["molecular_weight"], "data_type": "numeric"},
                {"name": "LogP", "value": molecule_data["properties"]["logp"], "data_type": "numeric"},
                {"name": "TPSA", "value": molecule_data["properties"].get("tpsa", 0), "data_type": "numeric"},
                {"name": "H-Bond Donors", "value": molecule_data["properties"]["h_bond_donors"], "data_type": "numeric"},
                {"name": "H-Bond Acceptors", "value": molecule_data["properties"]["h_bond_acceptors"], "data_type": "numeric"},
                {"name": "Melting Point", "value": molecule_data["properties"].get("melting_point"), "data_type": "numeric"},
                {"name": "Boiling Point", "value": molecule_data["properties"].get("boiling_point"), "data_type": "numeric"},
                {"name": "Glass Transition Temp", "value": molecule_data["properties"].get("glass_transition_temp"), "data_type": "numeric"},
                {"name": "Edge Case Type", "value": molecule_data.get("edge_case_type", ""), "data_type": "text"},
                
                # Cryoprotection scores (calculated based on properties)
                {"name": "Cryoprotection Hydrogen Bonding Score", "value": score_hydrogen_bonding(molecule_data), "data_type": "numeric"},
                {"name": "Cryoprotection Logp Score", "value": score_logp(molecule_data), "data_type": "numeric"},
                {"name": "Cryoprotection Molecular Size Score", "value": score_molecular_size(molecule_data), "data_type": "numeric"},
                {"name": "Cryoprotection Tpsa Score", "value": score_tpsa(molecule_data), "data_type": "numeric"},
                {"name": "Cryoprotection Score", "value": calculate_total_score(molecule_data), "data_type": "numeric"}
            ]
            
            for prop in properties_to_insert:
                if prop["value"] is not None:
                    property_id = str(uuid.uuid4())
                    
                    # Get property type ID
                    property_type = MolecularProperty.get_property_type_by_name(prop["name"])
                    if not property_type:
                        logger.warning(f"Property type {prop['name']} not found. Skipping.")
                        continue
                    
                    property_data = {
                        "id": property_id,
                        "molecule_id": molecule_id,
                        "property_type_id": property_type["id"],
                        "numeric_value": prop["value"] if prop["data_type"] == "numeric" else None,
                        "text_value": prop["value"] if prop["data_type"] == "text" else None,
                        "boolean_value": prop["value"] if prop["data_type"] == "boolean" else None,
                        "created_at": datetime.now().isoformat(),
                        "updated_at": datetime.now().isoformat(),
                        "created_by": DEFAULT_USER_ID
                    }
                    
                    MolecularProperty.insert(property_data)
            
            logger.info(f"Added edge case molecule: {molecule_data['name']}")
            
        except Exception as e:
            logger.error(f"Error adding edge case molecule {molecule_data['name']}: {str(e)}")
    
    logger.info(f"Successfully loaded {len(molecule_id_map)} edge case molecules.")
    return molecule_id_map

def load_mixtures(molecule_id_map: Dict[str, str]) -> None:
    """
    Load the mixtures dataset.
    
    Args:
        molecule_id_map: Dictionary mapping molecule names to their IDs
    """
    logger.info("Loading mixtures dataset...")
    data = load_json_file(MIXTURES_DATASET)
    
    if not data or "mixtures" not in data:
        logger.error("Invalid mixtures dataset format.")
        return
    
    for mixture_data in data["mixtures"]:
        try:
            # Create mixture record
            mixture_id = str(uuid.uuid4())
            
            mixture = {
                "id": mixture_id,
                "name": mixture_data["name"],
                "description": mixture_data["description"],
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat(),
                "created_by": DEFAULT_USER_ID
            }
            
            # Insert mixture
            Mixture.insert(mixture)
            
            # Insert mixture components
            for component in mixture_data["components"]:
                molecule_name = component["molecule_name"]
                
                if molecule_name not in molecule_id_map:
                    logger.warning(f"Molecule {molecule_name} not found in molecule ID map. Skipping component.")
                    continue
                
                component_id = str(uuid.uuid4())
                
                component_data = {
                    "id": component_id,
                    "mixture_id": mixture_id,
                    "molecule_id": molecule_id_map[molecule_name],
                    "concentration": component["concentration"],
                    "concentration_unit": component["concentration_unit"],
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat(),
                    "created_by": DEFAULT_USER_ID
                }
                
                MixtureComponent.insert(component_data)
            
            logger.info(f"Added mixture: {mixture_data['name']}")
            
        except Exception as e:
            logger.error(f"Error adding mixture {mixture_data['name']}: {str(e)}")
    
    logger.info(f"Successfully loaded {len(data['mixtures'])} mixtures.")

def score_hydrogen_bonding(molecule_data: Dict[str, Any]) -> float:
    """Calculate hydrogen bonding score for cryoprotection."""
    donors = molecule_data["properties"]["h_bond_donors"]
    acceptors = molecule_data["properties"]["h_bond_acceptors"]
    total = donors + acceptors
    
    # Ideal range is 3-8 total H-bonds
    if total < 3:
        return 3.0 + (total / 3.0) * 4.0
    elif total <= 8:
        return 7.0 + ((8 - total) / 5.0) * 3.0
    else:
        return 7.0 - ((total - 8) / 10.0) * 4.0

def score_logp(molecule_data: Dict[str, Any]) -> float:
    """Calculate LogP score for cryoprotection."""
    logp = molecule_data["properties"]["logp"]
    
    # Ideal range is -2 to 1
    if logp < -4:
        return 4.0 + ((logp + 4) / 2.0) * 3.0
    elif logp < -2:
        return 7.0 + ((logp + 2) / 2.0) * 3.0
    elif logp <= 1:
        return 10.0 - ((logp + 2) / 3.0) * 3.0
    else:
        return 7.0 - ((logp - 1) / 4.0) * 4.0

def score_molecular_size(molecule_data: Dict[str, Any]) -> float:
    """Calculate molecular size score for cryoprotection."""
    mw = molecule_data["properties"]["molecular_weight"]
    
    # Ideal range is 60-180 Da
    if mw < 30:
        return 3.0 + (mw / 30.0) * 4.0
    elif mw < 60:
        return 7.0 + ((mw - 30) / 30.0) * 3.0
    elif mw <= 180:
        return 10.0 - ((mw - 60) / 120.0) * 2.0
    elif mw <= 400:
        return 8.0 - ((mw - 180) / 220.0) * 4.0
    else:
        return 4.0 - ((mw - 400) / 600.0) * 4.0

def score_tpsa(molecule_data: Dict[str, Any]) -> float:
    """Calculate TPSA score for cryoprotection."""
    tpsa = molecule_data["properties"].get("tpsa", 0)
    
    # Ideal range is 40-90 Å²
    if tpsa < 20:
        return 5.0 + (tpsa / 20.0) * 3.0
    elif tpsa < 40:
        return 8.0 + ((tpsa - 20) / 20.0) * 2.0
    elif tpsa <= 90:
        return 10.0 - ((tpsa - 40) / 50.0) * 1.0
    elif tpsa <= 140:
        return 9.0 - ((tpsa - 90) / 50.0) * 3.0
    else:
        return 6.0 - ((tpsa - 140) / 60.0) * 3.0

def calculate_total_score(molecule_data: Dict[str, Any]) -> float:
    """Calculate overall cryoprotection score."""
    h_bond_score = score_hydrogen_bonding(molecule_data)
    logp_score = score_logp(molecule_data)
    size_score = score_molecular_size(molecule_data)
    tpsa_score = score_tpsa(molecule_data)
    
    # Weighted average
    weights = {
        "h_bond": 0.35,
        "logp": 0.25,
        "size": 0.20,
        "tpsa": 0.20
    }
    
    total_score = (
        h_bond_score * weights["h_bond"] +
        logp_score * weights["logp"] +
        size_score * weights["size"] +
        tpsa_score * weights["tpsa"]
    )
    
    return round(total_score, 1)

def main():
    """Main function to load test data."""
    parser = argparse.ArgumentParser(description="Load test data into CryoProtect database")
    parser.add_argument("--dataset", choices=["core", "mixtures", "edge", "all"], default="all",
                        help="Specify which dataset to load (default: all)")
    parser.add_argument("--clear", action="store_true", help="Clear existing data before loading")
    
    args = parser.parse_args()
    
    if args.clear:
        clear_existing_data()
    
    molecule_id_map = {}
    
    if args.dataset in ["core", "all"]:
        core_id_map = load_core_cryoprotectants()
        molecule_id_map.update(core_id_map)
    
    if args.dataset in ["edge", "all"]:
        edge_id_map = load_edge_cases()
        molecule_id_map.update(edge_id_map)
    
    if args.dataset in ["mixtures", "all"] and molecule_id_map:
        load_mixtures(molecule_id_map)
    
    logger.info("Test data loading complete!")

if __name__ == "__main__":
    main()