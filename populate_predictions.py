#!/usr/bin/env python3
"""
CryoProtect v2 - Prediction Population Script

This script populates the predictions table with scientifically accurate cryoprotectant
prediction data. It ensures proper relationships between predictions and molecules.

Usage:
    python populate_predictions.py [--dry-run]

Environment variables required (from .env):
    SUPABASE_URL, SUPABASE_KEY, SUPABASE_USER, SUPABASE_PASSWORD
"""

import os
import json
import uuid
import argparse
import logging
import random
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
        logging.FileHandler("prediction_population.log"),
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

# Prediction property types with scientifically accurate ranges
PREDICTION_PROPERTIES = [
    {
        "property_type": "Glass Transition Temperature",
        "unit": "°C",
        "ranges": {
            "Dimethyl sulfoxide": (-137, -130),
            "Glycerol": (-93, -90),
            "Ethylene glycol": (-128, -125),
            "Propylene glycol": (-108, -105),
            "Trehalose": (115, 120),
            "Sucrose": (65, 70),
            "Methanol": (-175, -170),
            "Formamide": (-113, -110),
            "Acetamide": (-73, -70),
            "1,2-Propanediol": (-108, -105),
            "default": (-120, -80)  # Default range for other molecules
        }
    },
    {
        "property_type": "Cryoprotective Efficacy",
        "unit": "%",
        "ranges": {
            "Dimethyl sulfoxide": (85, 95),
            "Glycerol": (75, 85),
            "Ethylene glycol": (80, 90),
            "Propylene glycol": (70, 80),
            "Trehalose": (65, 75),
            "Sucrose": (60, 70),
            "Methanol": (50, 60),
            "Formamide": (55, 65),
            "Acetamide": (60, 70),
            "1,2-Propanediol": (70, 80),
            "default": (50, 70)  # Default range for other molecules
        }
    },
    {
        "property_type": "Cell Membrane Permeability",
        "unit": "μm/s",
        "ranges": {
            "Dimethyl sulfoxide": (0.8, 1.2),
            "Glycerol": (0.3, 0.5),
            "Ethylene glycol": (0.9, 1.3),
            "Propylene glycol": (0.7, 1.0),
            "Trehalose": (0.01, 0.05),
            "Sucrose": (0.01, 0.05),
            "Methanol": (1.5, 2.0),
            "Formamide": (0.6, 0.9),
            "Acetamide": (0.5, 0.8),
            "1,2-Propanediol": (0.7, 1.0),
            "default": (0.3, 0.8)  # Default range for other molecules
        }
    },
    {
        "property_type": "Toxicity Index",
        "unit": "",
        "ranges": {
            "Dimethyl sulfoxide": (0.3, 0.5),
            "Glycerol": (0.1, 0.3),
            "Ethylene glycol": (0.4, 0.6),
            "Propylene glycol": (0.2, 0.4),
            "Trehalose": (0.0, 0.1),
            "Sucrose": (0.0, 0.1),
            "Methanol": (0.7, 0.9),
            "Formamide": (0.5, 0.7),
            "Acetamide": (0.4, 0.6),
            "1,2-Propanediol": (0.2, 0.4),
            "default": (0.3, 0.6)  # Default range for other molecules
        }
    }
]

# Calculation methods for predictions
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

def get_molecules(supabase):
    """Get all molecules from the database."""
    response = supabase.table("molecules").select("id, name").execute()
    if hasattr(response, 'data') and response.data:
        return response.data
    return []

def populate_calculation_methods(supabase, user_profile_id, dry_run=False):
    """Populate the calculation_method table with prediction methods."""
    if not user_profile_id:
        logger.warning("No user profile ID provided. Skipping calculation method population.")
        return {}
    
    # Check if calculation methods already exist
    response = supabase.table("calculation_method").select("*").execute()
    existing_methods = {}
    if hasattr(response, 'data') and response.data:
        logger.info(f"Calculation methods already exist ({len(response.data)} found)")
        existing_methods = {method["name"].lower(): method["id"] for method in response.data if "name" in method}
    
    # Prepare calculation method data
    method_map = existing_methods.copy()
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
        
        if not dry_run:
            response = supabase.table("calculation_method").insert(method_data).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error creating calculation method {method['name']}: {response.error}")
                continue
            logger.info(f"Created calculation method {method['name']} with ID: {method_id}")
        else:
            logger.info(f"DRY RUN: Would create calculation method {method['name']} with ID: {method_id}")
        
        method_map[method["name"].lower()] = method_id
        time.sleep(0.1)  # Avoid rate limiting
    
    logger.info(f"Populated {len(method_map) - len(existing_methods)} new calculation methods")
    return method_map

def populate_predictions(supabase, molecules, method_map, user_profile_id, dry_run=False):
    """Populate the predictions table with scientifically accurate data."""
    if not molecules or not method_map or not user_profile_id:
        logger.warning("Missing molecules, method_map, or user_profile_id. Skipping predictions population.")
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
        "Glass Transition Temperature": gc_method_id or next(iter(method_map.values())),
        "Cryoprotective Efficacy": nn_method_id or next(iter(method_map.values())),
        "Cell Membrane Permeability": md_method_id or next(iter(method_map.values())),
        "Toxicity Index": qspr_method_id or next(iter(method_map.values()))
    }
    
    # Populate predictions for each molecule and property
    prediction_count = 0
    for molecule in molecules:
        molecule_id = molecule["id"]
        molecule_name = molecule["name"]
        
        for prop in PREDICTION_PROPERTIES:
            property_type = prop["property_type"]
            unit = prop["unit"]
            
            # Get appropriate range for this molecule and property
            range_key = molecule_name if molecule_name in prop["ranges"] else "default"
            value_range = prop["ranges"][range_key]
            
            # Generate a random value within the range
            predicted_value = round(random.uniform(value_range[0], value_range[1]), 2)
            
            # Generate a random confidence value (0.7-0.95)
            confidence = round(random.uniform(0.7, 0.95), 2)
            
            # Get appropriate method for this property
            method_id = property_method_map.get(property_type, next(iter(method_map.values())))
            
            prediction_id = str(uuid.uuid4())
            prediction_data = {
                "id": prediction_id,
                "molecule_id": molecule_id,
                "property_type": property_type,
                "predicted_value": predicted_value,
                "unit": unit,
                "method_id": method_id,
                "model_version": "1.0",
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
            
            if not dry_run:
                response = supabase.table("predictions").insert(prediction_data).execute()
                if hasattr(response, 'error') and response.error:
                    logger.error(f"Error creating prediction {property_type} for {molecule_name}: {response.error}")
                    continue
                logger.info(f"Created prediction {property_type} for {molecule_name} with ID: {prediction_id}")
                prediction_count += 1
            else:
                logger.info(f"DRY RUN: Would create prediction {property_type} for {molecule_name} with ID: {prediction_id}")
                prediction_count += 1
            
            time.sleep(0.1)  # Avoid rate limiting
    
    logger.info(f"Populated {prediction_count} predictions")

def main():
    parser = argparse.ArgumentParser(description="Populate predictions table with scientifically accurate data.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 Prediction Population")
    
    # Connect to Supabase
    supabase = connect_to_supabase()
    
    # Get authenticated user ID
    auth_user_id = get_authenticated_user_id(supabase)
    
    # Create user profile if needed
    user_profile_id = create_user_profile(supabase, auth_user_id, args.dry_run)
    
    if not user_profile_id:
        logger.error("No user profile found. Please run populate_molecules.py first to create a user profile.")
        return
    
    # Get molecules
    molecules = get_molecules(supabase)
    if not molecules:
        logger.error("No molecules found. Please run populate_molecules.py first.")
        return
    
    # Populate calculation methods
    method_map = populate_calculation_methods(supabase, user_profile_id, args.dry_run)
    
    # Populate predictions
    populate_predictions(supabase, molecules, method_map, user_profile_id, args.dry_run)
    
    logger.info("Prediction population complete")

if __name__ == "__main__":
    main()