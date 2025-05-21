#!/usr/bin/env python3
"""
CryoProtect v2 - Experiment Population Script

This script populates the experiments and experiment_properties tables with scientifically accurate
cryoprotectant experiment data. It ensures proper relationships between experiments and mixtures.

Usage:
    python populate_experiments.py [--dry-run]

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
        logging.FileHandler("experiment_population.log"),
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

# Scientifically accurate experiments
EXPERIMENTS = [
    {
        "name": "Slow Freezing Protocol - DMSO/EG",
        "description": "Standard slow freezing protocol using DMSO/EG mixture",
        "mixture_name": "DMSO/EG Vitrification Solution",
        "preparation_protocol": "1. Equilibrate cells in 10% DMSO/EG solution for 15 min at 4°C\n2. Cool at 1°C/min to -80°C\n3. Plunge into liquid nitrogen",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Cell Viability", "value": 78.5, "unit": "%"},
            {"property_type": "Ice Crystal Formation", "value": 12.3, "unit": "%"},
            {"property_type": "Recovery Rate", "value": 65.2, "unit": "%"}
        ]
    },
    {
        "name": "Vitrification Protocol - Trehalose/Glycerol",
        "description": "Rapid vitrification protocol using Trehalose/Glycerol mixture",
        "mixture_name": "Trehalose/Glycerol Mixture",
        "preparation_protocol": "1. Equilibrate cells in 20% Trehalose/Glycerol solution for 5 min at room temperature\n2. Plunge directly into liquid nitrogen",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Cell Viability", "value": 92.1, "unit": "%"},
            {"property_type": "Ice Crystal Formation", "value": 0.5, "unit": "%"},
            {"property_type": "Recovery Rate", "value": 88.7, "unit": "%"}
        ]
    },
    {
        "name": "Controlled Rate Freezing - EAFS",
        "description": "Controlled rate freezing using EAFS solution",
        "mixture_name": "EAFS Solution",
        "preparation_protocol": "1. Equilibrate cells in EAFS solution for 10 min at 4°C\n2. Cool at 0.5°C/min to -40°C\n3. Cool at 10°C/min to -80°C\n4. Plunge into liquid nitrogen",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Cell Viability", "value": 85.7, "unit": "%"},
            {"property_type": "Ice Crystal Formation", "value": 5.2, "unit": "%"},
            {"property_type": "Recovery Rate", "value": 79.3, "unit": "%"}
        ]
    },
    {
        "name": "Oocyte Cryopreservation - PROH/Sucrose",
        "description": "Specialized protocol for oocyte cryopreservation using PROH/Sucrose",
        "mixture_name": "PROH/Sucrose Solution",
        "preparation_protocol": "1. Expose oocytes to 1.5M PROH for 10 min\n2. Transfer to 1.5M PROH + 0.2M sucrose for 5 min\n3. Plunge into liquid nitrogen",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Cell Viability", "value": 81.3, "unit": "%"},
            {"property_type": "Ice Crystal Formation", "value": 3.7, "unit": "%"},
            {"property_type": "Recovery Rate", "value": 72.8, "unit": "%"},
            {"property_type": "Fertilization Rate", "value": 68.5, "unit": "%"}
        ]
    },
    {
        "name": "Fish Sperm Cryopreservation - Methanol/Glycerol",
        "description": "Protocol optimized for fish sperm using Methanol/Glycerol mixture",
        "mixture_name": "Methanol/Glycerol Fish Cryoprotectant",
        "preparation_protocol": "1. Mix sperm with Methanol/Glycerol solution at 1:3 ratio\n2. Equilibrate for 2 min at 4°C\n3. Freeze in 0.5mL straws at 30°C/min to -80°C\n4. Store in liquid nitrogen",
        "temperature": -196,
        "temperature_unit": "°C",
        "pressure": 101.3,
        "pressure_unit": "kPa",
        "properties": [
            {"property_type": "Motility", "value": 73.2, "unit": "%"},
            {"property_type": "Membrane Integrity", "value": 68.9, "unit": "%"},
            {"property_type": "Fertilization Capacity", "value": 65.4, "unit": "%"}
        ]
    }
]

# Calculation methods for experiment properties
CALCULATION_METHODS = [
    {
        "name": "Differential Scanning Calorimetry",
        "description": "Experimental method to measure glass transition temperature and other thermal properties",
        "method_type": "experimental",
        "reference": "Wowk, B. Cryobiology 60.1 (2010): 11-22."
    },
    {
        "name": "Flow Cytometry",
        "description": "Method for measuring cell viability and membrane integrity",
        "method_type": "experimental",
        "reference": "Fuller, B. J. CryoLetters 25.6 (2004): 375-388."
    },
    {
        "name": "Light Microscopy",
        "description": "Visual assessment of ice crystal formation and cell morphology",
        "method_type": "experimental",
        "reference": "Seki, S. & Mazur, P. Biology of Reproduction 78.2 (2008): 186-195."
    },
    {
        "name": "Experimental Validation",
        "description": "Direct experimental measurement of cryoprotective efficacy",
        "method_type": "experimental",
        "reference": "Fahy, G. M. et al. Cryobiology 21.4 (1984): 407-426."
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

def get_project_id(supabase):
    """Get the first available project ID."""
    response = supabase.table("project").select("id").execute()
    if hasattr(response, 'data') and response.data:
        return response.data[0]["id"]
    return None

def get_mixture_map(supabase):
    """Get a mapping of mixture names to IDs."""
    response = supabase.table("mixtures").select("id, name").execute()
    if hasattr(response, 'data') and response.data:
        return {mix["name"].lower(): mix["id"] for mix in response.data if "name" in mix}
    return {}

def populate_calculation_methods(supabase, user_profile_id, dry_run=False):
    """Populate the calculation_method table."""
    if not user_profile_id:
        logger.warning("No user profile ID provided. Skipping calculation method population.")
        return {}
    
    # Check if calculation methods already exist
    response = supabase.table("calculation_method").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Calculation methods already exist ({len(response.data)} found)")
        return {method["name"].lower(): method["id"] for method in response.data if "name" in method}
    
    # Prepare calculation method data
    method_map = {}
    for method in CALCULATION_METHODS:
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
    
    logger.info(f"Populated {len(method_map)} calculation methods")
    return method_map

def populate_experiments(supabase, mixture_map, user_profile_id, project_id, dry_run=False):
    """Populate the experiments table with scientifically accurate data."""
    if not mixture_map or not user_profile_id:
        logger.warning("Missing mixture_map or user_profile_id. Skipping experiments population.")
        return {}
    
    # Check if experiments already exist
    response = supabase.table("experiments").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Experiments already exist ({len(response.data)} found)")
        return {exp["name"].lower(): exp["id"] for exp in response.data if "name" in exp}
    
    # Prepare experiment data
    experiment_map = {}
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
            "project_id": project_id,
            "preparation_protocol": exp["preparation_protocol"],
            "temperature": exp["temperature"],
            "temperature_unit": exp["temperature_unit"],
            "pressure": exp["pressure"],
            "pressure_unit": exp["pressure_unit"],
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
            response = supabase.table("experiments").insert(experiment_data).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error creating experiment {exp['name']}: {response.error}")
                continue
            logger.info(f"Created experiment {exp['name']} with ID: {experiment_id}")
        else:
            logger.info(f"DRY RUN: Would create experiment {exp['name']} with ID: {experiment_id}")
        
        experiment_map[exp["name"].lower()] = experiment_id
        time.sleep(0.1)  # Avoid rate limiting
    
    logger.info(f"Populated {len(experiment_map)} experiments")
    return experiment_map

def populate_experiment_properties(supabase, experiment_map, method_map, user_profile_id, dry_run=False):
    """Populate the experiment_properties table with properties for each experiment."""
    if not experiment_map or not method_map or not user_profile_id:
        logger.warning("Missing experiment_map, method_map, or user_profile_id. Skipping property population.")
        return
    
    # Check if properties already exist
    response = supabase.table("experiment_properties").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Experiment properties already exist ({len(response.data)} found)")
        return
    
    # Get a default method ID for properties
    default_method_id = next(iter(method_map.values())) if method_map else None
    
    # Populate properties for each experiment
    property_count = 0
    for exp in EXPERIMENTS:
        experiment_id = experiment_map.get(exp["name"].lower())
        if not experiment_id:
            logger.warning(f"Experiment ID not found for {exp['name']}. Skipping properties.")
            continue
        
        for prop in exp["properties"]:
            # Try to find an appropriate method for this property type
            method_id = None
            if "Cell Viability" in prop["property_type"] or "Membrane" in prop["property_type"]:
                method_id = method_map.get("flow cytometry", default_method_id)
            elif "Ice Crystal" in prop["property_type"]:
                method_id = method_map.get("light microscopy", default_method_id)
            else:
                method_id = method_map.get("experimental validation", default_method_id)
            
            property_id = str(uuid.uuid4())
            property_data = {
                "id": property_id,
                "experiment_id": experiment_id,
                "property_type": prop["property_type"],
                "value": prop["value"],
                "unit": prop["unit"],
                "method_id": method_id,
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
                response = supabase.table("experiment_properties").insert(property_data).execute()
                if hasattr(response, 'error') and response.error:
                    logger.error(f"Error creating property {prop['property_type']} for {exp['name']}: {response.error}")
                    continue
                logger.info(f"Created property {prop['property_type']} for {exp['name']} with ID: {property_id}")
                property_count += 1
            else:
                logger.info(f"DRY RUN: Would create property {prop['property_type']} for {exp['name']} with ID: {property_id}")
                property_count += 1
            
            time.sleep(0.1)  # Avoid rate limiting
    
    logger.info(f"Populated {property_count} experiment properties")

def main():
    parser = argparse.ArgumentParser(description="Populate experiments and experiment_properties tables with scientifically accurate data.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 Experiment Population")
    
    # Connect to Supabase
    supabase = connect_to_supabase()
    
    # Get authenticated user ID
    auth_user_id = get_authenticated_user_id(supabase)
    
    # Create user profile if needed
    user_profile_id = create_user_profile(supabase, auth_user_id, args.dry_run)
    
    if not user_profile_id:
        logger.error("No user profile found. Please run populate_molecules.py first to create a user profile.")
        return
    
    # Get project ID
    project_id = get_project_id(supabase)
    if not project_id:
        logger.warning("No project found. Experiments will be created without project association.")
    
    # Get mixture map
    mixture_map = get_mixture_map(supabase)
    if not mixture_map:
        logger.error("No mixtures found. Please run populate_mixtures.py first.")
        return
    
    # Populate calculation methods
    method_map = populate_calculation_methods(supabase, user_profile_id, args.dry_run)
    
    # Populate experiments
    experiment_map = populate_experiments(supabase, mixture_map, user_profile_id, project_id, args.dry_run)
    
    # Populate experiment properties
    populate_experiment_properties(supabase, experiment_map, method_map, user_profile_id, args.dry_run)
    
    logger.info("Experiment population complete")

if __name__ == "__main__":
    main()