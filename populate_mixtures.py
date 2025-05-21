#!/usr/bin/env python3
"""
CryoProtect v2 - Mixture Population Script

This script populates the mixtures and mixture_components tables with scientifically accurate
cryoprotectant mixture data. It ensures proper relationships between mixtures and molecules.

Usage:
    python populate_mixtures.py [--dry-run]

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
        logging.FileHandler("mixture_population.log"),
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

# Scientifically accurate mixtures
MIXTURES = [
    {
        "name": "DMSO/EG Vitrification Solution",
        "description": "1:1 molar mixture of DMSO and ethylene glycol, widely used for oocyte and embryo vitrification.",
        "reference": "Fahy et al., Cryobiology 21.4 (1984): 407-426; Rall & Fahy, Nature 313 (1985): 573-575.",
        "components": [
            {"name": "Dimethyl sulfoxide", "amount": 7.5, "amount_unit": "mol/L", "role": "cryoprotectant"},
            {"name": "Ethylene glycol", "amount": 7.5, "amount_unit": "mol/L", "role": "cryoprotectant"},
        ]
    },
    {
        "name": "Trehalose/Glycerol Mixture",
        "description": "2:1 mass ratio of trehalose to glycerol, used for slow freezing of cells.",
        "reference": "Best, B. P., Rejuvenation Research 18.5 (2015): 422-436.",
        "components": [
            {"name": "Trehalose", "amount": 0.2, "amount_unit": "g/mL", "role": "cryoprotectant"},
            {"name": "Glycerol", "amount": 0.1, "amount_unit": "g/mL", "role": "cryoprotectant"},
        ]
    },
    {
        "name": "EAFS Solution",
        "description": "EAFS (ethylene glycol, acetamide, Ficoll, sucrose) solution for embryo vitrification.",
        "reference": "Kasai et al., Biology of Reproduction 28.3 (1983): 687-693.",
        "components": [
            {"name": "Ethylene glycol", "amount": 3.4, "amount_unit": "mol/L", "role": "cryoprotectant"},
            {"name": "Acetamide", "amount": 2.0, "amount_unit": "mol/L", "role": "cryoprotectant"},
            {"name": "Sucrose", "amount": 0.3, "amount_unit": "mol/L", "role": "osmoprotectant"},
        ]
    },
    {
        "name": "PROH/Sucrose Solution",
        "description": "1.5 mol/L 1,2-propanediol (PROH) with 0.2 mol/L sucrose, used for oocyte cryopreservation.",
        "reference": "Fabbri et al., Human Reproduction 16.7 (2001): 1469-1476.",
        "components": [
            {"name": "1,2-Propanediol", "amount": 1.5, "amount_unit": "mol/L", "role": "cryoprotectant"},
            {"name": "Sucrose", "amount": 0.2, "amount_unit": "mol/L", "role": "osmoprotectant"},
        ]
    },
    {
        "name": "Methanol/Glycerol Fish Cryoprotectant",
        "description": "1.3 mol/L methanol and 0.5 mol/L glycerol, used for fish sperm cryopreservation.",
        "reference": "Cabrita et al., Aquaculture 249.1-4 (2005): 533-546.",
        "components": [
            {"name": "Methanol", "amount": 1.3, "amount_unit": "mol/L", "role": "cryoprotectant"},
            {"name": "Glycerol", "amount": 0.5, "amount_unit": "mol/L", "role": "cryoprotectant"},
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

def get_project_id(supabase):
    """Get the first available project ID."""
    response = supabase.table("project").select("id").execute()
    if hasattr(response, 'data') and response.data:
        return response.data[0]["id"]
    return None

def get_molecule_map(supabase):
    """Get a mapping of molecule names to IDs."""
    response = supabase.table("molecules").select("id, name").execute()
    if hasattr(response, 'data') and response.data:
        return {mol["name"].lower(): mol["id"] for mol in response.data if "name" in mol}
    return {}

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
        
        if not dry_run:
            response = supabase.table("mixtures").insert(mixture_data).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error creating mixture {mix['name']}: {response.error}")
                continue
            logger.info(f"Created mixture {mix['name']} with ID: {mixture_id}")
        else:
            logger.info(f"DRY RUN: Would create mixture {mix['name']} with ID: {mixture_id}")
        
        mixture_map[mix["name"].lower()] = mixture_id
        time.sleep(0.1)  # Avoid rate limiting
    
    logger.info(f"Populated {len(mixture_map)} mixtures")
    return mixture_map

def populate_mixture_components(supabase, mixture_map, molecule_map, dry_run=False):
    """Populate the mixture_components table with components for each mixture."""
    if not mixture_map or not molecule_map:
        logger.warning("Missing mixture_map or molecule_map. Skipping component population.")
        return
    
    # Check if components already exist
    response = supabase.table("mixture_components").select("*").execute()
    if hasattr(response, 'data') and response.data:
        logger.info(f"Mixture components already exist ({len(response.data)} found)")
        return
    
    # Populate components for each mixture
    component_count = 0
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
                "amount": comp["amount"],
                "amount_unit": comp["amount_unit"],
                "role": comp["role"],
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat()
            }
            
            if not dry_run:
                response = supabase.table("mixture_components").insert(component_data).execute()
                if hasattr(response, 'error') and response.error:
                    logger.error(f"Error creating component {comp['name']} for {mix['name']}: {response.error}")
                    continue
                logger.info(f"Created component {comp['name']} for {mix['name']} with ID: {component_id}")
                component_count += 1
            else:
                logger.info(f"DRY RUN: Would create component {comp['name']} for {mix['name']} with ID: {component_id}")
                component_count += 1
            
            time.sleep(0.1)  # Avoid rate limiting
    
    logger.info(f"Populated {component_count} mixture components")

def main():
    parser = argparse.ArgumentParser(description="Populate mixtures and mixture_components tables with scientifically accurate data.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 Mixture Population")
    
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
        logger.warning("No project found. Mixtures will be created without project association.")
    
    # Get molecule map
    molecule_map = get_molecule_map(supabase)
    if not molecule_map:
        logger.error("No molecules found. Please run populate_molecules.py first.")
        return
    
    # Populate mixtures
    mixture_map = populate_mixtures(supabase, molecule_map, user_profile_id, project_id, args.dry_run)
    
    # Populate mixture components
    populate_mixture_components(supabase, mixture_map, molecule_map, args.dry_run)
    
    logger.info("Mixture population complete")

if __name__ == "__main__":
    main()