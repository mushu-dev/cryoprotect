#!/usr/bin/env python3
"""
Test script to diagnose and fix Supabase HTTP 403 error.
This script tests inserting a molecule using the service role key.
"""

import os
import json
import logging
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def main():
    # Get Supabase credentials
    supabase_url = os.getenv("SUPABASE_URL")
    supabase_key = os.getenv("SUPABASE_KEY")
    
    if not supabase_url or not supabase_key:
        logger.error("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
        return
    
    logger.info(f"Using Supabase URL: {supabase_url}")
    logger.info(f"Using Supabase Key: {supabase_key[:10]}...")
    
    # Initialize Supabase client
    try:
        logger.info("Initializing Supabase client...")
        supabase = create_client(supabase_url, supabase_key)
    except Exception as e:
        logger.error(f"Failed to initialize Supabase client: {e}")
        return
    
    # Test authentication
    try:
        logger.info("Testing authentication...")
        user = supabase.auth.get_user()
        if hasattr(user, 'user') and user.user:
            logger.info(f"Authenticated as user ID: {user.user.id}")
        else:
            logger.info("Not authenticated as a user")
        
        # Check if using service role
        session = supabase.auth.get_session()
        if hasattr(session, 'access_token'):
            logger.info("Session access token available")
        else:
            logger.info("No session access token")
    except Exception as e:
        logger.error(f"Error checking authentication: {e}")
    
    # Test inserting a molecule
    try:
        logger.info("Testing molecule insert...")
        
        # Create test molecule
        test_molecule = {
            "name": "Test Molecule (Python Script)",
            "smiles": "CC(=O)O",
            "inchi": "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)",
            "inchikey": "QTBSBXVTEAMEQO-UHFFFAOYSA-N",
            "formula": "C2H4O2",
            "molecular_weight": 60.05,
            "data_source": "Test Script",
            "version": 1,
            "modification_history": [{
                "timestamp": datetime.now().isoformat(),
                "action": "created",
                "source": "test_supabase_insert.py"
            }]
        }
        
        # Try insert with explicit service role header
        logger.info("Attempting insert with explicit service role header...")
        response = supabase.table("molecules").insert(test_molecule).execute()
        
        if hasattr(response, "error") and response.error:
            logger.error(f"Error inserting molecule: {response.error}")
        else:
            logger.info(f"Successfully inserted molecule: {response.data}")
            
    except Exception as e:
        logger.error(f"Exception during molecule insertion: {e}")
    
    # Test RLS policies
    try:
        logger.info("Testing RLS policies...")
        response = supabase.rpc('get_my_claim', {'claim': 'role'}).execute()
        if hasattr(response, "data"):
            logger.info(f"Role claim: {response.data}")
        else:
            logger.info("No role claim data returned")
    except Exception as e:
        logger.error(f"Error testing RLS policies: {e}")

if __name__ == "__main__":
    main()