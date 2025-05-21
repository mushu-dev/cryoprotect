#!/usr/bin/env python3
"""
CryoProtect v2 - Create Tables and Populate

This script checks if the required tables exist with the correct plural names,
and then runs the population scripts to populate the database.
"""

import os
import sys
import subprocess
import logging
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("create_tables_and_populate.log"),
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
    logger.error("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
    sys.exit(1)

def connect_to_supabase():
    """Connect to Supabase and authenticate."""
    try:
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
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        raise

def check_table_exists(supabase, table_name):
    """Check if a table exists in the database."""
    try:
        # Try to select a single row from the table
        response = supabase.table(table_name).select("*").limit(1).execute()
        # If we get here without an error, the table exists
        logger.info(f"Table '{table_name}' exists")
        return True
    except Exception as e:
        if "relation" in str(e) and "does not exist" in str(e):
            logger.warning(f"Table '{table_name}' does not exist")
            return False
        # If it's some other error, log it and assume the table doesn't exist
        logger.error(f"Error checking if table {table_name} exists: {str(e)}")
        return False

def run_script(script_name):
    """Run a Python script and return the result."""
    logger.info(f"Running {script_name}...")
    try:
        result = subprocess.run(["python", script_name], capture_output=True, text=True)
        if result.returncode == 0:
            logger.info(f"Successfully ran {script_name}")
            logger.info(f"Output: {result.stdout}")
            return True
        else:
            logger.error(f"Error running {script_name}: {result.stderr}")
            return False
    except Exception as e:
        logger.error(f"Exception running {script_name}: {str(e)}")
        return False

def main():
    """Check if tables exist and run population scripts."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Create Tables and Populate")
    print("=" * 80)
    
    try:
        # Connect to Supabase
        print("\nConnecting to Supabase...")
        supabase = connect_to_supabase()
        
        # Check if tables exist
        print("\nChecking if tables exist...")
        tables = [
            "molecules",
            "molecular_properties",
            "mixtures",
            "mixture_components",
            "predictions",
            "experiments"
        ]
        
        missing_tables = []
        for table in tables:
            if not check_table_exists(supabase, table):
                missing_tables.append(table)
        
        if missing_tables:
            print(f"\nThe following tables are missing: {', '.join(missing_tables)}")
            print("Please create these tables before running the population scripts.")
            return 1
        
        # Run population scripts in order
        print("\nPopulating tables...")
        
        # Populate molecules
        print("\nPopulating molecules...")
        if not run_script("populate_molecules.py"):
            print("Failed to populate molecules. Check the logs for details.")
            return 1
        
        # Populate mixtures
        print("\nPopulating mixtures...")
        if not run_script("populate_mixtures.py"):
            print("Failed to populate mixtures. Check the logs for details.")
            return 1
        
        # Populate predictions
        print("\nPopulating predictions...")
        if not run_script("populate_predictions.py"):
            print("Failed to populate predictions. Check the logs for details.")
            return 1
        
        # Populate experiments
        print("\nPopulating experiments...")
        if not run_script("populate_experiments.py"):
            print("Failed to populate experiments. Check the logs for details.")
            return 1
        
        print("\n" + "=" * 60)
        print("Database Population Complete")
        print("=" * 60)
        print("\nSuccessfully populated the following tables:")
        print("1. molecules")
        print("2. molecular_properties")
        print("3. mixtures")
        print("4. mixture_components")
        print("5. predictions")
        print("6. experiments")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error checking tables and populating database: {str(e)}")
        print(f"\nError checking tables and populating database: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())