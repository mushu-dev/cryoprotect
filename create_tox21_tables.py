#!/usr/bin/env python3
"""
Create the essential Tox21 tables for testing.

This script creates the minimum necessary tables for Tox21 toxicity data integration.
"""

import os
import sys
import logging
from supabase import create_client
from config import Config

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def create_tables(supabase):
    """Create the essential Tox21 tables."""
    try:
        # Create toxicity_data_source table
        logger.info("Creating toxicity_data_source table...")
        response = supabase.table("toxicity_data_source").insert({
            "name": "Tox21",
            "description": "Toxicology in the 21st Century federal collaboration",
            "url": "https://ntp.niehs.nih.gov/whatwestudy/tox21/index.html"
        }).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error creating toxicity_data_source: {response.error}")
            return False
        
        # Get the Tox21 source ID
        tox21_source_id = response.data[0]['id']
        logger.info(f"Created Tox21 source with ID: {tox21_source_id}")
        
        # Create calculation_method for Tox21 Score
        logger.info("Creating Tox21 Score calculation method...")
        response = supabase.table("calculation_method").insert({
            "name": "Tox21 Score",
            "description": "Toxicity score based on Tox21 assay data",
            "method_type": "computational",
            "reference": "Tox21 program"
        }).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error creating calculation method: {response.error}")
            return False
        
        logger.info("Tables created successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error creating tables: {str(e)}")
        return False

def main():
    """Main function to create Tox21 tables."""
    # Initialize Supabase client
    config = Config()
    supabase_url = os.environ.get("SUPABASE_URL") or config.SUPABASE_URL
    supabase_key = os.environ.get("SUPABASE_KEY") or config.SUPABASE_KEY
    
    if not supabase_url or not supabase_key:
        logger.error("Supabase URL or key not found")
        return False
    
    supabase = create_client(supabase_url, supabase_key)
    
    # Create tables
    success = create_tables(supabase)
    
    if success:
        logger.info("Tox21 tables created successfully")
        return True
    else:
        logger.error("Failed to create Tox21 tables")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)