#!/usr/bin/env python3
"""
CryoProtect v2 - Create Missing Database Tables

This script creates the missing 'predictions' and 'experiments' tables
by importing and using functions from fix_database_modular.py.
"""

import sys
import logging
from fix_database_modular import connect_to_supabase, create_predictions_table, create_experiments_table

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("create_missing_tables.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def main():
    """Create the missing database tables."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Create Missing Database Tables")
    print("=" * 80)
    
    try:
        # Connect to Supabase
        print("\nConnecting to Supabase...")
        supabase = connect_to_supabase()
        
        # Create predictions table
        print("\nCreating 'predictions' table...")
        if create_predictions_table(supabase):
            print("Successfully created 'predictions' table.")
        else:
            print("Failed to create 'predictions' table. Check the logs for details.")
            return 1
        
        # Create experiments table
        print("\nCreating 'experiments' table...")
        if create_experiments_table(supabase):
            print("Successfully created 'experiments' table.")
        else:
            print("Failed to create 'experiments' table. Check the logs for details.")
            return 1
        
        print("\n" + "=" * 60)
        print("Database Tables Creation Complete")
        print("=" * 60)
        print("\nSuccessfully created the following tables:")
        print("1. predictions")
        print("2. experiments")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error creating database tables: {str(e)}")
        print(f"\nError creating database tables: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())