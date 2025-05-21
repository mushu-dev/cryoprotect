#!/usr/bin/env python3
"""
CryoProtect v2 - Create Mixture Components Table

This script creates the missing 'mixture_components' table to match the updated
plural table naming convention.
"""

import sys
import logging
from fix_database_modular import connect_to_supabase, check_table_exists

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("create_mixture_components_table.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def create_mixture_components_table(supabase):
    """Create the mixture_components table if it doesn't exist."""
    if check_table_exists(supabase, "mixture_components"):
        logger.info("Table 'mixture_components' already exists.")
        return True
    
    logger.info("Creating 'mixture_components' table...")
    
    # Execute SQL to create the mixture_components table
    sql = """
    CREATE TABLE IF NOT EXISTS public.mixture_components (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        mixture_id UUID NOT NULL REFERENCES public.mixtures(id) ON DELETE CASCADE,
        molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
        amount NUMERIC,
        amount_unit VARCHAR(50),
        role VARCHAR(100),
        created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
    );
    
    CREATE INDEX IF NOT EXISTS idx_mixture_components_mixture_id ON public.mixture_components(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_mixture_components_molecule_id ON public.mixture_components(molecule_id);
    
    CREATE TRIGGER set_timestamp_mixture_components
    BEFORE UPDATE ON public.mixture_components
    FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    """
    
    try:
        # Execute the SQL using a raw query
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        logger.info("Successfully created 'mixture_components' table.")
        return True
    except Exception as e:
        logger.error(f"Error creating 'mixture_components' table: {str(e)}")
        return False

def main():
    """Create the mixture_components table."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Create Mixture Components Table")
    print("=" * 80)
    
    try:
        # Connect to Supabase
        print("\nConnecting to Supabase...")
        supabase = connect_to_supabase()
        
        # Create mixture_components table
        print("\nCreating 'mixture_components' table...")
        if create_mixture_components_table(supabase):
            print("Successfully created 'mixture_components' table.")
        else:
            print("Failed to create 'mixture_components' table. Check the logs for details.")
            return 1
        
        print("\n" + "=" * 60)
        print("Database Table Creation Complete")
        print("=" * 60)
        print("\nSuccessfully created the following table:")
        print("1. mixture_components")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error creating database table: {str(e)}")
        print(f"\nError creating database table: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())