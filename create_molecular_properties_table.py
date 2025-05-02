#!/usr/bin/env python3
"""
CryoProtect v2 - Create Molecular Properties Table

This script creates the missing 'molecular_properties' table to match the updated
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
        logging.FileHandler("create_molecular_properties_table.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def create_molecular_properties_table(supabase):
    """Create the molecular_properties table if it doesn't exist."""
    if check_table_exists(supabase, "molecular_properties"):
        logger.info("Table 'molecular_properties' already exists.")
        return True
    
    logger.info("Creating 'molecular_properties' table...")
    
    # Execute SQL to create the molecular_properties table
    sql = """
    CREATE TABLE IF NOT EXISTS public.molecular_properties (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
        property_type VARCHAR(255) NOT NULL,
        value NUMERIC,
        unit VARCHAR(50),
        created_by UUID,
        data_source VARCHAR(255),
        version INTEGER,
        modification_history JSONB,
        created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
    );
    
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON public.molecular_properties(molecule_id);
    
    CREATE TRIGGER set_timestamp_molecular_properties
    BEFORE UPDATE ON public.molecular_properties
    FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    """
    
    try:
        # Execute the SQL using a raw query
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        logger.info("Successfully created 'molecular_properties' table.")
        return True
    except Exception as e:
        logger.error(f"Error creating 'molecular_properties' table: {str(e)}")
        return False

def main():
    """Create the molecular_properties table."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Create Molecular Properties Table")
    print("=" * 80)
    
    try:
        # Connect to Supabase
        print("\nConnecting to Supabase...")
        supabase = connect_to_supabase()
        
        # Create molecular_properties table
        print("\nCreating 'molecular_properties' table...")
        if create_molecular_properties_table(supabase):
            print("Successfully created 'molecular_properties' table.")
        else:
            print("Failed to create 'molecular_properties' table. Check the logs for details.")
            return 1
        
        print("\n" + "=" * 60)
        print("Database Table Creation Complete")
        print("=" * 60)
        print("\nSuccessfully created the following table:")
        print("1. molecular_properties")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error creating database table: {str(e)}")
        print(f"\nError creating database table: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())