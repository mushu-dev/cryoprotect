#!/usr/bin/env python3
"""
CryoProtect v2 - Create and Populate Tables

This script creates all necessary tables with the correct plural names and
then runs the population scripts in the correct order to populate the database.
"""

import os
import sys
import subprocess
import logging
from fix_database_modular import connect_to_supabase, check_table_exists

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("create_and_populate_tables.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def create_molecules_table(supabase):
    """Create the molecules table if it doesn't exist."""
    if check_table_exists(supabase, "molecules"):
        logger.info("Table 'molecules' already exists.")
        return True
    
    logger.info("Creating 'molecules' table...")
    
    # Execute SQL to create the molecules table
    sql = """
    CREATE TABLE IF NOT EXISTS public.molecules (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        name VARCHAR(255) NOT NULL,
        smiles VARCHAR(1000) UNIQUE,
        inchi TEXT,
        inchikey VARCHAR(27),
        formula VARCHAR(100),
        molecular_weight NUMERIC,
        created_by UUID,
        data_source VARCHAR(255),
        version INTEGER,
        modification_history JSONB,
        created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
    );
    
    CREATE INDEX IF NOT EXISTS idx_molecules_name ON public.molecules(name);
    CREATE INDEX IF NOT EXISTS idx_molecules_inchikey ON public.molecules(inchikey);
    
    CREATE TRIGGER set_timestamp_molecules
    BEFORE UPDATE ON public.molecules
    FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    """
    
    try:
        # Execute the SQL using a raw query
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        logger.info("Successfully created 'molecules' table.")
        return True
    except Exception as e:
        logger.error(f"Error creating 'molecules' table: {str(e)}")
        return False

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

def create_mixtures_table(supabase):
    """Create the mixtures table if it doesn't exist."""
    if check_table_exists(supabase, "mixtures"):
        logger.info("Table 'mixtures' already exists.")
        return True
    
    logger.info("Creating 'mixtures' table...")
    
    # Execute SQL to create the mixtures table
    sql = """
    CREATE TABLE IF NOT EXISTS public.mixtures (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        name VARCHAR(255) NOT NULL,
        description TEXT,
        project_id UUID,
        created_by UUID,
        data_source VARCHAR(255),
        version INTEGER,
        modification_history JSONB,
        created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
    );
    
    CREATE INDEX IF NOT EXISTS idx_mixtures_name ON public.mixtures(name);
    
    CREATE TRIGGER set_timestamp_mixtures
    BEFORE UPDATE ON public.mixtures
    FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    """
    
    try:
        # Execute the SQL using a raw query
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        logger.info("Successfully created 'mixtures' table.")
        return True
    except Exception as e:
        logger.error(f"Error creating 'mixtures' table: {str(e)}")
        return False

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

def create_trigger_function(supabase):
    """Create the trigger function for updating timestamps if it doesn't exist."""
    logger.info("Creating trigger function for updating timestamps...")
    
    sql = """
    CREATE OR REPLACE FUNCTION public.trigger_set_timestamp()
    RETURNS TRIGGER AS $$
    BEGIN
      NEW.updated_at = NOW();
      RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    """
    
    try:
        # Execute the SQL using a raw query
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        logger.info("Successfully created trigger function.")
        return True
    except Exception as e:
        logger.error(f"Error creating trigger function: {str(e)}")
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
    """Create all necessary tables and populate them."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Create and Populate Tables")
    print("=" * 80)
    
    try:
        # Connect to Supabase
        print("\nConnecting to Supabase...")
        supabase = connect_to_supabase()
        
        # Create trigger function
        print("\nCreating trigger function...")
        if not create_trigger_function(supabase):
            print("Failed to create trigger function. Check the logs for details.")
            return 1
        
        # Create molecules table
        print("\nCreating 'molecules' table...")
        if not create_molecules_table(supabase):
            print("Failed to create 'molecules' table. Check the logs for details.")
            return 1
        
        # Create molecular_properties table
        print("\nCreating 'molecular_properties' table...")
        if not create_molecular_properties_table(supabase):
            print("Failed to create 'molecular_properties' table. Check the logs for details.")
            return 1
        
        # Create mixtures table
        print("\nCreating 'mixtures' table...")
        if not create_mixtures_table(supabase):
            print("Failed to create 'mixtures' table. Check the logs for details.")
            return 1
        
        # Create mixture_components table
        print("\nCreating 'mixture_components' table...")
        if not create_mixture_components_table(supabase):
            print("Failed to create 'mixture_components' table. Check the logs for details.")
            return 1
        
        # Create predictions and experiments tables using existing script
        print("\nCreating 'predictions' and 'experiments' tables...")
        if not run_script("create_missing_tables.py"):
            print("Failed to create 'predictions' and 'experiments' tables. Check the logs for details.")
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
        print("Database Creation and Population Complete")
        print("=" * 60)
        print("\nSuccessfully created and populated the following tables:")
        print("1. molecules")
        print("2. molecular_properties")
        print("3. mixtures")
        print("4. mixture_components")
        print("5. predictions")
        print("6. experiments")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error creating and populating database: {str(e)}")
        print(f"\nError creating and populating database: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())