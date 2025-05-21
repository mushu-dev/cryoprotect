#!/usr/bin/env python3
"""
CryoProtect v2 - Create Missing Tables for Population

This script creates the missing tables needed for the population scripts:
- molecular_properties
- mixtures
- mixture_components

It attempts to use direct SQL statements and provides fallback options.
"""

import os
import sys
import logging
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("create_missing_tables_for_population.log"),
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
        # Try to execute the SQL directly
        try:
            response = supabase.postgrest.rpc('postgres_execute', {'query': sql}).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error executing SQL: {response.error}")
                return False
            logger.info("Successfully created trigger function.")
            return True
        except Exception as e:
            logger.warning(f"Error using postgres_execute: {str(e)}")
            
            # Try alternative approach
            try:
                response = supabase.functions.invoke("apply-migration", {
                    "name": "create_trigger_function",
                    "sql": sql
                })
                if hasattr(response, 'error') and response.error:
                    logger.error(f"Error applying migration: {response.error}")
                    return False
                logger.info("Successfully created trigger function using migration.")
                return True
            except Exception as e2:
                logger.warning(f"Error using functions.invoke: {str(e2)}")
                
                # Provide instructions for manual creation
                logger.info("Failed to create trigger function. Please use the Supabase dashboard SQL editor with the following SQL:")
                logger.info(f"\n{sql}")
                
                # Ask user if they want to continue
                print("\nFailed to create trigger function automatically.")
                print("Please execute the following SQL in the Supabase dashboard SQL editor:")
                print(f"\n{sql}")
                response = input("\nHave you executed the SQL manually? (y/n): ")
                return response.lower() == 'y'
    except Exception as e:
        logger.error(f"Error creating trigger function: {str(e)}")
        return False

def create_molecular_properties_table(supabase):
    """Create the molecular_properties table if it doesn't exist."""
    if check_table_exists(supabase, "molecular_properties"):
        logger.info("Table 'molecular_properties' already exists.")
        return True
    
    logger.info("Creating 'molecular_properties' table...")
    
    # SQL for creating the molecular_properties table
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
    
    DO $$
    BEGIN
        IF NOT EXISTS (
            SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_molecular_properties'
        ) THEN
            CREATE TRIGGER set_timestamp_molecular_properties
            BEFORE UPDATE ON public.molecular_properties
            FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
        END IF;
    END
    $$;
    """
    
    try:
        # Try to execute the SQL directly
        try:
            response = supabase.postgrest.rpc('postgres_execute', {'query': sql}).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error executing SQL: {response.error}")
                return False
            logger.info("Successfully created 'molecular_properties' table.")
            return True
        except Exception as e:
            logger.warning(f"Error using postgres_execute: {str(e)}")
            
            # Try alternative approach
            try:
                response = supabase.functions.invoke("apply-migration", {
                    "name": "create_molecular_properties_table",
                    "sql": sql
                })
                if hasattr(response, 'error') and response.error:
                    logger.error(f"Error applying migration: {response.error}")
                    return False
                logger.info("Successfully created 'molecular_properties' table using migration.")
                return True
            except Exception as e2:
                logger.warning(f"Error using functions.invoke: {str(e2)}")
                
                # Provide instructions for manual creation
                logger.info("Failed to create 'molecular_properties' table. Please use the Supabase dashboard SQL editor with the following SQL:")
                logger.info(f"\n{sql}")
                
                # Ask user if they want to continue
                print("\nFailed to create 'molecular_properties' table automatically.")
                print("Please execute the following SQL in the Supabase dashboard SQL editor:")
                print(f"\n{sql}")
                response = input("\nHave you executed the SQL manually? (y/n): ")
                return response.lower() == 'y'
    except Exception as e:
        logger.error(f"Error creating 'molecular_properties' table: {str(e)}")
        return False

def create_mixtures_table(supabase):
    """Create the mixtures table if it doesn't exist."""
    if check_table_exists(supabase, "mixtures"):
        logger.info("Table 'mixtures' already exists.")
        return True
    
    logger.info("Creating 'mixtures' table...")
    
    # SQL for creating the mixtures table
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
    
    DO $$
    BEGIN
        IF NOT EXISTS (
            SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_mixtures'
        ) THEN
            CREATE TRIGGER set_timestamp_mixtures
            BEFORE UPDATE ON public.mixtures
            FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
        END IF;
    END
    $$;
    """
    
    try:
        # Try to execute the SQL directly
        try:
            response = supabase.postgrest.rpc('postgres_execute', {'query': sql}).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error executing SQL: {response.error}")
                return False
            logger.info("Successfully created 'mixtures' table.")
            return True
        except Exception as e:
            logger.warning(f"Error using postgres_execute: {str(e)}")
            
            # Try alternative approach
            try:
                response = supabase.functions.invoke("apply-migration", {
                    "name": "create_mixtures_table",
                    "sql": sql
                })
                if hasattr(response, 'error') and response.error:
                    logger.error(f"Error applying migration: {response.error}")
                    return False
                logger.info("Successfully created 'mixtures' table using migration.")
                return True
            except Exception as e2:
                logger.warning(f"Error using functions.invoke: {str(e2)}")
                
                # Provide instructions for manual creation
                logger.info("Failed to create 'mixtures' table. Please use the Supabase dashboard SQL editor with the following SQL:")
                logger.info(f"\n{sql}")
                
                # Ask user if they want to continue
                print("\nFailed to create 'mixtures' table automatically.")
                print("Please execute the following SQL in the Supabase dashboard SQL editor:")
                print(f"\n{sql}")
                response = input("\nHave you executed the SQL manually? (y/n): ")
                return response.lower() == 'y'
    except Exception as e:
        logger.error(f"Error creating 'mixtures' table: {str(e)}")
        return False

def create_mixture_components_table(supabase):
    """Create the mixture_components table if it doesn't exist."""
    if check_table_exists(supabase, "mixture_components"):
        logger.info("Table 'mixture_components' already exists.")
        return True
    
    logger.info("Creating 'mixture_components' table...")
    
    # SQL for creating the mixture_components table
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
    
    DO $$
    BEGIN
        IF NOT EXISTS (
            SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_mixture_components'
        ) THEN
            CREATE TRIGGER set_timestamp_mixture_components
            BEFORE UPDATE ON public.mixture_components
            FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
        END IF;
    END
    $$;
    """
    
    try:
        # Try to execute the SQL directly
        try:
            response = supabase.postgrest.rpc('postgres_execute', {'query': sql}).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error executing SQL: {response.error}")
                return False
            logger.info("Successfully created 'mixture_components' table.")
            return True
        except Exception as e:
            logger.warning(f"Error using postgres_execute: {str(e)}")
            
            # Try alternative approach
            try:
                response = supabase.functions.invoke("apply-migration", {
                    "name": "create_mixture_components_table",
                    "sql": sql
                })
                if hasattr(response, 'error') and response.error:
                    logger.error(f"Error applying migration: {response.error}")
                    return False
                logger.info("Successfully created 'mixture_components' table using migration.")
                return True
            except Exception as e2:
                logger.warning(f"Error using functions.invoke: {str(e2)}")
                
                # Provide instructions for manual creation
                logger.info("Failed to create 'mixture_components' table. Please use the Supabase dashboard SQL editor with the following SQL:")
                logger.info(f"\n{sql}")
                
                # Ask user if they want to continue
                print("\nFailed to create 'mixture_components' table automatically.")
                print("Please execute the following SQL in the Supabase dashboard SQL editor:")
                print(f"\n{sql}")
                response = input("\nHave you executed the SQL manually? (y/n): ")
                return response.lower() == 'y'
    except Exception as e:
        logger.error(f"Error creating 'mixture_components' table: {str(e)}")
        return False

def main():
    """Create the missing database tables."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Create Missing Tables for Population")
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
        
        # Verify tables were created
        print("\nVerifying tables were created...")
        molecular_properties_exists = check_table_exists(supabase, "molecular_properties")
        mixtures_exists = check_table_exists(supabase, "mixtures")
        mixture_components_exists = check_table_exists(supabase, "mixture_components")
        
        if molecular_properties_exists and mixtures_exists and mixture_components_exists:
            print("\n" + "=" * 60)
            print("Database Tables Creation Complete")
            print("=" * 60)
            print("\nSuccessfully created the following tables:")
            print("1. molecular_properties")
            print("2. mixtures")
            print("3. mixture_components")
            
            print("\nYou can now run the population scripts to populate the database.")
            return 0
        else:
            print("\nFailed to verify table creation. Please check the logs for details.")
            return 1
    
    except Exception as e:
        logger.error(f"Error creating database tables: {str(e)}")
        print(f"\nError creating database tables: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())