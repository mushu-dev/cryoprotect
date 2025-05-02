#!/usr/bin/env python3
"""
CryoProtect v2 - Fix Missing Database Tables

This script creates the missing 'predictions' and 'experiments' tables
using direct SQL statements instead of relying on the exec_sql function.
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
        logging.FileHandler("fix_missing_tables.log"),
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

def create_predictions_table(supabase):
    """Create the predictions table if it doesn't exist using direct SQL."""
    if check_table_exists(supabase, "predictions"):
        logger.info("Table 'predictions' already exists.")
        return True
    
    logger.info("Creating 'predictions' table...")
    
    # SQL for creating the predictions table
    sql_statements = [
        """
        CREATE TABLE IF NOT EXISTS public.predictions (
            id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
            mixture_id UUID NOT NULL REFERENCES public.mixtures(id) ON DELETE CASCADE,
            property_type_id UUID NOT NULL REFERENCES public.property_types(id),
            calculation_method_id UUID NOT NULL REFERENCES public.calculation_methods(id),
            numeric_value NUMERIC,
            text_value TEXT,
            boolean_value BOOLEAN,
            confidence NUMERIC,
            created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            created_by UUID REFERENCES auth.users(id),
            UNIQUE(mixture_id, property_type_id, calculation_method_id)
        );
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_predictions_mixture_id ON public.predictions(mixture_id);
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_predictions_property_type_id ON public.predictions(property_type_id);
        """,
        """
        DO $$
        BEGIN
            IF NOT EXISTS (
                SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_predictions'
            ) THEN
                CREATE TRIGGER set_timestamp_predictions
                BEFORE UPDATE ON public.predictions
                FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
            END IF;
        END
        $$;
        """
    ]
    
    try:
        # Execute each SQL statement directly using the Postgres connection
        for sql in sql_statements:
            # Use the REST API to execute SQL directly
            response = supabase.postgrest.rpc('postgres_execute', {'query': sql}).execute()
            
            # Check for errors
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error executing SQL: {response.error}")
                return False
        
        logger.info("Successfully created 'predictions' table.")
        return True
    except Exception as e:
        logger.error(f"Error creating 'predictions' table: {str(e)}")
        
        # Alternative approach: Try using the Supabase dashboard SQL editor
        logger.info("Failed to create table using API. Please try using the Supabase dashboard SQL editor with the following SQL:")
        for sql in sql_statements:
            logger.info(f"\n{sql}")
        
        return False

def create_experiments_table(supabase):
    """Create the experiments table if it doesn't exist using direct SQL."""
    if check_table_exists(supabase, "experiments"):
        logger.info("Table 'experiments' already exists.")
        return True
    
    logger.info("Creating 'experiments' table...")
    
    # SQL for creating the experiments table
    sql_statements = [
        """
        CREATE TABLE IF NOT EXISTS public.experiments (
            id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
            mixture_id UUID NOT NULL REFERENCES public.mixtures(id) ON DELETE CASCADE,
            property_type_id UUID NOT NULL REFERENCES public.property_types(id),
            numeric_value NUMERIC,
            text_value TEXT,
            boolean_value BOOLEAN,
            experimental_conditions TEXT,
            date_performed DATE,
            created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            created_by UUID REFERENCES auth.users(id)
        );
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_experiments_mixture_id ON public.experiments(mixture_id);
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_experiments_property_type_id ON public.experiments(property_type_id);
        """,
        """
        DO $$
        BEGIN
            IF NOT EXISTS (
                SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_experiments'
            ) THEN
                CREATE TRIGGER set_timestamp_experiments
                BEFORE UPDATE ON public.experiments
                FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
            END IF;
        END
        $$;
        """
    ]
    
    try:
        # Execute each SQL statement directly using the Postgres connection
        for sql in sql_statements:
            # Use the REST API to execute SQL directly
            response = supabase.postgrest.rpc('postgres_execute', {'query': sql}).execute()
            
            # Check for errors
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error executing SQL: {response.error}")
                return False
        
        logger.info("Successfully created 'experiments' table.")
        return True
    except Exception as e:
        logger.error(f"Error creating 'experiments' table: {str(e)}")
        
        # Alternative approach: Try using the Supabase dashboard SQL editor
        logger.info("Failed to create table using API. Please try using the Supabase dashboard SQL editor with the following SQL:")
        for sql in sql_statements:
            logger.info(f"\n{sql}")
        
        return False

def create_tables_using_migration(supabase):
    """
    Alternative approach: Create tables using Supabase migrations.
    This is a fallback if direct SQL execution fails.
    """
    logger.info("Attempting to create tables using Supabase migration...")
    
    # SQL for creating both tables
    migration_sql = """
    -- Create predictions table
    CREATE TABLE IF NOT EXISTS public.predictions (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        mixture_id UUID NOT NULL REFERENCES public.mixtures(id) ON DELETE CASCADE,
        property_type_id UUID NOT NULL REFERENCES public.property_types(id),
        calculation_method_id UUID NOT NULL REFERENCES public.calculation_methods(id),
        numeric_value NUMERIC,
        text_value TEXT,
        boolean_value BOOLEAN,
        confidence NUMERIC,
        created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        created_by UUID REFERENCES auth.users(id),
        UNIQUE(mixture_id, property_type_id, calculation_method_id)
    );

    CREATE INDEX IF NOT EXISTS idx_predictions_mixture_id ON public.predictions(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_property_type_id ON public.predictions(property_type_id);

    DO $$
    BEGIN
        IF NOT EXISTS (
            SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_predictions'
        ) THEN
            CREATE TRIGGER set_timestamp_predictions
            BEFORE UPDATE ON public.predictions
            FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
        END IF;
    END
    $$;

    -- Create experiments table
    CREATE TABLE IF NOT EXISTS public.experiments (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        mixture_id UUID NOT NULL REFERENCES public.mixtures(id) ON DELETE CASCADE,
        property_type_id UUID NOT NULL REFERENCES public.property_types(id),
        numeric_value NUMERIC,
        text_value TEXT,
        boolean_value BOOLEAN,
        experimental_conditions TEXT,
        date_performed DATE,
        created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        created_by UUID REFERENCES auth.users(id)
    );

    CREATE INDEX IF NOT EXISTS idx_experiments_mixture_id ON public.experiments(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_experiments_property_type_id ON public.experiments(property_type_id);

    DO $$
    BEGIN
        IF NOT EXISTS (
            SELECT 1 FROM pg_trigger WHERE tgname = 'set_timestamp_experiments'
        ) THEN
            CREATE TRIGGER set_timestamp_experiments
            BEFORE UPDATE ON public.experiments
            FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
        END IF;
    END
    $$;
    """
    
    try:
        # Try to use the migrations API if available
        response = supabase.functions.invoke("apply-migration", {
            "name": "create_missing_tables",
            "sql": migration_sql
        })
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error applying migration: {response.error}")
            return False
        
        logger.info("Successfully created tables using migration.")
        return True
    except Exception as e:
        logger.error(f"Error creating tables using migration: {str(e)}")
        
        # Provide instructions for manual creation
        logger.info("Failed to create tables using migration. Please use the Supabase dashboard SQL editor with the following SQL:")
        logger.info(f"\n{migration_sql}")
        
        return False

def main():
    """Create the missing database tables."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Fix Missing Database Tables")
    print("=" * 80)
    
    try:
        # Connect to Supabase
        print("\nConnecting to Supabase...")
        supabase = connect_to_supabase()
        
        # Try direct SQL approach first
        print("\nAttempting to create tables using direct SQL...")
        predictions_success = create_predictions_table(supabase)
        experiments_success = create_experiments_table(supabase)
        
        # If direct SQL fails, try migration approach
        if not predictions_success or not experiments_success:
            print("\nDirect SQL approach failed. Trying migration approach...")
            migration_success = create_tables_using_migration(supabase)
            
            if not migration_success:
                print("\nAll approaches failed. Please check the logs for SQL statements to run manually.")
                return 1
        
        # Verify tables were created
        print("\nVerifying tables were created...")
        predictions_exists = check_table_exists(supabase, "predictions")
        experiments_exists = check_table_exists(supabase, "experiments")
        
        if predictions_exists and experiments_exists:
            print("\n" + "=" * 60)
            print("Database Tables Creation Complete")
            print("=" * 60)
            print("\nSuccessfully created the following tables:")
            print("1. predictions")
            print("2. experiments")
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