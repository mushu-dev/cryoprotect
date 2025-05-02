#!/usr/bin/env python3
"""
CryoProtect v2 - Database Table Fix Script

This script fixes the mismatch between table names in the database schema and the actual tables.
It creates the "predictions" and "experiments" tables if they don't exist and copies data from
the singular-named tables ("prediction", "experiment") to the plural-named tables.

Usage:
    python fix_database_tables.py
"""

import os
import logging
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("fix_database_tables.log"),
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

def connect_to_supabase():
    """Connect to Supabase and authenticate."""
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

def check_table_exists(supabase, table_name):
    """Check if a table exists in the database."""
    try:
        # Try to select a single row from the table
        response = supabase.table(table_name).select("*").limit(1).execute()
        # If we get here without an error, the table exists
        return True
    except Exception as e:
        if "relation" in str(e) and "does not exist" in str(e):
            return False
        # If it's some other error, log it and assume the table doesn't exist
        logger.error(f"Error checking if table {table_name} exists: {str(e)}")
        return False

def create_predictions_table(supabase):
    """Create the predictions table if it doesn't exist."""
    if check_table_exists(supabase, "predictions"):
        logger.info("Table 'predictions' already exists.")
        return True
    
    logger.info("Creating 'predictions' table...")
    
    # Execute SQL to create the predictions table
    sql = """
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
    
    CREATE TRIGGER set_timestamp_predictions
    BEFORE UPDATE ON public.predictions
    FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    """
    
    try:
        # Execute the SQL using a raw query
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        logger.info("Successfully created 'predictions' table.")
        return True
    except Exception as e:
        logger.error(f"Error creating 'predictions' table: {str(e)}")
        return False

def create_experiments_table(supabase):
    """Create the experiments table if it doesn't exist."""
    if check_table_exists(supabase, "experiments"):
        logger.info("Table 'experiments' already exists.")
        return True
    
    logger.info("Creating 'experiments' table...")
    
    # Execute SQL to create the experiments table
    sql = """
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
    
    CREATE TRIGGER set_timestamp_experiments
    BEFORE UPDATE ON public.experiments
    FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    """
    
    try:
        # Execute the SQL using a raw query
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        logger.info("Successfully created 'experiments' table.")
        return True
    except Exception as e:
        logger.error(f"Error creating 'experiments' table: {str(e)}")
        return False

def copy_data_from_prediction_to_predictions(supabase):
    """Copy data from the prediction table to the predictions table."""
    if not check_table_exists(supabase, "prediction"):
        logger.warning("Table 'prediction' does not exist. No data to copy.")
        return True
    
    logger.info("Copying data from 'prediction' to 'predictions'...")
    
    # First, get the data from the prediction table
    try:
        response = supabase.table("prediction").select("*").execute()
        if not hasattr(response, 'data') or not response.data:
            logger.info("No data found in 'prediction' table.")
            return True
        
        prediction_data = response.data
        logger.info(f"Found {len(prediction_data)} records in 'prediction' table.")
        
        # Now insert the data into the predictions table
        for record in prediction_data:
            # Map fields from prediction to predictions
            # Adjust this mapping based on the actual field names in both tables
            new_record = {
                "id": record.get("id"),
                "mixture_id": record.get("mixture_id"),
                "property_type_id": record.get("property_type_id"),
                "calculation_method_id": record.get("method_id"),
                "numeric_value": record.get("predicted_value"),
                "confidence": record.get("confidence"),
                "created_at": record.get("created_at"),
                "updated_at": record.get("updated_at"),
                "created_by": record.get("created_by")
            }
            
            # Insert the record into the predictions table
            try:
                insert_response = supabase.table("predictions").upsert(new_record).execute()
                if hasattr(insert_response, 'error') and insert_response.error:
                    logger.error(f"Error inserting record into 'predictions': {insert_response.error}")
            except Exception as e:
                logger.error(f"Error inserting record into 'predictions': {str(e)}")
        
        logger.info(f"Successfully copied data from 'prediction' to 'predictions'.")
        return True
    except Exception as e:
        logger.error(f"Error copying data from 'prediction' to 'predictions': {str(e)}")
        return False

def copy_data_from_experiment_to_experiments(supabase):
    """Copy data from the experiment table to the experiments table."""
    if not check_table_exists(supabase, "experiment"):
        logger.warning("Table 'experiment' does not exist. No data to copy.")
        return True
    
    logger.info("Copying data from 'experiment' to 'experiments'...")
    
    # First, get the data from the experiment table
    try:
        response = supabase.table("experiment").select("*").execute()
        if not hasattr(response, 'data') or not response.data:
            logger.info("No data found in 'experiment' table.")
            return True
        
        experiment_data = response.data
        logger.info(f"Found {len(experiment_data)} records in 'experiment' table.")
        
        # Now insert the data into the experiments table
        for record in experiment_data:
            # Map fields from experiment to experiments
            # Adjust this mapping based on the actual field names in both tables
            new_record = {
                "id": record.get("id"),
                "mixture_id": record.get("mixture_id"),
                "property_type_id": record.get("property_type_id"),
                "numeric_value": record.get("value"),
                "experimental_conditions": record.get("experimental_conditions"),
                "date_performed": record.get("date_performed"),
                "created_at": record.get("created_at"),
                "updated_at": record.get("updated_at"),
                "created_by": record.get("created_by")
            }
            
            # Insert the record into the experiments table
            try:
                insert_response = supabase.table("experiments").upsert(new_record).execute()
                if hasattr(insert_response, 'error') and insert_response.error:
                    logger.error(f"Error inserting record into 'experiments': {insert_response.error}")
            except Exception as e:
                logger.error(f"Error inserting record into 'experiments': {str(e)}")
        
        logger.info(f"Successfully copied data from 'experiment' to 'experiments'.")
        return True
    except Exception as e:
        logger.error(f"Error copying data from 'experiment' to 'experiments': {str(e)}")
        return False

def main():
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Database Table Fix Script")
    print("=" * 80)
    
    # Connect to Supabase
    supabase = connect_to_supabase()
    
    # Create the predictions table if it doesn't exist
    create_predictions_table(supabase)
    
    # Create the experiments table if it doesn't exist
    create_experiments_table(supabase)
    
    # Copy data from prediction to predictions
    copy_data_from_prediction_to_predictions(supabase)
    
    # Copy data from experiment to experiments
    copy_data_from_experiment_to_experiments(supabase)
    
    print("\n" + "=" * 60)
    print("Database Table Fix Complete")
    print("=" * 60)
    print("The 'predictions' and 'experiments' tables have been created and populated with data.")
    print("You may now need to update your code to use the correct table names.")
    print("=" * 60)

if __name__ == "__main__":
    main()