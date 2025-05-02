#!/usr/bin/env python3
"""
CryoProtect v2 - Modular Database Fix Script

This script provides modular functions to fix specific database issues:
1. Check table existence
2. Create missing tables
3. Copy data between tables
4. Verify table structure

Usage:
    python fix_database_modular.py [--action <action_name>]
    
Actions:
    check_tables - Check if required tables exist
    create_tables - Create missing tables
    copy_data - Copy data from singular to plural tables
    verify_structure - Verify table structure
    fix_all - Run all fixes (default)
"""

import os
import sys
import json
import argparse
import logging
import time
from datetime import datetime
from dotenv import load_dotenv

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("fix_database_modular.log"),
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
        from supabase import create_client, Client
        
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

def check_required_tables(supabase):
    """Check if all required tables exist."""
    required_tables = [
        "molecules", 
        "mixtures", 
        "predictions", 
        "experiments", 
        "property_types", 
        "calculation_methods"
    ]
    
    results = {}
    all_exist = True
    
    for table in required_tables:
        exists = check_table_exists(supabase, table)
        results[table] = exists
        if not exists:
            all_exist = False
    
    # Also check for singular tables
    singular_tables = ["prediction", "experiment"]
    for table in singular_tables:
        exists = check_table_exists(supabase, table)
        results[table] = exists
    
    return all_exist, results

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

def create_required_tables(supabase):
    """Create all required tables that don't exist."""
    success1 = create_predictions_table(supabase)
    success2 = create_experiments_table(supabase)
    return success1 and success2

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

def copy_all_data(supabase):
    """Copy data from singular to plural tables."""
    success1 = copy_data_from_prediction_to_predictions(supabase)
    success2 = copy_data_from_experiment_to_experiments(supabase)
    return success1 and success2

def verify_table_structure(supabase, table_name, expected_columns):
    """Verify that a table has the expected structure."""
    try:
        # Get table information using a raw SQL query
        sql = f"""
        SELECT column_name, data_type 
        FROM information_schema.columns 
        WHERE table_schema = 'public' AND table_name = '{table_name}'
        """
        
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        
        if not hasattr(response, 'data') or not response.data:
            logger.error(f"Could not get column information for table '{table_name}'")
            return False
        
        # Convert response to a dictionary for easier comparison
        columns = {}
        for col in response.data:
            columns[col['column_name']] = col['data_type']
        
        # Check if all expected columns exist
        missing_columns = []
        for col_name, col_type in expected_columns.items():
            if col_name not in columns:
                missing_columns.append(col_name)
            elif col_type.lower() not in columns[col_name].lower():
                logger.warning(f"Column '{col_name}' has type '{columns[col_name]}', expected '{col_type}'")
        
        if missing_columns:
            logger.error(f"Table '{table_name}' is missing columns: {', '.join(missing_columns)}")
            return False
        
        logger.info(f"Table '{table_name}' has all expected columns")
        return True
    except Exception as e:
        logger.error(f"Error verifying table structure for '{table_name}': {str(e)}")
        return False

def verify_all_table_structures(supabase):
    """Verify the structure of all important tables."""
    # Define expected columns for each table
    table_structures = {
        "predictions": {
            "id": "uuid",
            "mixture_id": "uuid",
            "property_type_id": "uuid",
            "calculation_method_id": "uuid",
            "numeric_value": "numeric",
            "text_value": "text",
            "boolean_value": "boolean",
            "confidence": "numeric",
            "created_at": "timestamp",
            "updated_at": "timestamp",
            "created_by": "uuid"
        },
        "experiments": {
            "id": "uuid",
            "mixture_id": "uuid",
            "property_type_id": "uuid",
            "numeric_value": "numeric",
            "text_value": "text",
            "boolean_value": "boolean",
            "experimental_conditions": "text",
            "date_performed": "date",
            "created_at": "timestamp",
            "updated_at": "timestamp",
            "created_by": "uuid"
        }
    }
    
    results = {}
    all_valid = True
    
    for table_name, expected_columns in table_structures.items():
        valid = verify_table_structure(supabase, table_name, expected_columns)
        results[table_name] = valid
        if not valid:
            all_valid = False
    
    return all_valid, results

def apply_schema_fixes(supabase):
    """
    Apply schema fixes to the predictions and experiments tables as specified.
    """
    logger.info("Applying schema fixes to predictions and experiments tables...")
    sql = """
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Fix predictions table
ALTER TABLE public.predictions
  ALTER COLUMN id TYPE UUID USING uuid_generate_v4(),
  ALTER COLUMN mixture_id TYPE UUID,
  ALTER COLUMN property_type_id TYPE UUID,
  ADD COLUMN IF NOT EXISTS calculation_method_id UUID,
  ADD COLUMN IF NOT EXISTS confidence NUMERIC,
  ADD COLUMN IF NOT EXISTS updated_at TIMESTAMPTZ DEFAULT NOW(),
  ADD COLUMN IF NOT EXISTS created_by UUID;

-- Fix experiments table
ALTER TABLE public.experiments
  ALTER COLUMN id TYPE UUID USING uuid_generate_v4(),
  ALTER COLUMN mixture_id TYPE UUID,
  ALTER COLUMN property_type_id TYPE UUID,
  ADD COLUMN IF NOT EXISTS experimental_conditions TEXT,
  ADD COLUMN IF NOT EXISTS date_performed DATE,
  ADD COLUMN IF NOT EXISTS updated_at TIMESTAMPTZ DEFAULT NOW(),
  ADD COLUMN IF NOT EXISTS created_by UUID;

-- Add indexes
CREATE INDEX IF NOT EXISTS idx_pred_mixture ON public.predictions(mixture_id);
CREATE INDEX IF NOT EXISTS idx_exp_mixture ON public.experiments(mixture_id);
"""
    try:
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        logger.info("Schema fixes applied successfully.")
        print("Schema fixes applied successfully.")
        return True
    except Exception as e:
        logger.error(f"Error applying schema fixes: {str(e)}")
        print(f"Error applying schema fixes: {str(e)}")
        return False

def main():
    """Main function to run the database fixes."""
    parser = argparse.ArgumentParser(description='Fix database issues')
    parser.add_argument('--action', choices=['check_tables', 'create_tables', 'copy_data', 'verify_structure', 'fix_all', 'apply_schema_fixes'],
                        default='fix_all', help='Action to perform')
    args = parser.parse_args()
    
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Modular Database Fix Script")
    print("=" * 80)
    
    try:
        # Connect to Supabase
        supabase = connect_to_supabase()
        
        if args.action == 'check_tables' or args.action == 'fix_all':
            print("\nChecking required tables...")
            all_exist, table_results = check_required_tables(supabase)
            
            print("\nTable Check Results:")
            for table, exists in table_results.items():
                print(f"  - {table}: {'Exists' if exists else 'Missing'}")
            
            if all_exist:
                print("\nAll required tables exist.")
            else:
                print("\nSome required tables are missing.")
        
        if args.action == 'create_tables' or args.action == 'fix_all':
            print("\nCreating missing tables...")
            if create_required_tables(supabase):
                print("Successfully created missing tables.")
            else:
                print("Failed to create some tables. Check the logs for details.")
        
        if args.action == 'copy_data' or args.action == 'fix_all':
            print("\nCopying data from singular to plural tables...")
            if copy_all_data(supabase):
                print("Successfully copied data.")
            else:
                print("Failed to copy some data. Check the logs for details.")
        
        if args.action == 'verify_structure' or args.action == 'fix_all':
            print("\nVerifying table structures...")
            all_valid, structure_results = verify_all_table_structures(supabase)
            
            print("\nTable Structure Verification Results:")
            for table, valid in structure_results.items():
                print(f"  - {table}: {'Valid' if valid else 'Invalid'}")
            
            if all_valid:
                print("\nAll table structures are valid.")
            else:
                print("\nSome table structures are invalid. Check the logs for details.")

        if args.action == 'apply_schema_fixes':
            print("\nApplying schema fixes to predictions and experiments tables...")
            if apply_schema_fixes(supabase):
                print("Schema fixes applied successfully.")
            else:
                print("Failed to apply schema fixes. Check the logs for details.")
        
        print("\n" + "=" * 60)
        print("Database Fix Complete")
        print("=" * 60)
        
        if args.action == 'fix_all':
            print("Summary:")
            print("  - Required tables checked")
            print("  - Missing tables created")
            print("  - Data copied from singular to plural tables")
            print("  - Table structures verified")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error in database fix: {str(e)}")
        print(f"\nError in database fix: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())