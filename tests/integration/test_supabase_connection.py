#!/usr/bin/env python3
"""
CryoProtect v2 - Verify Supabase Connection

This script verifies the connection to the Supabase project and displays
basic information about the database schema.

Usage:
    python verify_supabase_connection.py
"""

import os
import sys
import json
import logging
from dotenv import load_dotenv

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("verify_supabase_connection.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Supabase project ID
SUPABASE_PROJECT_ID = "tsdlmynydfuypiugmkev"

def connect_to_supabase():
    """Connect to Supabase and authenticate."""
    try:
        from supabase import create_client, Client
        
        # Get Supabase URL and key from environment variables
        supabase_url = os.getenv("SUPABASE_URL")
        supabase_key = os.getenv("SUPABASE_KEY")
        
        if not supabase_url or not supabase_key:
            logger.error("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
            sys.exit(1)
        
        # Connect to Supabase
        supabase = create_client(supabase_url, supabase_key)
        logger.info("Successfully connected to Supabase")
        return supabase
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        sys.exit(1)

def check_table_exists(supabase, table_name):
    """Check if a table exists in the database."""
    try:
        # Try to select a single row from the table
        response = supabase.table(table_name).select("*").limit(1).execute()
        # If we get here without an error, the table exists
        logger.info(f"Table '{table_name}' exists")
        return True, len(response.data) if hasattr(response, 'data') else 0
    except Exception as e:
        if "relation" in str(e) and "does not exist" in str(e):
            logger.warning(f"Table '{table_name}' does not exist")
            return False, 0
        # If it's some other error, log it and assume the table doesn't exist
        logger.error(f"Error checking if table {table_name} exists: {str(e)}")
        return False, 0

def get_table_schema(supabase, table_name):
    """Get the schema of a table."""
    try:
        sql = f"""
        SELECT column_name, data_type, is_nullable
        FROM information_schema.columns
        WHERE table_schema = 'public' AND table_name = '{table_name}'
        ORDER BY ordinal_position;
        """
        
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        
        if hasattr(response, 'data'):
            return response.data
        else:
            logger.warning(f"No schema information found for table '{table_name}'")
            return []
    except Exception as e:
        logger.error(f"Error getting schema for table '{table_name}': {str(e)}")
        return []

def get_foreign_keys(supabase, table_name):
    """Get the foreign keys for a table."""
    try:
        sql = f"""
        SELECT
            tc.constraint_name,
            kcu.column_name,
            ccu.table_name AS foreign_table_name,
            ccu.column_name AS foreign_column_name
        FROM
            information_schema.table_constraints AS tc
            JOIN information_schema.key_column_usage AS kcu
              ON tc.constraint_name = kcu.constraint_name
              AND tc.table_schema = kcu.table_schema
            JOIN information_schema.constraint_column_usage AS ccu
              ON ccu.constraint_name = tc.constraint_name
              AND ccu.table_schema = tc.table_schema
        WHERE tc.constraint_type = 'FOREIGN KEY' AND tc.table_name = '{table_name}';
        """
        
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        
        if hasattr(response, 'data'):
            return response.data
        else:
            logger.warning(f"No foreign key information found for table '{table_name}'")
            return []
    except Exception as e:
        logger.error(f"Error getting foreign keys for table '{table_name}': {str(e)}")
        return []

def main():
    """Main function to verify Supabase connection and display schema information."""
    print("\n" + "=" * 80)
    print(f"CryoProtect v2 - Verify Supabase Connection (Project ID: {SUPABASE_PROJECT_ID})")
    print("=" * 80)
    
    try:
        # Connect to Supabase
        print("\nConnecting to Supabase...")
        supabase = connect_to_supabase()
        print("Successfully connected to Supabase")
        
        # Check required tables
        print("\nChecking required tables...")
        required_tables = [
            "molecules", 
            "mixtures", 
            "mixture_components", 
            "predictions", 
            "experiments", 
            "property_types", 
            "calculation_methods",
            "molecular_properties"
        ]
        
        table_results = {}
        for table in required_tables:
            exists, count = check_table_exists(supabase, table)
            table_results[table] = {"exists": exists, "count": count}
        
        print("\nTable Check Results:")
        for table, result in table_results.items():
            status = "Exists" if result["exists"] else "Missing"
            count = result["count"] if result["exists"] else 0
            print(f"  - {table}: {status} ({count} records)")
        
        # Get schema information for existing tables
        print("\nSchema Information:")
        for table, result in table_results.items():
            if result["exists"]:
                print(f"\n{table} Schema:")
                schema = get_table_schema(supabase, table)
                for column in schema:
                    nullable = "NULL" if column["is_nullable"] == "YES" else "NOT NULL"
                    print(f"  - {column['column_name']}: {column['data_type']} {nullable}")
                
                # Get foreign keys
                foreign_keys = get_foreign_keys(supabase, table)
                if foreign_keys:
                    print(f"\n{table} Foreign Keys:")
                    for fk in foreign_keys:
                        print(f"  - {fk['column_name']} -> {fk['foreign_table_name']}.{fk['foreign_column_name']}")
        
        print("\n" + "=" * 60)
        print("Supabase Connection Verification Complete")
        print("=" * 60)
        
        if all(result["exists"] for result in table_results.values()):
            print("\nAll required tables exist. You can proceed with running fix_relationships.py.")
        else:
            print("\nSome required tables are missing. Please create them before running fix_relationships.py.")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error verifying Supabase connection: {str(e)}")
        print(f"\nError verifying Supabase connection: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())