#!/usr/bin/env python3
"""
CryoProtect v2 - Verify Performance Indexes

This script verifies that the performance indexes have been successfully created
in the CryoProtect Supabase database.

Usage:
    python verify_performance_indexes.py
"""

import os
import sys
import json
import logging
from dotenv import load_dotenv
from supabase import create_client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("verify_performance_indexes.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def connect_to_supabase():
    """Connect to Supabase using service role key."""
    supabase_url = os.getenv("SUPABASE_URL")
    supabase_key = os.getenv("SUPABASE_KEY")
    
    if not supabase_url or not supabase_key:
        logger.error("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
        return None
    
    try:
        supabase = create_client(supabase_url, supabase_key)
        logger.info("Connected to Supabase")
        return supabase
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        return None

def execute_sql(supabase, sql, description):
    """Execute SQL using Supabase."""
    try:
        logger.info(f"Executing SQL: {description}")
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False, None
        
        return True, response.data
    except Exception as e:
        logger.error(f"Error executing SQL ({description}): {str(e)}")
        return False, None

def verify_indexes(supabase):
    """Verify that the performance indexes exist."""
    # List of indexes to verify
    expected_indexes = [
        {
            "name": "idx_predictions_mixture_property",
            "table": "predictions",
            "columns": ["mixture_id", "property_type_id"],
            "using": None
        },
        {
            "name": "idx_molecule_name_trgm",
            "table": "molecules",
            "columns": ["name"],
            "using": "gin"
        },
        {
            "name": "idx_mixture_name_trgm",
            "table": "mixtures",
            "columns": ["name"],
            "using": "gin"
        }
    ]
    
    # Get all indexes
    sql = """
    SELECT
        i.relname AS index_name,
        t.relname AS table_name,
        array_agg(a.attname ORDER BY k.indnatts) AS column_names,
        am.amname AS index_method
    FROM
        pg_index k
        JOIN pg_class i ON i.oid = k.indexrelid
        JOIN pg_class t ON t.oid = k.indrelid
        JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(k.indkey)
        JOIN pg_am am ON am.oid = i.relam
    WHERE
        t.relkind = 'r'
        AND t.relnamespace = 'public'::regnamespace
        AND i.relnamespace = 'public'::regnamespace
    GROUP BY
        i.relname,
        t.relname,
        am.amname
    ORDER BY
        t.relname,
        i.relname;
    """
    
    success, indexes = execute_sql(supabase, sql, "Get all indexes")
    if not success:
        return False
    
    # Check if each expected index exists
    missing_indexes = []
    for expected in expected_indexes:
        found = False
        for idx in indexes:
            if idx["index_name"] == expected["name"]:
                # Check if it's on the right table
                if idx["table_name"] != expected["table"]:
                    logger.warning(f"Index {expected['name']} is on table {idx['table_name']} instead of {expected['table']}")
                    continue
                
                # Check if it has the right columns
                columns_match = True
                for col in expected["columns"]:
                    if col not in idx["column_names"]:
                        logger.warning(f"Index {expected['name']} is missing column {col}")
                        columns_match = False
                        break
                
                # Check if it uses the right method
                if expected["using"] and idx["index_method"] != expected["using"]:
                    logger.warning(f"Index {expected['name']} uses method {idx['index_method']} instead of {expected['using']}")
                    continue
                
                if columns_match:
                    found = True
                    logger.info(f"Index {expected['name']} exists on table {expected['table']}")
                    break
        
        if not found:
            missing_indexes.append(expected["name"])
            logger.error(f"Index {expected['name']} does not exist on table {expected['table']}")
    
    # Check if pg_trgm extension is installed
    sql = """
    SELECT EXISTS (
        SELECT 1 FROM pg_extension WHERE extname = 'pg_trgm'
    );
    """
    
    success, result = execute_sql(supabase, sql, "Check pg_trgm extension")
    if not success:
        return False
    
    pg_trgm_installed = result[0]["exists"]
    if not pg_trgm_installed:
        logger.error("pg_trgm extension is not installed")
        return False
    
    logger.info("pg_trgm extension is installed")
    
    # Return True if all indexes exist
    if missing_indexes:
        logger.error(f"Missing indexes: {', '.join(missing_indexes)}")
        return False
    
    logger.info("All performance indexes exist")
    return True

def main():
    """Main function to verify performance indexes."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Verify Performance Indexes")
    print("=" * 80 + "\n")
    
    # Connect to Supabase
    supabase = connect_to_supabase()
    if not supabase:
        return 1
    
    # Verify indexes
    success = verify_indexes(supabase)
    
    # Print summary
    print("\n" + "=" * 60)
    print("Performance Indexes Verification Summary")
    print("=" * 60)
    
    if success:
        print("\nStatus: PASS")
        print("All performance indexes are correctly installed.")
    else:
        print("\nStatus: FAIL")
        print("One or more performance indexes are missing.")
        print("Please run the apply_performance_indexes script.")
    
    print("\nExpected Indexes:")
    print("  - idx_predictions_mixture_property on predictions(mixture_id, property_type_id)")
    print("  - idx_molecule_name_trgm on molecules using gin(name gin_trgm_ops)")
    print("  - idx_mixture_name_trgm on mixtures using gin(name gin_trgm_ops)")
    
    print("\nFor detailed information, check the log file.")
    print("=" * 60 + "\n")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())