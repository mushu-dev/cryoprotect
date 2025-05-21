#!/usr/bin/env python3
"""
test_supabase_connection.py - Test connection to Supabase and MCP tools

This script verifies that:
1. The Supabase MCP tools are installed and working
2. The connection to the Supabase project is successful
3. The database can be queried

Usage:
    python test_supabase_connection.py
"""

import sys
import json
import logging
from supabase_mcp_tools import (
    get_project_details,
    list_tables,
    execute_sql_on_supabase
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Project ID
PROJECT_ID = "tsdlmynydfuypiugmkev"

def test_project_connection():
    """Test connection to the Supabase project"""
    try:
        logger.info("Testing connection to Supabase project...")
        project = get_project_details(PROJECT_ID)
        
        logger.info(f"✅ Successfully connected to project: {project['name']} (Region: {project['region']})")
        logger.info(f"   Project ID: {project['id']}")
        logger.info(f"   Organization ID: {project['organization_id']}")
        logger.info(f"   Status: {project['status']}")
        logger.info(f"   Database Version: {project['database']['version']}")
        
        return True
    except Exception as e:
        logger.error(f"❌ Failed to connect to Supabase project: {str(e)}")
        return False

def test_list_tables():
    """Test listing tables in the database"""
    try:
        logger.info("Testing ability to list tables...")
        tables = list_tables(PROJECT_ID, ["public"])
        
        # Count tables by schema
        schema_counts = {}
        for table in tables:
            schema = table.get("schema", "unknown")
            if schema not in schema_counts:
                schema_counts[schema] = 0
            schema_counts[schema] += 1
        
        logger.info(f"✅ Successfully listed tables. Found {len(tables)} tables total.")
        for schema, count in schema_counts.items():
            logger.info(f"   Schema '{schema}': {count} tables")
        
        # List public tables
        public_tables = [table["name"] for table in tables if table.get("schema") == "public"]
        logger.info(f"Public tables: {', '.join(public_tables)}")
        
        return True
    except Exception as e:
        logger.error(f"❌ Failed to list tables: {str(e)}")
        return False

def test_execute_query():
    """Test executing a simple query"""
    try:
        logger.info("Testing ability to execute SQL queries...")
        
        # Simple query to get database version
        result = execute_sql_on_supabase(
            PROJECT_ID,
            "SELECT version() as postgres_version;"
        )
        
        if result and len(result) > 0:
            postgres_version = result[0].get("postgres_version", "Unknown")
            logger.info(f"✅ Successfully executed query. PostgreSQL version: {postgres_version}")
            
            # Test a more complex query to count tables
            result = execute_sql_on_supabase(
                PROJECT_ID,
                """
                SELECT 
                    table_schema, 
                    COUNT(*) as table_count
                FROM 
                    information_schema.tables
                WHERE 
                    table_schema NOT IN ('pg_catalog', 'information_schema')
                GROUP BY 
                    table_schema
                ORDER BY 
                    table_schema;
                """
            )
            
            logger.info("Table counts by schema:")
            for row in result:
                logger.info(f"   {row['table_schema']}: {row['table_count']} tables")
            
            return True
        else:
            logger.error("❌ Query returned no results")
            return False
    except Exception as e:
        logger.error(f"❌ Failed to execute query: {str(e)}")
        return False

def main():
    """Main function to run all tests"""
    logger.info("Starting Supabase connection tests...")
    
    tests = [
        ("Project Connection", test_project_connection),
        ("List Tables", test_list_tables),
        ("Execute Query", test_execute_query)
    ]
    
    results = []
    
    for test_name, test_func in tests:
        logger.info(f"\n=== Testing {test_name} ===")
        success = test_func()
        results.append((test_name, success))
    
    # Print summary
    logger.info("\n=== Test Summary ===")
    all_passed = True
    for test_name, success in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        logger.info(f"{status} - {test_name}")
        if not success:
            all_passed = False
    
    if all_passed:
        logger.info("\n🎉 All tests passed! Your Supabase connection is working correctly.")
        logger.info("You can now run the schema standardization script.")
        return 0
    else:
        logger.error("\n❌ Some tests failed. Please fix the issues before running the schema standardization script.")
        return 1

if __name__ == "__main__":
    sys.exit(main())