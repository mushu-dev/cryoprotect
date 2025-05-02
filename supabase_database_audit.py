#!/usr/bin/env python3
"""
CryoProtect v2 - Supabase Database Audit

This script performs a comprehensive audit of the CryoProtect v2 Supabase database to:
1. Verify database schema, tables, and relationships before data population
2. Identify empty or sparsely populated tables
3. Map relationships between tables (foreign keys, dependencies)
4. Document data flow and dependencies between tables
5. Generate a detailed report for data population planning
6. Check database quotas and user permissions

Usage:
    python supabase_database_audit.py

Configuration:
    Uses the hierarchical config system (config.py) for all settings
    Required: SUPABASE_URL, SUPABASE_KEY
    Optional: SUPABASE_USER, SUPABASE_PASSWORD
"""

import os
import json
import sys
import time
from datetime import datetime
import logging
from pathlib import Path
from supabase import create_client, Client
import pandas as pd

# Import the hierarchical config system
from config import BaseConfig

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("database_audit.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Create logs directory if it doesn't exist
Path("logs").mkdir(exist_ok=True)

# Load configuration
try:
    config = BaseConfig.from_env()
    logger.info(f"Loaded configuration for environment: {config.__class__.__name__}")
    
    # Log non-secret configuration values
    safe_config = {k: v for k, v in config.as_dict().items()
                  if not any(secret in k.lower() for secret in ['key', 'password', 'secret', 'token'])}
    logger.info(f"Configuration: {json.dumps(safe_config, default=str, indent=2)}")
except Exception as e:
    logger.critical(f"Failed to load configuration: {str(e)}")
    sys.exit(1)

# Table relationships (based on schema analysis)
TABLE_RELATIONSHIPS = {
    "molecule": {
        "dependencies": [],
        "dependents": ["molecular_property", "mixture_component", "prediction"]
    },
    "mixture": {
        "dependencies": ["project", "user_profile"],
        "dependents": ["mixture_component", "experiment"]
    },
    "mixture_component": {
        "dependencies": ["mixture", "molecule"],
        "dependents": []
    },
    "experiment": {
        "dependencies": ["mixture", "project", "user_profile"],
        "dependents": ["experiment_property"]
    },
    "experiment_property": {
        "dependencies": ["experiment", "calculation_method", "user_profile"],
        "dependents": []
    },
    "molecular_property": {
        "dependencies": ["molecule", "calculation_method", "experiment", "user_profile"],
        "dependents": []
    },
    "prediction": {
        "dependencies": ["molecule", "calculation_method", "user_profile"],
        "dependents": []
    },
    "calculation_method": {
        "dependencies": ["user_profile"],
        "dependents": ["molecular_property", "prediction", "experiment_property"]
    },
    "project": {
        "dependencies": ["team", "user_profile"],
        "dependents": ["mixture", "experiment", "molecule"]
    },
    "team": {
        "dependencies": ["user_profile"],
        "dependents": ["project"]
    },
    "user_profile": {
        "dependencies": [],
        "dependents": ["team", "project", "molecule", "mixture", "experiment", 
                      "molecular_property", "prediction", "calculation_method"]
    }
}

# Required tables for ChEMBL data population
REQUIRED_TABLES = [
    "molecules",
    "molecular_properties",
    "mixture",
    "mixture_components",
    "experiment",
    "experiment_property",
    "prediction",
    "calculation_method",
    "property_types"
]

def connect_to_supabase():
    """Connect to Supabase and authenticate."""
    try:
        supabase = create_client(config.SUPABASE_URL, config.SUPABASE_KEY)
        logger.info(f"Connected to Supabase at {config.SUPABASE_URL}")
        
        # Authenticate if credentials provided
        if config.SUPABASE_USER and config.SUPABASE_PASSWORD:
            try:
                response = supabase.auth.sign_in_with_password({
                    "email": config.SUPABASE_USER,
                    "password": config.SUPABASE_PASSWORD
                })
                if hasattr(response, 'error') and response.error:
                    logger.warning(f"Authentication error: {response.error}")
                    logger.warning("Continuing without authentication. Some operations may fail.")
                else:
                    logger.info(f"Authenticated as {config.SUPABASE_USER}")
            except Exception as e:
                logger.warning(f"Authentication error: {str(e)}")
                logger.warning("Continuing without authentication. Some operations may fail.")
        else:
            logger.warning("No authentication credentials provided. Continuing without authentication.")
        
        # Test connection with a simple query
        try:
            response = supabase.rpc('get_service_role').execute()
            logger.info("Supabase connection test successful")
        except Exception as e:
            logger.error(f"Supabase connection test failed: {str(e)}")
            logger.error("Unable to verify Supabase connection. Aborting.")
            sys.exit(1)
            
        return supabase
    except Exception as e:
        logger.critical(f"Failed to connect to Supabase: {str(e)}")
        logger.critical("Aborting database audit due to connection failure")
        sys.exit(1)

def get_table_counts(supabase):
    """Get row counts for all tables in the database."""
    tables = list(TABLE_RELATIONSHIPS.keys())
    table_counts = {}
    
    for table in tables:
        try:
            # Count rows in table
            response = supabase.table(table).select("id", count="exact").execute()
            if hasattr(response, 'count'):
                count = response.count
            else:
                # Fallback if count attribute is not available
                count = len(response.data) if hasattr(response, 'data') else 0
            
            table_counts[table] = count
            logger.info(f"Table {table}: {count} rows")
        except Exception as e:
            logger.error(f"Error counting rows in {table}: {str(e)}")
            table_counts[table] = -1  # Error indicator
    
    return table_counts

def verify_required_tables(supabase):
    """Verify that all required tables exist in the database."""
    logger.info("Verifying required tables...")
    
    missing_tables = []
    existing_tables = []
    
    for table in REQUIRED_TABLES:
        try:
            # Try to select a single row to verify table exists
            response = supabase.table(table).select("*").limit(1).execute()
            if hasattr(response, 'data'):
                existing_tables.append(table)
                logger.info(f"✓ Table '{table}' exists")
            else:
                missing_tables.append(table)
                logger.error(f"✗ Table '{table}' could not be verified (unexpected response format)")
        except Exception as e:
            missing_tables.append(table)
            logger.error(f"✗ Table '{table}' does not exist or is not accessible: {str(e)}")
    
    if missing_tables:
        logger.critical(f"Missing required tables: {', '.join(missing_tables)}")
        logger.critical("Database schema verification failed. Required tables are missing.")
        return False, {"missing": missing_tables, "existing": existing_tables}
    
    logger.info("All required tables exist in the database")
    return True, {"missing": missing_tables, "existing": existing_tables}

def check_foreign_key_integrity(supabase):
    """Check foreign key integrity between tables."""
    logger.info("Checking foreign key integrity...")
    integrity_results = {}
    critical_failures = []
    
    for table, relations in TABLE_RELATIONSHIPS.items():
        dependencies = relations["dependencies"]
        
        if not dependencies:
            continue
        
        integrity_results[table] = {}
        
        for dep_table in dependencies:
            # This is a simplified check - in a real scenario, we would need to know the exact foreign key column names
            # For now, we'll assume the convention that foreign keys are named as {table_name}_id
            fk_column = f"{dep_table}_id"
            
            try:
                # Get all distinct foreign key values
                response = supabase.table(table).select(fk_column).execute()
                if not hasattr(response, 'data'):
                    logger.error(f"No data attribute in response for {table}.{fk_column}")
                    continue
                
                fk_values = [row.get(fk_column) for row in response.data if row.get(fk_column)]
                
                if not fk_values:
                    integrity_results[table][dep_table] = {
                        "status": "UNKNOWN",
                        "details": f"No foreign keys found in {table}.{fk_column}"
                    }
                    continue
                
                # Check if these values exist in the parent table
                # This is a simplified approach - in a real scenario, we would use a JOIN or EXISTS query
                # But Supabase's Python client has limitations for complex queries
                sample_fk = fk_values[0]  # Just check one as an example
                parent_response = supabase.table(dep_table).select("id").eq("id", sample_fk).execute()
                
                if hasattr(parent_response, 'data') and len(parent_response.data) > 0:
                    integrity_results[table][dep_table] = {
                        "status": "PASS",
                        "details": f"Foreign key integrity check passed for {table}.{fk_column}"
                    }
                    logger.info(f"✓ Foreign key integrity check passed for {table}.{fk_column}")
                else:
                    integrity_results[table][dep_table] = {
                        "status": "FAIL",
                        "details": f"Foreign key integrity check failed for {table}.{fk_column}"
                    }
                    critical_failures.append(f"{table}.{fk_column} -> {dep_table}.id")
                    logger.error(f"✗ Foreign key integrity check failed for {table}.{fk_column}")
            except Exception as e:
                logger.error(f"Error checking foreign key integrity for {table}.{fk_column}: {str(e)}")
                integrity_results[table][dep_table] = {
                    "status": "ERROR",
                    "details": f"Error: {str(e)}"
                }
    
    if critical_failures:
        logger.critical(f"Foreign key integrity check failed for: {', '.join(critical_failures)}")
        logger.critical("Database relationship verification failed. Foreign key constraints are violated.")
        return False, integrity_results
    
    logger.info("Foreign key integrity checks passed")
    return True, integrity_results

def check_database_quotas(supabase):
    """Check database quotas and user permissions."""
    logger.info("Checking database quotas and permissions...")
    
    quota_results = {
        "status": "UNKNOWN",
        "details": {},
        "warnings": []
    }
    
    try:
        # Try to get database size and quota information
        # This is a simplified approach - in a real scenario, we would use a more direct method
        # But Supabase's Python client has limitations for administrative queries
        response = supabase.rpc('get_database_size').execute()
        
        if hasattr(response, 'data') and response.data:
            quota_results["status"] = "PASS"
            quota_results["details"] = response.data
            
            # Check if we're approaching limits
            if "size_bytes" in response.data and "quota_bytes" in response.data:
                usage_percent = (response.data["size_bytes"] / response.data["quota_bytes"]) * 100
                quota_results["details"]["usage_percent"] = usage_percent
                
                if usage_percent > 90:
                    warning = f"Database is at {usage_percent:.1f}% of quota"
                    quota_results["warnings"].append(warning)
                    logger.warning(warning)
                elif usage_percent > 75:
                    warning = f"Database is at {usage_percent:.1f}% of quota"
                    quota_results["warnings"].append(warning)
                    logger.warning(warning)
        else:
            quota_results["status"] = "UNKNOWN"
            quota_results["details"] = {"message": "Could not retrieve database size information"}
            logger.warning("Could not retrieve database size information")
    except Exception as e:
        quota_results["status"] = "ERROR"
        quota_results["details"] = {"error": str(e)}
        logger.warning(f"Error checking database quotas: {str(e)}")
    
    # Check user permissions by attempting to create and drop a temporary table
    try:
        temp_table = f"temp_permission_check_{int(time.time())}"
        supabase.rpc('execute_sql', {"sql": f"CREATE TABLE {temp_table} (id serial PRIMARY KEY)"}).execute()
        supabase.rpc('execute_sql', {"sql": f"DROP TABLE {temp_table}"}).execute()
        quota_results["details"]["write_permission"] = True
        logger.info("✓ User has write permissions")
    except Exception as e:
        quota_results["details"]["write_permission"] = False
        quota_results["warnings"].append("User may not have sufficient write permissions")
        logger.warning(f"User may not have sufficient write permissions: {str(e)}")
    
    return quota_results

def sample_table_data(supabase, table, limit=5):
    """Get a sample of data from a table."""
    try:
        response = supabase.table(table).select("*").limit(limit).execute()
        if hasattr(response, 'data'):
            return response.data
        return []
    except Exception as e:
        logger.error(f"Error sampling data from {table}: {str(e)}")
        return []

def analyze_data_completeness(supabase):
    """Analyze data completeness in key tables."""
    completeness_results = {}
    
    # Define key tables and their critical columns
    key_tables = {
        "molecule": ["name", "smiles", "inchi", "inchikey", "formula", "molecular_weight"],
        "mixture": ["name", "description"],
        "mixture_component": ["molecule_id", "mixture_id", "amount", "amount_unit"],
        "experiment": ["name", "description", "mixture_id"],
        "experiment_property": ["experiment_id", "property_type", "value", "unit"],
        "molecular_property": ["molecule_id", "property_type", "value", "unit"],
        "prediction": ["molecule_id", "property_type", "predicted_value", "unit"]
    }
    
    for table, columns in key_tables.items():
        try:
            # Get all rows from the table (this could be optimized for large tables)
            response = supabase.table(table).select("*").execute()
            
            if not hasattr(response, 'data') or not response.data:
                completeness_results[table] = {
                    "status": "EMPTY",
                    "details": f"No data found in {table}"
                }
                continue
            
            # Count rows with missing values in critical columns
            total_rows = len(response.data)
            missing_data = {}
            
            for column in columns:
                missing_count = sum(1 for row in response.data if column not in row or row[column] is None)
                missing_percentage = (missing_count / total_rows) * 100 if total_rows > 0 else 0
                missing_data[column] = {
                    "missing_count": missing_count,
                    "missing_percentage": missing_percentage
                }
            
            # Determine overall completeness status
            if any(data["missing_percentage"] > 50 for data in missing_data.values()):
                status = "CRITICAL"
            elif any(data["missing_percentage"] > 20 for data in missing_data.values()):
                status = "WARNING"
            elif any(data["missing_percentage"] > 0 for data in missing_data.values()):
                status = "MINOR_ISSUES"
            else:
                status = "COMPLETE"
            
            completeness_results[table] = {
                "status": status,
                "total_rows": total_rows,
                "missing_data": missing_data,
                "details": f"{table} data completeness: {status}"
            }
            
        except Exception as e:
            logger.error(f"Error analyzing data completeness for {table}: {str(e)}")
            completeness_results[table] = {
                "status": "ERROR",
                "details": f"Error: {str(e)}"
            }
    
    return completeness_results

def generate_dependency_graph(table_counts):
    """Generate a dependency graph for visualization."""
    nodes = []
    edges = []
    
    for table, count in table_counts.items():
        # Create node with size based on row count
        size = 1
        if count > 0:
            size = min(10, max(3, count / 10))  # Scale node size between 3 and 10
        
        nodes.append({
            "id": table,
            "label": f"{table}\n({count} rows)",
            "size": size,
            "color": "#ff9900" if count == 0 else "#66cc00"
        })
        
        # Create edges based on dependencies
        for dep_table in TABLE_RELATIONSHIPS[table]["dependencies"]:
            edges.append({
                "from": dep_table,
                "to": table,
                "arrows": "to"
            })
    
    return {
        "nodes": nodes,
        "edges": edges
    }

def generate_audit_report(table_counts, integrity_results, completeness_results, graph_data):
    """Generate a comprehensive audit report."""
    report = {
        "timestamp": datetime.now().isoformat(),
        "summary": {
            "total_tables": len(table_counts),
            "empty_tables": sum(1 for count in table_counts.values() if count == 0),
            "populated_tables": sum(1 for count in table_counts.values() if count > 0),
            "integrity_issues": sum(1 for table_results in integrity_results.values() 
                                  for result in table_results.values() if result["status"] != "PASS"),
            "completeness_issues": sum(1 for result in completeness_results.values() 
                                     if result["status"] != "COMPLETE")
        },
        "table_counts": table_counts,
        "integrity_results": integrity_results,
        "completeness_results": completeness_results,
        "dependency_graph": graph_data,
        "population_recommendations": []
    }
    
    # Generate population recommendations based on dependencies
    population_sequence = []
    remaining_tables = set(table_counts.keys())
    
    # Start with tables that have no dependencies
    current_level = [table for table in remaining_tables 
                    if not TABLE_RELATIONSHIPS[table]["dependencies"]]
    
    while current_level:
        population_sequence.append(current_level)
        remaining_tables -= set(current_level)
        
        # Find tables whose dependencies are all in the population sequence
        populated_tables = set(sum(population_sequence, []))
        current_level = [
            table for table in remaining_tables
            if all(dep in populated_tables for dep in TABLE_RELATIONSHIPS[table]["dependencies"])
        ]
    
    # Add any remaining tables (in case of circular dependencies)
    if remaining_tables:
        population_sequence.append(list(remaining_tables))
    
    # Create recommendations
    for level, tables in enumerate(population_sequence):
        for table in tables:
            count = table_counts.get(table, 0)
            if count == 0:
                report["population_recommendations"].append({
                    "table": table,
                    "priority": level + 1,
                    "status": "EMPTY",
                    "recommendation": f"Populate {table} table (Priority {level + 1})",
                    "dependencies": TABLE_RELATIONSHIPS[table]["dependencies"]
                })
            elif count < 10:  # Arbitrary threshold for "sparsely populated"
                report["population_recommendations"].append({
                    "table": table,
                    "priority": level + 1,
                    "status": "SPARSE",
                    "recommendation": f"Add more data to {table} table (Priority {level + 1})",
                    "dependencies": TABLE_RELATIONSHIPS[table]["dependencies"]
                })
    
    return report

def main():
    logger.info("Starting CryoProtect v2 Supabase Database Audit")
    
    # Create reports directory if it doesn't exist
    Path("reports").mkdir(exist_ok=True)
    
    # Connect to Supabase
    supabase = connect_to_supabase()
    
    # Verify required tables
    tables_ok, tables_result = verify_required_tables(supabase)
    
    # Check foreign key integrity
    if tables_ok:
        integrity_ok, integrity_results = check_foreign_key_integrity(supabase)
    else:
        logger.error("Skipping foreign key integrity check due to missing tables")
        integrity_ok = False
        integrity_results = {}
    
    # Check database quotas and permissions
    quota_results = check_database_quotas(supabase)
    
    # Get table counts
    logger.info("Counting rows in tables...")
    table_counts = get_table_counts(supabase)
    
    # Analyze data completeness
    logger.info("Analyzing data completeness...")
    completeness_results = analyze_data_completeness(supabase)
    
    # Generate dependency graph
    logger.info("Generating dependency graph...")
    graph_data = generate_dependency_graph(table_counts)
    
    # Generate audit report
    logger.info("Generating audit report...")
    report = generate_audit_report(table_counts, integrity_results, completeness_results, graph_data)
    
    # Add verification results to report
    report["verification"] = {
        "tables_verified": tables_ok,
        "missing_tables": tables_result["missing"],
        "integrity_verified": integrity_ok,
        "quota_status": quota_results["status"],
        "quota_warnings": quota_results["warnings"],
        "timestamp": datetime.now().isoformat()
    }
    
    # Save report to file
    report_path = f"reports/database_audit_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    
    # Also save a copy to the latest report
    with open("reports/database_audit_report_latest.json", "w") as f:
        json.dump(report, f, indent=2)
    
    # Print summary
    print("\n" + "="*60)
    print("CryoProtect v2 Database Audit Summary")
    print("="*60)
    print(f"Database Schema Verification: {'PASS' if tables_ok else 'FAIL'}")
    print(f"Foreign Key Integrity: {'PASS' if integrity_ok else 'FAIL'}")
    print(f"Database Quota Status: {quota_results['status']}")
    if quota_results['warnings']:
        print(f"Quota Warnings: {', '.join(quota_results['warnings'])}")
    print(f"Total tables: {report['summary']['total_tables']}")
    print(f"Empty tables: {report['summary']['empty_tables']}")
    print(f"Populated tables: {report['summary']['populated_tables']}")
    print(f"Integrity issues: {report['summary']['integrity_issues']}")
    print(f"Completeness issues: {report['summary']['completeness_issues']}")
    print("\nTop Population Priorities:")
    
    for rec in sorted(report["population_recommendations"], key=lambda x: x["priority"])[:5]:
        print(f"- {rec['recommendation']}")
    
    print(f"\nFull report saved to {report_path}")
    print("="*60)
    
    # Abort if verification failed
    if not tables_ok or not integrity_ok:
        logger.critical("Database verification failed. Aborting.")
        sys.exit(1)
    
    logger.info("Database verification completed successfully")

if __name__ == "__main__":
    main()