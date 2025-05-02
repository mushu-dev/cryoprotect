#!/usr/bin/env python3
"""
verify_database_integrity.py

This script verifies the integrity of the CryoProtect v2 database by checking:
1. All tables have data
2. Foreign key relationships are valid
3. Data is scientifically consistent

Usage:
    python verify_database_integrity.py
"""

import os
import json
import logging
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# Initialize Supabase client
supabase: Client = create_client(SUPABASE_URL, SUPABASE_KEY)

def authenticate():
    """Attempt to authenticate with Supabase."""
    try:
        email = os.getenv("SUPABASE_USER")
        password = os.getenv("SUPABASE_PASSWORD")
        if email and password:
            supabase.auth.sign_in_with_password({"email": email, "password": password})
            logger.info(f"Successfully authenticated as {email}")
        else:
            logger.warning("No authentication credentials provided. Continuing anonymously.")
    except Exception as e:
        logger.warning(f"Authentication error: {e}")
        logger.warning("Continuing without authentication. Some operations may fail.")

def count_table_rows(table_name):
    """Count rows in a table."""
    try:
        response = supabase.table(table_name).select("id").execute()
        count = len(response.data)
        logger.info(f"Table {table_name}: {count} rows")
        return count
    except Exception as e:
        logger.error(f"Error counting rows in {table_name}: {e}")
        return 0

def verify_foreign_key(table_name, column_name, referenced_table, referenced_column="id"):
    """Verify that foreign keys reference existing records."""
    try:
        # Get all values of the foreign key column
        response = supabase.table(table_name).select(column_name).execute()
        if not response.data:
            logger.info(f"No data in {table_name}.{column_name} to verify")
            return True
        
        # Extract foreign key values (excluding nulls)
        fk_values = [row[column_name] for row in response.data if row.get(column_name)]
        if not fk_values:
            logger.info(f"No non-null values in {table_name}.{column_name} to verify")
            return True
        
        # Check if all foreign keys exist in the referenced table
        for fk_value in fk_values:
            ref_response = supabase.table(referenced_table).select("id").eq(referenced_column, fk_value).execute()
            if not ref_response.data:
                logger.error(f"Foreign key violation: {table_name}.{column_name} = {fk_value} does not exist in {referenced_table}.{referenced_column}")
                return False
        
        logger.info(f"All foreign keys in {table_name}.{column_name} are valid")
        return True
    except Exception as e:
        logger.error(f"Error verifying foreign key {table_name}.{column_name}: {e}")
        return False

def verify_scientific_consistency():
    """Verify that the data is scientifically consistent."""
    issues = []
    
    # Check that all mixtures have at least one component
    try:
        mixtures = supabase.table("mixture").select("id,name").execute().data
        for mixture in mixtures:
            components = supabase.table("mixture_component").select("*").eq("mixture_id", mixture["id"]).execute().data
            if not components:
                issues.append(f"Mixture '{mixture['name']}' has no components")
    except Exception as e:
        logger.error(f"Error checking mixture components: {e}")
    
    # Check that all experiments have a valid temperature
    try:
        experiments = supabase.table("experiment").select("id,name,temperature").execute().data
        for experiment in experiments:
            temp = experiment.get("temperature")
            if temp is None or temp > 100:  # Most cryopreservation happens below 100°C
                issues.append(f"Experiment '{experiment['name']}' has invalid temperature: {temp}")
    except Exception as e:
        logger.error(f"Error checking experiment temperatures: {e}")
    
    # Check that all predictions have a confidence value
    try:
        predictions = supabase.table("prediction").select("id,property_type,confidence").execute().data
        for prediction in predictions:
            confidence = prediction.get("confidence")
            if confidence is None or confidence < 0 or confidence > 1:
                issues.append(f"Prediction for '{prediction['property_type']}' has invalid confidence: {confidence}")
    except Exception as e:
        logger.error(f"Error checking prediction confidence: {e}")
    
    return issues

def main():
    """Main function to verify database integrity."""
    logger.info("Starting CryoProtect v2 Database Integrity Verification")
    
    # Authenticate with Supabase
    authenticate()
    
    # Check that all tables have data
    tables = [
        "molecule", "molecular_property", "mixture", "mixture_component",
        "experiment", "experiment_property", "prediction", "calculation_method",
        "project", "team", "user_profile"
    ]
    
    empty_tables = []
    for table in tables:
        count = count_table_rows(table)
        if count == 0:
            empty_tables.append(table)
    
    # Verify foreign key relationships
    fk_issues = []
    if not verify_foreign_key("mixture_component", "mixture_id", "mixture"):
        fk_issues.append("mixture_component.mixture_id → mixture.id")
    if not verify_foreign_key("mixture_component", "molecule_id", "molecule"):
        fk_issues.append("mixture_component.molecule_id → molecule.id")
    if not verify_foreign_key("experiment", "mixture_id", "mixture"):
        fk_issues.append("experiment.mixture_id → mixture.id")
    if not verify_foreign_key("experiment", "project_id", "project"):
        fk_issues.append("experiment.project_id → project.id")
    if not verify_foreign_key("experiment_property", "experiment_id", "experiment"):
        fk_issues.append("experiment_property.experiment_id → experiment.id")
    if not verify_foreign_key("molecular_property", "molecule_id", "molecule"):
        fk_issues.append("molecular_property.molecule_id → molecule.id")
    if not verify_foreign_key("prediction", "molecule_id", "molecule"):
        fk_issues.append("prediction.molecule_id → molecule.id")
    if not verify_foreign_key("prediction", "method_id", "calculation_method"):
        fk_issues.append("prediction.method_id → calculation_method.id")
    
    # Verify scientific consistency
    scientific_issues = verify_scientific_consistency()
    
    # Generate report
    report = {
        "timestamp": datetime.now().isoformat(),
        "empty_tables": empty_tables,
        "foreign_key_issues": fk_issues,
        "scientific_issues": scientific_issues,
        "status": "PASS" if not (empty_tables or fk_issues or scientific_issues) else "FAIL"
    }
    
    # Print summary
    print("\n" + "=" * 60)
    print("CryoProtect v2 Database Integrity Verification Summary")
    print("=" * 60)
    print(f"Status: {report['status']}")
    print(f"Empty Tables: {len(empty_tables)}")
    print(f"Foreign Key Issues: {len(fk_issues)}")
    print(f"Scientific Consistency Issues: {len(scientific_issues)}")
    print("=" * 60)
    
    # Save report to file
    with open("database_integrity_report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info("Verification complete. Report saved to database_integrity_report.json")
    
    return report["status"] == "PASS"

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)