#!/usr/bin/env python3
"""
CryoProtect v2 - Database Integrity Remediation Script

This script remediates critical database integrity issues by:
1. Repopulating the molecule table using populate_molecules.py
2. Removing orphaned records in mixture_component, molecular_property, and prediction tables
3. Verifying database integrity after remediation

Usage:
    python remediate_database_integrity.py [--dry-run]

Environment variables required (from .env):
    SUPABASE_URL, SUPABASE_KEY, SUPABASE_USER, SUPABASE_PASSWORD
"""

import os
import sys
import json
import logging
import argparse
import subprocess
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("database_remediation.log"),
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

def run_populate_molecules(dry_run=False):
    """Run the populate_molecules.py script to repopulate the molecule table."""
    logger.info("Running populate_molecules.py to repopulate the molecule table...")
    
    cmd = [sys.executable, "populate_molecules.py"]
    if dry_run:
        cmd.append("--dry-run")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info("populate_molecules.py completed successfully")
        logger.debug(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running populate_molecules.py: {e}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        return False

def identify_orphaned_records(supabase, table_name, fk_column, referenced_table, referenced_column="id"):
    """Identify orphaned records in a table where foreign keys don't exist in the referenced table."""
    logger.info(f"Identifying orphaned records in {table_name} where {fk_column} doesn't exist in {referenced_table}...")
    
    try:
        # Get all values of the foreign key column
        response = supabase.table(table_name).select(f"id,{fk_column}").execute()
        if not response.data:
            logger.info(f"No data in {table_name} to check for orphans")
            return []
        
        # Extract foreign key values (excluding nulls)
        records = [row for row in response.data if row.get(fk_column)]
        if not records:
            logger.info(f"No non-null values in {table_name}.{fk_column} to check for orphans")
            return []
        
        # Check each record to see if its foreign key exists in the referenced table
        orphaned_records = []
        for record in records:
            fk_value = record[fk_column]
            ref_response = supabase.table(referenced_table).select("id").eq(referenced_column, fk_value).execute()
            if not ref_response.data:
                orphaned_records.append(record)
                logger.debug(f"Orphaned record found: {table_name}.id = {record['id']} with {fk_column} = {fk_value}")
        
        logger.info(f"Found {len(orphaned_records)} orphaned records in {table_name}")
        return orphaned_records
    except Exception as e:
        logger.error(f"Error identifying orphaned records in {table_name}: {e}")
        return []

def remove_orphaned_records(supabase, table_name, orphaned_records, dry_run=False):
    """Remove orphaned records from a table."""
    if not orphaned_records:
        logger.info(f"No orphaned records to remove from {table_name}")
        return 0
    
    logger.info(f"Removing {len(orphaned_records)} orphaned records from {table_name}...")
    
    if dry_run:
        logger.info(f"DRY RUN: Would remove {len(orphaned_records)} orphaned records from {table_name}")
        return len(orphaned_records)
    
    try:
        # Extract IDs of orphaned records
        orphaned_ids = [record["id"] for record in orphaned_records]
        
        # Delete orphaned records in batches to avoid API limitations
        batch_size = 100
        deleted_count = 0
        
        for i in range(0, len(orphaned_ids), batch_size):
            batch_ids = orphaned_ids[i:i+batch_size]
            
            # Delete records in this batch
            response = supabase.table(table_name).delete().in_("id", batch_ids).execute()
            
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error deleting orphaned records from {table_name}: {response.error}")
            else:
                deleted_count += len(batch_ids)
                logger.info(f"Deleted batch of {len(batch_ids)} orphaned records from {table_name}")
        
        logger.info(f"Successfully removed {deleted_count} orphaned records from {table_name}")
        return deleted_count
    except Exception as e:
        logger.error(f"Error removing orphaned records from {table_name}: {e}")
        return 0

def run_verify_database_integrity():
    """Run the verify_database_integrity.py script to verify database integrity."""
    logger.info("Running verify_database_integrity.py to verify database integrity...")
    
    try:
        result = subprocess.run([sys.executable, "verify_database_integrity.py"], 
                               capture_output=True, text=True, check=True)
        logger.info("verify_database_integrity.py completed successfully")
        logger.debug(result.stdout)
        
        # Parse the database_integrity_report.json file
        try:
            with open("database_integrity_report.json", "r") as f:
                report = json.load(f)
            return report
        except Exception as e:
            logger.error(f"Error parsing database_integrity_report.json: {e}")
            return None
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running verify_database_integrity.py: {e}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Remediate database integrity issues.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 Database Integrity Remediation")
    
    # Step 1: Repopulate the molecule table
    if not run_populate_molecules(args.dry_run):
        logger.error("Failed to repopulate the molecule table. Aborting remediation.")
        return False
    
    # Step 2: Connect to Supabase
    supabase = connect_to_supabase()
    
    # Step 3: Identify and remove orphaned records in dependent tables
    remediation_summary = {
        "timestamp": datetime.now().isoformat(),
        "actions": [],
        "status": "UNKNOWN"
    }
    
    # Handle mixture_component table
    orphaned_mixture_components = identify_orphaned_records(supabase, "mixture_component", "molecule_id", "molecule")
    removed_mixture_components = remove_orphaned_records(supabase, "mixture_component", orphaned_mixture_components, args.dry_run)
    remediation_summary["actions"].append({
        "table": "mixture_component",
        "orphaned_records_found": len(orphaned_mixture_components),
        "orphaned_records_removed": removed_mixture_components
    })
    
    # Handle molecular_property table
    orphaned_molecular_properties = identify_orphaned_records(supabase, "molecular_property", "molecule_id", "molecule")
    removed_molecular_properties = remove_orphaned_records(supabase, "molecular_property", orphaned_molecular_properties, args.dry_run)
    remediation_summary["actions"].append({
        "table": "molecular_property",
        "orphaned_records_found": len(orphaned_molecular_properties),
        "orphaned_records_removed": removed_molecular_properties
    })
    
    # Handle prediction table
    orphaned_predictions = identify_orphaned_records(supabase, "prediction", "molecule_id", "molecule")
    removed_predictions = remove_orphaned_records(supabase, "prediction", orphaned_predictions, args.dry_run)
    remediation_summary["actions"].append({
        "table": "prediction",
        "orphaned_records_found": len(orphaned_predictions),
        "orphaned_records_removed": removed_predictions
    })
    
    # Step 4: Verify database integrity after remediation
    if not args.dry_run:
        integrity_report = run_verify_database_integrity()
        if integrity_report:
            remediation_summary["verification_result"] = integrity_report["status"]
            remediation_summary["status"] = "SUCCESS" if integrity_report["status"] == "PASS" else "FAILED"
        else:
            remediation_summary["verification_result"] = "ERROR"
            remediation_summary["status"] = "FAILED"
    else:
        remediation_summary["verification_result"] = "SKIPPED (DRY RUN)"
        remediation_summary["status"] = "DRY RUN"
    
    # Save remediation summary to file
    with open("database_remediation_summary.json", "w") as f:
        json.dump(remediation_summary, f, indent=2)
    
    # Print summary
    print("\n" + "=" * 60)
    print("CryoProtect v2 Database Integrity Remediation Summary")
    print("=" * 60)
    print(f"Status: {remediation_summary['status']}")
    for action in remediation_summary["actions"]:
        print(f"{action['table']}: Found {action['orphaned_records_found']} orphaned records, removed {action['orphaned_records_removed']}")
    if "verification_result" in remediation_summary:
        print(f"Verification Result: {remediation_summary['verification_result']}")
    print("=" * 60)
    
    logger.info("Remediation complete. Summary saved to database_remediation_summary.json")
    
    return remediation_summary["status"] == "SUCCESS"

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)