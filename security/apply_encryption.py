#!/usr/bin/env python3
"""
Apply Encryption to Existing Data

This script applies encryption to sensitive data in the database that was previously unencrypted.
It reads the encryption configuration to determine which fields need encryption and updates
the database accordingly.

Usage:
    python security/apply_encryption.py [--dry-run] [--table TABLE_NAME]

Options:
    --dry-run       Show what would be encrypted without making changes
    --table TABLE   Only encrypt data in the specified table
"""

import os
import sys
import json
import argparse
import logging
from datetime import datetime
from typing import Dict, List, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from security.encryption import get_encryption_service
from security.encryption_config import get_encrypted_fields, ENCRYPTED_FIELDS

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/encryption_migration.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("encryption_migration")

def get_supabase_client():
    """Get the Supabase client."""
    try:
        from supabase_adapter import get_supabase_client as get_client
        return get_client()
    except ImportError:
        logger.error("Could not import supabase_adapter. Make sure it's in your PYTHONPATH.")
        sys.exit(1)

def get_table_data(supabase, table_name: str) -> List[Dict[str, Any]]:
    """
    Get all data from a table.
    
    Args:
        supabase: Supabase client
        table_name: Name of the table
        
    Returns:
        List of records from the table
    """
    try:
        response = supabase.table(table_name).select("*").execute()
        if response.error:
            logger.error(f"Error fetching data from {table_name}: {response.error}")
            return []
        return response.data
    except Exception as e:
        logger.error(f"Exception fetching data from {table_name}: {str(e)}")
        return []

def update_record(supabase, table_name: str, record_id: str, updates: Dict[str, Any], id_field: str = "id") -> bool:
    """
    Update a record in the database.
    
    Args:
        supabase: Supabase client
        table_name: Name of the table
        record_id: ID of the record to update
        updates: Dictionary of field updates
        id_field: Name of the ID field (default: "id")
        
    Returns:
        True if successful, False otherwise
    """
    try:
        response = supabase.table(table_name).update(updates).eq(id_field, record_id).execute()
        if response.error:
            logger.error(f"Error updating record {record_id} in {table_name}: {response.error}")
            return False
        return True
    except Exception as e:
        logger.error(f"Exception updating record {record_id} in {table_name}: {str(e)}")
        return False

def encrypt_table_data(supabase, table_name: str, dry_run: bool = False) -> Dict[str, Any]:
    """
    Encrypt sensitive fields in a table.
    
    Args:
        supabase: Supabase client
        table_name: Name of the table
        dry_run: If True, don't actually update the database
        
    Returns:
        Dictionary with encryption statistics
    """
    encrypted_fields = get_encrypted_fields(table_name)
    if not encrypted_fields:
        logger.info(f"No fields to encrypt in table {table_name}")
        return {"table": table_name, "records_processed": 0, "records_updated": 0, "fields_encrypted": 0}
    
    logger.info(f"Encrypting fields {', '.join(encrypted_fields)} in table {table_name}")
    
    # Get encryption service
    encryption_service = get_encryption_service()
    
    # Get all records from the table
    records = get_table_data(supabase, table_name)
    if not records:
        logger.warning(f"No records found in table {table_name}")
        return {"table": table_name, "records_processed": 0, "records_updated": 0, "fields_encrypted": 0}
    
    records_processed = 0
    records_updated = 0
    fields_encrypted = 0
    
    # Process each record
    for record in records:
        records_processed += 1
        record_id = record.get("id")
        if not record_id:
            logger.warning(f"Record without ID in table {table_name}, skipping")
            continue
        
        # Check which fields need encryption
        updates = {}
        for field in encrypted_fields:
            if field in record and record[field] is not None:
                # Skip if already encrypted (check for base64 pattern)
                if isinstance(record[field], str) and record[field].startswith("key_"):
                    continue
                
                # Encrypt the field
                encrypted_value = encryption_service.encrypt(record[field])
                updates[field] = encrypted_value
                fields_encrypted += 1
        
        # Update the record if there are changes
        if updates:
            if dry_run:
                logger.info(f"Would update record {record_id} in {table_name} with encrypted values for {', '.join(updates.keys())}")
            else:
                success = update_record(supabase, table_name, record_id, updates)
                if success:
                    records_updated += 1
                    logger.info(f"Updated record {record_id} in {table_name} with encrypted values for {', '.join(updates.keys())}")
    
    logger.info(f"Processed {records_processed} records in {table_name}, updated {records_updated} records, encrypted {fields_encrypted} fields")
    return {
        "table": table_name,
        "records_processed": records_processed,
        "records_updated": records_updated,
        "fields_encrypted": fields_encrypted
    }

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Apply encryption to sensitive data in the database")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be encrypted without making changes")
    parser.add_argument("--table", help="Only encrypt data in the specified table")
    args = parser.parse_args()
    
    # Ensure logs directory exists
    os.makedirs("logs", exist_ok=True)
    
    # Get Supabase client
    supabase = get_supabase_client()
    
    # Determine which tables to process
    tables_to_process = [args.table] if args.table else list(ENCRYPTED_FIELDS.keys())
    
    # Process each table
    results = []
    for table_name in tables_to_process:
        logger.info(f"Processing table {table_name}")
        result = encrypt_table_data(supabase, table_name, args.dry_run)
        results.append(result)
    
    # Generate summary report
    total_records_processed = sum(r["records_processed"] for r in results)
    total_records_updated = sum(r["records_updated"] for r in results)
    total_fields_encrypted = sum(r["fields_encrypted"] for r in results)
    
    summary = {
        "timestamp": datetime.utcnow().isoformat(),
        "dry_run": args.dry_run,
        "tables_processed": len(results),
        "total_records_processed": total_records_processed,
        "total_records_updated": total_records_updated,
        "total_fields_encrypted": total_fields_encrypted,
        "details": results
    }
    
    # Save summary report
    report_path = f"reports/encryption_migration_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.json"
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w") as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Encryption migration complete. Processed {total_records_processed} records, updated {total_records_updated} records, encrypted {total_fields_encrypted} fields.")
    logger.info(f"Summary report saved to {report_path}")

if __name__ == "__main__":
    main()