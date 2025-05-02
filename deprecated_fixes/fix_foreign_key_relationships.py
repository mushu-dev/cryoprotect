#!/usr/bin/env python3
"""
CryoProtect v2 - Fix Foreign Key Relationships

This script fixes foreign key relationship issues in the CryoProtect database:
1. Identifies and removes duplicate foreign key constraints
2. Updates all foreign key references to use the new plural table names
3. Ensures all foreign keys have proper indexes and ON DELETE behavior
4. Adds missing foreign key constraints
5. Includes verification and rollback mechanisms

Usage:
    python fix_foreign_key_relationships.py [--dry-run] [--verify-only] [--rollback]
"""

import os
import sys
import json
import logging
import argparse
import time
from datetime import datetime
from dotenv import load_dotenv

# Set up logging
log_filename = f"fix_foreign_keys_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Track changes for potential rollback
applied_changes = []

# Table mapping (singular to plural)
TABLE_MAPPING = {
    "molecule": "molecules",
    "mixture": "mixtures",
    "prediction": "predictions",
    "experiment": "experiments",
    "experiment_property": "experiment_properties",
    "mixture_component": "mixture_components",
    "calculation_method": "calculation_methods",
    "property_type": "property_types",
    "project": "projects",
    "team": "teams"
}

# Define deletion behavior for each table relationship
# Format: (parent_table, child_table, child_column): "CASCADE" or "SET NULL"
DELETE_BEHAVIOR = {
    # Default is CASCADE for most relationships
    ("molecules", "molecular_properties", "molecule_id"): "CASCADE",
    ("molecules", "mixture_components", "molecule_id"): "CASCADE",
    ("molecules", "predictions", "molecule_id"): "CASCADE",
    ("mixtures", "mixture_components", "mixture_id"): "CASCADE",
    ("mixtures", "experiments", "mixture_id"): "CASCADE",
    ("mixtures", "predictions", "mixture_id"): "CASCADE",
    ("property_types", "molecular_properties", "property_type_id"): "CASCADE",
    ("property_types", "experiment_properties", "property_type_id"): "CASCADE",
    ("property_types", "predictions", "property_type_id"): "CASCADE",
    ("calculation_methods", "predictions", "calculation_method_id"): "CASCADE",
    ("experiments", "experiment_properties", "experiment_id"): "CASCADE",
    ("teams", "projects", "team_id"): "CASCADE",
    # Use SET NULL for optional relationships
    ("teams", "user_profile", "team_id"): "SET NULL",
    # Add more relationships as needed
}

def get_supabase_client():
    """Get a Supabase client with service role key."""
    try:
        from supabase import create_client, Client
        
        SUPABASE_URL = os.getenv("SUPABASE_URL")
        SUPABASE_KEY = os.getenv("SUPABASE_KEY")
        
        if not SUPABASE_URL or not SUPABASE_KEY:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
        
        return create_client(SUPABASE_URL, SUPABASE_KEY)
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        sys.exit(1)

def execute_sql(supabase, sql, description, dry_run=False):
    """Execute SQL using the Supabase client."""
    if dry_run:
        logger.info(f"DRY RUN: Would execute SQL for: {description}")
        logger.info(f"SQL: {sql}")
        return True, None
    
    try:
        logger.info(f"Executing SQL: {description}")
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False, response.error
        
        logger.info(f"SQL executed successfully: {description}")
        
        # Track the change for potential rollback
        applied_changes.append({
            "description": description,
            "timestamp": datetime.now().isoformat(),
            "sql": sql
        })
        
        return True, response.data
    except Exception as e:
        logger.error(f"Error executing SQL ({description}): {str(e)}")
        return False, str(e)

def backup_database(supabase, dry_run=False):
    """Create a backup of the current database state."""
    if dry_run:
        logger.info("DRY RUN: Would create database backup")
        return "dry-run-backup.json"
    
    logger.info("Creating database backup before making changes...")
    
    # Get all foreign key constraints
    fk_sql = """
    SELECT
        tc.constraint_name,
        tc.table_name,
        kcu.column_name,
        ccu.table_name AS referenced_table,
        ccu.column_name AS referenced_column,
        rc.delete_rule
    FROM
        information_schema.table_constraints AS tc
        JOIN information_schema.key_column_usage AS kcu
          ON tc.constraint_name = kcu.constraint_name
          AND tc.table_schema = kcu.table_schema
        JOIN information_schema.constraint_column_usage AS ccu
          ON ccu.constraint_name = tc.constraint_name
          AND ccu.table_schema = tc.table_schema
        JOIN information_schema.referential_constraints AS rc
          ON rc.constraint_name = tc.constraint_name
    WHERE tc.constraint_type = 'FOREIGN KEY'
        AND tc.table_schema = 'public';
    """
    
    success, fk_data = execute_sql(supabase, fk_sql, "Getting all foreign key constraints for backup")
    
    if not success:
        logger.error(f"Failed to get foreign key constraints for backup: {fk_data}")
        return None
    
    # Get all indexes
    idx_sql = """
    SELECT
        tablename,
        indexname,
        indexdef
    FROM
        pg_indexes
    WHERE
        schemaname = 'public';
    """
    
    success, idx_data = execute_sql(supabase, idx_sql, "Getting all indexes for backup")
    
    if not success:
        logger.error(f"Failed to get indexes for backup: {idx_data}")
        return None
    
    # Create backup data
    backup_data = {
        "foreign_keys": fk_data,
        "indexes": idx_data,
        "timestamp": datetime.now().isoformat()
    }
    
    # Save backup to file
    backup_filename = f"fk_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    try:
        with open(backup_filename, 'w') as f:
            json.dump(backup_data, f, indent=2)
        logger.info(f"Database backup saved to {backup_filename}")
        return backup_filename
    except Exception as e:
        logger.error(f"Error saving backup to file: {str(e)}")
        return None

def get_all_foreign_keys(supabase):
    """Get all foreign key constraints in the database with deletion behavior."""
    sql = """
    SELECT
        tc.constraint_name,
        tc.table_name,
        kcu.column_name,
        ccu.table_name AS referenced_table,
        ccu.column_name AS referenced_column,
        rc.delete_rule
    FROM
        information_schema.table_constraints AS tc
        JOIN information_schema.key_column_usage AS kcu
          ON tc.constraint_name = kcu.constraint_name
          AND tc.table_schema = kcu.table_schema
        JOIN information_schema.constraint_column_usage AS ccu
          ON ccu.constraint_name = tc.constraint_name
          AND ccu.table_schema = tc.table_schema
        JOIN information_schema.referential_constraints AS rc
          ON rc.constraint_name = tc.constraint_name
    WHERE tc.constraint_type = 'FOREIGN KEY'
        AND tc.table_schema = 'public';
    """
    
    success, result = execute_sql(supabase, sql, "Getting all foreign key constraints")
    
    if not success:
        logger.error(f"Failed to get foreign key constraints: {result}")
        return None
    
    # Handle different possible result formats
    if result:
        if isinstance(result, list):
            if len(result) > 0:
                # Check if the first item has the expected keys
                first_item = result[0]
                if isinstance(first_item, dict) and 'constraint_name' in first_item:
                    return result
                else:
                    logger.error(f"Unexpected result format from get_all_foreign_keys: {type(first_item)}")
                    logger.error(f"First item content: {first_item}")
                    return []
            else:
                logger.warning("No foreign key constraints found")
                return []
        elif isinstance(result, dict):
            # Check if it's a success response
            if 'success' in result and result['success']:
                logger.warning("Received success response but no foreign key data")
                return []
            else:
                logger.error(f"Unexpected dictionary result: {result}")
                return []
        else:
            logger.error(f"Unexpected result type: {type(result)}")
            return []
    else:
        logger.warning("No result from foreign key query")
        return []

def get_all_tables(supabase):
    """Get all tables in the database."""
    # Let's try a direct approach using a hardcoded list of expected tables
    # This is a fallback for when the database queries don't work as expected
    expected_tables = [
        "molecules", "mixtures", "mixture_components", "predictions",
        "experiments", "experiment_properties", "property_types",
        "calculation_methods", "projects", "teams", "user_profile",
        "molecular_properties"
    ]
    
    # First try the standard method
    sql = """
    SELECT table_name
    FROM information_schema.tables
    WHERE table_schema = 'public'
    AND table_type = 'BASE TABLE';
    """
    
    success, result = execute_sql(supabase, sql, "Getting all tables")
    
    if success:
        # Handle different possible result formats
        if result:
            if isinstance(result, list):
                if len(result) > 0:
                    if all(isinstance(item, dict) for item in result):
                        # Result is a list of dictionaries
                        return [table['table_name'] for table in result]
                    elif all(isinstance(item, str) for item in result):
                        # Result is a list of strings
                        return result
    
    # If the first method didn't work, try an alternative method
    alt_sql = """
    SELECT tablename FROM pg_catalog.pg_tables
    WHERE schemaname = 'public';
    """
    
    try:
        alt_success, alt_result = execute_sql(supabase, alt_sql, "Getting all tables (alternative method)")
        if alt_success:
            if isinstance(alt_result, list) and len(alt_result) > 0:
                if all(isinstance(item, dict) for item in alt_result):
                    return [table['tablename'] for table in alt_result]
                elif all(isinstance(item, str) for item in alt_result):
                    return alt_result
    except Exception as e:
        logger.error(f"Error in alternative table query: {str(e)}")
    
    # If all else fails, use the hardcoded list
    logger.warning("Using hardcoded table list as fallback")
    return expected_tables

def get_table_columns(supabase, table_name):
    """Get all columns for a specific table."""
    sql = f"""
    SELECT column_name, data_type, is_nullable
    FROM information_schema.columns
    WHERE table_schema = 'public'
    AND table_name = '{table_name}'
    ORDER BY ordinal_position;
    """
    
    success, result = execute_sql(supabase, sql, f"Getting columns for table {table_name}")
    
    if not success:
        logger.error(f"Failed to get columns for table {table_name}: {result}")
        return None
    
    # Handle different possible result formats
    if result:
        if isinstance(result, list) and len(result) > 0:
            # Check if the first item has the expected keys
            first_item = result[0]
            if isinstance(first_item, dict) and 'column_name' in first_item:
                return result
            else:
                logger.error(f"Unexpected result format from get_table_columns: {type(first_item)}")
                logger.error(f"First item content: {first_item}")
                return []
        else:
            logger.warning(f"No columns found for table {table_name}")
            return []
    else:
        return []

def get_duplicate_foreign_keys(foreign_keys):
    """Identify duplicate foreign key constraints."""
    # Group foreign keys by table, column, referenced table, and referenced column
    grouped_fks = {}
    for fk in foreign_keys:
        key = (fk['table_name'], fk['column_name'], fk['referenced_table'], fk['referenced_column'])
        if key not in grouped_fks:
            grouped_fks[key] = []
        grouped_fks[key].append(fk)
    
    # Find groups with more than one foreign key
    duplicates = []
    for key, fks in grouped_fks.items():
        if len(fks) > 1:
            duplicates.extend(fks[1:])  # Keep the first one, mark the rest as duplicates
    
    return duplicates

def get_outdated_foreign_keys(foreign_keys, table_mapping):
    """Identify foreign keys that reference old singular table names."""
    outdated_fks = []
    
    for fk in foreign_keys:
        referenced_table = fk['referenced_table']
        if referenced_table in table_mapping:
            # This foreign key references an old singular table name
            outdated_fks.append(fk)
    
    return outdated_fks

def get_incorrect_delete_behavior_fks(foreign_keys, delete_behavior):
    """Identify foreign keys with incorrect ON DELETE behavior."""
    incorrect_fks = []
    
    for fk in foreign_keys:
        table_name = fk['table_name']
        column_name = fk['column_name']
        referenced_table = fk['referenced_table']
        current_delete_rule = fk['delete_rule']
        
        # Check if this relationship has a specific delete behavior defined
        key = (referenced_table, table_name, column_name)
        expected_delete_rule = None
        
        # Try to find the relationship in DELETE_BEHAVIOR
        if key in delete_behavior:
            expected_delete_rule = delete_behavior[key]
        else:
            # Check if the referenced table is in singular form
            for singular, plural in TABLE_MAPPING.items():
                if referenced_table == singular:
                    key = (plural, table_name, column_name)
                    if key in delete_behavior:
                        expected_delete_rule = delete_behavior[key]
                        break
        
        # If no specific rule is defined, default to CASCADE
        if expected_delete_rule is None:
            expected_delete_rule = "CASCADE"
        
        # Compare current rule with expected rule
        if current_delete_rule != expected_delete_rule:
            fk['expected_delete_rule'] = expected_delete_rule
            incorrect_fks.append(fk)
    
    return incorrect_fks

def get_missing_indexes(supabase, foreign_keys):
    """Identify foreign keys without corresponding indexes."""
    # Get all indexes
    sql = """
    SELECT
        tablename,
        indexname,
        indexdef
    FROM
        pg_indexes
    WHERE
        schemaname = 'public';
    """
    
    success, indexes = execute_sql(supabase, sql, "Getting all indexes")
    
    if not success:
        logger.error(f"Failed to get indexes: {indexes}")
        return None
    
    # Handle different possible result formats
    if not indexes or not isinstance(indexes, list):
        logger.warning("No indexes found or unexpected result format")
        return foreign_keys  # Assume all foreign keys need indexes
    
    # Check if the result has the expected format
    if len(indexes) > 0 and not isinstance(indexes[0], dict):
        logger.error(f"Unexpected index result format: {type(indexes[0])}")
        return foreign_keys  # Assume all foreign keys need indexes
    
    # Create a dictionary of table -> list of indexed columns
    table_indexes = {}
    for idx in indexes:
        # Skip if the index doesn't have the expected keys
        if not all(key in idx for key in ['tablename', 'indexname', 'indexdef']):
            continue
            
        table_name = idx['tablename']
        if table_name not in table_indexes:
            table_indexes[table_name] = []
        
        # Extract column names from indexdef
        # This is a simplified approach and might not work for all index types
        indexdef = idx['indexdef']
        if '(' in indexdef and ')' in indexdef:
            columns_str = indexdef.split('(')[1].split(')')[0]
            columns = [col.strip() for col in columns_str.split(',')]
            table_indexes[table_name].extend(columns)
    
    # Find foreign keys without indexes
    missing_indexes = []
    for fk in foreign_keys:
        # Skip if the foreign key doesn't have the expected keys
        if not all(key in fk for key in ['table_name', 'column_name']):
            continue
            
        table_name = fk['table_name']
        column_name = fk['column_name']
        
        if table_name not in table_indexes or column_name not in table_indexes[table_name]:
            missing_indexes.append(fk)
    
    return missing_indexes

def get_missing_foreign_keys(supabase, tables, foreign_keys):
    """Identify potential missing foreign key constraints."""
    missing_fks = []
    
    # Get existing foreign key relationships
    existing_fks = set()
    for fk in foreign_keys:
        key = (fk['table_name'], fk['column_name'], fk['referenced_table'], fk['referenced_column'])
        existing_fks.add(key)
    
    # Check each table for potential foreign key columns
    for table_name in tables:
        columns = get_table_columns(supabase, table_name)
        if not columns:
            continue
        
        for column in columns:
            column_name = column['column_name']
            
            # Skip primary key columns
            if column_name == 'id':
                continue
            
            # Check if column name suggests a foreign key
            if column_name.endswith('_id'):
                # Determine potential referenced table
                potential_table = column_name[:-3]  # Remove '_id'
                
                # Check if it's a plural table name
                if potential_table in TABLE_MAPPING.values():
                    referenced_table = potential_table
                # Check if it's a singular table name
                elif potential_table in TABLE_MAPPING:
                    referenced_table = TABLE_MAPPING[potential_table]
                else:
                    # Not a recognized table name pattern
                    continue
                
                # Check if the referenced table exists
                if referenced_table in tables:
                    # Check if this relationship already exists
                    key = (table_name, column_name, referenced_table, 'id')
                    if key not in existing_fks:
                        missing_fks.append({
                            'table_name': table_name,
                            'column_name': column_name,
                            'referenced_table': referenced_table,
                            'referenced_column': 'id'
                        })
    
    return missing_fks

def begin_transaction(supabase):
    """Begin a database transaction."""
    sql = "BEGIN;"
    success, _ = execute_sql(supabase, sql, "Beginning transaction")
    return success

def commit_transaction(supabase):
    """Commit the current transaction."""
    sql = "COMMIT;"
    success, _ = execute_sql(supabase, sql, "Committing transaction")
    return success

def rollback_transaction(supabase):
    """Rollback the current transaction."""
    sql = "ROLLBACK;"
    success, _ = execute_sql(supabase, sql, "Rolling back transaction")
    return success

def remove_duplicate_foreign_keys(supabase, duplicate_fks, dry_run=False):
    """Remove duplicate foreign key constraints."""
    if not duplicate_fks:
        logger.info("No duplicate foreign key constraints found")
        return True
    
    logger.info(f"Found {len(duplicate_fks)} duplicate foreign key constraints to remove")
    
    success = True
    for fk in duplicate_fks:
        constraint_name = fk['constraint_name']
        table_name = fk['table_name']
        
        sql = f"""
        ALTER TABLE public.{table_name}
        DROP CONSTRAINT IF EXISTS {constraint_name};
        """
        
        result, _ = execute_sql(supabase, sql, f"Removing duplicate foreign key constraint {constraint_name}", dry_run)
        if not result and not dry_run:
            success = False
    
    return success

def update_outdated_foreign_keys(supabase, outdated_fks, table_mapping, dry_run=False):
    """Update foreign keys that reference old singular table names."""
    if not outdated_fks:
        logger.info("No outdated foreign key references found")
        return True
    
    logger.info(f"Found {len(outdated_fks)} outdated foreign key references to update")
    
    success = True
    for fk in outdated_fks:
        constraint_name = fk['constraint_name']
        table_name = fk['table_name']
        column_name = fk['column_name']
        referenced_table = fk['referenced_table']
        referenced_column = fk['referenced_column']
        
        # Get the new plural table name
        new_referenced_table = table_mapping.get(referenced_table, referenced_table)
        
        # Determine appropriate ON DELETE behavior
        key = (new_referenced_table, table_name, column_name)
        delete_behavior = DELETE_BEHAVIOR.get(key, "CASCADE")
        
        # Drop the old constraint
        drop_sql = f"""
        ALTER TABLE public.{table_name}
        DROP CONSTRAINT IF EXISTS {constraint_name};
        """
        
        # Add the new constraint
        add_sql = f"""
        ALTER TABLE public.{table_name}
        ADD CONSTRAINT {constraint_name}
        FOREIGN KEY ({column_name})
        REFERENCES public.{new_referenced_table} ({referenced_column})
        ON DELETE {delete_behavior};
        """
        
        combined_sql = drop_sql + "\n" + add_sql
        
        result, _ = execute_sql(supabase, combined_sql, f"Updating foreign key constraint {constraint_name}", dry_run)
        if not result and not dry_run:
            success = False
    
    return success

def fix_delete_behavior(supabase, incorrect_fks, dry_run=False):
    """Fix foreign keys with incorrect ON DELETE behavior."""
    if not incorrect_fks:
        logger.info("No foreign keys with incorrect ON DELETE behavior found")
        return True
    
    logger.info(f"Found {len(incorrect_fks)} foreign keys with incorrect ON DELETE behavior")
    
    success = True
    for fk in incorrect_fks:
        constraint_name = fk['constraint_name']
        table_name = fk['table_name']
        column_name = fk['column_name']
        referenced_table = fk['referenced_table']
        referenced_column = fk['referenced_column']
        expected_delete_rule = fk['expected_delete_rule']
        
        # Drop the old constraint
        drop_sql = f"""
        ALTER TABLE public.{table_name}
        DROP CONSTRAINT IF EXISTS {constraint_name};
        """
        
        # Add the new constraint with correct ON DELETE behavior
        add_sql = f"""
        ALTER TABLE public.{table_name}
        ADD CONSTRAINT {constraint_name}
        FOREIGN KEY ({column_name})
        REFERENCES public.{referenced_table} ({referenced_column})
        ON DELETE {expected_delete_rule};
        """
        
        combined_sql = drop_sql + "\n" + add_sql
        
        result, _ = execute_sql(supabase, combined_sql, f"Fixing ON DELETE behavior for {constraint_name}", dry_run)
        if not result and not dry_run:
            success = False
    
    return success

def add_missing_indexes(supabase, missing_indexes, dry_run=False):
    """Add indexes for foreign keys that don't have them."""
    if not missing_indexes:
        logger.info("No missing indexes found")
        return True
    
    logger.info(f"Found {len(missing_indexes)} foreign keys without indexes")
    
    success = True
    for fk in missing_indexes:
        table_name = fk['table_name']
        column_name = fk['column_name']
        
        # Create index name
        index_name = f"idx_{table_name}_{column_name}"
        
        sql = f"""
        CREATE INDEX IF NOT EXISTS {index_name}
        ON public.{table_name} ({column_name});
        """
        
        result, _ = execute_sql(supabase, sql, f"Creating index {index_name}", dry_run)
        if not result and not dry_run:
            success = False
    
    return success

def add_missing_foreign_keys(supabase, missing_fks, dry_run=False):
    """Add missing foreign key constraints."""
    if not missing_fks:
        logger.info("No missing foreign key constraints found")
        return True
    
    logger.info(f"Found {len(missing_fks)} missing foreign key constraints")
    
    success = True
    for fk in missing_fks:
        table_name = fk['table_name']
        column_name = fk['column_name']
        referenced_table = fk['referenced_table']
        referenced_column = fk['referenced_column']
        
        # Determine appropriate constraint name
        constraint_name = f"{table_name}_{column_name}_fkey"
        
        # Determine appropriate ON DELETE behavior
        key = (referenced_table, table_name, column_name)
        delete_behavior = DELETE_BEHAVIOR.get(key, "CASCADE")
        
        sql = f"""
        ALTER TABLE public.{table_name}
        ADD CONSTRAINT {constraint_name}
        FOREIGN KEY ({column_name})
        REFERENCES public.{referenced_table} ({referenced_column})
        ON DELETE {delete_behavior};
        """
        
        result, _ = execute_sql(supabase, sql, f"Adding foreign key constraint {constraint_name}", dry_run)
        if not result and not dry_run:
            success = False
    
    return success

def verify_foreign_keys(supabase):
    """Verify that all foreign key constraints are correctly set up."""
    logger.info("Verifying foreign key constraints...")
    
    # Get all foreign keys
    foreign_keys = get_all_foreign_keys(supabase)
    if not foreign_keys:
        logger.error("Failed to get foreign key constraints for verification")
        return False
    
    # Get all tables
    tables = get_all_tables(supabase)
    if not tables:
        logger.error("Failed to get tables for verification")
        return False
    
    # Check for issues
    duplicate_fks = get_duplicate_foreign_keys(foreign_keys)
    outdated_fks = get_outdated_foreign_keys(foreign_keys, TABLE_MAPPING)
    incorrect_delete_fks = get_incorrect_delete_behavior_fks(foreign_keys, DELETE_BEHAVIOR)
    missing_indexes = get_missing_indexes(supabase, foreign_keys)
    missing_fks = get_missing_foreign_keys(supabase, tables, foreign_keys)
    
    if missing_indexes is None:
        logger.error("Failed to check for missing indexes during verification")
        return False
    
    # Check if there are any issues
    has_issues = (
        len(duplicate_fks) > 0 or
        len(outdated_fks) > 0 or
        len(incorrect_delete_fks) > 0 or
        len(missing_indexes) > 0 or
        len(missing_fks) > 0
    )
    
    # Print verification results
    logger.info("\n" + "=" * 80)
    logger.info("CryoProtect v2 Foreign Key Verification Results")
    logger.info("=" * 80)
    
    logger.info(f"\nIssues found:")
    logger.info(f"- Duplicate foreign key constraints: {len(duplicate_fks)}")
    logger.info(f"- Outdated foreign key references: {len(outdated_fks)}")
    logger.info(f"- Foreign keys with incorrect ON DELETE behavior: {len(incorrect_delete_fks)}")
    logger.info(f"- Foreign keys without indexes: {len(missing_indexes)}")
    logger.info(f"- Missing foreign key constraints: {len(missing_fks)}")
    
    if has_issues:
        logger.info("\nStatus: VERIFICATION FAILED")
        logger.info("There are still foreign key issues that need to be fixed.")
    else:
        logger.info("\nStatus: VERIFICATION PASSED")
        logger.info("All foreign key constraints are correctly set up.")
    
    logger.info("=" * 80 + "\n")
    
    return not has_issues

def rollback_changes_from_backup(supabase, backup_filename, dry_run=False):
    """Rollback changes using a backup file."""
    if dry_run:
        logger.info("DRY RUN: Would rollback changes from backup")
        return True
    
    logger.info(f"Rolling back changes from backup file: {backup_filename}")
    
    try:
        # Load backup data
        with open(backup_filename, 'r') as f:
            backup_data = json.load(f)
        
        if 'foreign_keys' not in backup_data or 'indexes' not in backup_data:
            logger.error("Invalid backup file format")
            return False
        
        # Begin transaction
        if not begin_transaction(supabase):
            logger.error("Failed to begin transaction for rollback")
            return False
        
        # Get current foreign keys
        current_fks = get_all_foreign_keys(supabase)
        if not current_fks:
            logger.error("Failed to get current foreign key constraints for rollback")
            rollback_transaction(supabase)
            return False
        
        # Drop all current foreign key constraints
        for fk in current_fks:
            constraint_name = fk['constraint_name']
            table_name = fk['table_name']
            
            sql = f"""
            ALTER TABLE public.{table_name}
            DROP CONSTRAINT IF EXISTS {constraint_name};
            """
            
            success, _ = execute_sql(supabase, sql, f"Dropping constraint {constraint_name} for rollback")
            if not success:
                logger.error(f"Failed to drop constraint {constraint_name} for rollback")
                rollback_transaction(supabase)
                return False
        
        # Recreate foreign key constraints from backup
        for fk in backup_data['foreign_keys']:
            constraint_name = fk['constraint_name']
            table_name = fk['table_name']
            column_name = fk['column_name']
            referenced_table = fk['referenced_table']
            referenced_column = fk['referenced_column']
            delete_rule = fk['delete_rule']
            
            sql = f"""
            ALTER TABLE public.{table_name}
            ADD CONSTRAINT {constraint_name}
            FOREIGN KEY ({column_name})
            REFERENCES public.{referenced_table} ({referenced_column})
            ON DELETE {delete_rule};
            """
            
            success, _ = execute_sql(supabase, sql, f"Recreating constraint {constraint_name} from backup")
            if not success:
                logger.error(f"Failed to recreate constraint {constraint_name} from backup")
                rollback_transaction(supabase)
                return False
        
        # Commit transaction
        if not commit_transaction(supabase):
            logger.error("Failed to commit transaction for rollback")
            rollback_transaction(supabase)
            return False
        
        logger.info("Successfully rolled back changes from backup")
        return True
    except Exception as e:
        logger.error(f"Error rolling back changes from backup: {str(e)}")
        rollback_transaction(supabase)
        return False

def main():
    """Main function to fix foreign key relationships."""
    parser = argparse.ArgumentParser(description="Fix foreign key relationships in CryoProtect database")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without making changes")
    parser.add_argument("--verify-only", action="store_true", help="Only verify the foreign key constraints without making changes")
    parser.add_argument("--rollback", action="store_true", help="Rollback changes using the most recent backup")
    parser.add_argument("--backup-file", help="Specify a backup file to use for rollback")
    args = parser.parse_args()
    
    logger.info("Starting foreign key relationship fixes")
    
    # Connect to Supabase
    supabase = get_supabase_client()
    
    # Handle rollback if requested
    if args.rollback:
        backup_file = args.backup_file
        if not backup_file:
            # Find the most recent backup file
            backup_files = [f for f in os.listdir('.') if f.startswith('fk_backup_') and f.endswith('.json')]
            if not backup_files:
                logger.error("No backup files found for rollback")
                return 1
            backup_file = sorted(backup_files)[-1]  # Get the most recent backup
        
        success = rollback_changes_from_backup(supabase, backup_file, args.dry_run)
        return 0 if success else 1
    
    # Handle verify-only if requested
    if args.verify_only:
        success = verify_foreign_keys(supabase)
        return 0 if success else 1
    
    # Create a backup before making changes
    backup_file = backup_database(supabase, args.dry_run)
    if not backup_file and not args.dry_run:
        logger.error("Failed to create database backup")
        return 1
    
    # Begin transaction
    if not args.dry_run:
        if not begin_transaction(supabase):
            logger.error("Failed to begin transaction")
            return 1
    
    # Get all foreign key constraints
    foreign_keys = get_all_foreign_keys(supabase)
    if not foreign_keys and not args.dry_run:
        logger.error("Failed to get foreign key constraints")
        rollback_transaction(supabase)
        return 1
    
    # Get all tables
    tables = get_all_tables(supabase)
    if not tables and not args.dry_run:
        logger.error("Failed to get tables")
        rollback_transaction(supabase)
        return 1
    
    logger.info(f"Found {len(foreign_keys)} foreign key constraints")
    
    # Identify issues
    duplicate_fks = get_duplicate_foreign_keys(foreign_keys)
    outdated_fks = get_outdated_foreign_keys(foreign_keys, TABLE_MAPPING)
    incorrect_delete_fks = get_incorrect_delete_behavior_fks(foreign_keys, DELETE_BEHAVIOR)
    missing_indexes = get_missing_indexes(supabase, foreign_keys)
    missing_fks = get_missing_foreign_keys(supabase, tables, foreign_keys)
    
    if missing_indexes is None and not args.dry_run:
        logger.error("Failed to check for missing indexes")
        rollback_transaction(supabase)
        return 1
    
    # Fix issues
    success = True
    
    # 1. Remove duplicate foreign key constraints
    if not remove_duplicate_foreign_keys(supabase, duplicate_fks, args.dry_run):
        logger.error("Failed to remove duplicate foreign key constraints")
        success = False
    
    # 2. Update outdated foreign key references
    if not update_outdated_foreign_keys(supabase, outdated_fks, TABLE_MAPPING, args.dry_run):
        logger.error("Failed to update outdated foreign key references")
        success = False
    
    # 3. Fix foreign keys with incorrect ON DELETE behavior
    if not fix_delete_behavior(supabase, incorrect_delete_fks, args.dry_run):
        logger.error("Failed to fix foreign keys with incorrect ON DELETE behavior")
        success = False
    
    # 4. Add missing indexes
    if not add_missing_indexes(supabase, missing_indexes, args.dry_run):
        logger.error("Failed to add missing indexes")
        success = False
    
    # 5. Add missing foreign key constraints
    if not add_missing_foreign_keys(supabase, missing_fks, args.dry_run):
        logger.error("Failed to add missing foreign key constraints")
        success = False
    
    # Commit or rollback transaction
    if not args.dry_run:
        if success:
            if not commit_transaction(supabase):
                logger.error("Failed to commit transaction")
                rollback_transaction(supabase)
                success = False
        else:
            logger.warning("Errors occurred, rolling back changes")
            rollback_transaction(supabase)
    
    # Verify the fixes if not in dry-run mode
    verification_success = True
    if not args.dry_run and success:
        logger.info("Verifying fixes...")
        verification_success = verify_foreign_keys(supabase)
        if not verification_success:
            logger.warning("Verification failed after applying fixes")
    
    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("CryoProtect v2 Foreign Key Relationship Fixes Summary")
    logger.info("=" * 80)
    
    if args.dry_run:
        logger.info("\nDRY RUN - No changes were made")
    
    logger.info(f"\nIssues found:")
    logger.info(f"- Duplicate foreign key constraints: {len(duplicate_fks)}")
    logger.info(f"- Outdated foreign key references: {len(outdated_fks)}")
    logger.info(f"- Foreign keys with incorrect ON DELETE behavior: {len(incorrect_delete_fks)}")
    logger.info(f"- Foreign keys without indexes: {len(missing_indexes)}")
    logger.info(f"- Missing foreign key constraints: {len(missing_fks)}")
    
    if args.dry_run:
        logger.info("\nStatus: DRY RUN COMPLETED")
    elif success and verification_success:
        logger.info("\nStatus: SUCCESS")
        logger.info("All foreign key relationship issues were fixed successfully.")
    elif success and not verification_success:
        logger.info("\nStatus: PARTIAL SUCCESS")
        logger.info("Changes were applied but verification found remaining issues.")
    else:
        logger.info("\nStatus: FAILED")
        logger.info("Some issues were encountered. Changes were rolled back.")
    
    logger.info("\nFor detailed information, check the log file.")
    logger.info("=" * 80 + "\n")
    
    return 0 if (success and (verification_success or args.dry_run)) else 1

if __name__ == "__main__":
    sys.exit(main())
