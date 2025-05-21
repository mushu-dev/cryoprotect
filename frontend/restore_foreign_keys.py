#!/usr/bin/env python3
"""
Script to restore foreign key constraints in the CryoProtect database.
This script applies the SQL migration to restore the foreign key constraints
that were previously dropped.
"""

import os
import sys
import argparse
import json
import logging
from datetime import datetime
import db_utils

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def execute_migration():
    """Execute the migration to restore foreign key constraints."""
    sql_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
        'migrations', 
        '026_restore_foreign_keys.sql'
    )
    
    try:
        if db_utils.execute_sql_file(sql_file):
            logger.info("Successfully executed migration to restore foreign key constraints")
            return True
        else:
            logger.error("Failed to execute migration")
            return False
    except Exception as e:
        logger.error(f"Error executing migration: {e}")
        return False

def get_foreign_key_stats():
    """Get statistics on foreign key constraints before and after the migration."""
    query = """
    SELECT
        tc.table_schema, 
        tc.constraint_name, 
        tc.table_name as source_table, 
        kcu.column_name as source_column, 
        ccu.table_name AS target_table,
        ccu.column_name AS target_column 
    FROM 
        information_schema.table_constraints AS tc 
    JOIN 
        information_schema.key_column_usage AS kcu ON tc.constraint_name = kcu.constraint_name
        AND tc.table_schema = kcu.table_schema
    JOIN 
        information_schema.constraint_column_usage AS ccu ON ccu.constraint_name = tc.constraint_name
        AND ccu.table_schema = tc.table_schema
    WHERE 
        tc.constraint_type = 'FOREIGN KEY' 
        AND tc.table_schema = 'public'
    ORDER BY 
        tc.table_name, 
        kcu.column_name;
    """
    
    try:
        foreign_keys = db_utils.execute_query(query, cursor_factory=db_utils.RealDictCursor)
        
        # Group by source table
        tables = {}
        for fk in foreign_keys:
            source_table = fk['source_table']
            if source_table not in tables:
                tables[source_table] = []
            tables[source_table].append(fk)
        
        return {
            'foreign_keys': foreign_keys,
            'total_foreign_keys': len(foreign_keys),
            'tables_with_foreign_keys': len(tables),
            'foreign_keys_by_table': tables
        }
    except Exception as e:
        logger.error(f"Error getting foreign key stats: {e}")
        return {
            'foreign_keys': [],
            'total_foreign_keys': 0,
            'tables_with_foreign_keys': 0,
            'foreign_keys_by_table': {}
        }

def check_orphaned_records():
    """Check for orphaned records that would violate foreign key constraints."""
    orphaned_records = []
    
    # Define relationships to check
    relationships = [
        {'source_table': 'molecular_properties', 'source_column': 'molecule_id', 'target_table': 'molecules', 'target_column': 'id'},
        {'source_table': 'molecular_properties', 'source_column': 'property_type_id', 'target_table': 'property_types', 'target_column': 'id'},
        {'source_table': 'mixture_components', 'source_column': 'mixture_id', 'target_table': 'mixtures', 'target_column': 'id'},
        {'source_table': 'mixture_components', 'source_column': 'molecule_id', 'target_table': 'molecules', 'target_column': 'id'},
        {'source_table': 'molecule_experiments', 'source_column': 'molecule_id', 'target_table': 'molecules', 'target_column': 'id'},
        {'source_table': 'consolidated_molecules', 'source_column': 'primary_molecule_id', 'target_table': 'molecules', 'target_column': 'id'},
    ]
    
    for rel in relationships:
        query = f"""
        SELECT COUNT(*) as orphaned_count 
        FROM {rel['source_table']} 
        WHERE {rel['source_column']} IS NOT NULL 
        AND NOT EXISTS (
            SELECT 1 
            FROM {rel['target_table']} 
            WHERE {rel['target_column']} = {rel['source_table']}.{rel['source_column']}
        )
        """
        
        try:
            result = db_utils.execute_query(query, cursor_factory=db_utils.RealDictCursor)[0]
            orphaned_count = result['orphaned_count']
            
            if orphaned_count > 0:
                orphaned_records.append({
                    'source_table': rel['source_table'],
                    'source_column': rel['source_column'],
                    'target_table': rel['target_table'],
                    'target_column': rel['target_column'],
                    'orphaned_count': orphaned_count
                })
        except Exception as e:
            logger.error(f"Error checking orphaned records for {rel['source_table']}.{rel['source_column']}: {e}")
    
    return orphaned_records

def save_report(stats, filename_prefix="foreign_key_restoration_report"):
    """Save the statistics report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    try:
        with open(filename, 'w') as f:
            json.dump(stats, f, indent=2)
        logger.info(f"Report saved to {filename}")
    except Exception as e:
        logger.error(f"Error saving report: {e}")

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Restore foreign key constraints in the CryoProtect database.")
    parser.add_argument("--check-only", action="store_true", help="Check for orphaned records but don't execute the migration")
    parser.add_argument("--fix-orphaned", action="store_true", help="Fix orphaned records by removing them (USE WITH CAUTION)")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Get foreign key stats before migration
    logger.info("Getting foreign key stats before migration...")
    before_stats = get_foreign_key_stats()
    
    # Check for orphaned records
    logger.info("Checking for orphaned records...")
    orphaned_records = check_orphaned_records()
    
    if orphaned_records:
        logger.warning(f"Found {len(orphaned_records)} relationships with orphaned records")
        for rel in orphaned_records:
            logger.warning(f"  {rel['source_table']}.{rel['source_column']} -> {rel['target_table']}.{rel['target_column']}: {rel['orphaned_count']} orphaned records")
        
        if args.fix_orphaned:
            logger.info("Fixing orphaned records...")
            for rel in orphaned_records:
                query = f"""
                DELETE FROM {rel['source_table']} 
                WHERE {rel['source_column']} IS NOT NULL 
                AND NOT EXISTS (
                    SELECT 1 
                    FROM {rel['target_table']} 
                    WHERE {rel['target_column']} = {rel['source_table']}.{rel['source_column']}
                )
                """
                
                try:
                    db_utils.execute_query(query, fetch=False)
                    logger.info(f"  Removed orphaned records from {rel['source_table']}")
                except Exception as e:
                    logger.error(f"  Error fixing orphaned records in {rel['source_table']}: {e}")
    else:
        logger.info("No orphaned records found")
    
    if args.check_only:
        logger.info("Check-only mode. Not executing migration.")
        sys.exit(0)
    
    # Execute the migration
    logger.info("Executing migration to restore foreign key constraints...")
    success = execute_migration()
    
    if not success:
        logger.error("Migration failed. Exiting.")
        sys.exit(1)
    
    # Get foreign key stats after migration
    logger.info("Getting foreign key stats after migration...")
    after_stats = get_foreign_key_stats()
    
    # Save report
    stats = {
        'timestamp': datetime.now().isoformat(),
        'before': before_stats,
        'after': after_stats,
        'orphaned_records': orphaned_records,
        'fixed_orphaned': args.fix_orphaned and len(orphaned_records) > 0
    }
    save_report(stats)
    
    logger.info(f"Foreign key constraints before: {before_stats['total_foreign_keys']}")
    logger.info(f"Foreign key constraints after: {after_stats['total_foreign_keys']}")
    logger.info(f"Increase: {after_stats['total_foreign_keys'] - before_stats['total_foreign_keys']}")

if __name__ == "__main__":
    main()