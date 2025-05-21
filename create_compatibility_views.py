#!/usr/bin/env python3
"""
CryoProtect v2 - Create Compatibility Views

This script creates database views that map singular table names to their plural counterparts,
providing backward compatibility during the migration process.

Tables:
1. molecule → molecules
2. mixture → mixtures
3. experiment → experiments
4. prediction → predictions
5. project → projects

Author: CryoProtect Team
Date: 2025-04-18
"""

import os
import sys
import logging
import argparse
from datetime import datetime
import psycopg2
from psycopg2 import sql
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"compatibility_views_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Database connection parameters
DB_URL = os.getenv("DATABASE_URL") or os.getenv("SUPABASE_DB_URL")
if not DB_URL:
    logger.error("Database URL not found in environment variables")
    sys.exit(1)

# Table mappings
TABLE_MAPPINGS = [
    {"singular": "molecule", "plural": "molecules"},
    {"singular": "mixture", "plural": "mixtures"},
    {"singular": "experiment", "plural": "experiments"},
    {"singular": "prediction", "plural": "predictions"},
    {"singular": "project", "plural": "projects"}
]

def connect_to_db():
    """
    Establish a connection to the database.
    
    Returns:
        tuple: (connection, cursor) or (None, None) if connection fails
    """
    try:
        conn = psycopg2.connect(DB_URL)
        conn.autocommit = False  # We'll manage transactions explicitly
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        logger.info("Successfully connected to the database")
        return conn, cursor
    except Exception as e:
        logger.error(f"Error connecting to the database: {str(e)}")
        return None, None

def get_table_columns(cursor, table_name):
    """
    Get column names for a table.
    
    Args:
        cursor: Database cursor
        table_name (str): Name of the table
        
    Returns:
        list: List of column names
    """
    try:
        cursor.execute("""
            SELECT column_name
            FROM information_schema.columns
            WHERE table_name = %s
            ORDER BY ordinal_position
        """, (table_name,))
        return [row['column_name'] for row in cursor.fetchall()]
    except Exception as e:
        logger.error(f"Error getting columns for table {table_name}: {str(e)}")
        return []

def check_view_exists(cursor, view_name):
    """
    Check if a view already exists.
    
    Args:
        cursor: Database cursor
        view_name (str): Name of the view
        
    Returns:
        bool: True if view exists, False otherwise
    """
    try:
        cursor.execute("""
            SELECT 1
            FROM information_schema.views
            WHERE table_schema = 'public' AND table_name = %s
        """, (view_name,))
        return cursor.fetchone() is not None
    except Exception as e:
        logger.error(f"Error checking if view {view_name} exists: {str(e)}")
        return False

def create_compatibility_view(cursor, singular_name, plural_name):
    """
    Create a compatibility view that maps a singular table name to its plural counterpart.
    
    Args:
        cursor: Database cursor
        singular_name (str): Singular table name
        plural_name (str): Plural table name
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Check if the plural table exists
        cursor.execute("""
            SELECT 1
            FROM information_schema.tables
            WHERE table_schema = 'public' AND table_name = %s
        """, (plural_name,))
        if not cursor.fetchone():
            logger.error(f"Plural table {plural_name} does not exist")
            return False
        
        # Get columns from the plural table
        columns = get_table_columns(cursor, plural_name)
        if not columns:
            logger.error(f"Could not get columns for table {plural_name}")
            return False
        
        # Check if view already exists
        if check_view_exists(cursor, singular_name):
            logger.info(f"View {singular_name} already exists, dropping it")
            cursor.execute(
                sql.SQL("DROP VIEW IF EXISTS {} CASCADE").format(
                    sql.Identifier(singular_name)
                )
            )
        
        # Create column identifiers
        column_identifiers = [sql.Identifier(col) for col in columns]
        
        # Create the view
        view_sql = sql.SQL("""
            CREATE OR REPLACE VIEW {} AS
            SELECT {}
            FROM {}
        """).format(
            sql.Identifier(singular_name),
            sql.SQL(", ").join(column_identifiers),
            sql.Identifier(plural_name)
        )
        
        cursor.execute(view_sql)
        logger.info(f"Created compatibility view: {singular_name} → {plural_name}")
        
        # Add a comment to the view
        comment_sql = sql.SQL("""
            COMMENT ON VIEW {} IS 'Compatibility view for {}, created on {}'
        """).format(
            sql.Identifier(singular_name),
            sql.Literal(plural_name),
            sql.Literal(datetime.now().isoformat())
        )
        cursor.execute(comment_sql)
        
        return True
    except Exception as e:
        logger.error(f"Error creating compatibility view {singular_name}: {str(e)}")
        return False

def create_all_compatibility_views(dry_run=False):
    """
    Create compatibility views for all table mappings.
    
    Args:
        dry_run (bool): If True, don't commit changes
        
    Returns:
        dict: Report of created views
    """
    conn, cursor = connect_to_db()
    if not conn or not cursor:
        logger.error("Failed to connect to the database")
        return {"success": False, "message": "Database connection failed"}
    
    try:
        logger.info(f"Starting compatibility view creation (dry_run={dry_run})")
        
        created_views = []
        failed_views = []
        
        for mapping in TABLE_MAPPINGS:
            singular = mapping["singular"]
            plural = mapping["plural"]
            
            success = create_compatibility_view(cursor, singular, plural)
            
            if success:
                created_views.append({"singular": singular, "plural": plural})
            else:
                failed_views.append({"singular": singular, "plural": plural})
        
        if dry_run:
            logger.info("Dry run completed, rolling back changes")
            conn.rollback()
        else:
            logger.info("Committing changes")
            conn.commit()
        
        report = {
            "success": len(failed_views) == 0,
            "timestamp": datetime.now().isoformat(),
            "created_views": created_views,
            "failed_views": failed_views,
            "dry_run": dry_run
        }
        
        # Save report to file
        if not dry_run:
            report_file = f"compatibility_views_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            import json
            with open(report_file, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Report saved to {report_file}")
        
        return report
    except Exception as e:
        logger.error(f"Error creating compatibility views: {str(e)}")
        conn.rollback()
        return {"success": False, "message": str(e)}
    finally:
        cursor.close()
        conn.close()

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Create compatibility views for CryoProtect v2")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without committing changes")
    args = parser.parse_args()
    
    report = create_all_compatibility_views(dry_run=args.dry_run)
    
    if report["success"]:
        logger.info("Successfully created compatibility views")
        if args.dry_run:
            logger.info("Dry run completed successfully, no changes were committed")
    else:
        logger.error("Failed to create some compatibility views")
        if "message" in report:
            logger.error(f"Error: {report['message']}")
    
    return 0 if report["success"] else 1

if __name__ == "__main__":
    sys.exit(main())