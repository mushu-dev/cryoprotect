#!/usr/bin/env python3
"""
CryoProtect v2 Database Migration Script

This script safely migrates data from singular-named tables to plural-named tables
using the "expand and contract" approach. It creates backup tables, validates data integrity,
and supports rollback in case of failure.

Tables to migrate:
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
import time
import json
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
        logging.FileHandler(f"migration_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
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

# Migration configuration
TABLES_TO_MIGRATE = [
    {"singular": "molecule", "plural": "molecules"},
    {"singular": "mixture", "plural": "mixtures"},
    {"singular": "experiment", "plural": "experiments"},
    {"singular": "prediction", "plural": "predictions"},
    {"singular": "project", "plural": "projects"}
]

# Foreign key relationships to update
FK_RELATIONSHIPS = [
    {
        "table": "mixture_component",
        "columns": [
            {"name": "mixture_id", "references": "mixture", "new_references": "mixtures"},
            {"name": "molecule_id", "references": "molecule", "new_references": "molecules"}
        ]
    },
    {
        "table": "molecular_property",
        "columns": [
            {"name": "molecule_id", "references": "molecule", "new_references": "molecules"}
        ]
    },
    {
        "table": "prediction",
        "columns": [
            {"name": "molecule_id", "references": "molecule", "new_references": "molecules"},
            {"name": "mixture_id", "references": "mixture", "new_references": "mixtures"}
        ]
    },
    {
        "table": "experiment",
        "columns": [
            {"name": "molecule_id", "references": "molecule", "new_references": "molecules"},
            {"name": "mixture_id", "references": "mixture", "new_references": "mixtures"}
        ]
    },
    {
        "table": "project_membership",
        "columns": [
            {"name": "project_id", "references": "project", "new_references": "projects"}
        ]
    }
]

# Views that need to be updated
VIEWS_TO_UPDATE = [
    "experiment_with_results",
    "mixture_with_components",
    "mixtures_with_components",
    "molecule_with_properties",
    "molecules_with_properties"
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

def get_table_info(cursor, table_name):
    """
    Get table schema information including columns and constraints.
    
    Args:
        cursor: Database cursor
        table_name (str): Name of the table
        
    Returns:
        dict: Table information including columns and constraints
    """
    # Get column information
    cursor.execute("""
        SELECT column_name, data_type, is_nullable, column_default
        FROM information_schema.columns
        WHERE table_name = %s
        ORDER BY ordinal_position
    """, (table_name,))
    columns = cursor.fetchall()
    
    # Get primary key information
    cursor.execute("""
        SELECT c.column_name
        FROM information_schema.table_constraints tc
        JOIN information_schema.constraint_column_usage AS ccu USING (constraint_schema, constraint_name)
        JOIN information_schema.columns AS c ON c.table_schema = tc.constraint_schema
            AND tc.table_name = c.table_name AND ccu.column_name = c.column_name
        WHERE tc.constraint_type = 'PRIMARY KEY' AND tc.table_name = %s
    """, (table_name,))
    primary_keys = [row['column_name'] for row in cursor.fetchall()]
    
    # Get foreign key information
    cursor.execute("""
        SELECT
            kcu.column_name,
            ccu.table_name AS foreign_table_name,
            ccu.column_name AS foreign_column_name,
            tc.constraint_name
        FROM information_schema.table_constraints AS tc
        JOIN information_schema.key_column_usage AS kcu
            ON tc.constraint_name = kcu.constraint_name
            AND tc.table_schema = kcu.table_schema
        JOIN information_schema.constraint_column_usage AS ccu
            ON ccu.constraint_name = tc.constraint_name
            AND ccu.table_schema = tc.table_schema
        WHERE tc.constraint_type = 'FOREIGN KEY' AND tc.table_name = %s
    """, (table_name,))
    foreign_keys = cursor.fetchall()
    
    # Get unique constraints
    cursor.execute("""
        SELECT
            tc.constraint_name,
            kcu.column_name
        FROM information_schema.table_constraints tc
        JOIN information_schema.key_column_usage kcu
            ON tc.constraint_name = kcu.constraint_name
        WHERE tc.constraint_type = 'UNIQUE'
            AND tc.table_name = %s
        ORDER BY tc.constraint_name, kcu.ordinal_position
    """, (table_name,))
    unique_constraints = {}
    for row in cursor.fetchall():
        if row['constraint_name'] not in unique_constraints:
            unique_constraints[row['constraint_name']] = []
        unique_constraints[row['constraint_name']].append(row['column_name'])
    
    # Get indexes
    cursor.execute("""
        SELECT
            i.relname AS index_name,
            a.attname AS column_name,
            ix.indisunique AS is_unique
        FROM
            pg_class t,
            pg_class i,
            pg_index ix,
            pg_attribute a
        WHERE
            t.oid = ix.indrelid
            AND i.oid = ix.indexrelid
            AND a.attrelid = t.oid
            AND a.attnum = ANY(ix.indkey)
            AND t.relkind = 'r'
            AND t.relname = %s
        ORDER BY
            i.relname, a.attnum
    """, (table_name,))
    indexes = {}
    for row in cursor.fetchall():
        if row['index_name'] not in indexes:
            indexes[row['index_name']] = {
                'columns': [],
                'is_unique': row['is_unique']
            }
        indexes[row['index_name']]['columns'].append(row['column_name'])
    
    return {
        'columns': columns,
        'primary_keys': primary_keys,
        'foreign_keys': foreign_keys,
        'unique_constraints': unique_constraints,
        'indexes': indexes
    }

def get_view_definition(cursor, view_name):
    """
    Get the SQL definition of a view.
    
    Args:
        cursor: Database cursor
        view_name (str): Name of the view
        
    Returns:
        str: SQL definition of the view
    """
    cursor.execute("""
        SELECT pg_get_viewdef(c.oid, true) as view_definition
        FROM pg_class c
        JOIN pg_namespace n ON n.oid = c.relnamespace
        WHERE c.relkind = 'v'
        AND n.nspname = 'public'
        AND c.relname = %s
    """, (view_name,))
    result = cursor.fetchone()
    return result['view_definition'] if result else None

def create_backup_table(cursor, table_name):
    """
    Create a backup of a table.
    
    Args:
        cursor: Database cursor
        table_name (str): Name of the table to backup
        
    Returns:
        str: Name of the backup table
    """
    backup_table = f"{table_name}_backup_{datetime.now().strftime('%Y%m%d')}"
    logger.info(f"Creating backup table: {backup_table}")
    
    cursor.execute(
        sql.SQL("CREATE TABLE {} AS TABLE {}").format(
            sql.Identifier(backup_table),
            sql.Identifier(table_name)
        )
    )
    
    # Copy constraints and indexes (simplified approach)
    table_info = get_table_info(cursor, table_name)
    
    # Create primary key
    if table_info['primary_keys']:
        pk_columns = ', '.join(table_info['primary_keys'])
        cursor.execute(
            sql.SQL("ALTER TABLE {} ADD PRIMARY KEY ({})").format(
                sql.Identifier(backup_table),
                sql.SQL(pk_columns)
            )
        )
    
    return backup_table

def create_plural_table(cursor, singular_table, plural_table):
    """
    Create a new plural-named table with the same structure as the singular table.
    
    Args:
        cursor: Database cursor
        singular_table (str): Name of the source table
        plural_table (str): Name of the target table
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        logger.info(f"Creating plural table: {plural_table}")
        
        # Get table information
        table_info = get_table_info(cursor, singular_table)
        
        # Create table with the same structure
        columns_sql = []
        for col in table_info['columns']:
            nullable = "NULL" if col['is_nullable'] == 'YES' else "NOT NULL"
            default = f"DEFAULT {col['column_default']}" if col['column_default'] else ""
            columns_sql.append(
                sql.SQL("{} {} {} {}").format(
                    sql.Identifier(col['column_name']),
                    sql.SQL(col['data_type']),
                    sql.SQL(nullable),
                    sql.SQL(default)
                )
            )
        
        create_table_sql = sql.SQL("CREATE TABLE {} ({})").format(
            sql.Identifier(plural_table),
            sql.SQL(", ").join(columns_sql)
        )
        cursor.execute(create_table_sql)
        
        # Add primary key
        if table_info['primary_keys']:
            pk_columns = [sql.Identifier(col) for col in table_info['primary_keys']]
            cursor.execute(
                sql.SQL("ALTER TABLE {} ADD PRIMARY KEY ({})").format(
                    sql.Identifier(plural_table),
                    sql.SQL(", ").join(pk_columns)
                )
            )
        
        # Add unique constraints
        for constraint_name, columns in table_info['unique_constraints'].items():
            # Create a new constraint name for the plural table
            new_constraint_name = constraint_name.replace(singular_table, plural_table)
            col_identifiers = [sql.Identifier(col) for col in columns]
            cursor.execute(
                sql.SQL("ALTER TABLE {} ADD CONSTRAINT {} UNIQUE ({})").format(
                    sql.Identifier(plural_table),
                    sql.Identifier(new_constraint_name),
                    sql.SQL(", ").join(col_identifiers)
                )
            )
        
        # Add indexes (excluding those created for constraints)
        for index_name, index_info in table_info['indexes'].items():
            # Skip indexes that are part of constraints
            if any(index_name.endswith(f"_{constraint}") for constraint in table_info['unique_constraints']):
                continue
            
            # Create a new index name for the plural table
            new_index_name = index_name.replace(singular_table, plural_table)
            col_identifiers = [sql.Identifier(col) for col in index_info['columns']]
            
            index_sql = sql.SQL("CREATE {} INDEX {} ON {} ({})").format(
                sql.SQL("UNIQUE" if index_info['is_unique'] else ""),
                sql.Identifier(new_index_name),
                sql.Identifier(plural_table),
                sql.SQL(", ").join(col_identifiers)
            )
            cursor.execute(index_sql)
        
        return True
    except Exception as e:
        logger.error(f"Error creating plural table {plural_table}: {str(e)}")
        return False

def copy_data(cursor, source_table, target_table):
    """
    Copy data from source table to target table.
    
    Args:
        cursor: Database cursor
        source_table (str): Name of the source table
        target_table (str): Name of the target table
        
    Returns:
        int: Number of rows copied
    """
    try:
        logger.info(f"Copying data from {source_table} to {target_table}")
        
        # Get column names
        cursor.execute("""
            SELECT column_name
            FROM information_schema.columns
            WHERE table_name = %s
            ORDER BY ordinal_position
        """, (source_table,))
        columns = [row['column_name'] for row in cursor.fetchall()]
        
        # Prepare column identifiers
        column_identifiers = [sql.Identifier(col) for col in columns]
        
        # Copy data
        cursor.execute(
            sql.SQL("INSERT INTO {} ({}) SELECT {} FROM {}").format(
                sql.Identifier(target_table),
                sql.SQL(", ").join(column_identifiers),
                sql.SQL(", ").join(column_identifiers),
                sql.Identifier(source_table)
            )
        )
        
        # Get number of rows copied
        cursor.execute(
            sql.SQL("SELECT COUNT(*) as count FROM {}").format(
                sql.Identifier(target_table)
            )
        )
        result = cursor.fetchone()
        row_count = result['count'] if result else 0
        
        logger.info(f"Copied {row_count} rows to {target_table}")
        return row_count
    except Exception as e:
        logger.error(f"Error copying data to {target_table}: {str(e)}")
        return 0

def update_foreign_keys(cursor, relationships):
    """
    Update foreign key relationships to reference the new plural tables.
    
    Args:
        cursor: Database cursor
        relationships (list): List of foreign key relationships to update
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        logger.info("Updating foreign key relationships")
        
        for rel in relationships:
            table = rel['table']
            for column in rel['columns']:
                col_name = column['name']
                old_ref = column['references']
                new_ref = column['new_references']
                
                # Get constraint name
                cursor.execute("""
                    SELECT tc.constraint_name
                    FROM information_schema.table_constraints AS tc
                    JOIN information_schema.key_column_usage AS kcu
                        ON tc.constraint_name = kcu.constraint_name
                    JOIN information_schema.constraint_column_usage AS ccu
                        ON ccu.constraint_name = tc.constraint_name
                    WHERE tc.constraint_type = 'FOREIGN KEY'
                        AND tc.table_name = %s
                        AND kcu.column_name = %s
                        AND ccu.table_name = %s
                """, (table, col_name, old_ref))
                
                constraint = cursor.fetchone()
                if constraint:
                    constraint_name = constraint['constraint_name']
                    
                    # Drop the old constraint
                    cursor.execute(
                        sql.SQL("ALTER TABLE {} DROP CONSTRAINT {}").format(
                            sql.Identifier(table),
                            sql.Identifier(constraint_name)
                        )
                    )
                    
                    # Create new constraint
                    cursor.execute(
                        sql.SQL("ALTER TABLE {} ADD CONSTRAINT {} FOREIGN KEY ({}) REFERENCES {} (id)").format(
                            sql.Identifier(table),
                            sql.Identifier(f"{table}_{col_name}_fkey_new"),
                            sql.Identifier(col_name),
                            sql.Identifier(new_ref)
                        )
                    )
                    
                    logger.info(f"Updated foreign key {table}.{col_name} to reference {new_ref}")
        
        return True
    except Exception as e:
        logger.error(f"Error updating foreign keys: {str(e)}")
        return False

def update_views(cursor, views):
    """
    Update views to reference the new plural tables.
    
    Args:
        cursor: Database cursor
        views (list): List of view names to update
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        logger.info("Updating views")
        
        for view_name in views:
            # Get view definition
            view_def = get_view_definition(cursor, view_name)
            if not view_def:
                logger.warning(f"View {view_name} not found, skipping")
                continue
            
            # Create a backup of the view definition
            backup_file = f"view_{view_name}_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}.sql"
            with open(backup_file, 'w') as f:
                f.write(view_def)
            
            # Update the view definition to use plural table names
            new_view_def = view_def
            for table_info in TABLES_TO_MIGRATE:
                singular = table_info['singular']
                plural = table_info['plural']
                # Replace table names but be careful with column names
                new_view_def = new_view_def.replace(f" {singular} ", f" {plural} ")
                new_view_def = new_view_def.replace(f" {singular},", f" {plural},")
                new_view_def = new_view_def.replace(f" {singular}.", f" {plural}.")
                new_view_def = new_view_def.replace(f"FROM {singular}", f"FROM {plural}")
                new_view_def = new_view_def.replace(f"JOIN {singular}", f"JOIN {plural}")
            
            # Drop and recreate the view
            cursor.execute(
                sql.SQL("DROP VIEW IF EXISTS {}").format(
                    sql.Identifier(view_name)
                )
            )
            
            cursor.execute(new_view_def)
            logger.info(f"Updated view: {view_name}")
        
        return True
    except Exception as e:
        logger.error(f"Error updating views: {str(e)}")
        return False

def validate_data_integrity(cursor, singular_table, plural_table):
    """
    Validate that data was correctly migrated between tables.
    
    Args:
        cursor: Database cursor
        singular_table (str): Name of the source table
        plural_table (str): Name of the target table
        
    Returns:
        bool: True if validation passes, False otherwise
    """
    try:
        logger.info(f"Validating data integrity between {singular_table} and {plural_table}")
        
        # Count rows in both tables
        cursor.execute(
            sql.SQL("SELECT COUNT(*) as count FROM {}").format(
                sql.Identifier(singular_table)
            )
        )
        singular_count = cursor.fetchone()['count']
        
        cursor.execute(
            sql.SQL("SELECT COUNT(*) as count FROM {}").format(
                sql.Identifier(plural_table)
            )
        )
        plural_count = cursor.fetchone()['count']
        
        if singular_count != plural_count:
            logger.error(f"Row count mismatch: {singular_table}={singular_count}, {plural_table}={plural_count}")
            return False
        
        # Get primary key column
        cursor.execute("""
            SELECT c.column_name
            FROM information_schema.table_constraints tc
            JOIN information_schema.constraint_column_usage AS ccu USING (constraint_schema, constraint_name)
            JOIN information_schema.columns AS c ON c.table_schema = tc.constraint_schema
                AND tc.table_name = c.table_name AND ccu.column_name = c.column_name
            WHERE tc.constraint_type = 'PRIMARY KEY' AND tc.table_name = %s
        """, (singular_table,))
        pk_columns = [row['column_name'] for row in cursor.fetchall()]
        
        if not pk_columns:
            logger.warning(f"No primary key found for {singular_table}, using all columns for validation")
            # Get all columns
            cursor.execute("""
                SELECT column_name
                FROM information_schema.columns
                WHERE table_name = %s
            """, (singular_table,))
            columns = [row['column_name'] for row in cursor.fetchall()]
        else:
            columns = pk_columns
        
        # Compare a sample of rows using primary key or all columns
        column_identifiers = [sql.Identifier(col) for col in columns]
        
        # Get a sample of IDs from the singular table
        cursor.execute(
            sql.SQL("SELECT {} FROM {} LIMIT 100").format(
                sql.SQL(", ").join(column_identifiers),
                sql.Identifier(singular_table)
            )
        )
        sample_rows = cursor.fetchall()
        
        # For each sample row, check if it exists in the plural table
        for row in sample_rows:
            where_conditions = []
            for col in columns:
                where_conditions.append(
                    sql.SQL("{} = %s").format(sql.Identifier(col))
                )
            
            cursor.execute(
                sql.SQL("SELECT COUNT(*) as count FROM {} WHERE {}").format(
                    sql.Identifier(plural_table),
                    sql.SQL(" AND ").join(where_conditions)
                ),
                [row[col] for col in columns]
            )
            
            match_count = cursor.fetchone()['count']
            if match_count != 1:
                logger.error(f"Data integrity check failed for {singular_table} row: {row}")
                return False
        
        logger.info(f"Data integrity validation passed for {singular_table} → {plural_table}")
        return True
    except Exception as e:
        logger.error(f"Error validating data integrity: {str(e)}")
        return False

def perform_migration(conn, cursor, dry_run=False):
    """
    Perform the migration from singular to plural table names.
    
    Args:
        conn: Database connection
        cursor: Database cursor
        dry_run (bool): If True, don't commit changes
        
    Returns:
        bool: True if migration was successful, False otherwise
    """
    try:
        logger.info(f"Starting migration (dry_run={dry_run})")
        
        # Create backup tables
        backup_tables = {}
        for table_info in TABLES_TO_MIGRATE:
            singular = table_info['singular']
            backup_tables[singular] = create_backup_table(cursor, singular)
        
        # Create plural tables and copy data
        for table_info in TABLES_TO_MIGRATE:
            singular = table_info['singular']
            plural = table_info['plural']
            
            # Check if plural table already exists
            cursor.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_name = %s
                ) as exists
            """, (plural,))
            table_exists = cursor.fetchone()['exists']
            
            if table_exists:
                logger.info(f"Plural table {plural} already exists")
                
                # Check if it has data
                cursor.execute(
                    sql.SQL("SELECT COUNT(*) as count FROM {}").format(
                        sql.Identifier(plural)
                    )
                )
                count = cursor.fetchone()['count']
                
                if count > 0:
                    logger.warning(f"Plural table {plural} already has {count} rows, skipping data copy")
                    continue
            else:
                # Create the plural table
                if not create_plural_table(cursor, singular, plural):
                    logger.error(f"Failed to create plural table {plural}")
                    if not dry_run:
                        conn.rollback()
                    return False
            
            # Copy data from singular to plural
            row_count = copy_data(cursor, singular, plural)
            if row_count == 0:
                logger.error(f"Failed to copy data to {plural}")
                if not dry_run:
                    conn.rollback()
                return False
            
            # Validate data integrity
            if not validate_data_integrity(cursor, singular, plural):
                logger.error(f"Data integrity validation failed for {singular} → {plural}")
                if not dry_run:
                    conn.rollback()
                return False
        
        # Update foreign key relationships
        if not update_foreign_keys(cursor, FK_RELATIONSHIPS):
            logger.error("Failed to update foreign key relationships")
            if not dry_run:
                conn.rollback()
            return False
        
        # Update views
        if not update_views(cursor, VIEWS_TO_UPDATE):
            logger.warning("Some views could not be updated")
            # Continue with migration, as view updates are not critical
        
        # If this is a dry run, roll back all changes
        if dry_run:
            logger.info("Dry run completed, rolling back changes")
            conn.rollback()
            return True
        
        # Commit the transaction
        conn.commit()
        logger.info("Migration completed successfully")
        
        # Create migration report
        migration_report = {
            "timestamp": datetime.now().isoformat(),
            "tables_migrated": TABLES_TO_MIGRATE,
            "backup_tables": backup_tables,
            "foreign_keys_updated": FK_RELATIONSHIPS,
            "views_updated": VIEWS_TO_UPDATE
        }
        
        with open(f"migration_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json", 'w') as f:
            json.dump(migration_report, f, indent=2)
        
        return True
    except Exception as e:
        logger.error(f"Error during migration: {str(e)}")
        if not dry_run:
            conn.rollback()
        return False

def rollback_migration(conn, cursor, backup_tables):
    """
    Rollback the migration using backup tables.
    
    Args:
        conn: Database connection
        cursor: Database cursor
        backup_tables (dict): Dictionary mapping original table names to backup table names
        
    Returns:
        bool: True if rollback was successful, False otherwise
    """
    try:
        logger.info("Starting migration rollback")
        
        # Restore data from backup tables
        for original_table, backup_table in backup_tables.items():
            # Truncate the original table
            cursor.execute(
                sql.SQL("TRUNCATE TABLE {} CASCADE").format(
                    sql.Identifier(original_table)
                )
            )
            
            # Copy data from backup
            row_count = copy_data(cursor, backup_table, original_table)
            if row_count == 0:
                logger.error(f"Failed to restore data from {backup_table} to {original_table}")
                conn.rollback()
                return False
            
            logger.info(f"Restored {row_count} rows from {backup_table} to {original_table}")
        
        # Commit the transaction
        conn.commit()
        logger.info("Rollback completed successfully")
        return True
    except Exception as e:
        logger.error(f"Error during rollback: {str(e)}")
        conn.rollback()
        return False

def main():
    """Main function to run the migration script."""
    parser = argparse.ArgumentParser(description='CryoProtect v2 Database Migration Script')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run without committing changes')
    parser.add_argument('--rollback', action='store_true', help='Rollback migration using backup tables')
    parser.add_argument('--report', type=str, help='Path to migration report for rollback')
    args = parser.parse_args()
    
    # Connect to the database
    conn, cursor = connect_to_db()
    if not conn or not cursor:
        sys.exit(1)
    
    try:
        if args.rollback:
            if not args.report:
                logger.error("Migration report is required for rollback")
                sys.exit(1)
            
            # Load migration report
            with open(args.report, 'r') as f:
                report = json.load(f)
            
            # Perform rollback
            if rollback_migration(conn, cursor, report['backup_tables']):
                logger.info("Migration rollback completed successfully")
            else:
                logger.error("Migration rollback failed")
                sys.exit(1)
        else:
            # Perform migration
            if perform_migration(conn, cursor, dry_run=args.dry_run):
                if args.dry_run:
                    logger.info("Dry run completed successfully")
                else:
                    logger.info("Migration completed successfully")
            else:
                logger.error("Migration failed")
                sys.exit(1)
    finally:
        # Close database connection
        if conn:
            conn.close()

if __name__ == "__main__":
    main()