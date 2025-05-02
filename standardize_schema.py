#!/usr/bin/env python3
"""
standardize_schema.py - Standardize the CryoProtect Supabase database schema

This script:
1. Converts all singular-named tables to plural
2. Ensures all tables use UUID as primary key with DEFAULT gen_random_uuid()
3. Adds proper REFERENCES constraints with indexes for all foreign keys
4. Includes rollback mechanism in case of failure

Usage:
    python standardize_schema.py
"""

import os
import sys
import json
import time
import logging
import tempfile
import traceback
from datetime import datetime
from pathlib import Path

# Create logs directory if it doesn't exist
logs_dir = Path("logs")
logs_dir.mkdir(exist_ok=True)

# Configure logging
log_file = logs_dir / f"schema_migration_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Table mapping (singular to plural)
TABLE_MAPPING = {
    "molecule": "molecules",
    "mixture": "mixtures",
    "prediction": "predictions",
    "experiment": "experiments",
    "experiment_property": "experiment_properties",
    "mixture_component": "mixture_components",
    "calculation_method": "calculation_methods",
    "property_type": "property_types",  # This is already plural in the DB
    "project": "projects",
    "team": "teams"
}

# Project ID
PROJECT_ID = "tsdlmynydfuypiugmkev"

# SQL templates
CREATE_TABLE_TEMPLATE = """
CREATE TABLE IF NOT EXISTS public.{plural_name} (
    {columns}
);
"""

COPY_DATA_TEMPLATE = """
INSERT INTO public.{plural_name} 
SELECT * FROM public.{singular_name};
"""

ADD_INDEX_TEMPLATE = """
CREATE INDEX IF NOT EXISTS idx_{plural_name}_{column_name} ON public.{plural_name} ({column_name});
"""

ADD_FK_CONSTRAINT_TEMPLATE = """
ALTER TABLE public.{plural_name} 
ADD CONSTRAINT fk_{plural_name}_{column_name} 
FOREIGN KEY ({column_name}) 
REFERENCES public.{referenced_table} ({referenced_column})
ON DELETE CASCADE;
"""

# Create a backup directory for database backups
def create_backup_directory():
    """Create a backup directory for database backups"""
    backup_dir = Path("database_backups")
    try:
        backup_dir.mkdir(exist_ok=True)
        logger.info(f"Created backup directory: {backup_dir}")
        return backup_dir
    except Exception as e:
        logger.error(f"Failed to create backup directory: {str(e)}")
        return None

# Function to execute SQL with Supabase MCP
def execute_sql(query, description=None, in_transaction=False):
    """Execute SQL query using Supabase MCP"""
    if description:
        logger.info(f"Executing: {description}")
    
    try:
        # Format the query for logging (remove extra whitespace)
        formatted_query = " ".join(line.strip() for line in query.split("\n") if line.strip())
        logger.info(f"SQL: {formatted_query[:100]}...")
        
        # Use Supabase MCP to execute the query
        import subprocess
        import json
        
        # Create a temporary file with the query
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sql', delete=False) as temp_file:
            temp_path = temp_file.name
            temp_file.write(query)
        
        try:
            # Add transaction control if requested
            if in_transaction and not query.strip().upper().startswith(("BEGIN", "COMMIT", "ROLLBACK")):
                transaction_query = f"BEGIN;\n{query}\nCOMMIT;"
            else:
                transaction_query = query
                
            # Execute the query using Supabase MCP
            result = subprocess.run(
                [
                    "python", "-c",
                    f"""
import json
import sys
try:
    from supabase_mcp_tools import execute_sql_on_supabase
    result = execute_sql_on_supabase('{PROJECT_ID}', open(r'{temp_path}').read())
    print(json.dumps(result))
except Exception as e:
    error_details = {{'error': str(e), 'traceback': traceback.format_exc()}}
    print(json.dumps(error_details))
                    """
                ],
                capture_output=True,
                text=True
            )
            
            # Parse the result
            try:
                result_json = json.loads(result.stdout)
                if "error" in result_json:
                    error_msg = result_json['error']
                    if 'traceback' in result_json:
                        logger.error(f"Error executing query: {error_msg}\nTraceback: {result_json['traceback']}")
                    else:
                        logger.error(f"Error executing query: {error_msg}")
                    return False, error_msg
                return True, result_json
            except json.JSONDecodeError:
                logger.error(f"Error parsing result: {result.stdout}")
                return False, result.stdout
        finally:
            # Clean up the temporary file
            try:
                os.unlink(temp_path)
            except Exception as e:
                logger.warning(f"Failed to remove temporary file {temp_path}: {str(e)}")
    
    except Exception as e:
        error_details = traceback.format_exc()
        logger.error(f"Error executing query: {str(e)}\n{error_details}")
        return False, str(e)

def get_table_schema(table_name):
    """Get the schema of a table"""
    query = f"""
    SELECT column_name, data_type, column_default, is_nullable
    FROM information_schema.columns
    WHERE table_schema = 'public' AND table_name = '{table_name}'
    ORDER BY ordinal_position;
    """
    
    success, result = execute_sql(query, f"Getting schema for {table_name}")
    if not success:
        logger.error(f"Failed to get schema for {table_name}: {result}")
        return None
    
    return result

def get_foreign_keys(table_name):
    """Get foreign keys for a table"""
    query = f"""
    SELECT
        tc.constraint_name,
        kcu.column_name,
        ccu.table_name AS referenced_table,
        ccu.column_name AS referenced_column
    FROM
        information_schema.table_constraints AS tc
        JOIN information_schema.key_column_usage AS kcu
          ON tc.constraint_name = kcu.constraint_name
          AND tc.table_schema = kcu.table_schema
        JOIN information_schema.constraint_column_usage AS ccu
          ON ccu.constraint_name = tc.constraint_name
          AND ccu.table_schema = tc.table_schema
    WHERE tc.constraint_type = 'FOREIGN KEY'
        AND tc.table_schema = 'public'
        AND tc.table_name = '{table_name}';
    """
    
    success, result = execute_sql(query, f"Getting foreign keys for {table_name}")
    if not success:
        logger.error(f"Failed to get foreign keys for {table_name}: {result}")
        return None
    
    return result

def table_exists(table_name):
    """Check if a table exists"""
    query = f"""
    SELECT EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' 
        AND table_name = '{table_name}'
    );
    """
    
    success, result = execute_sql(query, f"Checking if {table_name} exists")
    if not success:
        logger.error(f"Failed to check if {table_name} exists: {result}")
        return False
    
    return result[0]["exists"]

def create_migration_tables():
    """Create tables to track migration status"""
    query = """
    CREATE TABLE IF NOT EXISTS public.schema_migrations (
        id SERIAL PRIMARY KEY,
        migration_name TEXT NOT NULL,
        applied_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
        status TEXT NOT NULL,
        details JSONB
    );
    
    CREATE TABLE IF NOT EXISTS public.schema_migration_operations (
        id SERIAL PRIMARY KEY,
        migration_id INTEGER REFERENCES public.schema_migrations(id),
        operation_type TEXT NOT NULL,
        table_name TEXT NOT NULL,
        sql_up TEXT NOT NULL,
        sql_down TEXT NOT NULL,
        applied BOOLEAN DEFAULT FALSE,
        applied_at TIMESTAMP WITH TIME ZONE,
        status TEXT,
        error_message TEXT
    );
    """
    
    success, result = execute_sql(query, "Creating migration tracking tables", in_transaction=True)
    if not success:
        logger.error(f"Failed to create migration tracking tables: {result}")
        return False
    
    return True

def start_migration():
    """Start a new migration and return its ID"""
    query = """
    INSERT INTO public.schema_migrations
    (migration_name, status, details)
    VALUES
    ('standardize_schema', 'started', '{"tables": []}')
    RETURNING id;
    """
    
    success, result = execute_sql(query, "Starting migration", in_transaction=True)
    if not success:
        logger.error(f"Failed to start migration: {result}")
        return None
    
    return result[0]["id"]

def record_operation(migration_id, operation_type, table_name, sql_up, sql_down):
    """Record a migration operation"""
    query = f"""
    INSERT INTO public.schema_migration_operations 
    (migration_id, operation_type, table_name, sql_up, sql_down) 
    VALUES 
    ({migration_id}, '{operation_type}', '{table_name}', $${sql_up}$$, $${sql_down}$$)
    RETURNING id;
    """
    
    success, result = execute_sql(query, f"Recording {operation_type} operation for {table_name}")
    if not success:
        logger.error(f"Failed to record operation: {result}")
        return None
    
    return result[0]["id"]

def mark_operation_applied(operation_id, status="success", error_message=None):
    """Mark a migration operation as applied"""
    error_clause = f", error_message = '{error_message}'" if error_message else ""
    query = f"""
    UPDATE public.schema_migration_operations 
    SET applied = TRUE, 
        applied_at = NOW(), 
        status = '{status}'{error_clause}
    WHERE id = {operation_id};
    """
    
    success, result = execute_sql(query, f"Marking operation {operation_id} as applied")
    if not success:
        logger.error(f"Failed to mark operation as applied: {result}")
        return False
    
    return True

def update_migration_status(migration_id, status, details=None):
    """Update the status of a migration"""
    details_clause = f", details = '{json.dumps(details)}'" if details else ""
    query = f"""
    UPDATE public.schema_migrations 
    SET status = '{status}'{details_clause}
    WHERE id = {migration_id};
    """
    
    success, result = execute_sql(query, f"Updating migration {migration_id} status to {status}")
    if not success:
        logger.error(f"Failed to update migration status: {result}")
        return False
    
    return True

def create_plural_table(migration_id, singular_name, plural_name):
    """Create a new plural table based on the singular table schema"""
    # Get the schema of the singular table
    schema = get_table_schema(singular_name)
    if not schema:
        return False
    
    # Build the column definitions
    columns = []
    for column in schema:
        name = column["column_name"]
        data_type = column["data_type"]
        default = f"DEFAULT {column['column_default']}" if column["column_default"] else ""
        nullable = "NOT NULL" if column["is_nullable"] == "NO" else ""
        
        # Ensure primary key has DEFAULT gen_random_uuid()
        if name == "id" and data_type == "uuid" and not column["column_default"]:
            default = "DEFAULT gen_random_uuid()"
        
        columns.append(f"{name} {data_type} {default} {nullable}")
    
    # Create the table
    create_sql = CREATE_TABLE_TEMPLATE.format(
        plural_name=plural_name,
        columns=",\n    ".join(columns)
    )
    
    # Add primary key constraint
    create_sql += f"""
    ALTER TABLE public.{plural_name} ADD PRIMARY KEY (id);
    """
    
    # Record the operation
    drop_sql = f"DROP TABLE IF EXISTS public.{plural_name};"
    operation_id = record_operation(
        migration_id,
        "create_table",
        plural_name,
        create_sql,
        drop_sql
    )
    
    if not operation_id:
        return False
    
    # Execute the SQL with transaction support
    success, result = execute_sql(create_sql, f"Creating {plural_name} table", in_transaction=True)
    
    # Mark the operation as applied
    if success:
        mark_operation_applied(operation_id)
    else:
        mark_operation_applied(operation_id, "error", str(result))
    
    return success

def copy_data(migration_id, singular_name, plural_name):
    """Copy data from singular table to plural table"""
    copy_sql = COPY_DATA_TEMPLATE.format(
        singular_name=singular_name,
        plural_name=plural_name
    )
    
    # Record the operation
    truncate_sql = f"TRUNCATE TABLE public.{plural_name};"
    operation_id = record_operation(
        migration_id,
        "copy_data",
        plural_name,
        copy_sql,
        truncate_sql
    )
    
    if not operation_id:
        return False
    
    # Execute the SQL with transaction support
    success, result = execute_sql(copy_sql, f"Copying data from {singular_name} to {plural_name}", in_transaction=True)
    
    # Mark the operation as applied
    if success:
        mark_operation_applied(operation_id)
    else:
        mark_operation_applied(operation_id, "error", str(result))
    
    return success

def add_foreign_keys(migration_id, singular_name, plural_name):
    """Add foreign key constraints to the plural table"""
    # Get foreign keys from the singular table
    foreign_keys = get_foreign_keys(singular_name)
    if not foreign_keys:
        logger.info(f"No foreign keys found for {singular_name}")
        return True
    
    success = True
    
    for fk in foreign_keys:
        column_name = fk["column_name"]
        referenced_table = fk["referenced_table"]
        referenced_column = fk["referenced_column"]
        
        # If the referenced table is in our mapping, use the plural name
        if referenced_table in TABLE_MAPPING:
            referenced_table = TABLE_MAPPING[referenced_table]
        
        # Create index for the foreign key
        index_sql = ADD_INDEX_TEMPLATE.format(
            plural_name=plural_name,
            column_name=column_name
        )
        
        # Add foreign key constraint
        fk_sql = ADD_FK_CONSTRAINT_TEMPLATE.format(
            plural_name=plural_name,
            column_name=column_name,
            referenced_table=referenced_table,
            referenced_column=referenced_column
        )
        
        combined_sql = index_sql + "\n" + fk_sql
        
        # Record the operation
        drop_fk_sql = f"""
        ALTER TABLE public.{plural_name}
        DROP CONSTRAINT IF EXISTS fk_{plural_name}_{column_name};
        
        DROP INDEX IF EXISTS idx_{plural_name}_{column_name};
        """
        
        operation_id = record_operation(
            migration_id,
            "add_foreign_key",
            f"{plural_name}.{column_name}",
            combined_sql,
            drop_fk_sql
        )
        
        if not operation_id:
            success = False
            continue
        
        # Execute the SQL with transaction support
        op_success, result = execute_sql(combined_sql, f"Adding foreign key for {plural_name}.{column_name}", in_transaction=True)
        
        # Mark the operation as applied
        if op_success:
            mark_operation_applied(operation_id)
        else:
            mark_operation_applied(operation_id, "error", str(result))
            success = False
    
    return success

def rollback_migration(migration_id):
    """Rollback a migration by executing the down SQL for each operation in reverse order"""
    query = f"""
    SELECT id, operation_type, table_name, sql_down
    FROM public.schema_migration_operations
    WHERE migration_id = {migration_id} AND applied = TRUE
    ORDER BY id DESC;
    """
    
    success, operations = execute_sql(query, "Getting operations to rollback")
    if not success:
        logger.error(f"Failed to get operations to rollback: {operations}")
        return False
    
    rollback_success = True
    
    # Begin transaction for rollback
    begin_transaction = "BEGIN;"
    success, result = execute_sql(begin_transaction, "Beginning rollback transaction")
    if not success:
        logger.error(f"Failed to begin rollback transaction: {result}")
        return False
    
    try:
        for op in operations:
            logger.info(f"Rolling back {op['operation_type']} for {op['table_name']}")
            
            op_success, result = execute_sql(op["sql_down"], f"Rolling back operation {op['id']}")
            if not op_success:
                logger.error(f"Failed to rollback operation {op['id']}: {result}")
                rollback_success = False
                raise Exception(f"Rollback failed for operation {op['id']}: {result}")
            
            # Mark the operation as not applied
            update_query = f"""
            UPDATE public.schema_migration_operations
            SET applied = FALSE, applied_at = NULL, status = 'rolled_back'
            WHERE id = {op['id']};
            """
            
            update_success, update_result = execute_sql(update_query, f"Marking operation {op['id']} as rolled back")
            if not update_success:
                logger.error(f"Failed to mark operation {op['id']} as rolled back: {update_result}")
                rollback_success = False
                raise Exception(f"Failed to mark operation {op['id']} as rolled back: {update_result}")
        
        # Commit transaction if all operations succeeded
        commit_transaction = "COMMIT;"
        success, result = execute_sql(commit_transaction, "Committing rollback transaction")
        if not success:
            logger.error(f"Failed to commit rollback transaction: {result}")
            rollback_success = False
            raise Exception(f"Failed to commit rollback transaction: {result}")
        
    except Exception as e:
        # Rollback transaction if any operation failed
        logger.error(f"Error during rollback, rolling back transaction: {str(e)}")
        rollback_transaction = "ROLLBACK;"
        execute_sql(rollback_transaction, "Rolling back transaction due to error")
        rollback_success = False
    
    # Update migration status
    update_migration_status(migration_id, "rolled_back")
    
    return rollback_success

def create_database_backup():
    """Create a backup of the database before migration"""
    backup_dir = create_backup_directory()
    if not backup_dir:
        return False
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup_file = backup_dir / f"db_backup_before_schema_standardization_{timestamp}.sql"
    
    # Create a backup using pg_dump via Supabase MCP
    backup_query = f"""
    SELECT pg_dump_to_text();
    """
    
    logger.info(f"Creating database backup to {backup_file}")
    success, result = execute_sql(backup_query, "Creating database backup")
    
    if not success:
        logger.error(f"Failed to create database backup: {result}")
        return False
    
    try:
        # Write the backup to a file
        with open(backup_file, 'w') as f:
            f.write(result[0]['pg_dump_to_text'])
        logger.info(f"Database backup created successfully: {backup_file}")
        return True
    except Exception as e:
        logger.error(f"Failed to write backup to file: {str(e)}")
        return False

def main():
    """Main function to standardize the schema"""
    logger.info("Starting schema standardization")
    
    # Create a database backup before migration
    if not create_database_backup():
        logger.warning("Failed to create database backup, proceeding with caution")
        if input("Continue without backup? (y/n): ").lower() != 'y':
            logger.info("Migration aborted by user")
            return 1
    
    # Create migration tracking tables
    if not create_migration_tables():
        logger.error("Failed to create migration tracking tables")
        return 1
    
    # Start a new migration
    migration_id = start_migration()
    if not migration_id:
        logger.error("Failed to start migration")
        return 1
    
    logger.info(f"Started migration with ID {migration_id}")
    
    # Process each table
    migration_details = {"tables": []}
    overall_success = True
    
    for singular_name, plural_name in TABLE_MAPPING.items():
        table_result = {"singular": singular_name, "plural": plural_name, "status": "pending"}
        
        # Check if the singular table exists
        if not table_exists(singular_name):
            logger.warning(f"Table {singular_name} does not exist, skipping")
            table_result["status"] = "skipped"
            migration_details["tables"].append(table_result)
            continue
        
        # Check if the plural table already exists
        if table_exists(plural_name):
            logger.warning(f"Table {plural_name} already exists, skipping")
            table_result["status"] = "exists"
            migration_details["tables"].append(table_result)
            continue
        
        logger.info(f"Processing {singular_name} -> {plural_name}")
        
        # Create the plural table
        if not create_plural_table(migration_id, singular_name, plural_name):
            logger.error(f"Failed to create {plural_name} table")
            table_result["status"] = "error_create"
            migration_details["tables"].append(table_result)
            overall_success = False
            continue
        
        # Copy data from singular to plural
        if not copy_data(migration_id, singular_name, plural_name):
            logger.error(f"Failed to copy data from {singular_name} to {plural_name}")
            table_result["status"] = "error_copy"
            migration_details["tables"].append(table_result)
            overall_success = False
            continue
        
        # Add foreign keys
        if not add_foreign_keys(migration_id, singular_name, plural_name):
            logger.error(f"Failed to add foreign keys to {plural_name}")
            table_result["status"] = "error_fk"
            migration_details["tables"].append(table_result)
            overall_success = False
            continue
        
        table_result["status"] = "success"
        migration_details["tables"].append(table_result)
        logger.info(f"Successfully processed {singular_name} -> {plural_name}")
    
    # Update migration status
    if overall_success:
        update_migration_status(migration_id, "completed", migration_details)
        logger.info("Schema standardization completed successfully")
    else:
        logger.error("Schema standardization completed with errors")
        
        # Ask if we should rollback
        if input("Do you want to rollback the migration? (y/n): ").lower() == 'y':
            logger.info("Rolling back migration")
            if rollback_migration(migration_id):
                logger.info("Migration rolled back successfully")
                update_migration_status(migration_id, "rolled_back", migration_details)
            else:
                logger.error("Failed to rollback migration")
                update_migration_status(migration_id, "rollback_failed", migration_details)
        else:
            update_migration_status(migration_id, "completed_with_errors", migration_details)
    
    return 0 if overall_success else 1

if __name__ == "__main__":
    sys.exit(main())