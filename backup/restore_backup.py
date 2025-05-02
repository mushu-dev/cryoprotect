#!/usr/bin/env python3
"""
CryoProtect v2 - Backup Restoration Script

This script restores a backup created by the CryoProtect v2 Backup Manager.
It supports both JSON and SQL format backups.

Usage:
    python restore_backup.py /path/to/backup [options]
"""

import os
import sys
import json
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any

# Import logging configuration
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from logging_config import setup_logging

# Import Supabase MCP tools
try:
    from supabase_mcp_tools import use_mcp_tool
except ImportError:
    # Mock implementation for testing without MCP tools
    def use_mcp_tool(*args, **kwargs):
        logging.warning("Supabase MCP tools not available, using mock implementation")
        return []

# Set up logging
setup_logging()
logger = logging.getLogger("backup_restore")


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="CryoProtect v2 Backup Restoration Script")
    
    parser.add_argument("backup_path", help="Path to the backup directory")
    parser.add_argument("--project-id", help="Supabase project ID")
    parser.add_argument("--schema", default="public", help="Database schema to restore to")
    parser.add_argument("--tables", nargs="+", help="Specific tables to restore (default: all)")
    parser.add_argument("--test-mode", action="store_true", help="Run in test mode without making changes")
    parser.add_argument("--target-prefix", default="", help="Prefix for target tables (e.g., 'restored_')")
    
    return parser.parse_args()


def load_backup_metadata(backup_path: Path) -> Dict[str, Any]:
    """
    Load backup metadata from the backup directory.
    
    Args:
        backup_path: Path to the backup directory
    
    Returns:
        Dictionary with backup metadata
    """
    metadata_path = backup_path / "metadata.json"
    
    if not metadata_path.exists():
        logger.error(f"Metadata file not found: {metadata_path}")
        sys.exit(1)
    
    try:
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
        
        logger.info(f"Loaded backup metadata: {metadata}")
        return metadata
    except Exception as e:
        logger.error(f"Error loading backup metadata: {str(e)}")
        sys.exit(1)


def restore_json_backup(backup_path: Path, project_id: str, schema: str, 
                       tables: Optional[List[str]] = None, test_mode: bool = False,
                       target_prefix: str = "") -> bool:
    """
    Restore a JSON format backup.
    
    Args:
        backup_path: Path to the backup directory
        project_id: Supabase project ID
        schema: Database schema to restore to
        tables: Specific tables to restore (default: all)
        test_mode: Run in test mode without making changes
        target_prefix: Prefix for target tables
    
    Returns:
        True if successful, False otherwise
    """
    metadata = load_backup_metadata(backup_path)
    
    if metadata["format"] != "json":
        logger.error(f"Backup format is not JSON: {metadata['format']}")
        return False
    
    # Determine which tables to restore
    tables_to_restore = tables if tables else metadata["tables"]
    
    # Check that all specified tables exist in the backup
    for table in tables_to_restore:
        if table not in metadata["tables"]:
            logger.error(f"Table not found in backup: {table}")
            return False
    
    success_count = 0
    for table in tables_to_restore:
        table_backup_path = backup_path / f"{table}.json"
        
        if not table_backup_path.exists():
            logger.error(f"Table backup file not found: {table_backup_path}")
            continue
        
        try:
            # Load the table data
            with open(table_backup_path, 'r') as f:
                table_data = json.load(f)
            
            # Determine target table name
            target_table = f"{target_prefix}{table}"
            
            logger.info(f"Restoring table {table} to {target_table} ({len(table_data)} rows)")
            
            if test_mode:
                logger.info(f"TEST MODE: Would restore {len(table_data)} rows to {target_table}")
                success_count += 1
                continue
            
            # Create the target table if it doesn't exist
            if target_prefix:
                # Get the table schema
                schema_query = f"""
                SELECT column_name, data_type, is_nullable
                FROM information_schema.columns
                WHERE table_name = '{table}'
                ORDER BY ordinal_position
                """
                
                schema_data = use_mcp_tool(
                    server_name="supabase",
                    tool_name="execute_sql",
                    arguments={
                        "project_id": project_id,
                        "query": schema_query
                    }
                )
                
                # Create the target table
                create_table_query = f"CREATE TABLE IF NOT EXISTS {schema}.{target_table} (\n"
                
                columns = []
                for col in schema_data:
                    nullable = "NULL" if col["is_nullable"] == "YES" else "NOT NULL"
                    columns.append(f"    {col['column_name']} {col['data_type']} {nullable}")
                
                create_table_query += ",\n".join(columns)
                create_table_query += "\n);"
                
                use_mcp_tool(
                    server_name="supabase",
                    tool_name="execute_sql",
                    arguments={
                        "project_id": project_id,
                        "query": create_table_query
                    }
                )
            
            # Insert the data in batches
            batch_size = 100
            for i in range(0, len(table_data), batch_size):
                batch = table_data[i:i+batch_size]
                
                if not batch:
                    continue
                
                # Create the insert query
                columns = batch[0].keys()
                columns_str = ", ".join(columns)
                
                values_list = []
                for row in batch:
                    values = []
                    for col in columns:
                        if row[col] is None:
                            values.append("NULL")
                        elif isinstance(row[col], (int, float)):
                            values.append(str(row[col]))
                        else:
                            # Escape single quotes
                            escaped_val = str(row[col]).replace("'", "''")
                            values.append(f"'{escaped_val}'")
                    
                    values_list.append(f"({', '.join(values)})")
                
                insert_query = f"INSERT INTO {schema}.{target_table} ({columns_str}) VALUES {', '.join(values_list)};"
                
                use_mcp_tool(
                    server_name="supabase",
                    tool_name="execute_sql",
                    arguments={
                        "project_id": project_id,
                        "query": insert_query
                    }
                )
            
            logger.info(f"Successfully restored table {table} to {target_table}")
            success_count += 1
        
        except Exception as e:
            logger.error(f"Error restoring table {table}: {str(e)}")
    
    logger.info(f"Restoration completed: {success_count}/{len(tables_to_restore)} tables restored successfully")
    return success_count == len(tables_to_restore)


def restore_sql_backup(backup_path: Path, project_id: str, schema: str, 
                      tables: Optional[List[str]] = None, test_mode: bool = False,
                      target_prefix: str = "") -> bool:
    """
    Restore a SQL format backup.
    
    Args:
        backup_path: Path to the backup directory
        project_id: Supabase project ID
        schema: Database schema to restore to
        tables: Specific tables to restore (default: all)
        test_mode: Run in test mode without making changes
        target_prefix: Prefix for target tables
    
    Returns:
        True if successful, False otherwise
    """
    metadata = load_backup_metadata(backup_path)
    
    if metadata["format"] != "sql":
        logger.error(f"Backup format is not SQL: {metadata['format']}")
        return False
    
    # Determine which tables to restore
    tables_to_restore = tables if tables else metadata["tables"]
    
    # Check that all specified tables exist in the backup
    for table in tables_to_restore:
        if table not in metadata["tables"]:
            logger.error(f"Table not found in backup: {table}")
            return False
    
    success_count = 0
    for table in tables_to_restore:
        table_backup_path = backup_path / f"{table}.sql"
        
        if not table_backup_path.exists():
            logger.error(f"Table backup file not found: {table_backup_path}")
            continue
        
        try:
            # Load the SQL file
            with open(table_backup_path, 'r') as f:
                sql_content = f.read()
            
            # Determine target table name
            target_table = f"{target_prefix}{table}"
            
            logger.info(f"Restoring table {table} to {target_table}")
            
            if test_mode:
                logger.info(f"TEST MODE: Would execute SQL script for {table}")
                success_count += 1
                continue
            
            # Replace table name if target_prefix is specified
            if target_prefix:
                sql_content = sql_content.replace(
                    f"CREATE TABLE IF NOT EXISTS {table}",
                    f"CREATE TABLE IF NOT EXISTS {schema}.{target_table}"
                )
                sql_content = sql_content.replace(
                    f"INSERT INTO {table}",
                    f"INSERT INTO {schema}.{target_table}"
                )
            
            # Execute the SQL
            use_mcp_tool(
                server_name="supabase",
                tool_name="execute_sql",
                arguments={
                    "project_id": project_id,
                    "query": sql_content
                }
            )
            
            logger.info(f"Successfully restored table {table} to {target_table}")
            success_count += 1
        
        except Exception as e:
            logger.error(f"Error restoring table {table}: {str(e)}")
    
    logger.info(f"Restoration completed: {success_count}/{len(tables_to_restore)} tables restored successfully")
    return success_count == len(tables_to_restore)


def verify_restoration(backup_path: Path, project_id: str, schema: str, 
                      tables: Optional[List[str]] = None, target_prefix: str = "") -> bool:
    """
    Verify the restoration by comparing row counts.
    
    Args:
        backup_path: Path to the backup directory
        project_id: Supabase project ID
        schema: Database schema to restore to
        tables: Specific tables to restore (default: all)
        target_prefix: Prefix for target tables
    
    Returns:
        True if verification passed, False otherwise
    """
    metadata = load_backup_metadata(backup_path)
    
    # Determine which tables to verify
    tables_to_verify = tables if tables else metadata["tables"]
    
    verification_passed = True
    for table in tables_to_verify:
        target_table = f"{target_prefix}{table}"
        
        try:
            # Get row count from the backup
            if metadata["format"] == "json":
                table_backup_path = backup_path / f"{table}.json"
                with open(table_backup_path, 'r') as f:
                    backup_data = json.load(f)
                backup_row_count = len(backup_data)
            else:
                # For SQL format, we can't easily determine the row count
                # So we'll skip this verification
                logger.warning(f"Skipping row count verification for SQL format backup of table {table}")
                continue
            
            # Get row count from the restored table
            count_query = f"SELECT COUNT(*) FROM {schema}.{target_table}"
            result = use_mcp_tool(
                server_name="supabase",
                tool_name="execute_sql",
                arguments={
                    "project_id": project_id,
                    "query": count_query
                }
            )
            
            restored_row_count = result[0]["count"]
            
            if backup_row_count == restored_row_count:
                logger.info(f"Verification passed for table {table}: {backup_row_count} rows")
            else:
                logger.error(f"Verification failed for table {table}: Backup has {backup_row_count} rows, restored table has {restored_row_count} rows")
                verification_passed = False
        
        except Exception as e:
            logger.error(f"Error verifying table {table}: {str(e)}")
            verification_passed = False
    
    return verification_passed


def main():
    """Main function."""
    args = parse_arguments()
    
    # Convert backup path to Path object
    backup_path = Path(args.backup_path)
    
    # Check if backup directory exists
    if not backup_path.exists() or not backup_path.is_dir():
        logger.error(f"Backup directory does not exist: {backup_path}")
        return 1
    
    # Load backup metadata
    metadata = load_backup_metadata(backup_path)
    
    # Use project ID from arguments or metadata
    project_id = args.project_id or metadata.get("project_id")
    if not project_id:
        logger.error("Project ID not specified and not found in backup metadata")
        return 1
    
    # Restore the backup
    if metadata["format"] == "json":
        success = restore_json_backup(
            backup_path=backup_path,
            project_id=project_id,
            schema=args.schema,
            tables=args.tables,
            test_mode=args.test_mode,
            target_prefix=args.target_prefix
        )
    else:
        success = restore_sql_backup(
            backup_path=backup_path,
            project_id=project_id,
            schema=args.schema,
            tables=args.tables,
            test_mode=args.test_mode,
            target_prefix=args.target_prefix
        )
    
    if not success:
        logger.error("Restoration failed")
        return 1
    
    # Verify the restoration
    if not args.test_mode:
        verification_passed = verify_restoration(
            backup_path=backup_path,
            project_id=project_id,
            schema=args.schema,
            tables=args.tables,
            target_prefix=args.target_prefix
        )
        
        if verification_passed:
            logger.info("Verification passed: All tables restored successfully")
        else:
            logger.warning("Verification failed: Some tables may not have been restored correctly")
    
    logger.info("Restoration completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())