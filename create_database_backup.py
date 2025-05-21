#!/usr/bin/env python3
"""
CryoProtect v2 - Supabase Database Backup Script

This script creates a timestamped backup of all tables in the public schema
of the CryoProtect Supabase project. It uses the Supabase MCP tools for
database operations.

Usage:
    python create_database_backup.py

The script will:
1. Connect to the Supabase project
2. List all tables in the public schema
3. Export data from each table
4. Save the data in JSON format in the 'backups' directory
"""

import os
import sys
import json
import time
import logging
from datetime import datetime
import argparse
from pathlib import Path

# Import logging configuration
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from logging_config import setup_logging

# Import Supabase MCP tools
from supabase_mcp_tools import use_mcp_tool

# Set up logging
setup_logging()
logger = logging.getLogger("database_backup")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Create a backup of the CryoProtect Supabase database")
    parser.add_argument("--format", choices=["json", "sql"], default="json",
                        help="Output format for the backup (default: json)")
    parser.add_argument("--project-id", default="tsdlmynydfuypiugmkev",
                        help="Supabase project ID (default: tsdlmynydfuypiugmkev)")
    parser.add_argument("--schema", default="public",
                        help="Database schema to backup (default: public)")
    return parser.parse_args()

def create_backup_directory():
    """Create the backups directory if it doesn't exist."""
    backup_dir = Path("backups")
    if not backup_dir.exists():
        logger.info(f"Creating backup directory: {backup_dir}")
        backup_dir.mkdir(parents=True)
    return backup_dir

def get_timestamp():
    """Generate a timestamp string for the backup filename."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def backup_table_to_json(project_id, table_name, backup_path):
    """
    Backup a table to a JSON file using Supabase MCP tools.
    
    Args:
        project_id (str): The Supabase project ID
        table_name (str): The name of the table to backup
        backup_path (Path): The path to save the backup file
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Log the operation
        logger.info(f"Backing up table: {table_name}")
        
        # Use the MCP tool to execute SQL
        try:
            data = use_mcp_tool(
                server_name="supabase",
                tool_name="execute_sql",
                arguments={
                    "project_id": project_id,
                    "query": f"SELECT * FROM {table_name}"
                }
            )
            
            # Write the data to the backup file
            with open(backup_path, 'w') as f:
                json.dump(data, f, indent=2)
            
            logger.info(f"Successfully backed up table {table_name} to {backup_path}")
            return True
        except json.JSONDecodeError:
            logger.error(f"Failed to parse JSON output for table {table_name}")
            return False
    except Exception as e:
        logger.error(f"Error backing up table {table_name}: {str(e)}")
        return False

def backup_table_to_sql(project_id, table_name, backup_path):
    """
    Backup a table to a SQL file using Supabase MCP tools.
    
    Args:
        project_id (str): The Supabase project ID
        table_name (str): The name of the table to backup
        backup_path (Path): The path to save the backup file
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Log the operation
        logger.info(f"Backing up table: {table_name}")
        
        # First, get the table schema
        try:
            schema_data = use_mcp_tool(
                server_name="supabase",
                tool_name="execute_sql",
                arguments={
                    "project_id": project_id,
                    "query": f"SELECT column_name, data_type, is_nullable FROM information_schema.columns WHERE table_name = '{table_name}' ORDER BY ordinal_position"
                }
            )
            
            # Then, get the table data
            table_data = use_mcp_tool(
                server_name="supabase",
                tool_name="execute_sql",
                arguments={
                    "project_id": project_id,
                    "query": f"SELECT * FROM {table_name}"
                }
            )
            
            # Generate SQL statements
            with open(backup_path, 'w') as f:
                # Write table creation statement
                f.write(f"-- Table: {table_name}\n")
                f.write(f"CREATE TABLE IF NOT EXISTS {table_name} (\n")
                
                # Write column definitions
                columns = []
                for col in schema_data:
                    nullable = "NULL" if col["is_nullable"] == "YES" else "NOT NULL"
                    columns.append(f"    {col['column_name']} {col['data_type']} {nullable}")
                
                f.write(",\n".join(columns))
                f.write("\n);\n\n")
                
                # Write data insertion statements
                f.write(f"-- Data for table: {table_name}\n")
                for row in table_data:
                    columns = ", ".join(row.keys())
                    # Create a list of formatted values
                    formatted_values = []
                    for val in row.values():
                        if val is not None:
                            # Escape single quotes in SQL
                            escaped_val = str(val).replace("'", "''")
                            formatted_values.append(f"'{escaped_val}'")
                        else:
                            formatted_values.append("NULL")
                    
                    values = ", ".join(formatted_values)
                    f.write(f"INSERT INTO {table_name} ({columns}) VALUES ({values});\n")
            
            logger.info(f"Successfully backed up table {table_name} to {backup_path}")
            return True
        except json.JSONDecodeError:
            logger.error(f"Failed to parse JSON output for table {table_name}")
            return False
    except Exception as e:
        logger.error(f"Error backing up table {table_name}: {str(e)}")
        return False

def list_tables(project_id, schema="public"):
    """
    List all tables in the specified schema using Supabase MCP tools.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema to query
    
    Returns:
        list: A list of table names
    """
    try:
        logger.info(f"Listing tables in schema: {schema}")
        
        # Use the MCP tool to execute SQL
        try:
            data = use_mcp_tool(
                server_name="supabase",
                tool_name="execute_sql",
                arguments={
                    "project_id": project_id,
                    "query": f"SELECT table_name FROM information_schema.tables WHERE table_schema = '{schema}' AND table_type = 'BASE TABLE'"
                }
            )
            
            tables = [item["table_name"] for item in data]
            logger.info(f"Found {len(tables)} tables in schema {schema}")
            return tables
        except json.JSONDecodeError:
            logger.error("Failed to parse JSON output for table list")
            return []
    except Exception as e:
        logger.error(f"Error listing tables: {str(e)}")
        return []

def main():
    """Main function to create a database backup."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Create backup directory
        backup_dir = create_backup_directory()
        
        # Generate timestamp for the backup
        timestamp = get_timestamp()
        
        # Create a subdirectory for this backup
        backup_subdir = backup_dir / f"backup_{timestamp}"
        backup_subdir.mkdir(parents=True)
        
        # Log the start of the backup process
        logger.info(f"Starting backup of Supabase project {args.project_id}")
        logger.info(f"Backup format: {args.format}")
        logger.info(f"Backup directory: {backup_subdir}")
        
        # List all tables in the schema
        tables = list_tables(args.project_id, args.schema)
        
        if not tables:
            logger.error("No tables found or failed to list tables")
            return 1
        
        # Create a metadata file with backup information
        metadata = {
            "timestamp": timestamp,
            "project_id": args.project_id,
            "schema": args.schema,
            "format": args.format,
            "tables": tables
        }
        
        with open(backup_subdir / "metadata.json", "w") as f:
            json.dump(metadata, f, indent=2)
        
        # Backup each table
        success_count = 0
        for table in tables:
            file_extension = ".json" if args.format == "json" else ".sql"
            backup_path = backup_subdir / f"{table}{file_extension}"
            
            if args.format == "json":
                success = backup_table_to_json(args.project_id, table, backup_path)
            else:
                success = backup_table_to_sql(args.project_id, table, backup_path)
            
            if success:
                success_count += 1
        
        # Log the completion of the backup process
        logger.info(f"Backup completed: {success_count}/{len(tables)} tables backed up successfully")
        logger.info(f"Backup saved to: {backup_subdir}")
        
        # Create a summary file
        summary = {
            "timestamp": timestamp,
            "project_id": args.project_id,
            "schema": args.schema,
            "format": args.format,
            "total_tables": len(tables),
            "successful_tables": success_count,
            "backup_path": str(backup_subdir)
        }
        
        with open(backup_subdir / "summary.json", "w") as f:
            json.dump(summary, f, indent=2)
        
        return 0
    except Exception as e:
        logger.error(f"Backup failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())