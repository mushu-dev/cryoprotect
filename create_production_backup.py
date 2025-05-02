#!/usr/bin/env python3
"""
CryoProtect v2 - Production Database Backup Script

This script creates a comprehensive backup of the production database before
deploying performance improvements. It includes:
1. Full schema backup
2. Table data backup
3. Index definitions backup
4. RLS policy backup
5. Function backup

The backup is saved in both JSON and SQL formats for maximum recoverability.
"""

import os
import sys
import json
import time
import logging
import argparse
from datetime import datetime
from pathlib import Path

# Import logging configuration
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from logging_config import setup_logging
except ImportError:
    # Fallback logging setup if the import fails
    def setup_logging():
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers=[
                logging.FileHandler("production_backup.log"),
                logging.StreamHandler()
            ]
        )

# Import Supabase MCP tools
try:
    from supabase_mcp_tools import execute_sql_on_supabase
except ImportError:
    print("Error: supabase_mcp_tools.py not found. Make sure it's in the same directory.")
    sys.exit(1)

# Set up logging
setup_logging()
logger = logging.getLogger("production_backup")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Create a comprehensive backup of the production database")
    parser.add_argument("--project-id", default="tsdlmynydfuypiugmkev",
                        help="Supabase project ID (default: tsdlmynydfuypiugmkev)")
    parser.add_argument("--schema", default="public",
                        help="Database schema to backup (default: public)")
    parser.add_argument("--format", choices=["json", "sql", "both"], default="both",
                        help="Output format for the backup (default: both)")
    parser.add_argument("--output-dir", default="production_backups",
                        help="Directory to save backups (default: production_backups)")
    return parser.parse_args()

def create_backup_directory(base_dir):
    """
    Create the backup directory with timestamp.
    
    Args:
        base_dir (str): Base directory for backups
        
    Returns:
        Path: Path to the created backup directory
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = Path(base_dir) / f"backup_{timestamp}"
    
    if not backup_dir.exists():
        logger.info(f"Creating backup directory: {backup_dir}")
        backup_dir.mkdir(parents=True)
    
    return backup_dir

def backup_schema(project_id, schema, backup_dir):
    """
    Backup the database schema.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema to backup
        backup_dir (Path): Directory to save the backup
        
    Returns:
        dict: Schema information
    """
    logger.info(f"Backing up schema: {schema}")
    
    try:
        # Get schema information
        query = f"""
        SELECT 
            table_name,
            column_name,
            data_type,
            column_default,
            is_nullable,
            character_maximum_length
        FROM 
            information_schema.columns
        WHERE 
            table_schema = '{schema}'
        ORDER BY 
            table_name, ordinal_position;
        """
        
        schema_data = execute_sql_on_supabase(project_id, query)
        
        # Save schema information
        schema_file = backup_dir / "schema.json"
        with open(schema_file, "w") as f:
            json.dump(schema_data, f, indent=2)
        
        logger.info(f"Schema backup saved to: {schema_file}")
        return schema_data
    
    except Exception as e:
        logger.error(f"Error backing up schema: {str(e)}")
        return None

def backup_tables(project_id, schema, backup_dir, format="both"):
    """
    Backup all tables in the schema.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema to backup
        backup_dir (Path): Directory to save the backup
        format (str): Output format (json, sql, or both)
        
    Returns:
        dict: Information about backed up tables
    """
    logger.info(f"Backing up tables in schema: {schema}")
    
    try:
        # Get list of tables
        query = f"""
        SELECT 
            table_name
        FROM 
            information_schema.tables
        WHERE 
            table_schema = '{schema}'
            AND table_type = 'BASE TABLE';
        """
        
        tables_result = execute_sql_on_supabase(project_id, query)
        tables = [table["table_name"] for table in tables_result]
        
        # Create tables directory
        tables_dir = backup_dir / "tables"
        if not tables_dir.exists():
            tables_dir.mkdir()
        
        # Backup each table
        table_info = {}
        for table in tables:
            logger.info(f"Backing up table: {table}")
            
            # Get table data
            data_query = f"SELECT * FROM {schema}.{table};"
            table_data = execute_sql_on_supabase(project_id, data_query)
            
            # Save as JSON if requested
            if format in ["json", "both"]:
                json_file = tables_dir / f"{table}.json"
                with open(json_file, "w") as f:
                    json.dump(table_data, f, indent=2)
            
            # Save as SQL if requested
            if format in ["sql", "both"]:
                sql_file = tables_dir / f"{table}.sql"
                with open(sql_file, "w") as f:
                    # Write table creation statement
                    f.write(f"-- Table: {schema}.{table}\n")
                    
                    # Get table schema
                    schema_query = f"""
                    SELECT 
                        column_name, 
                        data_type, 
                        column_default,
                        is_nullable
                    FROM 
                        information_schema.columns
                    WHERE 
                        table_schema = '{schema}'
                        AND table_name = '{table}'
                    ORDER BY 
                        ordinal_position;
                    """
                    
                    columns = execute_sql_on_supabase(project_id, schema_query)
                    
                    # Write CREATE TABLE statement
                    f.write(f"CREATE TABLE IF NOT EXISTS {schema}.{table} (\n")
                    column_defs = []
                    for col in columns:
                        nullable = "NULL" if col["is_nullable"] == "YES" else "NOT NULL"
                        default = f"DEFAULT {col['column_default']}" if col["column_default"] else ""
                        column_defs.append(f"    {col['column_name']} {col['data_type']} {nullable} {default}".strip())
                    
                    f.write(",\n".join(column_defs))
                    f.write("\n);\n\n")
                    
                    # Write INSERT statements
                    f.write(f"-- Data for table: {schema}.{table}\n")
                    for row in table_data:
                        columns = ", ".join(row.keys())
                        # Format values properly
                        values = []
                        for val in row.values():
                            if val is None:
                                values.append("NULL")
                            elif isinstance(val, (int, float)):
                                values.append(str(val))
                            else:
                                # Escape single quotes
                                val_str = str(val).replace("'", "''")
                                values.append(f"'{val_str}'")
                        
                        values_str = ", ".join(values)
                        f.write(f"INSERT INTO {schema}.{table} ({columns}) VALUES ({values_str});\n")
            
            # Record table info
            table_info[table] = {
                "row_count": len(table_data),
                "json_backup": str(tables_dir / f"{table}.json") if format in ["json", "both"] else None,
                "sql_backup": str(tables_dir / f"{table}.sql") if format in ["sql", "both"] else None
            }
        
        # Save table info
        table_info_file = backup_dir / "table_info.json"
        with open(table_info_file, "w") as f:
            json.dump(table_info, f, indent=2)
        
        logger.info(f"Table backups saved to: {tables_dir}")
        return table_info
    
    except Exception as e:
        logger.error(f"Error backing up tables: {str(e)}")
        return None

def backup_indexes(project_id, schema, backup_dir):
    """
    Backup all indexes in the schema.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema to backup
        backup_dir (Path): Directory to save the backup
        
    Returns:
        list: Index definitions
    """
    logger.info(f"Backing up indexes in schema: {schema}")
    
    try:
        # Get index definitions
        query = f"""
        SELECT
            t.relname AS table_name,
            i.relname AS index_name,
            a.attname AS column_name,
            ix.indisunique AS is_unique,
            ix.indisprimary AS is_primary,
            pg_get_indexdef(ix.indexrelid) AS index_definition
        FROM
            pg_index ix
            JOIN pg_class i ON i.oid = ix.indexrelid
            JOIN pg_class t ON t.oid = ix.indrelid
            JOIN pg_namespace n ON n.oid = t.relnamespace
            JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(ix.indkey)
        WHERE
            n.nspname = '{schema}'
        ORDER BY
            t.relname, i.relname, a.attnum;
        """
        
        indexes = execute_sql_on_supabase(project_id, query)
        
        # Save index definitions
        indexes_file = backup_dir / "indexes.json"
        with open(indexes_file, "w") as f:
            json.dump(indexes, f, indent=2)
        
        # Save as SQL
        indexes_sql_file = backup_dir / "indexes.sql"
        with open(indexes_sql_file, "w") as f:
            f.write(f"-- Indexes for schema: {schema}\n\n")
            
            # Group by index name
            index_defs = {}
            for idx in indexes:
                if idx["index_definition"] not in index_defs:
                    index_defs[idx["index_definition"]] = True
            
            # Write index definitions
            for index_def in index_defs:
                f.write(f"{index_def};\n\n")
        
        logger.info(f"Index backup saved to: {indexes_file} and {indexes_sql_file}")
        return indexes
    
    except Exception as e:
        logger.error(f"Error backing up indexes: {str(e)}")
        return None

def backup_rls_policies(project_id, schema, backup_dir):
    """
    Backup all RLS policies in the schema.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema to backup
        backup_dir (Path): Directory to save the backup
        
    Returns:
        list: RLS policy definitions
    """
    logger.info(f"Backing up RLS policies in schema: {schema}")
    
    try:
        # Get RLS policy definitions
        query = f"""
        SELECT
            schemaname,
            tablename,
            policyname,
            permissive,
            roles,
            cmd,
            qual,
            with_check
        FROM
            pg_policies
        WHERE
            schemaname = '{schema}';
        """
        
        policies = execute_sql_on_supabase(project_id, query)
        
        # Save policy definitions
        policies_file = backup_dir / "rls_policies.json"
        with open(policies_file, "w") as f:
            json.dump(policies, f, indent=2)
        
        # Save as SQL
        policies_sql_file = backup_dir / "rls_policies.sql"
        with open(policies_sql_file, "w") as f:
            f.write(f"-- RLS policies for schema: {schema}\n\n")
            
            # Enable RLS on tables
            for policy in policies:
                table = policy["tablename"]
                f.write(f"ALTER TABLE {schema}.{table} ENABLE ROW LEVEL SECURITY;\n")
            
            f.write("\n")
            
            # Create policies
            for policy in policies:
                table = policy["tablename"]
                policy_name = policy["policyname"]
                cmd = policy["cmd"]
                roles = policy["roles"]
                using = policy["qual"]
                with_check = policy["with_check"]
                
                f.write(f"CREATE POLICY {policy_name} ON {schema}.{table}\n")
                f.write(f"    FOR {cmd}\n")
                
                if roles and roles != "{public}":
                    roles_str = ", ".join(roles[1:-1].split(","))
                    f.write(f"    TO {roles_str}\n")
                
                if using:
                    f.write(f"    USING ({using})\n")
                
                if with_check:
                    f.write(f"    WITH CHECK ({with_check})\n")
                
                f.write(";\n\n")
        
        logger.info(f"RLS policy backup saved to: {policies_file} and {policies_sql_file}")
        return policies
    
    except Exception as e:
        logger.error(f"Error backing up RLS policies: {str(e)}")
        return None

def backup_functions(project_id, schema, backup_dir):
    """
    Backup all functions in the schema.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema to backup
        backup_dir (Path): Directory to save the backup
        
    Returns:
        list: Function definitions
    """
    logger.info(f"Backing up functions in schema: {schema}")
    
    try:
        # Get function definitions
        query = f"""
        SELECT
            p.proname AS function_name,
            pg_get_functiondef(p.oid) AS function_definition
        FROM
            pg_proc p
            JOIN pg_namespace n ON p.pronamespace = n.oid
        WHERE
            n.nspname = '{schema}';
        """
        
        functions = execute_sql_on_supabase(project_id, query)
        
        # Save function definitions
        functions_file = backup_dir / "functions.json"
        with open(functions_file, "w") as f:
            json.dump(functions, f, indent=2)
        
        # Save as SQL
        functions_sql_file = backup_dir / "functions.sql"
        with open(functions_sql_file, "w") as f:
            f.write(f"-- Functions for schema: {schema}\n\n")
            
            for func in functions:
                f.write(f"{func['function_definition']};\n\n")
        
        logger.info(f"Function backup saved to: {functions_file} and {functions_sql_file}")
        return functions
    
    except Exception as e:
        logger.error(f"Error backing up functions: {str(e)}")
        return None

def create_full_backup_script(backup_dir, schema):
    """
    Create a single SQL script that can restore the entire database.
    
    Args:
        backup_dir (Path): Directory with backups
        schema (str): The database schema
        
    Returns:
        str: Path to the full backup script
    """
    logger.info("Creating full backup script")
    
    try:
        # Create full backup script
        full_backup_file = backup_dir / "full_backup.sql"
        
        with open(full_backup_file, "w") as f:
            f.write(f"-- CryoProtect v2 Full Database Backup\n")
            f.write(f"-- Schema: {schema}\n")
            f.write(f"-- Created: {datetime.now().isoformat()}\n\n")
            
            # Create schema
            f.write(f"-- Create schema\n")
            f.write(f"CREATE SCHEMA IF NOT EXISTS {schema};\n\n")
            
            # Add table definitions and data
            f.write(f"-- Tables and data\n")
            tables_dir = backup_dir / "tables"
            if tables_dir.exists():
                for sql_file in sorted(tables_dir.glob("*.sql")):
                    with open(sql_file, "r") as table_file:
                        f.write(table_file.read())
                    f.write("\n\n")
            
            # Add indexes
            indexes_sql_file = backup_dir / "indexes.sql"
            if indexes_sql_file.exists():
                with open(indexes_sql_file, "r") as idx_file:
                    f.write(idx_file.read())
                f.write("\n\n")
            
            # Add functions
            functions_sql_file = backup_dir / "functions.sql"
            if functions_sql_file.exists():
                with open(functions_sql_file, "r") as func_file:
                    f.write(func_file.read())
                f.write("\n\n")
            
            # Add RLS policies
            policies_sql_file = backup_dir / "rls_policies.sql"
            if policies_sql_file.exists():
                with open(policies_sql_file, "r") as policy_file:
                    f.write(policy_file.read())
                f.write("\n\n")
        
        logger.info(f"Full backup script saved to: {full_backup_file}")
        return str(full_backup_file)
    
    except Exception as e:
        logger.error(f"Error creating full backup script: {str(e)}")
        return None

def create_backup_summary(backup_dir, project_id, schema, format, components):
    """
    Create a summary of the backup.
    
    Args:
        backup_dir (Path): Directory with backups
        project_id (str): The Supabase project ID
        schema (str): The database schema
        format (str): Backup format
        components (dict): Backup components
        
    Returns:
        dict: Backup summary
    """
    logger.info("Creating backup summary")
    
    try:
        # Create summary
        summary = {
            "timestamp": datetime.now().isoformat(),
            "project_id": project_id,
            "schema": schema,
            "format": format,
            "backup_directory": str(backup_dir),
            "components": {}
        }
        
        # Add component information
        for component, data in components.items():
            if data:
                if isinstance(data, list):
                    summary["components"][component] = {
                        "count": len(data),
                        "status": "SUCCESS"
                    }
                elif isinstance(data, dict):
                    summary["components"][component] = {
                        "count": len(data),
                        "status": "SUCCESS"
                    }
                else:
                    summary["components"][component] = {
                        "status": "SUCCESS" if data else "FAILED"
                    }
            else:
                summary["components"][component] = {
                    "status": "FAILED"
                }
        
        # Determine overall status
        failed_components = [c for c, info in summary["components"].items() if info["status"] == "FAILED"]
        if not failed_components:
            summary["status"] = "SUCCESS"
        elif len(failed_components) < len(components):
            summary["status"] = "COMPLETED_WITH_WARNINGS"
            summary["warnings"] = f"Failed to backup components: {', '.join(failed_components)}"
        else:
            summary["status"] = "ERROR"
            summary["error"] = "All backup components failed"
        
        # Save summary
        summary_file = backup_dir / "backup_summary.json"
        with open(summary_file, "w") as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Backup summary saved to: {summary_file}")
        return summary
    
    except Exception as e:
        logger.error(f"Error creating backup summary: {str(e)}")
        return None

def main():
    """Main function to create a production database backup."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Create backup directory
        backup_dir = create_backup_directory(args.output_dir)
        
        logger.info(f"Starting production backup of Supabase project {args.project_id}")
        logger.info(f"Schema: {args.schema}")
        logger.info(f"Format: {args.format}")
        logger.info(f"Backup directory: {backup_dir}")
        
        # Perform backups
        backup_components = {}
        
        # Backup schema
        backup_components["schema"] = backup_schema(args.project_id, args.schema, backup_dir)
        
        # Backup tables
        backup_components["tables"] = backup_tables(args.project_id, args.schema, backup_dir, args.format)
        
        # Backup indexes
        backup_components["indexes"] = backup_indexes(args.project_id, args.schema, backup_dir)
        
        # Backup RLS policies
        backup_components["rls_policies"] = backup_rls_policies(args.project_id, args.schema, backup_dir)
        
        # Backup functions
        backup_components["functions"] = backup_functions(args.project_id, args.schema, backup_dir)
        
        # Create full backup script
        backup_components["full_backup_script"] = create_full_backup_script(backup_dir, args.schema)
        
        # Create backup summary
        summary = create_backup_summary(backup_dir, args.project_id, args.schema, args.format, backup_components)
        
        # Print summary
        if summary:
            logger.info(f"Backup completed with status: {summary['status']}")
            
            if summary["status"] == "COMPLETED_WITH_WARNINGS":
                logger.warning(summary["warnings"])
            elif summary["status"] == "ERROR":
                logger.error(summary["error"])
            
            logger.info(f"Backup saved to: {backup_dir}")
            
            # Return appropriate exit code
            if summary["status"] == "SUCCESS":
                return 0
            elif summary["status"] == "COMPLETED_WITH_WARNINGS":
                return 1
            else:
                return 2
        else:
            logger.error("Failed to create backup summary")
            return 3
    
    except Exception as e:
        logger.error(f"Backup failed: {str(e)}")
        return 4

if __name__ == "__main__":
    sys.exit(main())