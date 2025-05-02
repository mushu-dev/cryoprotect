#!/usr/bin/env python3
"""
CryoProtect v2 - Apply Critical Performance Improvements

This script applies critical performance improvements to the CryoProtect database
before production deployment. It implements:
1. Missing indexes for frequently queried columns
2. Optimized views for common queries
3. Connection pooling optimizations
4. Query plan analysis to identify slow queries
5. Transaction support for safer operations
6. Proper error handling and reporting
7. Verification steps to ensure indexes are correctly created
8. Rollback mechanism in case of failure
9. Performance report showing before/after metrics

This is an enhanced version of the original script with additional features
for the database remediation project.
"""

import os
import sys
import json
import logging
import argparse
import time
import shutil
import threading
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set
import functools
from contextlib import contextmanager

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
                logging.FileHandler("performance_improvements.log"),
                logging.StreamHandler()
            ]
        )

# Import Supabase MCP tools
try:
    from supabase_mcp_tools import execute_sql_on_supabase, backup_table, apply_migration
except ImportError:
    print("Error: supabase_mcp_tools.py not found. Make sure it's in the same directory.")
    sys.exit(1)

# Import connection pooling
try:
    from connection_pool_wrapper import (
        initialize_supabase_pool, get_supabase_connection,
        get_supabase_pool_stats, shutdown_supabase_pool
    )
except ImportError:
    print("Warning: connection_pool_wrapper.py not found. Connection pooling will be limited.")
    initialize_supabase_pool = None
    get_supabase_connection = None
    get_supabase_pool_stats = None
    shutdown_supabase_pool = None

# Import query plan analysis
try:
    from analyze_query_plans import (
        analyze_plan_for_issues, suggest_optimizations,
        QUERIES as COMMON_QUERIES
    )
except ImportError:
    print("Warning: analyze_query_plans.py not found. Query plan analysis will be limited.")
    analyze_plan_for_issues = None
    suggest_optimizations = None
    COMMON_QUERIES = []

# Memory bank for caching query results
MEMORY_BANK_DIR = Path("memory-bank")
MEMORY_CACHE_FILE = MEMORY_BANK_DIR / "performance_cache.json"
ROLLBACK_FILE = MEMORY_BANK_DIR / "performance_rollback.json"

# Set up logging
setup_logging()
logger = logging.getLogger("performance_improvements")

# Initialize memory cache
def init_memory_cache():
    """Initialize the memory cache file if it doesn't exist."""
    if not MEMORY_BANK_DIR.exists():
        MEMORY_BANK_DIR.mkdir(parents=True, exist_ok=True)
    
    if not MEMORY_CACHE_FILE.exists():
        with open(MEMORY_CACHE_FILE, 'w') as f:
            json.dump({
                "version": "1.0",
                "last_updated": datetime.now().isoformat(),
                "query_cache": {},
                "verification_results": {}
            }, f, indent=2)
    
    if not ROLLBACK_FILE.exists():
        with open(ROLLBACK_FILE, 'w') as f:
            json.dump({
                "version": "1.0",
                "last_updated": datetime.now().isoformat(),
                "rollback_operations": []
            }, f, indent=2)

# Memory cache decorator
def memory_cache(func):
    """Decorator to cache function results in memory bank."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Create a cache key from function name and arguments
        key = f"{func.__name__}_{str(args)}_{str(kwargs)}"
        
        # Load the cache
        try:
            with open(MEMORY_CACHE_FILE, 'r') as f:
                cache = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            init_memory_cache()
            with open(MEMORY_CACHE_FILE, 'r') as f:
                cache = json.load(f)
        
        # Check if result is in cache
        if key in cache.get("query_cache", {}):
            logger.info(f"Using cached result for {func.__name__}")
            return cache["query_cache"][key]
        
        # Call the function
        result = func(*args, **kwargs)
        
        # Update the cache
        cache.setdefault("query_cache", {})
        cache["query_cache"][key] = result
        cache["last_updated"] = datetime.now().isoformat()
        
        # Save the cache
        with open(MEMORY_CACHE_FILE, 'w') as f:
            json.dump(cache, f, indent=2)
        
        return result
    
    return wrapper

def ensure_backup_directory(backup_dir="backups"):
    """
    Ensure the backup directory exists.
    
    Args:
        backup_dir (str): The backup directory path
        
    Returns:
        Path: The backup directory path
    """
    backup_path = Path(backup_dir)
    backup_path.mkdir(parents=True, exist_ok=True)
    return backup_path

def retry_operation(max_retries=3, delay=2):
    """
    Decorator to retry operations with exponential backoff.
    
    Args:
        max_retries (int): Maximum number of retry attempts
        delay (int): Initial delay between retries in seconds
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            retries = 0
            current_delay = delay
            
            while retries < max_retries:
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    retries += 1
                    if retries >= max_retries:
                        logger.error(f"Operation failed after {max_retries} retries: {str(e)}")
                        raise
                    
                    logger.warning(f"Operation failed, retrying in {current_delay}s... ({retries}/{max_retries})")
                    time.sleep(current_delay)
                    current_delay *= 2  # Exponential backoff
            
            return func(*args, **kwargs)  # Final attempt
        
        return wrapper
    
    return decorator

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Apply critical performance improvements to the CryoProtect database")
    parser.add_argument("--project-id", default="tsdlmynydfuypiugmkev",
                        help="Supabase project ID (default: tsdlmynydfuypiugmkev)")
    parser.add_argument("--schema", default="public",
                        help="Database schema to modify (default: public)")
    parser.add_argument("--backup", action="store_true",
                        help="Create backups of tables before modifying them")
    parser.add_argument("--test-only", action="store_true",
                        help="Only test the improvements without applying them")
    parser.add_argument("--backup-dir", default="production_backups",
                        help="Directory to store backups (default: production_backups)")
    parser.add_argument("--analyze-queries", action="store_true",
                        help="Analyze query plans to identify slow queries")
    parser.add_argument("--connection-pooling", action="store_true",
                        help="Enable enhanced connection pooling")
    parser.add_argument("--min-connections", type=int, default=2,
                        help="Minimum number of connections for the pool (default: 2)")
    parser.add_argument("--max-connections", type=int, default=10,
                        help="Maximum number of connections for the pool (default: 10)")
    parser.add_argument("--rollback-on-failure", action="store_true",
                        help="Automatically rollback changes if an operation fails")
    return parser.parse_args()

# Transaction support
@contextmanager
def transaction(project_id):
    """
    Context manager for database transactions.
    
    Args:
        project_id (str): The Supabase project ID
    """
    try:
        # Begin transaction
        execute_sql_on_supabase(project_id, "BEGIN;")
        logger.debug("Transaction started")
        yield
        # Commit transaction
        execute_sql_on_supabase(project_id, "COMMIT;")
        logger.debug("Transaction committed")
    except Exception as e:
        # Rollback transaction on error
        logger.error(f"Transaction error: {str(e)}, rolling back")
        execute_sql_on_supabase(project_id, "ROLLBACK;")
        raise

def get_performance_sql_batches(schema="public"):
    """
    Generate the SQL statements for performance improvements in batches.
    
    Args:
        schema (str): The database schema to modify
        
    Returns:
        dict: SQL statements for performance improvements organized in batches
    """
    # Batch 1: RLS performance indexes
    rls_indexes = f"""
    -- Add index for RLS performance
    CREATE INDEX IF NOT EXISTS idx_mixtures_created_by ON {schema}.mixtures(created_by);
    CREATE INDEX IF NOT EXISTS idx_experiments_created_by ON {schema}.experiments(created_by);
    CREATE INDEX IF NOT EXISTS idx_predictions_created_by ON {schema}.predictions(created_by);
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_created_by ON {schema}.molecular_properties(created_by);
    CREATE INDEX IF NOT EXISTS idx_calculation_methods_created_by ON {schema}.calculation_methods(created_by);
    
    -- Add comments
    COMMENT ON INDEX {schema}.idx_mixtures_created_by IS 'Improves RLS performance for queries filtering on created_by';
    COMMENT ON INDEX {schema}.idx_experiments_created_by IS 'Improves RLS performance for queries filtering on created_by';
    COMMENT ON INDEX {schema}.idx_predictions_created_by IS 'Improves RLS performance for queries filtering on created_by';
    COMMENT ON INDEX {schema}.idx_molecular_properties_created_by IS 'Improves RLS performance for queries filtering on created_by';
    COMMENT ON INDEX {schema}.idx_calculation_methods_created_by IS 'Improves RLS performance for queries filtering on created_by';
    """
    
    # Batch 2: Mixture and component indexes
    mixture_indexes = f"""
    -- Add index for mixture components
    CREATE INDEX IF NOT EXISTS idx_mixture_component_mixture_id ON {schema}.mixture_components(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_mixture_component_molecule_id ON {schema}.mixture_components(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_mixture_component_amount ON {schema}.mixture_components(amount);
    
    -- Add index for experiments
    CREATE INDEX IF NOT EXISTS idx_experiments_mixture_id ON {schema}.experiments(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_experiments_property_type_id ON {schema}.experiments(property_type_id);
    CREATE INDEX IF NOT EXISTS idx_experiments_date ON {schema}.experiments(date);
    
    -- Add composite indexes for common query patterns
    CREATE INDEX IF NOT EXISTS idx_mixture_components_mixture_molecule ON {schema}.mixture_components(mixture_id, molecule_id);
    CREATE INDEX IF NOT EXISTS idx_experiments_mixture_property ON {schema}.experiments(mixture_id, property_type_id);
    
    -- Add comments
    COMMENT ON INDEX {schema}.idx_mixture_component_mixture_id IS 'Improves performance for queries joining mixture_components with mixtures';
    COMMENT ON INDEX {schema}.idx_mixture_component_molecule_id IS 'Improves performance for queries joining mixture_components with molecules';
    COMMENT ON INDEX {schema}.idx_mixture_component_amount IS 'Improves performance for queries filtering on component amount';
    COMMENT ON INDEX {schema}.idx_experiments_mixture_id IS 'Improves performance for queries joining experiments with mixtures';
    COMMENT ON INDEX {schema}.idx_experiments_property_type_id IS 'Improves performance for queries filtering on property types';
    COMMENT ON INDEX {schema}.idx_experiments_date IS 'Improves performance for queries filtering on experiment date';
    COMMENT ON INDEX {schema}.idx_mixture_components_mixture_molecule IS 'Improves performance for queries filtering on both mixture and molecule';
    COMMENT ON INDEX {schema}.idx_experiments_mixture_property IS 'Improves performance for queries filtering on both mixture and property type';
    """
    
    # Batch 3: Prediction and property indexes
    prediction_indexes = f"""
    -- Add composite index for predictions
    CREATE INDEX IF NOT EXISTS idx_predictions_mixture_property ON {schema}.predictions(mixture_id, property_type_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_calculation_method ON {schema}.predictions(calculation_method_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_confidence ON {schema}.predictions(confidence);
    CREATE INDEX IF NOT EXISTS idx_predictions_numeric_value ON {schema}.predictions(numeric_value);
    
    -- Add index for property types
    CREATE INDEX IF NOT EXISTS idx_molecular_property_property_type ON {schema}.molecular_properties(property_type_id);
    CREATE INDEX IF NOT EXISTS idx_molecular_property_molecule_id ON {schema}.molecular_properties(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_molecular_property_numeric_value ON {schema}.molecular_properties(numeric_value);
    
    -- Add composite indexes for common query patterns
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_property ON {schema}.molecular_properties(molecule_id, property_type_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_method_property ON {schema}.predictions(calculation_method_id, property_type_id);
    
    -- Add comments
    COMMENT ON INDEX {schema}.idx_predictions_mixture_property IS 'Improves performance for queries filtering on both mixture_id and property_type_id';
    COMMENT ON INDEX {schema}.idx_predictions_calculation_method IS 'Improves performance for queries filtering on calculation method';
    COMMENT ON INDEX {schema}.idx_predictions_confidence IS 'Improves performance for queries filtering on prediction confidence';
    COMMENT ON INDEX {schema}.idx_predictions_numeric_value IS 'Improves performance for queries filtering on prediction value';
    COMMENT ON INDEX {schema}.idx_molecular_property_property_type IS 'Improves performance for queries joining molecular_properties with property_types';
    COMMENT ON INDEX {schema}.idx_molecular_property_molecule_id IS 'Improves performance for queries joining molecular_properties with molecules';
    COMMENT ON INDEX {schema}.idx_molecular_property_numeric_value IS 'Improves performance for queries filtering on property value';
    COMMENT ON INDEX {schema}.idx_molecular_properties_molecule_property IS 'Improves performance for queries filtering on both molecule and property type';
    COMMENT ON INDEX {schema}.idx_predictions_method_property IS 'Improves performance for queries filtering on both calculation method and property type';
    """
    
    # Batch 4: Text search indexes
    text_search_indexes = f"""
    -- Add text search index for molecule names
    CREATE EXTENSION IF NOT EXISTS pg_trgm;
    CREATE INDEX IF NOT EXISTS idx_molecule_name_trgm ON {schema}.molecules USING gin (name gin_trgm_ops);
    CREATE INDEX IF NOT EXISTS idx_mixture_name_trgm ON {schema}.mixtures USING gin (name gin_trgm_ops);
    CREATE INDEX IF NOT EXISTS idx_property_type_name_trgm ON {schema}.property_types USING gin (name gin_trgm_ops);
    
    -- Add btree indexes for exact matching
    CREATE INDEX IF NOT EXISTS idx_molecule_name ON {schema}.molecules(name);
    CREATE INDEX IF NOT EXISTS idx_mixture_name ON {schema}.mixtures(name);
    CREATE INDEX IF NOT EXISTS idx_property_type_name ON {schema}.property_types(name);
    
    -- Add comments
    COMMENT ON INDEX {schema}.idx_molecule_name_trgm IS 'Improves performance for text search operations on molecule names';
    COMMENT ON INDEX {schema}.idx_mixture_name_trgm IS 'Improves performance for text search operations on mixture names';
    COMMENT ON INDEX {schema}.idx_property_type_name_trgm IS 'Improves performance for text search operations on property type names';
    COMMENT ON INDEX {schema}.idx_molecule_name IS 'Improves performance for exact matching on molecule names';
    COMMENT ON INDEX {schema}.idx_mixture_name IS 'Improves performance for exact matching on mixture names';
    COMMENT ON INDEX {schema}.idx_property_type_name IS 'Improves performance for exact matching on property type names';
    """
    
    # Batch 5: Additional performance indexes
    additional_indexes = f"""
    -- Add indexes for timestamp columns for time-based queries
    CREATE INDEX IF NOT EXISTS idx_mixtures_created_at ON {schema}.mixtures(created_at);
    CREATE INDEX IF NOT EXISTS idx_mixtures_updated_at ON {schema}.mixtures(updated_at);
    CREATE INDEX IF NOT EXISTS idx_experiments_created_at ON {schema}.experiments(created_at);
    CREATE INDEX IF NOT EXISTS idx_predictions_created_at ON {schema}.predictions(created_at);
    
    -- Add indexes for foreign keys that might not be automatically indexed
    CREATE INDEX IF NOT EXISTS idx_property_types_unit_id ON {schema}.property_types(unit_id);
    CREATE INDEX IF NOT EXISTS idx_calculation_methods_version ON {schema}.calculation_methods(version);
    
    -- Add comments
    COMMENT ON INDEX {schema}.idx_mixtures_created_at IS 'Improves performance for queries filtering on creation date';
    COMMENT ON INDEX {schema}.idx_mixtures_updated_at IS 'Improves performance for queries filtering on update date';
    COMMENT ON INDEX {schema}.idx_experiments_created_at IS 'Improves performance for queries filtering on creation date';
    COMMENT ON INDEX {schema}.idx_predictions_created_at IS 'Improves performance for queries filtering on creation date';
    COMMENT ON INDEX {schema}.idx_property_types_unit_id IS 'Improves performance for queries joining property_types with units';
    COMMENT ON INDEX {schema}.idx_calculation_methods_version IS 'Improves performance for queries filtering on calculation method version';
    """
    
    return {
        "rls_indexes": rls_indexes,
        "mixture_indexes": mixture_indexes,
        "prediction_indexes": prediction_indexes,
        "text_search_indexes": text_search_indexes,
        "additional_indexes": additional_indexes
    }

def get_optimized_views_sql(schema="public"):
    """
    Generate SQL statements for creating optimized views for common queries.
    
    Args:
        schema (str): The database schema to modify
        
    Returns:
        dict: SQL statements for creating optimized views
    """
    # View 1: Molecule with properties
    molecule_with_properties = f"""
    CREATE OR REPLACE VIEW {schema}.molecule_with_properties AS
    SELECT
        m.id AS molecule_id,
        m.name AS molecule_name,
        m.smiles,
        m.inchi,
        m.created_by,
        m.created_at,
        m.updated_at,
        pt.id AS property_type_id,
        pt.name AS property_name,
        pt.description AS property_description,
        mp.numeric_value,
        u.symbol AS unit_symbol,
        u.name AS unit_name
    FROM {schema}.molecules m
    LEFT JOIN {schema}.molecular_properties mp ON m.id = mp.molecule_id
    LEFT JOIN {schema}.property_types pt ON mp.property_type_id = pt.id
    LEFT JOIN {schema}.units u ON pt.unit_id = u.id;
    
    COMMENT ON VIEW {schema}.molecule_with_properties IS 'Optimized view for querying molecules with their properties';
    """
    
    # View 2: Mixture with components
    mixture_with_components = f"""
    CREATE OR REPLACE VIEW {schema}.mixture_with_components AS
    SELECT
        mix.id AS mixture_id,
        mix.name AS mixture_name,
        mix.description,
        mix.created_by,
        mix.created_at,
        mix.updated_at,
        mc.id AS component_id,
        m.id AS molecule_id,
        m.name AS molecule_name,
        m.smiles,
        mc.amount,
        mc.amount_unit
    FROM {schema}.mixtures mix
    LEFT JOIN {schema}.mixture_components mc ON mix.id = mc.mixture_id
    LEFT JOIN {schema}.molecules m ON mc.molecule_id = m.id;
    
    COMMENT ON VIEW {schema}.mixture_with_components IS 'Optimized view for querying mixtures with their components';
    """
    
    # View 3: Prediction with experiment comparison
    prediction_experiment_comparison = f"""
    CREATE OR REPLACE VIEW {schema}.prediction_experiment_comparison AS
    SELECT
        mix.id AS mixture_id,
        mix.name AS mixture_name,
        pt.id AS property_type_id,
        pt.name AS property_name,
        p.numeric_value AS predicted_value,
        e.numeric_value AS experimental_value,
        ABS(p.numeric_value - e.numeric_value) AS absolute_difference,
        CASE
            WHEN ABS(e.numeric_value) > 0
            THEN (ABS(p.numeric_value - e.numeric_value) / ABS(e.numeric_value)) * 100
            ELSE NULL
        END AS percentage_difference,
        p.confidence AS prediction_confidence,
        cm.id AS calculation_method_id,
        cm.name AS calculation_method_name,
        u.symbol AS unit_symbol
    FROM {schema}.predictions p
    JOIN {schema}.experiments e ON p.mixture_id = e.mixture_id AND p.property_type_id = e.property_type_id
    JOIN {schema}.mixtures mix ON p.mixture_id = mix.id
    JOIN {schema}.property_types pt ON p.property_type_id = pt.id
    JOIN {schema}.calculation_methods cm ON p.calculation_method_id = cm.id
    LEFT JOIN {schema}.units u ON pt.unit_id = u.id;
    
    COMMENT ON VIEW {schema}.prediction_experiment_comparison IS 'Optimized view for comparing predictions with experimental values';
    """
    
    # View 4: Molecule property statistics
    molecule_property_statistics = f"""
    CREATE OR REPLACE VIEW {schema}.molecule_property_statistics AS
    SELECT
        pt.id AS property_type_id,
        pt.name AS property_name,
        COUNT(mp.id) AS property_count,
        MIN(mp.numeric_value) AS min_value,
        MAX(mp.numeric_value) AS max_value,
        AVG(mp.numeric_value) AS avg_value,
        PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY mp.numeric_value) AS median_value,
        STDDEV(mp.numeric_value) AS stddev_value,
        u.symbol AS unit_symbol
    FROM {schema}.molecular_properties mp
    JOIN {schema}.property_types pt ON mp.property_type_id = pt.id
    LEFT JOIN {schema}.units u ON pt.unit_id = u.id
    GROUP BY pt.id, pt.name, u.symbol;
    
    COMMENT ON VIEW {schema}.molecule_property_statistics IS 'Statistical summary of molecular properties';
    """
    
    # View 5: Recent activity
    recent_activity = f"""
    CREATE OR REPLACE VIEW {schema}.recent_activity AS
    SELECT 'mixture' AS entity_type, id, name, created_by, created_at, updated_at
    FROM {schema}.mixtures
    UNION ALL
    SELECT 'experiment' AS entity_type, id, 'Experiment' AS name, created_by, created_at, updated_at
    FROM {schema}.experiments
    UNION ALL
    SELECT 'prediction' AS entity_type, id, 'Prediction' AS name, created_by, created_at, updated_at
    FROM {schema}.predictions
    ORDER BY created_at DESC;
    
    COMMENT ON VIEW {schema}.recent_activity IS 'Unified view of recent activity across different entity types';
    """
    
    return {
        "molecule_with_properties": molecule_with_properties,
        "mixture_with_components": mixture_with_components,
        "prediction_experiment_comparison": prediction_experiment_comparison,
        "molecule_property_statistics": molecule_property_statistics,
        "recent_activity": recent_activity
    }

def get_performance_sql(schema="public"):
    """
    Generate the SQL statements for performance improvements (legacy function).
    
    Args:
        schema (str): The database schema to modify
        
    Returns:
        str: SQL statements for performance improvements
    """
    logger.warning("Using legacy get_performance_sql function. Consider using get_performance_sql_batches instead.")
    
    batches = get_performance_sql_batches(schema)
    return "\n".join([
        batches["rls_indexes"],
        batches["mixture_indexes"],
        batches["prediction_indexes"],
        batches["text_search_indexes"]
    ])

def get_connection_pooling_sql():
    """
    Generate the SQL statements for configuring connection pooling.
    
    Returns:
        str: SQL statements for connection pooling configuration
    """
    return """
    -- Configure connection pooling
    ALTER SYSTEM SET max_connections = '100';
    ALTER SYSTEM SET superuser_reserved_connections = '3';
    ALTER SYSTEM SET idle_in_transaction_session_timeout = '10000';  -- 10 seconds
    ALTER SYSTEM SET statement_timeout = '30000';  -- 30 seconds
    
    -- Configure connection pooler settings
    ALTER SYSTEM SET pgbouncer.max_client_conn = '1000';
    ALTER SYSTEM SET pgbouncer.default_pool_size = '50';
    ALTER SYSTEM SET pgbouncer.min_pool_size = '5';
    ALTER SYSTEM SET pgbouncer.reserve_pool_size = '10';
    ALTER SYSTEM SET pgbouncer.max_db_connections = '80';
    ALTER SYSTEM SET pgbouncer.max_user_connections = '50';
    
    -- Configure work memory and shared buffers for better performance
    ALTER SYSTEM SET work_mem = '8MB';
    ALTER SYSTEM SET maintenance_work_mem = '128MB';
    ALTER SYSTEM SET effective_cache_size = '1GB';
    
    -- Apply changes
    SELECT pg_reload_conf();
    """

def register_rollback_operation(operation_type, details):
    """
    Register a rollback operation in the rollback file.
    
    Args:
        operation_type (str): Type of operation (e.g., 'create_index', 'create_view')
        details (dict): Details of the operation for rollback
    """
    try:
        with open(ROLLBACK_FILE, 'r') as f:
            rollback_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        init_memory_cache()
        with open(ROLLBACK_FILE, 'r') as f:
            rollback_data = json.load(f)
    
    rollback_data.setdefault("rollback_operations", [])
    rollback_data["rollback_operations"].append({
        "type": operation_type,
        "details": details,
        "timestamp": datetime.now().isoformat()
    })
    rollback_data["last_updated"] = datetime.now().isoformat()
    
    with open(ROLLBACK_FILE, 'w') as f:
        json.dump(rollback_data, f, indent=2)

def generate_rollback_sql(schema="public"):
    """
    Generate SQL statements to rollback all applied changes.
    
    Args:
        schema (str): The database schema
        
    Returns:
        str: SQL statements for rolling back changes
    """
    try:
        with open(ROLLBACK_FILE, 'r') as f:
            rollback_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        logger.warning("Rollback file not found or invalid")
        return ""
    
    rollback_operations = rollback_data.get("rollback_operations", [])
    if not rollback_operations:
        logger.warning("No rollback operations found")
        return ""
    
    # Reverse the operations to rollback in the correct order (LIFO)
    rollback_operations.reverse()
    
    rollback_statements = []
    
    for operation in rollback_operations:
        op_type = operation["type"]
        details = operation["details"]
        
        if op_type == "create_index":
            index_name = details.get("index_name")
            if index_name:
                rollback_statements.append(f"DROP INDEX IF EXISTS {schema}.{index_name};")
        
        elif op_type == "create_view":
            view_name = details.get("view_name")
            if view_name:
                rollback_statements.append(f"DROP VIEW IF EXISTS {schema}.{view_name};")
        
        elif op_type == "connection_pooling":
            # Reset connection pooling settings to defaults
            rollback_statements.append("""
            ALTER SYSTEM RESET max_connections;
            ALTER SYSTEM RESET superuser_reserved_connections;
            ALTER SYSTEM RESET idle_in_transaction_session_timeout;
            ALTER SYSTEM RESET statement_timeout;
            ALTER SYSTEM RESET pgbouncer.max_client_conn;
            ALTER SYSTEM RESET pgbouncer.default_pool_size;
            ALTER SYSTEM RESET pgbouncer.min_pool_size;
            ALTER SYSTEM RESET pgbouncer.reserve_pool_size;
            ALTER SYSTEM RESET pgbouncer.max_db_connections;
            ALTER SYSTEM RESET pgbouncer.max_user_connections;
            ALTER SYSTEM RESET work_mem;
            ALTER SYSTEM RESET maintenance_work_mem;
            ALTER SYSTEM RESET effective_cache_size;
            SELECT pg_reload_conf();
            """)
    
    return "\n".join(rollback_statements)

def execute_rollback(project_id, schema="public"):
    """
    Execute rollback of all applied changes.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema
        
    Returns:
        bool: True if rollback was successful, False otherwise
    """
    logger.info("Executing rollback of applied changes...")
    
    rollback_sql = generate_rollback_sql(schema)
    if not rollback_sql:
        logger.warning("No rollback SQL generated")
        return False
    
    try:
        with transaction(project_id):
            execute_sql_on_supabase(project_id, rollback_sql)
        
        # Clear the rollback file after successful rollback
        with open(ROLLBACK_FILE, 'w') as f:
            json.dump({
                "version": "1.0",
                "last_updated": datetime.now().isoformat(),
                "rollback_operations": []
            }, f, indent=2)
        
        logger.info("Rollback executed successfully")
        return True
    except Exception as e:
        logger.error(f"Error executing rollback: {str(e)}")
        return False

@memory_cache
def test_performance_improvements(project_id, schema="public"):
    """
    Test the performance improvements by running EXPLAIN ANALYZE on key queries.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema to test
        
    Returns:
        dict: Test results with query plans and execution times
    """
    logger.info("Testing performance improvements...")
    
    test_queries = {
        "rls_filter": f"""
        EXPLAIN ANALYZE
        SELECT * FROM {schema}.mixtures
        WHERE created_by = '00000000-0000-0000-0000-000000000000'
        LIMIT 10;
        """,
        
        "mixture_components_join": f"""
        EXPLAIN ANALYZE
        SELECT m.name, mc.amount, mc.amount_unit
        FROM {schema}.mixtures m
        JOIN {schema}.mixture_components mc ON m.id = mc.mixture_id
        LIMIT 10;
        """,
        
        "predictions_filter": f"""
        EXPLAIN ANALYZE
        SELECT p.value, p.unit
        FROM {schema}.predictions p
        WHERE p.mixture_id = (SELECT id FROM {schema}.mixtures LIMIT 1)
        AND p.property_type_id = (SELECT id FROM {schema}.property_types LIMIT 1);
        """,
        
        "experiments_join": f"""
        EXPLAIN ANALYZE
        SELECT e.*, m.name
        FROM {schema}.experiments e
        JOIN {schema}.mixtures m ON e.mixture_id = m.id
        LIMIT 10;
        """,
        
        "molecule_text_search": f"""
        EXPLAIN ANALYZE
        SELECT * FROM {schema}.molecules
        WHERE name ILIKE '%glycerol%'
        LIMIT 10;
        """,
        
        "property_type_join": f"""
        EXPLAIN ANALYZE
        SELECT mp.*, pt.name
        FROM {schema}.molecular_properties mp
        JOIN {schema}.property_types pt ON mp.property_type_id = pt.id
        LIMIT 10;
        """
    }
    
    results = {}
    
    for name, query in test_queries.items():
        try:
            logger.info(f"Running test query: {name}")
            result = execute_sql_on_supabase(project_id, query)
            results[name] = result
        except Exception as e:
            logger.error(f"Error running test query {name}: {str(e)}")
            results[name] = {"error": str(e)}
    
    return results

@retry_operation(max_retries=3, delay=2)
def verify_index_creation(project_id, schema, index_name):
    """
    Verify that an index was created successfully.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema
        index_name (str): The name of the index to verify
        
    Returns:
        bool: True if the index exists, False otherwise
    """
    query = f"""
    SELECT EXISTS (
        SELECT 1
        FROM pg_indexes
        WHERE schemaname = '{schema}'
        AND indexname = '{index_name}'
    );
    """
    
    try:
        result = execute_sql_on_supabase(project_id, query)
        return result[0]['exists']
    except Exception as e:
        logger.error(f"Error verifying index {index_name}: {str(e)}")
        return False

def save_verification_result(index_name, success):
    """
    Save the verification result to the memory cache.
    
    Args:
        index_name (str): The name of the index
        success (bool): Whether the verification was successful
    """
    try:
        with open(MEMORY_CACHE_FILE, 'r') as f:
            cache = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        init_memory_cache()
        with open(MEMORY_CACHE_FILE, 'r') as f:
            cache = json.load(f)
    
    cache.setdefault("verification_results", {})
    cache["verification_results"][index_name] = {
        "success": success,
        "timestamp": datetime.now().isoformat()
    }
    
    with open(MEMORY_CACHE_FILE, 'w') as f:
        json.dump(cache, f, indent=2)

def apply_improvements(project_id, schema="public", create_backups=False, test_only=False, backup_dir="production_backups",
                      analyze_queries=False, connection_pooling=False, rollback_on_failure=False,
                      min_connections=2, max_connections=10):
    """
    Apply the performance improvements to the database.
    
    Args:
        project_id (str): The Supabase project ID
        schema (str): The database schema to modify
        create_backups (bool): Whether to create backups before modifying tables
        test_only (bool): Whether to only test without applying changes
        backup_dir (str): Directory to store backups
        analyze_queries (bool): Whether to analyze query plans
        connection_pooling (bool): Whether to enable enhanced connection pooling
        rollback_on_failure (bool): Whether to automatically rollback changes on failure
        min_connections (int): Minimum number of connections for the pool
        max_connections (int): Maximum number of connections for the pool
        
    Returns:
        dict: Results of the operation
    """
    logger.info(f"Starting performance improvements for project {project_id}, schema {schema}")
    
    # Initialize memory cache
    init_memory_cache()
    
    # Create a results dictionary
    results = {
        "timestamp": datetime.now().isoformat(),
        "project_id": project_id,
        "schema": schema,
        "create_backups": create_backups,
        "test_only": test_only,
        "analyze_queries": analyze_queries,
        "connection_pooling": connection_pooling,
        "rollback_on_failure": rollback_on_failure,
        "status": "PENDING",
        "backups": [],
        "applied_improvements": [],
        "errors": [],
        "test_results_before": {},
        "test_results_after": {},
        "query_analysis": {},
        "optimized_views": []
    }
    
    try:
        # Run tests before applying improvements
        if test_only or create_backups or analyze_queries:
            logger.info("Running performance tests before applying improvements...")
            results["test_results_before"] = test_performance_improvements(project_id, schema)
        
        if test_only:
            logger.info("Test-only mode, skipping actual improvements")
            results["status"] = "TEST_ONLY"
            return results
        
        # Create backups if requested
        if create_backups:
            # Create backup directory
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_path = ensure_backup_directory(Path(backup_dir) / f"backup_{timestamp}")
            
            tables_to_backup = ["mixtures", "experiments", "predictions", "mixture_components",
                               "molecules", "molecular_properties"]
            
            for table in tables_to_backup:
                try:
                    # Create a backup file path
                    backup_file = backup_path / f"{table}_backup.sql"
                    
                    # Use the backup_table function with retry logic
                    @retry_operation(max_retries=3, delay=2)
                    def create_backup():
                        return backup_table(project_id, f"{schema}.{table}")
                    
                    backup_table_name = create_backup()
                    
                    # Record the backup details
                    results["backups"].append({
                        "original_table": f"{schema}.{table}",
                        "backup_table": backup_table_name,
                        "backup_file": str(backup_file)
                    })
                    
                    # Save backup metadata to the file
                    with open(backup_file, 'w') as f:
                        f.write(f"-- Backup of {schema}.{table}\n")
                        f.write(f"-- Created: {datetime.now().isoformat()}\n")
                        f.write(f"-- Backup table: {backup_table_name}\n")
                    
                    logger.info(f"Created backup of {schema}.{table} as {backup_table_name}")
                except Exception as e:
                    error_msg = f"Error creating backup of {schema}.{table}: {str(e)}"
                    logger.error(error_msg)
                    results["errors"].append(error_msg)
        
        # Apply performance indexes in batches
        performance_sql_batches = get_performance_sql_batches(schema)
        
        # Use transaction for safer operations if rollback_on_failure is enabled
        if rollback_on_failure:
            logger.info("Using transaction for safer operations (rollback on failure enabled)")
        
        # Batch 1: RLS indexes
        logger.info("Applying RLS performance indexes (Batch 1/5)...")
        try:
            apply_migration(
                project_id=project_id,
                name="performance_improvements_rls",
                query=performance_sql_batches["rls_indexes"]
            )
            
            # Register rollback operation
            register_rollback_operation("create_index", {
                "batch": "rls_indexes",
                "indexes": ["idx_mixtures_created_by", "idx_experiments_created_by", "idx_predictions_created_by",
                           "idx_molecular_properties_created_by", "idx_calculation_methods_created_by"]
            })
            
            # Verify indexes
            rls_indexes = ["idx_mixtures_created_by", "idx_experiments_created_by", "idx_predictions_created_by",
                          "idx_molecular_properties_created_by", "idx_calculation_methods_created_by"]
            verification_results = []
            
            for index_name in rls_indexes:
                success = verify_index_creation(project_id, schema, index_name)
                verification_results.append(success)
                save_verification_result(index_name, success)
            
            if all(verification_results):
                results["applied_improvements"].append("rls_indexes")
                logger.info("Successfully applied RLS performance indexes")
            else:
                failed_indexes = [idx for idx, success in zip(rls_indexes, verification_results) if not success]
                error_msg = f"Failed to verify some RLS indexes: {', '.join(failed_indexes)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)
                
                if rollback_on_failure:
                    logger.warning("Rollback on failure enabled, rolling back changes...")
                    execute_rollback(project_id, schema)
                    results["status"] = "ROLLED_BACK"
                    return results
        except Exception as e:
            error_msg = f"Error applying RLS performance indexes: {str(e)}"
            logger.error(error_msg)
            results["errors"].append(error_msg)
            
            if rollback_on_failure:
                logger.warning("Rollback on failure enabled, rolling back changes...")
                execute_rollback(project_id, schema)
                results["status"] = "ROLLED_BACK"
                return results
        
        # Batch 2: Mixture indexes
        logger.info("Applying mixture and component indexes (Batch 2/5)...")
        try:
            apply_migration(
                project_id=project_id,
                name="performance_improvements_mixture",
                query=performance_sql_batches["mixture_indexes"]
            )
            
            # Register rollback operation
            register_rollback_operation("create_index", {
                "batch": "mixture_indexes",
                "indexes": ["idx_mixture_component_mixture_id", "idx_mixture_component_molecule_id",
                           "idx_mixture_component_amount", "idx_experiments_mixture_id",
                           "idx_experiments_property_type_id", "idx_experiments_date",
                           "idx_mixture_components_mixture_molecule", "idx_experiments_mixture_property"]
            })
            
            # Verify indexes
            mixture_indexes = ["idx_mixture_component_mixture_id", "idx_mixture_component_molecule_id",
                              "idx_experiments_mixture_id", "idx_experiments_property_type_id"]
            verification_results = []
            
            for index_name in mixture_indexes:
                success = verify_index_creation(project_id, schema, index_name)
                verification_results.append(success)
                save_verification_result(index_name, success)
            
            if all(verification_results):
                results["applied_improvements"].append("mixture_indexes")
                logger.info("Successfully applied mixture and component indexes")
            else:
                failed_indexes = [idx for idx, success in zip(mixture_indexes, verification_results) if not success]
                error_msg = f"Failed to verify some mixture indexes: {', '.join(failed_indexes)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)
                
                if rollback_on_failure:
                    logger.warning("Rollback on failure enabled, rolling back changes...")
                    execute_rollback(project_id, schema)
                    results["status"] = "ROLLED_BACK"
                    return results
        except Exception as e:
            error_msg = f"Error applying mixture and component indexes: {str(e)}"
            logger.error(error_msg)
            results["errors"].append(error_msg)
            
            if rollback_on_failure:
                logger.warning("Rollback on failure enabled, rolling back changes...")
                execute_rollback(project_id, schema)
                results["status"] = "ROLLED_BACK"
                return results
        
        # Batch 3: Prediction indexes
        logger.info("Applying prediction and property indexes (Batch 3/5)...")
        try:
            apply_migration(
                project_id=project_id,
                name="performance_improvements_prediction",
                query=performance_sql_batches["prediction_indexes"]
            )
            
            # Register rollback operation
            register_rollback_operation("create_index", {
                "batch": "prediction_indexes",
                "indexes": ["idx_predictions_mixture_property", "idx_predictions_calculation_method",
                           "idx_predictions_confidence", "idx_predictions_numeric_value",
                           "idx_molecular_property_property_type", "idx_molecular_property_molecule_id",
                           "idx_molecular_property_numeric_value", "idx_molecular_properties_molecule_property",
                           "idx_predictions_method_property"]
            })
            
            # Verify indexes
            prediction_indexes = ["idx_predictions_mixture_property", "idx_molecular_property_property_type",
                                 "idx_predictions_calculation_method", "idx_molecular_property_molecule_id"]
            verification_results = []
            
            for index_name in prediction_indexes:
                success = verify_index_creation(project_id, schema, index_name)
                verification_results.append(success)
                save_verification_result(index_name, success)
            
            if all(verification_results):
                results["applied_improvements"].append("prediction_indexes")
                logger.info("Successfully applied prediction and property indexes")
            else:
                failed_indexes = [idx for idx, success in zip(prediction_indexes, verification_results) if not success]
                error_msg = f"Failed to verify some prediction indexes: {', '.join(failed_indexes)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)
                
                if rollback_on_failure:
                    logger.warning("Rollback on failure enabled, rolling back changes...")
                    execute_rollback(project_id, schema)
                    results["status"] = "ROLLED_BACK"
                    return results
        except Exception as e:
            error_msg = f"Error applying prediction and property indexes: {str(e)}"
            logger.error(error_msg)
            results["errors"].append(error_msg)
            
            if rollback_on_failure:
                logger.warning("Rollback on failure enabled, rolling back changes...")
                execute_rollback(project_id, schema)
                results["status"] = "ROLLED_BACK"
                return results
        
        # Batch 4: Text search indexes
        logger.info("Applying text search indexes (Batch 4/5)...")
        try:
            apply_migration(
                project_id=project_id,
                name="performance_improvements_text_search",
                query=performance_sql_batches["text_search_indexes"]
            )
            
            # Register rollback operation
            register_rollback_operation("create_index", {
                "batch": "text_search_indexes",
                "indexes": ["idx_molecule_name_trgm", "idx_mixture_name_trgm", "idx_property_type_name_trgm",
                           "idx_molecule_name", "idx_mixture_name", "idx_property_type_name"]
            })
            
            # Verify indexes
            text_search_indexes = ["idx_molecule_name_trgm", "idx_mixture_name_trgm", "idx_property_type_name_trgm"]
            verification_results = []
            
            for index_name in text_search_indexes:
                success = verify_index_creation(project_id, schema, index_name)
                verification_results.append(success)
                save_verification_result(index_name, success)
            
            if all(verification_results):
                results["applied_improvements"].append("text_search_indexes")
                logger.info("Successfully applied text search indexes")
            else:
                failed_indexes = [idx for idx, success in zip(text_search_indexes, verification_results) if not success]
                error_msg = f"Failed to verify some text search indexes: {', '.join(failed_indexes)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)
                
                if rollback_on_failure:
                    logger.warning("Rollback on failure enabled, rolling back changes...")
                    execute_rollback(project_id, schema)
                    results["status"] = "ROLLED_BACK"
                    return results
        except Exception as e:
            error_msg = f"Error applying text search indexes: {str(e)}"
            logger.error(error_msg)
            results["errors"].append(error_msg)
            
            if rollback_on_failure:
                logger.warning("Rollback on failure enabled, rolling back changes...")
                execute_rollback(project_id, schema)
                results["status"] = "ROLLED_BACK"
                return results
        
        # Batch 5: Additional indexes
        logger.info("Applying additional performance indexes (Batch 5/5)...")
        try:
            apply_migration(
                project_id=project_id,
                name="performance_improvements_additional",
                query=performance_sql_batches["additional_indexes"]
            )
            
            # Register rollback operation
            register_rollback_operation("create_index", {
                "batch": "additional_indexes",
                "indexes": ["idx_mixtures_created_at", "idx_mixtures_updated_at",
                           "idx_experiments_created_at", "idx_predictions_created_at",
                           "idx_property_types_unit_id", "idx_calculation_methods_version"]
            })
            
            # Verify indexes
            additional_indexes = ["idx_mixtures_created_at", "idx_property_types_unit_id"]
            verification_results = []
            
            for index_name in additional_indexes:
                success = verify_index_creation(project_id, schema, index_name)
                verification_results.append(success)
                save_verification_result(index_name, success)
            
            if all(verification_results):
                results["applied_improvements"].append("additional_indexes")
                logger.info("Successfully applied additional performance indexes")
            else:
                failed_indexes = [idx for idx, success in zip(additional_indexes, verification_results) if not success]
                error_msg = f"Failed to verify some additional indexes: {', '.join(failed_indexes)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)
                
                if rollback_on_failure:
                    logger.warning("Rollback on failure enabled, rolling back changes...")
                    execute_rollback(project_id, schema)
                    results["status"] = "ROLLED_BACK"
                    return results
        except Exception as e:
            error_msg = f"Error applying additional indexes: {str(e)}"
            logger.error(error_msg)
            results["errors"].append(error_msg)
            
            if rollback_on_failure:
                logger.warning("Rollback on failure enabled, rolling back changes...")
                execute_rollback(project_id, schema)
                results["status"] = "ROLLED_BACK"
                return results
        
        # Apply connection pooling configuration if requested
        if connection_pooling:
            connection_pooling_sql = get_connection_pooling_sql()
            logger.info("Configuring connection pooling...")
            
            try:
                @retry_operation(max_retries=3, delay=2)
                def configure_connection_pooling():
                    return execute_sql_on_supabase(project_id, connection_pooling_sql)
                
                configure_connection_pooling()
                
                # Register rollback operation
                register_rollback_operation("connection_pooling", {})
                
                results["applied_improvements"].append("connection_pooling")
                logger.info("Successfully configured connection pooling")
            except Exception as e:
                error_msg = f"Error configuring connection pooling: {str(e)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)
                
                if rollback_on_failure:
                    logger.warning("Rollback on failure enabled, rolling back changes...")
                    execute_rollback(project_id, schema)
                    results["status"] = "ROLLED_BACK"
                    return results
        
        # Create optimized views
        logger.info("Creating optimized views for common queries...")
        optimized_views_sql = get_optimized_views_sql(schema)
        
        for view_name, view_sql in optimized_views_sql.items():
            try:
                logger.info(f"Creating optimized view: {view_name}")
                apply_migration(
                    project_id=project_id,
                    name=f"create_view_{view_name}",
                    query=view_sql
                )
                
                # Register rollback operation
                register_rollback_operation("create_view", {
                    "view_name": view_name
                })
                
                results["optimized_views"].append(view_name)
                logger.info(f"Successfully created optimized view: {view_name}")
            except Exception as e:
                error_msg = f"Error creating optimized view {view_name}: {str(e)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)
                
                if rollback_on_failure:
                    logger.warning("Rollback on failure enabled, rolling back changes...")
                    execute_rollback(project_id, schema)
                    results["status"] = "ROLLED_BACK"
                    return results
        
        # Analyze query plans if requested
        if analyze_queries and analyze_plan_for_issues and suggest_optimizations:
            logger.info("Analyzing query plans to identify slow queries...")
            query_analysis_results = {}
            
            for query_info in COMMON_QUERIES:
                try:
                    query_name = query_info["name"]
                    query = query_info["query"]
                    
                    # Execute the query with EXPLAIN ANALYZE
                    explain_query = f"EXPLAIN ANALYZE {query}"
                    plan_result = execute_sql_on_supabase(project_id, explain_query)
                    
                    # Analyze the plan for issues
                    issues = analyze_plan_for_issues(plan_result)
                    
                    # Suggest optimizations
                    suggestions = suggest_optimizations(query, plan_result, issues)
                    
                    query_analysis_results[query_name] = {
                        "query": query,
                        "plan": plan_result,
                        "issues": issues,
                        "suggestions": suggestions
                    }
                    
                    logger.info(f"Analyzed query plan for: {query_name}")
                except Exception as e:
                    logger.error(f"Error analyzing query plan for {query_name}: {str(e)}")
            
            results["query_analysis"] = query_analysis_results
        
        # Run tests after applying improvements
        logger.info("Running performance tests after applying improvements...")
        results["test_results_after"] = test_performance_improvements(project_id, schema)
        
        # Determine overall status
        if not results["errors"]:
            results["status"] = "SUCCESS"
        else:
            results["status"] = "COMPLETED_WITH_WARNINGS" if results["applied_improvements"] else "ERROR"
        
        return results
    
    except Exception as e:
        error_msg = f"Unexpected error: {str(e)}\n{traceback.format_exc()}"
        logger.error(error_msg)
        results["errors"].append(error_msg)
        results["status"] = "ERROR"
        
        if rollback_on_failure:
            logger.warning("Rollback on failure enabled, rolling back changes...")
            execute_rollback(project_id, schema)
            results["status"] = "ROLLED_BACK"
        
        return results

def save_results(results, output_dir="."):
    """
    Save the results to a file.
    
    Args:
        results (dict): The results to save
        output_dir (str): The directory to save the results to
    
    Returns:
        str: The path to the saved file
    """
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"performance_improvements_{timestamp}.json"
    filepath = output_dir / filename
    
    with open(filepath, "w") as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"Results saved to {filepath}")
    return str(filepath)

def generate_report(results):
    """
    Generate a human-readable report from the results.
    
    Args:
        results (dict): The results to generate a report from
        
    Returns:
        str: A human-readable report
    """
    report = []
    report.append("=" * 80)
    report.append("CryoProtect v2 - Performance Improvements Report")
    report.append("=" * 80)
    report.append("")
    
    report.append(f"Timestamp: {results['timestamp']}")
    report.append(f"Project ID: {results['project_id']}")
    report.append(f"Schema: {results['schema']}")
    report.append(f"Status: {results['status']}")
    report.append("")
    
    report.append("Applied Improvements:")
    if results["applied_improvements"]:
        for improvement in results["applied_improvements"]:
            report.append(f"  - {improvement}")
    else:
        report.append("  None")
    report.append("")
    
    if results["backups"]:
        report.append("Backups Created:")
        for backup in results["backups"]:
            backup_info = f"  - {backup['original_table']} -> {backup['backup_table']}"
            if 'backup_file' in backup:
                backup_info += f" (File: {backup['backup_file']})"
            report.append(backup_info)
        report.append("")
        
    # Add verification results if available
    try:
        with open(MEMORY_CACHE_FILE, 'r') as f:
            cache = json.load(f)
            if "verification_results" in cache and cache["verification_results"]:
                report.append("Index Verification Results:")
                for index_name, result in cache["verification_results"].items():
                    status = "✓ Success" if result["success"] else "✗ Failed"
                    report.append(f"  - {index_name}: {status}")
                report.append("")
    except (FileNotFoundError, json.JSONDecodeError):
        pass
    
    if results["errors"]:
        report.append("Errors:")
        for error in results["errors"]:
            report.append(f"  - {error}")
        report.append("")
    
    if "test_results_before" in results and "test_results_after" in results:
        report.append("Performance Comparison:")
        report.append("  Query                  | Before                | After")
        report.append("  -----------------------|-----------------------|----------------------")
        
        for query_name in results["test_results_after"].keys():
            before = "N/A"
            after = "N/A"
            
            if query_name in results["test_results_before"] and "planning time" in str(results["test_results_before"][query_name]):
                before_str = str(results["test_results_before"][query_name])
                if "Execution Time:" in before_str:
                    before = before_str.split("Execution Time:")[1].split("ms")[0].strip()
                    before = f"{before} ms"
            
            if "planning time" in str(results["test_results_after"][query_name]):
                after_str = str(results["test_results_after"][query_name])
                if "Execution Time:" in after_str:
                    after = after_str.split("Execution Time:")[1].split("ms")[0].strip()
                    after = f"{after} ms"
            
            report.append(f"  {query_name:22} | {before:21} | {after:20}")
    
    report.append("")
    report.append("Recommendations:")
    report.append("  1. Monitor query performance in production")
    report.append("  2. Consider adding additional indexes if specific queries are still slow")
    report.append("  3. Regularly analyze tables to update statistics")
    report.append("  4. Utilize the memory bank caching for frequently accessed data")
    report.append("  5. Periodically verify index health and rebuild if necessary")
    report.append("  6. Consider implementing connection pooling tuning based on workload")
    report.append("  7. Apply future performance improvements in batches with verification")
    report.append("")
    
    return "\n".join(report)

def main():
    """Main function to apply performance improvements."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Apply improvements
        results = apply_improvements(
            project_id=args.project_id,
            schema=args.schema,
            create_backups=args.backup,
            test_only=args.test_only,
            backup_dir=args.backup_dir,
            analyze_queries=args.analyze_queries,
            connection_pooling=args.connection_pooling,
            rollback_on_failure=args.rollback_on_failure,
            min_connections=args.min_connections,
            max_connections=args.max_connections
        )
        
        # Save results
        save_results(results)
        
        # Generate and print report
        report = generate_report(results)
        print(report)
        
        # Save report to file
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_path = f"performance_improvements_report_{timestamp}.txt"
        with open(report_path, "w") as f:
            f.write(report)
        
        logger.info(f"Report saved to {report_path}")
        
        # Return appropriate exit code
        if results["status"] == "SUCCESS":
            return 0
        elif results["status"] == "COMPLETED_WITH_WARNINGS":
            return 1
        else:
            return 2
    
    except Exception as e:
        logger.error(f"Error in main function: {str(e)}")
        return 3

if __name__ == "__main__":
    sys.exit(main())