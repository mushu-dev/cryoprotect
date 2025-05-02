#!/usr/bin/env python3
"""
CryoProtect v2 - Patch Connection Pooling

This script patches the api/utils.py file to use the connection pool wrapper.
It creates a backup of the original file and applies the changes with minimal downtime.
"""

import os
import sys
import shutil
import logging
import argparse
import re
from datetime import datetime
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("patch_connection_pooling.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("patch_connection_pooling")

def create_backup(file_path):
    """
    Create a backup of the file.
    
    Args:
        file_path: Path to the file to backup
        
    Returns:
        Path to the backup file
    """
    backup_path = f"{file_path}.bak.{datetime.now().strftime('%Y%m%d%H%M%S')}"
    shutil.copy2(file_path, backup_path)
    logger.info(f"Created backup of {file_path} at {backup_path}")
    return backup_path

def patch_utils_file(utils_path, dry_run=False):
    """
    Patch the api/utils.py file to use the connection pool.
    
    Args:
        utils_path: Path to the api/utils.py file
        dry_run: If True, only show what would be changed without making changes
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Read the original file
        with open(utils_path, 'r') as f:
            content = f.read()
        
        # Create a backup
        if not dry_run:
            backup_path = create_backup(utils_path)
        
        # Modify the imports
        import_pattern = r"from supabase import create_client, Client"
        import_replacement = (
            "from supabase import create_client, Client\n"
            "# Import connection pool wrapper\n"
            "try:\n"
            "    from connection_pool_wrapper import (\n"
            "        initialize_supabase_pool, get_pooled_supabase_client,\n"
            "        release_supabase_client, get_supabase_pool_stats,\n"
            "        shutdown_supabase_pool, ConnectionPoolError\n"
            "    )\n"
            "    USE_CONNECTION_POOL = True\n"
            "    logger = logging.getLogger(__name__)\n"
            "    logger.info(\"Using connection pool for Supabase client\")\n"
            "except ImportError:\n"
            "    USE_CONNECTION_POOL = False\n"
            "    logger = logging.getLogger(__name__)\n"
            "    logger.warning(\"Connection pool wrapper not found, using direct connections\")"
        )
        
        content = re.sub(import_pattern, import_replacement, content)
        
        # Modify the get_supabase_client function
        client_pattern = r"def get_supabase_client\(\) -> Client:(.*?)return g\.supabase"
        client_replacement = (
            "def get_supabase_client() -> Client:\n"
            "    \"\"\"\n"
            "    Get or create a Supabase client.\n"
            "    \n"
            "    Returns:\n"
            "        Client: A Supabase client instance\n"
            "    \"\"\"\n"
            "    # Check if we're using the connection pool\n"
            "    if USE_CONNECTION_POOL:\n"
            "        # Initialize the pool if it's not already initialized\n"
            "        if not hasattr(g, 'supabase_pool_initialized'):\n"
            "            supabase_url = current_app.config['SUPABASE_URL']\n"
            "            supabase_key = current_app.config['SUPABASE_KEY']\n"
            "            \n"
            "            if not supabase_url or not supabase_key:\n"
            "                raise ValueError(\"SUPABASE_URL and SUPABASE_KEY must be set in configuration\")\n"
            "            \n"
            "            # Get pool settings from config or use defaults\n"
            "            min_connections = current_app.config.get('SUPABASE_MIN_CONNECTIONS', 2)\n"
            "            max_connections = current_app.config.get('SUPABASE_MAX_CONNECTIONS', 10)\n"
            "            connection_timeout = current_app.config.get('SUPABASE_CONNECTION_TIMEOUT', 30)\n"
            "            \n"
            "            # Initialize the pool\n"
            "            try:\n"
            "                initialize_supabase_pool(\n"
            "                    supabase_url=supabase_url,\n"
            "                    supabase_key=supabase_key,\n"
            "                    min_connections=min_connections,\n"
            "                    max_connections=max_connections,\n"
            "                    connection_timeout=connection_timeout\n"
            "                )\n"
            "                g.supabase_pool_initialized = True\n"
            "                logger.info(\"Supabase connection pool initialized\")\n"
            "            except Exception as e:\n"
            "                logger.error(f\"Error initializing Supabase connection pool: {str(e)}\")\n"
            "                # Fall back to direct connections\n"
            "                USE_CONNECTION_POOL = False\n"
            "        \n"
            "        # Get a client from the pool\n"
            "        if not hasattr(g, 'supabase'):\n"
            "            try:\n"
            "                g.supabase = get_pooled_supabase_client()\n"
            "            except ConnectionPoolError as e:\n"
            "                logger.error(f\"Error getting Supabase client from pool: {str(e)}\")\n"
            "                # Fall back to direct connection\n"
            "                supabase_url = current_app.config['SUPABASE_URL']\n"
            "                supabase_key = current_app.config['SUPABASE_KEY']\n"
            "                g.supabase = create_client(supabase_url, supabase_key)\n"
            "        \n"
            "        return g.supabase\n"
            "    \n"
            "    # Original implementation (direct connections)\n"
            "    if not hasattr(g, 'supabase'):\n"
            "        supabase_url = current_app.config['SUPABASE_URL']\n"
            "        supabase_key = current_app.config['SUPABASE_KEY']\n"
            "        \n"
            "        if not supabase_url or not supabase_key:\n"
            "            raise ValueError(\"SUPABASE_URL and SUPABASE_KEY must be set in configuration\")\n"
            "        \n"
            "        g.supabase = create_client(supabase_url, supabase_key)\n"
            "    \n"
            "    return g.supabase"
        )
        
        content = re.sub(client_pattern, client_replacement, content, flags=re.DOTALL)
        
        # Add teardown function to release connections back to the pool
        teardown_pattern = r"def get_user_id\(\):(.*?)return None"
        teardown_replacement = (
            "def get_user_id():\n"
            "    \"\"\"\n"
            "    Get the current user ID from the authenticated user.\n"
            "    \n"
            "    Returns:\n"
            "        str: User ID or None if not authenticated\n"
            "    \"\"\"\n"
            "    if hasattr(g, 'user') and g.user:\n"
            "        return g.user.id\n"
            "    \n"
            "    return None\n"
            "\n"
            "\n"
            "def release_supabase_connection():\n"
            "    \"\"\"\n"
            "    Release the Supabase client connection back to the pool.\n"
            "    This function should be called in the teardown_request handler.\n"
            "    \"\"\"\n"
            "    if not USE_CONNECTION_POOL:\n"
            "        return\n"
            "    \n"
            "    if hasattr(g, 'supabase'):\n"
            "        try:\n"
            "            release_supabase_client(g.supabase)\n"
            "            delattr(g, 'supabase')\n"
            "        except Exception as e:\n"
            "            logger.error(f\"Error releasing Supabase client: {str(e)}\")"
        )
        
        content = re.sub(teardown_pattern, teardown_replacement, content, flags=re.DOTALL)
        
        # Write the modified content back to the file
        if not dry_run:
            with open(utils_path, 'w') as f:
                f.write(content)
            logger.info(f"Successfully patched {utils_path}")
        else:
            logger.info(f"Dry run: would patch {utils_path}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error patching {utils_path}: {str(e)}")
        return False

def patch_app_file(app_path, dry_run=False):
    """
    Patch the app.py file to release connections in the teardown_request handler.
    
    Args:
        app_path: Path to the app.py file
        dry_run: If True, only show what would be changed without making changes
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Read the original file
        with open(app_path, 'r') as f:
            content = f.read()
        
        # Create a backup
        if not dry_run:
            backup_path = create_backup(app_path)
        
        # Add import for release_supabase_connection
        import_pattern = r"from api.utils import get_supabase_client, authenticate_user, token_required"
        import_replacement = (
            "from api.utils import get_supabase_client, authenticate_user, token_required, release_supabase_connection"
        )
        
        content = re.sub(import_pattern, import_replacement, content)
        
        # Modify the teardown_request handler
        teardown_pattern = r"@app\.teardown_request\s+def teardown_request\(exception=None\):(.*?)if hasattr\(g, 'supabase'\):\s+del g\.supabase"
        teardown_replacement = (
            "@app.teardown_request\n"
            "    def teardown_request(exception=None):\n"
            "        # Clean up resources\n"
            "        if hasattr(g, 'supabase'):\n"
            "            # Release connection back to the pool if using connection pooling\n"
            "            release_supabase_connection()\n"
            "            # For non-pooled connections, just delete the reference\n"
            "            if hasattr(g, 'supabase'):\n"
            "                del g.supabase"
        )
        
        content = re.sub(teardown_pattern, teardown_replacement, content, flags=re.DOTALL)
        
        # Add connection pool configuration to the app config
        config_pattern = r"app\.config\.update\(\{\s+'APISPEC_SPEC': APISpec\("
        config_replacement = (
            "# Configure connection pooling\n"
            "    app.config.update({\n"
            "        'SUPABASE_MIN_CONNECTIONS': 2,\n"
            "        'SUPABASE_MAX_CONNECTIONS': 10,\n"
            "        'SUPABASE_CONNECTION_TIMEOUT': 30\n"
            "    })\n"
            "    \n"
            "    app.config.update({\n"
            "        'APISPEC_SPEC': APISpec("
        )
        
        content = re.sub(config_pattern, config_replacement, content)
        
        # Write the modified content back to the file
        if not dry_run:
            with open(app_path, 'w') as f:
                f.write(content)
            logger.info(f"Successfully patched {app_path}")
        else:
            logger.info(f"Dry run: would patch {app_path}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error patching {app_path}: {str(e)}")
        return False

def update_config_file(config_path, dry_run=False):
    """
    Update the config.py file to add connection pool settings.
    
    Args:
        config_path: Path to the config.py file
        dry_run: If True, only show what would be changed without making changes
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Read the original file
        with open(config_path, 'r') as f:
            content = f.read()
        
        # Create a backup
        if not dry_run:
            backup_path = create_backup(config_path)
        
        # Add connection pool settings to the Config class
        config_pattern = r"class Config:\s+\"\"\"Base configuration.\"\"\"\s+(.*?)# API settings"
        config_replacement = (
            "class Config:\n"
            "    \"\"\"Base configuration.\"\"\"\n"
            "    SECRET_KEY = os.getenv('SECRET_KEY', 'dev-key-please-change-in-production')\n"
            "    DEBUG = False\n"
            "    TESTING = False\n"
            "    \n"
            "    # Supabase connection\n"
            "    SUPABASE_URL = os.getenv('SUPABASE_URL')\n"
            "    SUPABASE_KEY = os.getenv('SUPABASE_KEY')\n"
            "    \n"
            "    # Authentication\n"
            "    SUPABASE_USER = os.getenv('SUPABASE_USER')\n"
            "    SUPABASE_PASSWORD = os.getenv('SUPABASE_PASSWORD')\n"
            "    \n"
            "    # Connection pool settings\n"
            "    SUPABASE_MIN_CONNECTIONS = int(os.getenv('SUPABASE_MIN_CONNECTIONS', '2'))\n"
            "    SUPABASE_MAX_CONNECTIONS = int(os.getenv('SUPABASE_MAX_CONNECTIONS', '10'))\n"
            "    SUPABASE_CONNECTION_TIMEOUT = int(os.getenv('SUPABASE_CONNECTION_TIMEOUT', '30'))\n"
            "    SUPABASE_CONNECTION_LIFETIME = int(os.getenv('SUPABASE_CONNECTION_LIFETIME', '3600'))\n"
            "    SUPABASE_IDLE_TIMEOUT = int(os.getenv('SUPABASE_IDLE_TIMEOUT', '300'))\n"
            "    \n"
            "    # API settings"
        )
        
        content = re.sub(config_pattern, config_replacement, content, flags=re.DOTALL)
        
        # Write the modified content back to the file
        if not dry_run:
            with open(config_path, 'w') as f:
                f.write(content)
            logger.info(f"Successfully updated {config_path}")
        else:
            logger.info(f"Dry run: would update {config_path}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error updating {config_path}: {str(e)}")
        return False

def create_test_script(output_path, dry_run=False):
    """
    Create a test script to verify the connection pool works under load.
    
    Args:
        output_path: Path to save the test script
        dry_run: If True, only show what would be created without creating it
        
    Returns:
        True if successful, False otherwise
    """
    try:
        test_script_content = """#!/usr/bin/env python3
"""
        test_script_content += '''"""
CryoProtect v2 - Test Connection Pooling

This script tests the connection pool under load conditions.
It simulates concurrent requests to verify the connection pool works correctly.
"""

import os
import sys
import time
import json
import logging
import argparse
import threading
import random
import statistics
from datetime import datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("test_connection_pooling.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("test_connection_pooling")

# Import connection pooling
try:
    from implement_connection_pooling import (
        initialize_connection_pool, execute_query,
        get_connection_pool_stats, shutdown_connection_pool
    )
except ImportError:
    logger.error("Error: implement_connection_pooling.py not found. Make sure it's in the same directory.")
    sys.exit(1)

# Import Supabase MCP tools
try:
    from supabase_mcp_tools import execute_sql_on_supabase
except ImportError:
    logger.error("Error: supabase_mcp_tools.py not found. Make sure it's in the same directory.")
    sys.exit(1)

def run_query(query, project_id=None, use_pool=True):
    """
    Run a query using either the connection pool or direct execution.
    
    Args:
        query: The SQL query to execute
        project_id: The Supabase project ID (required if use_pool=False)
        use_pool: Whether to use the connection pool
        
    Returns:
        Tuple of (result, execution_time)
    """
    start_time = time.time()
    
    try:
        if use_pool:
            result = execute_query(query)
        else:
            result = execute_sql_on_supabase(project_id, query)
        
        execution_time = time.time() - start_time
        return result, execution_time
    
    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"Error executing query: {str(e)}")
        return None, execution_time

def worker(worker_id, num_queries, query, project_id, use_pool):
    """
    Worker function to execute queries.
    
    Args:
        worker_id: ID of the worker
        num_queries: Number of queries to execute
        query: The SQL query to execute
        project_id: The Supabase project ID
        use_pool: Whether to use the connection pool
        
    Returns:
        List of execution times
    """
    execution_times = []
    errors = 0
    
    for i in range(num_queries):
        try:
            _, execution_time = run_query(query, project_id, use_pool)
            execution_times.append(execution_time)
            logger.debug(f"Worker {worker_id}, Query {i}: {execution_time:.4f}s")
        except Exception as e:
            errors += 1
            logger.error(f"Worker {worker_id}, Query {i} failed: {str(e)}")
    
    return {
        "worker_id": worker_id,
        "execution_times": execution_times,
        "errors": errors,
        "avg_time": statistics.mean(execution_times) if execution_times else 0,
        "min_time": min(execution_times) if execution_times else 0,
        "max_time": max(execution_times) if execution_times else 0
    }

def run_load_test(project_id, num_workers, queries_per_worker, use_pool=True):
    """
    Run a load test with multiple workers.
    
    Args:
        project_id: The Supabase project ID
        num_workers: Number of concurrent workers
        queries_per_worker: Number of queries per worker
        use_pool: Whether to use the connection pool
        
    Returns:
        Dictionary with test results
    """
    # Test query
    query = "SELECT COUNT(*) FROM public.molecules;"
    
    # Initialize connection pool if using it
    if use_pool:
        initialize_connection_pool(
            project_id=project_id,
            min_connections=2,
            max_connections=10,
            connection_timeout=30
        )
    
    # Run the load test
    start_time = time.time()
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i in range(num_workers):
            future = executor.submit(
                worker, i, queries_per_worker, query, project_id, use_pool
            )
            futures.append(future)
        
        # Wait for all workers to complete
        worker_results = [future.result() for future in futures]
    
    end_time = time.time()
    
    # Calculate statistics
    all_times = []
    total_errors = 0
    
    for result in worker_results:
        all_times.extend(result["execution_times"])
        total_errors += result["errors"]
    
    # Get pool stats if using the pool
    pool_stats = get_connection_pool_stats() if use_pool else None
    
    # Shutdown the pool if using it
    if use_pool:
        shutdown_connection_pool()
    
    # Calculate overall statistics
    total_queries = num_workers * queries_per_worker
    successful_queries = total_queries - total_errors
    
    stats = {
        "timestamp": datetime.now().isoformat(),
        "use_pool": use_pool,
        "num_workers": num_workers,
        "queries_per_worker": queries_per_worker,
        "total_queries": total_queries,
        "successful_queries": successful_queries,
        "total_time": end_time - start_time,
        "queries_per_second": successful_queries / (end_time - start_time),
        "success_rate": successful_queries / total_queries if total_queries > 0 else 0,
        "avg_time": statistics.mean(all_times) if all_times else 0,
        "min_time": min(all_times) if all_times else 0,
        "max_time": max(all_times) if all_times else 0,
        "median_time": statistics.median(all_times) if all_times else 0,
        "stdev_time": statistics.stdev(all_times) if len(all_times) > 1 else 0,
        "total_errors": total_errors,
        "worker_results": worker_results,
        "pool_stats": pool_stats
    }
    
    return stats

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Test connection pooling for CryoProtect database")
    parser.add_argument("--project-id", default="tsdlmynydfuypiugmkev",
                        help="Supabase project ID (default: tsdlmynydfuypiugmkev)")
    parser.add_argument("--workers", type=int, default=5,
                        help="Number of concurrent workers (default: 5)")
    parser.add_argument("--queries", type=int, default=20,
                        help="Number of queries per worker (default: 20)")
    parser.add_argument("--compare", action="store_true",
                        help="Compare performance with and without connection pooling")
    parser.add_argument("--output", default="connection_pool_test_results.json",
                        help="Output file for test results (default: connection_pool_test_results.json)")
    return parser.parse_args()

def main():
    """Main function to test connection pooling."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        logger.info(f"Testing connection pooling for project {args.project_id}")
        logger.info(f"Workers: {args.workers}, Queries per worker: {args.queries}")
        
        results = {}
        
        # Run test with connection pooling
        logger.info("Running test with connection pooling...")
        results["with_pool"] = run_load_test(
            project_id=args.project_id,
            num_workers=args.workers,
            queries_per_worker=args.queries,
            use_pool=True
        )
        
        logger.info(f"Test with pool: {results['with_pool']['queries_per_second']:.2f} queries/sec, " +
                   f"{results['with_pool']['success_rate'] * 100:.2f}% success rate")
        
        # Run test without connection pooling if requested
        if args.compare:
            logger.info("Running test without connection pooling...")
            results["without_pool"] = run_load_test(
                project_id=args.project_id,
                num_workers=args.workers,
                queries_per_worker=args.queries,
                use_pool=False
            )
            
            logger.info(f"Test without pool: {results['without_pool']['queries_per_second']:.2f} queries/sec, " +
                       f"{results['without_pool']['success_rate'] * 100:.2f}% success rate")
            
            # Calculate improvement
            with_pool_qps = results["with_pool"]["queries_per_second"]
            without_pool_qps = results["without_pool"]["queries_per_second"]
            
            if without_pool_qps > 0:
                improvement = (with_pool_qps - without_pool_qps) / without_pool_qps * 100
                logger.info(f"Performance improvement: {improvement:.2f}%")
                
                results["comparison"] = {
                    "improvement_percent": improvement,
                    "with_pool_qps": with_pool_qps,
                    "without_pool_qps": without_pool_qps
                }
        
        # Save results to file
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"Test results saved to {args.output}")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error testing connection pooling: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
'''
        
        if not dry_run:
            with open(output_path, 'w') as f:
                f.write(test_script_content)
            # Make the script executable
            os.chmod(output_path, 0o755)
            logger.info(f"Created test script at {output_path}")
        else:
            logger.info(f"Dry run: would create test script at {output_path}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error creating test script: {str(e)}")
        return False

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Patch connection pooling for CryoProtect")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be changed without making changes")
    parser.add_argument("--utils-path", default="api/utils.py",
                        help="Path to the api/utils.py file (default: api/utils.py)")
    parser.add_argument("--app-path", default="app.py",
                        help="Path to the app.py file (default: app.py)")
    parser.add_argument("--config-path", default="config.py",
                        help="Path to the config.py file (default: config.py)")
    parser.add_argument("--test-script", default="test_connection_pooling.py",
                        help="Path to save the test script (default: test_connection_pooling.py)")
    return parser.parse_args()

def main():
    """Main function to patch connection pooling."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        logger.info("Patching connection pooling for CryoProtect")
        logger.info(f"Dry run: {args.dry_run}")
        
        # Patch the api/utils.py file
        utils_success = patch_utils_file(args.utils_path, args.dry_run)
        
        # Patch the app.py file
        app_success = patch_app_file(args.app_path, args.dry_run)
        
        # Update the config.py file
        config_success = update_config_file(args.config_path, args.dry_run)
        
        # Create the test script
        test_script_success = create_test_script(args.test_script, args.dry_run)
        
        # Check if all operations were successful
        if utils_success and app_success and config_success and test_script_success:
            logger.info("Successfully patched connection pooling")
            return 0
        else:
            logger.error("Some operations failed")
            return 1
    
    except Exception as e:
        logger.error(f"Error patching connection pooling: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())