#!/usr/bin/env python3
"""
CryoProtect v2 - Comprehensive Integration Test Script

This script verifies that all implemented fixes work together correctly:
1. Performance improvements (indexes are applied correctly)
2. Connection pooling (connections are pooled and reused)
3. API integration (endpoints work with new database structure)

It generates a detailed report of test results and uses memory bank for caching.
"""

import os
import sys
import json
import time
import logging
import argparse
import requests
import functools
import threading
import statistics
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Union
from concurrent.futures import ThreadPoolExecutor

# Import logging configuration
try:
    from logging_config import setup_logging
except ImportError:
    # Fallback logging setup if the import fails
    def setup_logging():
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers=[
                logging.FileHandler("verify_all_fixes.log"),
                logging.StreamHandler()
            ]
        )

# Try to import Supabase client
try:
    from supabase import create_client, Client
except ImportError:
    print("Error: Supabase Python client not found. Please install it with 'pip install supabase'.")
    sys.exit(1)

# Try to import connection pool wrapper
try:
    from connection_pool_wrapper import (
        get_supabase_connection, initialize_supabase_pool,
        get_supabase_pool_stats, shutdown_supabase_pool
    )
except ImportError:
    print("Error: connection_pool_wrapper.py not found. Make sure it's in the same directory.")
    sys.exit(1)

# Memory bank for caching test results
MEMORY_BANK_DIR = Path("memory-bank")
MEMORY_CACHE_FILE = MEMORY_BANK_DIR / "integration_test_cache.json"
API_ENDPOINTS_FILE = MEMORY_BANK_DIR / "api_endpoints.json"

# Set up logging
setup_logging()
logger = logging.getLogger("verify_all_fixes")

# Constants
DEFAULT_API_URL = "http://localhost:5000"
DEFAULT_SUPABASE_URL = os.getenv("SUPABASE_URL", "")
DEFAULT_SUPABASE_KEY = os.getenv("SUPABASE_KEY", "")
DEFAULT_CONCURRENT_USERS = 5
DEFAULT_TEST_DURATION = 30  # seconds
DEFAULT_REQUEST_TIMEOUT = 10  # seconds

# Status levels
class Status:
    SUCCESS = "SUCCESS"
    COMPLETED_WITH_WARNINGS = "COMPLETED_WITH_WARNINGS"
    ERROR = "ERROR"

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
                "test_results": {},
                "performance_metrics": {},
                "connection_pool_metrics": {},
                "api_integration_metrics": {}
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
        
        # Check if result is in cache and not older than 1 hour
        if "test_results" in cache and key in cache["test_results"]:
            cached_time = datetime.fromisoformat(cache["test_results"][key]["timestamp"])
            current_time = datetime.now()
            time_diff = (current_time - cached_time).total_seconds()
            
            # Use cached result if it's less than 1 hour old
            if time_diff < 3600:  # 1 hour in seconds
                logger.info(f"Using cached result for {func.__name__} (cached {time_diff:.1f} seconds ago)")
                return cache["test_results"][key]["result"]
        
        # Call the function
        result = func(*args, **kwargs)
        
        # Update the cache
        cache.setdefault("test_results", {})
        cache["test_results"][key] = {
            "result": result,
            "timestamp": datetime.now().isoformat()
        }
        
        # Save the cache
        with open(MEMORY_CACHE_FILE, 'w') as f:
            json.dump(cache, f, indent=2)
        
        return result
    
    return wrapper

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Verify all CryoProtect fixes work together correctly")
    parser.add_argument("--supabase-url", default=DEFAULT_SUPABASE_URL,
                        help=f"Supabase URL (default: {DEFAULT_SUPABASE_URL or 'from env'})")
    parser.add_argument("--supabase-key", default=DEFAULT_SUPABASE_KEY,
                        help=f"Supabase API key (default: {DEFAULT_SUPABASE_KEY or 'from env'})")
    parser.add_argument("--api-url", default=DEFAULT_API_URL,
                        help=f"API URL (default: {DEFAULT_API_URL})")
    parser.add_argument("--concurrent-users", type=int, default=DEFAULT_CONCURRENT_USERS,
                        help=f"Number of concurrent users for load testing (default: {DEFAULT_CONCURRENT_USERS})")
    parser.add_argument("--test-duration", type=int, default=DEFAULT_TEST_DURATION,
                        help=f"Duration of load tests in seconds (default: {DEFAULT_TEST_DURATION})")
    parser.add_argument("--skip-performance", action="store_true",
                        help="Skip performance tests")
    parser.add_argument("--skip-connection-pooling", action="store_true",
                        help="Skip connection pooling tests")
    parser.add_argument("--skip-api", action="store_true",
                        help="Skip API integration tests")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose output")
    return parser.parse_args()

def load_api_endpoints():
    """Load API endpoints from the memory bank."""
    try:
        with open(API_ENDPOINTS_FILE, 'r') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        logger.warning(f"API endpoints file not found or invalid: {API_ENDPOINTS_FILE}")
        return {
            "base_url": DEFAULT_API_URL,
            "endpoints": {
                "health": {"path": "/health", "method": "GET"},
                "molecules": {"path": "/api/molecules", "method": "GET"},
                "mixtures": {"path": "/api/mixtures", "method": "GET"}
            }
        }

#
# Performance Tests
#

@memory_cache
def verify_indexes_exist(supabase_url, supabase_key):
    """
    Verify that all required indexes exist in the database.
    
    Args:
        supabase_url: Supabase URL
        supabase_key: Supabase API key
        
    Returns:
        dict: Results of the verification
    """
    logger.info("Verifying indexes exist...")
    
    # List of indexes that should exist
    required_indexes = [
        # RLS performance indexes
        "idx_mixtures_created_by",
        "idx_experiments_created_by",
        "idx_predictions_created_by",
        
        # Mixture and component indexes
        "idx_mixture_component_mixture_id",
        "idx_experiments_mixture_id",
        
        # Prediction and property indexes
        "idx_predictions_mixture_property",
        "idx_molecular_property_property_type",
        
        # Text search indexes
        "idx_molecule_name_trgm"
    ]
    
    # Initialize Supabase client
    supabase = create_client(supabase_url, supabase_key)
    
    # Query to get all indexes
    query = """
    SELECT indexname, tablename
    FROM pg_indexes
    WHERE schemaname = 'public'
    ORDER BY tablename, indexname;
    """
    
    try:
        # Execute the query
        response = supabase.rpc('exec_sql', {'sql': query}).execute()
        
        if hasattr(response, 'data'):
            indexes = response.data
        else:
            # Handle different response format
            indexes = response[1] if isinstance(response, list) and len(response) > 1 else []
        
        # Check which required indexes exist
        existing_indexes = [idx['indexname'] for idx in indexes]
        missing_indexes = [idx for idx in required_indexes if idx not in existing_indexes]
        
        # Calculate percentage of indexes that exist
        total_required = len(required_indexes)
        total_existing = total_required - len(missing_indexes)
        percentage = (total_existing / total_required) * 100 if total_required > 0 else 0
        
        # Determine status
        status = Status.SUCCESS if percentage == 100 else (
            Status.COMPLETED_WITH_WARNINGS if percentage >= 75 else Status.ERROR
        )
        
        return {
            "status": status,
            "total_required": total_required,
            "total_existing": total_existing,
            "percentage": percentage,
            "existing_indexes": existing_indexes,
            "missing_indexes": missing_indexes,
            "all_indexes": indexes
        }
    except Exception as e:
        logger.error(f"Error verifying indexes: {str(e)}")
        return {
            "status": Status.ERROR,
            "error": str(e),
            "total_required": len(required_indexes),
            "total_existing": 0,
            "percentage": 0,
            "existing_indexes": [],
            "missing_indexes": required_indexes
        }

@memory_cache
def test_query_performance(supabase_url, supabase_key):
    """
    Test the performance of key database queries.
    
    Args:
        supabase_url: Supabase URL
        supabase_key: Supabase API key
        
    Returns:
        dict: Results of the performance tests
    """
    logger.info("Testing query performance...")
    
    # Initialize Supabase client
    supabase = create_client(supabase_url, supabase_key)
    
    # Test queries
    test_queries = {
        "rls_filter": """
        EXPLAIN ANALYZE
        SELECT * FROM public.mixtures
        WHERE created_by = '00000000-0000-0000-0000-000000000000'
        LIMIT 10;
        """,
        
        "mixture_components_join": """
        EXPLAIN ANALYZE
        SELECT m.name, mc.amount, mc.amount_unit
        FROM public.mixtures m
        JOIN public.mixture_components mc ON m.id = mc.mixture_id
        LIMIT 10;
        """,
        
        "predictions_filter": """
        EXPLAIN ANALYZE
        SELECT p.value, p.unit
        FROM public.predictions p
        WHERE p.mixture_id = (SELECT id FROM public.mixtures LIMIT 1)
        AND p.property_type_id = (SELECT id FROM public.property_types LIMIT 1);
        """,
        
        "experiments_join": """
        EXPLAIN ANALYZE
        SELECT e.*, m.name
        FROM public.experiments e
        JOIN public.mixtures m ON e.mixture_id = m.id
        LIMIT 10;
        """,
        
        "molecule_text_search": """
        EXPLAIN ANALYZE
        SELECT * FROM public.molecules
        WHERE name ILIKE '%glycerol%'
        LIMIT 10;
        """
    }
    
    results = {}
    
    for name, query in test_queries.items():
        try:
            logger.info(f"Running test query: {name}")
            
            # Measure execution time
            start_time = time.time()
            response = supabase.rpc('exec_sql', {'sql': query}).execute()
            execution_time = time.time() - start_time
            
            # Extract query plan
            if hasattr(response, 'data'):
                query_plan = response.data
            else:
                # Handle different response format
                query_plan = response[1] if isinstance(response, list) and len(response) > 1 else []
            
            # Extract execution time from query plan
            plan_execution_time = None
            for line in query_plan:
                if isinstance(line, dict) and 'QUERY PLAN' in line:
                    plan_text = line['QUERY PLAN']
                    if 'Execution Time:' in plan_text:
                        time_str = plan_text.split('Execution Time:')[1].strip()
                        if 'ms' in time_str:
                            plan_execution_time = float(time_str.split('ms')[0].strip())
            
            # Determine status based on execution time
            status = Status.SUCCESS
            if plan_execution_time is not None:
                if plan_execution_time > 500:
                    status = Status.ERROR
                elif plan_execution_time > 200:
                    status = Status.COMPLETED_WITH_WARNINGS
            
            results[name] = {
                "status": status,
                "client_execution_time_ms": execution_time * 1000,
                "plan_execution_time_ms": plan_execution_time,
                "query_plan": query_plan
            }
        except Exception as e:
            logger.error(f"Error running test query {name}: {str(e)}")
            results[name] = {
                "status": Status.ERROR,
                "error": str(e)
            }
    
    # Calculate overall status
    statuses = [result["status"] for result in results.values()]
    overall_status = Status.SUCCESS
    if Status.ERROR in statuses:
        overall_status = Status.ERROR
    elif Status.COMPLETED_WITH_WARNINGS in statuses:
        overall_status = Status.COMPLETED_WITH_WARNINGS
    
    return {
        "status": overall_status,
        "query_results": results
    }

def run_performance_tests(supabase_url, supabase_key):
    """
    Run all performance tests.
    
    Args:
        supabase_url: Supabase URL
        supabase_key: Supabase API key
        
    Returns:
        dict: Results of all performance tests
    """
    logger.info("Running performance tests...")
    
    # Verify indexes exist
    index_results = verify_indexes_exist(supabase_url, supabase_key)
    
    # Only test query performance if indexes exist
    query_results = {"status": Status.ERROR, "error": "Skipped due to missing indexes"}
    if index_results["status"] != Status.ERROR:
        query_results = test_query_performance(supabase_url, supabase_key)
    
    # Determine overall status
    overall_status = Status.SUCCESS
    if index_results["status"] == Status.ERROR or query_results["status"] == Status.ERROR:
        overall_status = Status.ERROR
    elif index_results["status"] == Status.COMPLETED_WITH_WARNINGS or query_results["status"] == Status.COMPLETED_WITH_WARNINGS:
        overall_status = Status.COMPLETED_WITH_WARNINGS
    
    return {
        "status": overall_status,
        "index_verification": index_results,
        "query_performance": query_results,
        "timestamp": datetime.now().isoformat()
    }

#
# Connection Pooling Tests
#

def test_connection_pooling(supabase_url, supabase_key, concurrent_users=5, test_duration=30):
    """
    Test connection pooling by simulating concurrent users.
    
    Args:
        supabase_url: Supabase URL
        supabase_key: Supabase API key
        concurrent_users: Number of concurrent users to simulate
        test_duration: Duration of the test in seconds
        
    Returns:
        dict: Results of the connection pooling tests
    """
    logger.info(f"Testing connection pooling with {concurrent_users} concurrent users for {test_duration} seconds...")
    
    # Initialize connection pool
    try:
        initialize_supabase_pool(
            supabase_url=supabase_url,
            supabase_key=supabase_key,
            min_connections=2,
            max_connections=concurrent_users + 2  # Add buffer
        )
    except Exception as e:
        logger.error(f"Error initializing connection pool: {str(e)}")
        return {
            "status": Status.ERROR,
            "error": str(e)
        }
    
    # Get initial pool stats
    initial_stats = get_supabase_pool_stats()
    
    # Define a worker function for each simulated user
    def worker(user_id):
        results = []
        start_time = time.time()
        query_count = 0
        
        # Run queries until the test duration is reached
        while time.time() - start_time < test_duration:
            try:
                # Use the connection pool to execute a simple query
                with get_supabase_connection() as supabase:
                    query_start = time.time()
                    response = supabase.from_('molecules').select('*').limit(5).execute()
                    query_end = time.time()
                    
                    # Record the result
                    results.append({
                        "user_id": user_id,
                        "query_time_ms": (query_end - query_start) * 1000,
                        "success": True,
                        "timestamp": datetime.now().isoformat()
                    })
                    
                    query_count += 1
                    
                    # Small delay to simulate user think time
                    time.sleep(0.1)
            except Exception as e:
                # Record the error
                results.append({
                    "user_id": user_id,
                    "error": str(e),
                    "success": False,
                    "timestamp": datetime.now().isoformat()
                })
                
                # Longer delay after an error
                time.sleep(0.5)
        
        return {
            "user_id": user_id,
            "query_count": query_count,
            "results": results
        }
    
    # Start worker threads
    with ThreadPoolExecutor(max_workers=concurrent_users) as executor:
        futures = [executor.submit(worker, i) for i in range(concurrent_users)]
        user_results = [future.result() for future in futures]
    
    # Get final pool stats
    final_stats = get_supabase_pool_stats()
    
    # Shutdown the connection pool
    shutdown_supabase_pool()
    
    # Analyze results
    all_query_times = []
    total_queries = 0
    total_errors = 0
    
    for user_result in user_results:
        total_queries += user_result["query_count"]
        for result in user_result["results"]:
            if result["success"] and "query_time_ms" in result:
                all_query_times.append(result["query_time_ms"])
            else:
                total_errors += 1
    
    # Calculate statistics
    avg_query_time = statistics.mean(all_query_times) if all_query_times else 0
    p95_query_time = statistics.quantiles(all_query_times, n=20)[19] if len(all_query_times) >= 20 else (max(all_query_times) if all_query_times else 0)
    error_rate = (total_errors / total_queries) * 100 if total_queries > 0 else 0
    
    # Determine status
    status = Status.SUCCESS
    if error_rate > 5 or avg_query_time > 500:
        status = Status.ERROR
    elif error_rate > 1 or avg_query_time > 200:
        status = Status.COMPLETED_WITH_WARNINGS
    
    return {
        "status": status,
        "concurrent_users": concurrent_users,
        "test_duration_seconds": test_duration,
        "total_queries": total_queries,
        "total_errors": total_errors,
        "error_rate_percent": error_rate,
        "average_query_time_ms": avg_query_time,
        "p95_query_time_ms": p95_query_time,
        "initial_pool_stats": initial_stats,
        "final_pool_stats": final_stats,
        "user_results": user_results,
        "timestamp": datetime.now().isoformat()
    }

#
# API Integration Tests
#

def test_api_endpoints(api_url, timeout=DEFAULT_REQUEST_TIMEOUT):
    """
    Test API endpoints to verify they work with the new database structure.
    
    Args:
        api_url: Base URL of the API
        timeout: Request timeout in seconds
        
    Returns:
        dict: Results of the API endpoint tests
    """
    logger.info(f"Testing API endpoints at {api_url}...")
    
    # Load API endpoints from memory bank
    api_config = load_api_endpoints()
    base_url = api_config.get("base_url", api_url)
    endpoints = api_config.get("endpoints", {})
    
    # If no endpoints are defined, use default endpoints
    if not endpoints:
        endpoints = {
            "health": {"path": "/health", "method": "GET"},
            "molecules": {"path": "/api/molecules", "method": "GET"},
            "mixtures": {"path": "/api/mixtures", "method": "GET"}
        }
    
    results = {}
    
    # Test each endpoint
    for name, endpoint in endpoints.items():
        path = endpoint.get("path", "")
        method = endpoint.get("method", "GET")
        auth_required = endpoint.get("authentication_required", False)
        
        # Skip endpoints that require authentication for now
        if auth_required:
            results[name] = {
                "status": "SKIPPED",
                "reason": "Authentication required"
            }
            continue
        
        # Construct the full URL
        url = f"{api_url}{path}"
        
        try:
            logger.info(f"Testing endpoint: {method} {url}")
            
            # Make the request
            start_time = time.time()
            if method == "GET":
                response = requests.get(url, timeout=timeout)
            elif method == "POST":
                response = requests.post(url, json={}, timeout=timeout)
            else:
                results[name] = {
                    "status": "SKIPPED",
                    "reason": f"Unsupported method: {method}"
                }
                continue
            
            execution_time = time.time() - start_time
            
            # Check the response
            if response.status_code < 400:
                status = Status.SUCCESS
            elif response.status_code < 500:
                status = Status.COMPLETED_WITH_WARNINGS
            else:
                status = Status.ERROR
            
            # Try to parse JSON response
            try:
                response_data = response.json()
            except:
                response_data = {"text": response.text[:500] + ("..." if len(response.text) > 500 else "")}
            
            results[name] = {
                "status": status,
                "url": url,
                "method": method,
                "response_code": response.status_code,
                "execution_time_ms": execution_time * 1000,
                "response_data": response_data
            }
        except Exception as e:
            logger.error(f"Error testing endpoint {name}: {str(e)}")
            results[name] = {
                "status": Status.ERROR,
                "url": url,
                "method": method,
                "error": str(e)
            }
    
    # Calculate overall status
    statuses = [result["status"] for result in results.values() if result["status"] != "SKIPPED"]
    overall_status = Status.SUCCESS
    if Status.ERROR in statuses:
        overall_status = Status.ERROR
    elif Status.COMPLETED_WITH_WARNINGS in statuses:
        overall_status = Status.COMPLETED_WITH_WARNINGS
    
    return {
        "status": overall_status,
        "endpoint_results": results,
        "timestamp": datetime.now().isoformat()
    }

def test_api_with_different_roles(api_url, timeout=DEFAULT_REQUEST_TIMEOUT):
    """
    Test API endpoints with different user roles to verify RLS policies.
    
    Args:
        api_url: Base URL of the API
        timeout: Request timeout in seconds
        
    Returns:
        dict: Results of the API role tests
    """
    logger.info(f"Testing API with different user roles...")
    
    # Define test users with different roles
    test_users = {
        "anonymous": {
            "description": "Anonymous user (no authentication)",
            "token": None
        },
        "regular_user": {
            "description": "Regular authenticated user",
            "token": "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJyb2xlIjoidXNlciIsInVzZXJfaWQiOiIxMjM0NTY3OC05MDEyLTM0NTYtNzg5MC0xMjM0NTY3ODkwMTIifQ.token"
        },
        "admin_user": {
            "description": "Admin user",
            "token": "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJyb2xlIjoiYWRtaW4iLCJ1c2VyX2lkIjoiOTg3NjU0MzItMTIzNC01Njc4LTkwMTItMzQ1Njc4OTAxMjM0In0.token"
        }
    }
    
    # Define test endpoints that should have different behavior based on role
    test_endpoints = [
        {
            "name": "list_mixtures",
            "path": "/api/mixtures",
            "method": "GET",
            "expected_results": {
                "anonymous": {"status_code": 401},  # Unauthorized
                "regular_user": {"status_code": 200},  # OK, but should only see own mixtures
                "admin_user": {"status_code": 200}  # OK, should see all mixtures
            }
        },
        {
            "name": "list_experiments",
            "path": "/api/experiments",
            "method": "GET",
            "expected_results": {
                "anonymous": {"status_code": 401},  # Unauthorized
                "regular_user": {"status_code": 200},  # OK, but should only see own experiments
                "admin_user": {"status_code": 200}  # OK, should see all experiments
            }
        }
    ]
    
    results = {}
    
    # Test each endpoint with each user role
    for endpoint in test_endpoints:
        endpoint_results = {}
        
        for role, user in test_users.items():
            try:
                logger.info(f"Testing {endpoint['name']} as {role}...")
                
                # Prepare headers
                headers = {}
                if user["token"]:
                    headers["Authorization"] = f"Bearer {user['token']}"
                
                # Make the request
                url = f"{api_url}{endpoint['path']}"
                start_time = time.time()
                
                if endpoint["method"] == "GET":
                    response = requests.get(url, headers=headers, timeout=timeout)
                elif endpoint["method"] == "POST":
                    response = requests.post(url, headers=headers, json={}, timeout=timeout)
                else:
                    endpoint_results[role] = {
                        "status": "SKIPPED",
                        "reason": f"Unsupported method: {endpoint['method']}"
                    }
                    continue
                
                execution_time = time.time() - start_time
                
                # Check if the response matches expectations
                expected_status = endpoint["expected_results"].get(role, {}).get("status_code")
                if expected_status and response.status_code == expected_status:
                    status = Status.SUCCESS
                else:
                    status = Status.ERROR
                
                # Try to parse JSON response
                try:
                    response_data = response.json()
                except:
                    response_data = {"text": response.text[:500] + ("..." if len(response.text) > 500 else "")}
                
                endpoint_results[role] = {
                    "status": status,
                    "url": url,
                    "method": endpoint["method"],
                    "response_code": response.status_code,
                    "expected_code": expected_status,
                    "execution_time_ms": execution_time * 1000,
                    "response_data": response_data
                }
            except Exception as e:
                logger.error(f"Error testing {endpoint['name']} as {role}: {str(e)}")
                endpoint_results[role] = {
                    "status": Status.ERROR,
                    "url": url,
                    "method": endpoint["method"],
                    "error": str(e)
                }
        
        # Calculate overall status for this endpoint
        statuses = [result["status"] for result in endpoint_results.values() if result["status"] != "SKIPPED"]
        overall_status = Status.SUCCESS
        if Status.ERROR in statuses:
            overall_status = Status.ERROR
        elif Status.COMPLETED_WITH_WARNINGS in statuses:
            overall_status = Status.COMPLETED_WITH_WARNINGS
        
        results[endpoint["name"]] = {
            "status": overall_status,
            "role_results": endpoint_results
        }
    
    # Calculate overall status
    statuses = [result["status"] for result in results.values()]
    overall_status = Status.SUCCESS
    if Status.ERROR in statuses:
        overall_status = Status.ERROR
    elif Status.COMPLETED_WITH_WARNINGS in statuses:
        overall_status = Status.COMPLETED_WITH_WARNINGS
    
    return {
        "status": overall_status,
        "endpoint_results": results,
        "timestamp": datetime.now().isoformat()
    }

def run_api_tests(api_url):
    """
    Run all API integration tests.
    
    Args:
        api_url: Base URL of the API
        
    Returns:
        dict: Results of all API tests
    """
    logger.info("Running API integration tests...")
    
    # Test API endpoints
    endpoint_results = test_api_endpoints(api_url)
    
    # Test API with different roles
    role_results = test_api_with_different_roles(api_url)
    
    # Determine overall status
    overall_status = Status.SUCCESS
    if endpoint_results["status"] == Status.ERROR or role_results["status"] == Status.ERROR:
        overall_status = Status.ERROR
    elif endpoint_results["status"] == Status.COMPLETED_WITH_WARNINGS or role_results["status"] == Status.COMPLETED_WITH_WARNINGS:
        overall_status = Status.COMPLETED_WITH_WARNINGS
    
    return {
        "status": overall_status,
        "endpoint_tests": endpoint_results,
        "role_tests": role_results,
        "timestamp": datetime.now().isoformat()
    }

#
# Report Generation
#

def generate_report(results):
    """
    Generate a detailed report of test results.
    
    Args:
        results: Dictionary containing all test results
        
    Returns:
        str: Formatted report text
    """
    report = []
    report.append("=" * 80)
    report.append("CryoProtect v2 - Comprehensive Integration Test Report")
    report.append("=" * 80)
    report.append(f"Generated: {datetime.now().isoformat()}")
    report.append(f"Overall Status: {results['status']}")
    report.append("")
    
    # Summary
    report.append("-" * 80)
    report.append("Summary")
    report.append("-" * 80)
    
    # Performance tests summary
    if "performance_tests" in results:
        perf = results["performance_tests"]
        report.append(f"Performance Tests: {perf['status']}")
        
        if "index_verification" in perf:
            idx = perf["index_verification"]
            report.append(f"  - Indexes: {idx['status']} ({idx.get('percentage', 0):.1f}% complete)")
            if "missing_indexes" in idx and idx["missing_indexes"]:
                report.append(f"    - Missing indexes: {', '.join(idx['missing_indexes'])}")
        
        if "query_performance" in perf:
            qp = perf["query_performance"]
            report.append(f"  - Query Performance: {qp['status']}")
            
            # Add details for slow queries
            if "query_results" in qp:
                slow_queries = []
                for name, result in qp["query_results"].items():
                    if result["status"] != Status.SUCCESS:
                        exec_time = result.get("plan_execution_time_ms", "N/A")
                        slow_queries.append(f"{name} ({exec_time} ms)")
                
                if slow_queries:
                    report.append(f"    - Slow queries: {', '.join(slow_queries)}")
    
    # Connection pooling tests summary
    if "connection_pooling_tests" in results:
        cp = results["connection_pooling_tests"]
        report.append(f"Connection Pooling Tests: {cp['status']}")
        report.append(f"  - Concurrent Users: {cp.get('concurrent_users', 'N/A')}")
        report.append(f"  - Total Queries: {cp.get('total_queries', 'N/A')}")
        report.append(f"  - Error Rate: {cp.get('error_rate_percent', 'N/A'):.2f}%")
        report.append(f"  - Avg Query Time: {cp.get('average_query_time_ms', 'N/A'):.2f} ms")
        report.append(f"  - P95 Query Time: {cp.get('p95_query_time_ms', 'N/A'):.2f} ms")
    
    # API tests summary
    if "api_tests" in results:
        api = results["api_tests"]
        report.append(f"API Integration Tests: {api['status']}")
        
        if "endpoint_tests" in api:
            et = api["endpoint_tests"]
            total = len(et.get("endpoint_results", {}))
            success = sum(1 for r in et.get("endpoint_results", {}).values() 
                         if r.get("status") == Status.SUCCESS)
            report.append(f"  - Endpoints: {success}/{total} successful")
        
        if "role_tests" in api:
            rt = api["role_tests"]
            total = len(rt.get("endpoint_results", {}))
            success = sum(1 for r in rt.get("endpoint_results", {}).values() 
                         if r.get("status") == Status.SUCCESS)
            report.append(f"  - Role-based Access: {success}/{total} successful")
    
    # Detailed results
    report.append("")
    report.append("-" * 80)
    report.append("Detailed Results")
    report.append("-" * 80)
    
    # Performance tests details
    if "performance_tests" in results:
        perf = results["performance_tests"]
        report.append("Performance Tests:")
        
        if "index_verification" in perf:
            idx = perf["index_verification"]
            report.append(f"  Index Verification: {idx['status']}")
            report.append(f"  - Required Indexes: {idx.get('total_required', 0)}")
            report.append(f"  - Existing Indexes: {idx.get('total_existing', 0)}")
            report.append(f"  - Percentage Complete: {idx.get('percentage', 0):.1f}%")
            
            if "missing_indexes" in idx and idx["missing_indexes"]:
                report.append("  - Missing Indexes:")
                for missing in idx["missing_indexes"]:
                    report.append(f"    - {missing}")
        
        if "query_performance" in perf:
            qp = perf["query_performance"]
            report.append(f"  Query Performance: {qp['status']}")
            
            if "query_results" in qp:
                report.append("  - Query Results:")
                for name, result in qp["query_results"].items():
                    status = result["status"]
                    exec_time = result.get("plan_execution_time_ms", "N/A")
                    report.append(f"    - {name}: {status} ({exec_time} ms)")
    
    # Connection pooling tests details
    if "connection_pooling_tests" in results:
        cp = results["connection_pooling_tests"]
        report.append("Connection Pooling Tests:")
        report.append(f"  Status: {cp['status']}")
        report.append(f"  - Concurrent Users: {cp.get('concurrent_users', 'N/A')}")
        report.append(f"  - Test Duration: {cp.get('test_duration_seconds', 'N/A')} seconds")
        report.append(f"  - Total Queries: {cp.get('total_queries', 'N/A')}")
        report.append(f"  - Total Errors: {cp.get('total_errors', 'N/A')}")
        report.append(f"  - Error Rate: {cp.get('error_rate_percent', 'N/A'):.2f}%")
        report.append(f"  - Average Query Time: {cp.get('average_query_time_ms', 'N/A'):.2f} ms")
        report.append(f"  - P95 Query Time: {cp.get('p95_query_time_ms', 'N/A'):.2f} ms")
        
        if "initial_pool_stats" in cp and "final_pool_stats" in cp:
            report.append("  - Pool Statistics:")
            report.append("    Initial:")
            for k, v in cp["initial_pool_stats"].items():
                report.append(f"      {k}: {v}")
            report.append("    Final:")
            for k, v in cp["final_pool_stats"].items():
                report.append(f"      {k}: {v}")
    
    # API tests details
    if "api_tests" in results:
        api = results["api_tests"]
        report.append("API Integration Tests:")
        report.append(f"  Status: {api['status']}")
        
        if "endpoint_tests" in api:
            et = api["endpoint_tests"]
            report.append("  Endpoint Tests:")
            report.append(f"    Status: {et['status']}")
            
            if "endpoint_results" in et:
                report.append("    Results:")
                for name, result in et["endpoint_results"].items():
                    status = result["status"]
                    code = result.get("response_code", "N/A")
                    time_ms = result.get("execution_time_ms", "N/A")
                    if isinstance(time_ms, (int, float)):
                        time_ms = f"{time_ms:.2f} ms"
                    report.append(f"      {name}: {status} (Code: {code}, Time: {time_ms})")
        
        if "role_tests" in api:
            rt = api["role_tests"]
            report.append("  Role-based Access Tests:")
            report.append(f"    Status: {rt['status']}")
            
            if "endpoint_results" in rt:
                report.append("    Results:")
                for name, result in rt["endpoint_results"].items():
                    status = result["status"]
                    report.append(f"      {name}: {status}")
                    
                    if "role_results" in result:
                        for role, role_result in result["role_results"].items():
                            r_status = role_result["status"]
                            expected = role_result.get("expected_code", "N/A")
                            actual = role_result.get("response_code", "N/A")
                            report.append(f"        {role}: {r_status} (Expected: {expected}, Actual: {actual})")
    
    # Recommendations
    report.append("")
    report.append("-" * 80)
    report.append("Recommendations")
    report.append("-" * 80)
    
    recommendations = []
    
    # Performance recommendations
    if "performance_tests" in results:
        perf = results["performance_tests"]
        
        if "index_verification" in perf:
            idx = perf["index_verification"]
            if "missing_indexes" in idx and idx["missing_indexes"]:
                recommendations.append("Apply missing indexes to improve query performance:")
                for missing in idx["missing_indexes"]:
                    recommendations.append(f"  - Create index {missing}")
        
        if "query_performance" in perf:
            qp = perf["query_performance"]
            if "query_results" in qp:
                slow_queries = []
                for name, result in qp["query_results"].items():
                    if result["status"] != Status.SUCCESS:
                        slow_queries.append(name)
                
                if slow_queries:
                    recommendations.append("Optimize slow queries:")
                    for query in slow_queries:
                        recommendations.append(f"  - {query}")
    
    # Connection pooling recommendations
    if "connection_pooling_tests" in results:
        cp = results["connection_pooling_tests"]
        
        if cp.get("error_rate_percent", 0) > 1:
            recommendations.append("Improve connection pooling reliability:")
            recommendations.append("  - Increase max_connections parameter")
            recommendations.append("  - Implement better error handling and retry logic")
        
        if cp.get("average_query_time_ms", 0) > 200:
            recommendations.append("Optimize connection pooling performance:")
            recommendations.append("  - Tune pool size parameters")
            recommendations.append("  - Consider implementing query caching")
    
    # API recommendations
    if "api_tests" in results:
        api = results["api_tests"]
        
        if "endpoint_tests" in api:
            et = api["endpoint_tests"]
            if et["status"] != Status.SUCCESS:
                recommendations.append("Fix API endpoint issues:")
                for name, result in et.get("endpoint_results", {}).items():
                    if result["status"] != Status.SUCCESS:
                        recommendations.append(f"  - {name}: Check response code {result.get('response_code', 'N/A')}")
        
        if "role_tests" in api:
            rt = api["role_tests"]
            if rt["status"] != Status.SUCCESS:
                recommendations.append("Fix role-based access issues:")
                for name, result in rt.get("endpoint_results", {}).items():
                    if result["status"] != Status.SUCCESS:
                        recommendations.append(f"  - {name}: Verify RLS policies")
    
    # Add recommendations to report
    if recommendations:
        for rec in recommendations:
            report.append(rec)
    else:
        report.append("No specific recommendations. All tests passed successfully.")
    
    return "\n".join(report)

def save_report(results, report_text, output_dir="."):
    """
    Save test results and report to files.
    
    Args:
        results: Dictionary containing all test results
        report_text: Formatted report text
        output_dir: Directory to save the files
        
    Returns:
        tuple: Paths to the saved files
    """
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Generate timestamp for filenames
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save JSON results
    json_path = output_path / f"integration_test_results_{timestamp}.json"
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save text report
    report_path = output_path / f"integration_test_report_{timestamp}.txt"
    with open(report_path, 'w') as f:
        f.write(report_text)
    
    logger.info(f"Results saved to {json_path}")
    logger.info(f"Report saved to {report_path}")
    
    return json_path, report_path

def main():
    """Main function to run all tests."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Configure logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize memory cache
    init_memory_cache()
    
    # Initialize results dictionary
    results = {
        "timestamp": datetime.now().isoformat(),
        "status": Status.SUCCESS
    }
    
    # Run performance tests
    if not args.skip_performance:
        logger.info("Starting performance tests...")
        performance_results = run_performance_tests(args.supabase_url, args.supabase_key)
        results["performance_tests"] = performance_results
        
        if performance_results["status"] == Status.ERROR:
            results["status"] = Status.ERROR
        elif performance_results["status"] == Status.COMPLETED_WITH_WARNINGS and results["status"] != Status.ERROR:
            results["status"] = Status.COMPLETED_WITH_WARNINGS
    else:
        logger.info("Skipping performance tests...")
    
    # Run connection pooling tests
    if not args.skip_connection_pooling:
        logger.info("Starting connection pooling tests...")
        pooling_results = test_connection_pooling(
            args.supabase_url, 
            args.supabase_key,
            args.concurrent_users,
            args.test_duration
        )
        results["connection_pooling_tests"] = pooling_results
        
        if pooling_results["status"] == Status.ERROR:
            results["status"] = Status.ERROR
        elif pooling_results["status"] == Status.COMPLETED_WITH_WARNINGS and results["status"] != Status.ERROR:
            results["status"] = Status.COMPLETED_WITH_WARNINGS
    else:
        logger.info("Skipping connection pooling tests...")
    
    # Run API integration tests
    if not args.skip_api:
        logger.info("Starting API integration tests...")
        api_results = run_api_tests(args.api_url)
        results["api_tests"] = api_results
        
        if api_results["status"] == Status.ERROR:
            results["status"] = Status.ERROR
        elif api_results["status"] == Status.COMPLETED_WITH_WARNINGS and results["status"] != Status.ERROR:
            results["status"] = Status.COMPLETED_WITH_WARNINGS
    else:
        logger.info("Skipping API integration tests...")
    
    # Generate report
    logger.info("Generating report...")
    report_text = generate_report(results)
    
    # Save results and report
    json_path, report_path = save_report(results, report_text)
    
    # Print summary
    print("\n" + "=" * 80)
    print(f"Integration Test Results: {results['status']}")
    print("=" * 80)
    print(f"Detailed report saved to: {report_path}")
    print(f"JSON results saved to: {json_path}")
    
    # Return exit code based on status
    if results["status"] == Status.ERROR:
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
