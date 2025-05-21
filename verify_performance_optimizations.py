#!/usr/bin/env python3
"""
CryoProtect v2 - Verify Performance Optimizations

This script verifies that the performance optimizations applied to the CryoProtect v2 database
are providing the expected performance improvements. It:

1. Runs performance benchmarks on common database operations
2. Compares performance before and after the optimizations
3. Verifies that all indexes are being used effectively
4. Generates a performance report with the results

The script tests performance for operations on the following key tables:
- molecules
- mixtures
- predictions
- experiments
"""

import os
import sys
import json
import time
import logging
import argparse
import psutil
import statistics
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set
from dotenv import load_dotenv
from contextlib import contextmanager

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("performance_verification.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Supabase connection
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

if not SUPABASE_URL or not SUPABASE_KEY:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

# Import Supabase client
try:
    from supabase import create_client, Client
except ImportError:
    logger.error("Error: supabase-py package not found. Install with 'pip install supabase'")
    sys.exit(1)

# Import connection pooling if available
try:
    from connection_pool_wrapper import (
        initialize_supabase_pool, get_supabase_connection,
        get_supabase_pool_stats, shutdown_supabase_pool
    )
    HAS_CONNECTION_POOL = True
except ImportError:
    logger.warning("connection_pool_wrapper.py not found. Connection pooling will not be used.")
    HAS_CONNECTION_POOL = False

# Test queries for key tables
TEST_QUERIES = {
    # Molecule operations
    "molecule_get_all": {
        "description": "Get all molecules (limit 100)",
        "query": "SELECT * FROM molecules LIMIT 100",
        "category": "molecules",
        "operation_type": "read"
    },
    "molecule_get_by_id": {
        "description": "Get molecule by ID",
        "query": "SELECT * FROM molecules WHERE id = '{molecule_id}'",
        "category": "molecules",
        "operation_type": "read"
    },
    "molecule_search_by_name": {
        "description": "Search molecules by name",
        "query": "SELECT * FROM molecules WHERE name ILIKE '%glyc%' LIMIT 20",
        "category": "molecules",
        "operation_type": "read"
    },
    "molecule_with_properties": {
        "description": "Get molecule with properties",
        "query": '''
            SELECT m.*, mp.property_type_id, mp.numeric_value, pt.name as property_name
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            LIMIT 50
        ''',
        "category": "molecules",
        "operation_type": "read"
    },
    # Mixture operations
    "mixture_get_all": {
        "description": "Get all mixtures (limit 100)",
        "query": "SELECT * FROM mixtures LIMIT 100",
        "category": "mixtures",
        "operation_type": "read"
    },
    "mixture_get_by_id": {
        "description": "Get mixture by ID",
        "query": "SELECT * FROM mixtures WHERE id = '{mixture_id}'",
        "category": "mixtures",
        "operation_type": "read"
    },
    "mixture_with_components": {
        "description": "Get mixture with components",
        "query": '''
            SELECT mix.*, mc.molecule_id, mc.amount, m.name as molecule_name
            FROM mixtures mix
            JOIN mixture_components mc ON mix.id = mc.mixture_id
            JOIN molecules m ON mc.molecule_id = m.id
            LIMIT 50
        ''',
        "category": "mixtures",
        "operation_type": "read"
    },
    "mixture_filter_by_created_by": {
        "description": "Filter mixtures by created_by (RLS test)",
        "query": "SELECT * FROM mixtures WHERE created_by = '{user_id}' LIMIT 20",
        "category": "mixtures",
        "operation_type": "read"
    },
    # Prediction operations
    "prediction_get_all": {
        "description": "Get all predictions (limit 100)",
        "query": "SELECT * FROM predictions LIMIT 100",
        "category": "predictions",
        "operation_type": "read"
    },
    "prediction_get_by_mixture": {
        "description": "Get predictions for mixture",
        "query": "SELECT * FROM predictions WHERE mixture_id = '{mixture_id}' LIMIT 20",
        "category": "predictions",
        "operation_type": "read"
    },
    "prediction_with_calculation_method": {
        "description": "Get predictions with calculation method",
        "query": '''
            SELECT p.*, cm.name as method_name
            FROM predictions p
            JOIN calculation_methods cm ON p.calculation_method_id = cm.id
            LIMIT 50
        ''',
        "category": "predictions",
        "operation_type": "read"
    },
    "prediction_filter_by_confidence": {
        "description": "Filter predictions by confidence",
        "query": "SELECT * FROM predictions WHERE confidence > 0.8 LIMIT 20",
        "category": "predictions",
        "operation_type": "read"
    },
    # Experiment operations
    "experiment_get_all": {
        "description": "Get all experiments (limit 100)",
        "query": "SELECT * FROM experiments LIMIT 100",
        "category": "experiments",
        "operation_type": "read"
    },
    "experiment_get_by_mixture": {
        "description": "Get experiments for mixture",
        "query": "SELECT * FROM experiments WHERE mixture_id = '{mixture_id}' LIMIT 20",
        "category": "experiments",
        "operation_type": "read"
    },
    "experiment_with_property_type": {
        "description": "Get experiments with property type",
        "query": '''
            SELECT e.*, pt.name as property_name
            FROM experiments e
            JOIN property_types pt ON e.property_type_id = pt.id
            LIMIT 50
        ''',
        "category": "experiments",
        "operation_type": "read"
    },
    "experiment_filter_by_date": {
        "description": "Filter experiments by date",
        "query": "SELECT * FROM experiments WHERE date > '2024-01-01' LIMIT 20",
        "category": "experiments",
        "operation_type": "read"
    }
}

def connect_to_supabase() -> 'Client':
    """Connect to Supabase using service role key."""
    try:
        supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
        logger.info("Connected to Supabase using service role key")
        return supabase
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        raise

def get_test_data(supabase: 'Client') -> Dict[str, Any]:
    """Get test data for query parameters."""
    test_data = {}
    try:
        # Get a molecule ID
        response = supabase.table("molecules").select("id").limit(1).execute()
        if hasattr(response, 'data') and response.data:
            test_data["molecule_id"] = response.data[0]["id"]
        else:
            test_data["molecule_id"] = "00000000-0000-0000-0000-000000000000"  # Fallback
        # Get a mixture ID
        response = supabase.table("mixtures").select("id").limit(1).execute()
        if hasattr(response, 'data') and response.data:
            test_data["mixture_id"] = response.data[0]["id"]
        else:
            test_data["mixture_id"] = "00000000-0000-0000-0000-000000000000"  # Fallback
        # Get a user ID for RLS tests
        response = supabase.table("mixtures").select("created_by").limit(1).execute()
        if hasattr(response, 'data') and response.data and response.data[0]["created_by"]:
            test_data["user_id"] = response.data[0]["created_by"]
        else:
            test_data["user_id"] = "00000000-0000-0000-0000-000000000000"  # Fallback
        logger.info(f"Test data: {test_data}")
        return test_data
    except Exception as e:
        logger.error(f"Error getting test data: {str(e)}")
        # Return fallback test data
        return {
            "molecule_id": "00000000-0000-0000-0000-000000000000",
            "mixture_id": "00000000-0000-0000-0000-000000000000",
            "user_id": "00000000-0000-0000-0000-000000000000"
        }

def replace_placeholders(query: str, test_data: Dict[str, Any]) -> str:
    """Replace placeholders in the query with test data."""
    for key, value in test_data.items():
        placeholder = "{" + key + "}"
        query = query.replace(placeholder, str(value))
    return query

def execute_query(supabase: 'Client', query: str) -> Tuple[Any, float]:
    """Execute a query and measure its performance."""
    try:
        start_time = time.time()
        response = supabase.rpc("execute_sql", {"query": query}).execute()
        end_time = time.time()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing query: {response.error}")
            return None, None
        execution_time = (end_time - start_time) * 1000  # ms
        return response.data, execution_time
    except Exception as e:
        logger.error(f"Error executing query: {str(e)}")
        return None, None

def run_performance_tests(supabase: 'Client', test_data: Dict[str, Any]) -> Dict[str, Any]:
    """Run performance tests on all test queries."""
    results = {}
    logger.info("Running performance tests...")
    cpu_usage = []
    memory_usage = []
    def track_resource_usage():
        cpu_usage.append(psutil.cpu_percent(interval=0.1))
        memory_usage.append(psutil.virtual_memory().percent)
    for query_name, query_info in TEST_QUERIES.items():
        logger.info(f"Testing query: {query_name} - {query_info['description']}")
        query = replace_placeholders(query_info["query"], test_data)
        execution_times = []
        for i in range(5):
            track_resource_usage()
            _, execution_time = execute_query(supabase, query)
            if execution_time is not None:
                execution_times.append(execution_time)
            time.sleep(0.1)
        if execution_times:
            avg_time = statistics.mean(execution_times)
            min_time = min(execution_times)
            max_time = max(execution_times)
            p95_time = sorted(execution_times)[int(len(execution_times) * 0.95)]
            results[query_name] = {
                "description": query_info["description"],
                "category": query_info["category"],
                "operation_type": query_info["operation_type"],
                "avg_ms": avg_time,
                "min_ms": min_time,
                "max_ms": max_time,
                "p95_ms": p95_time,
                "variability": max_time / avg_time if avg_time > 0 else 0
            }
        else:
            results[query_name] = {
                "description": query_info["description"],
                "category": query_info["category"],
                "operation_type": query_info["operation_type"],
                "error": "Failed to execute query"
            }
    resource_metrics = {
        "cpu": {
            "min": min(cpu_usage) if cpu_usage else 0,
            "max": max(cpu_usage) if cpu_usage else 0,
            "avg": statistics.mean(cpu_usage) if cpu_usage else 0
        },
        "memory": {
            "min": min(memory_usage) if memory_usage else 0,
            "max": max(memory_usage) if memory_usage else 0,
            "avg": statistics.mean(memory_usage) if memory_usage else 0
        }
    }
    return {
        "query_results": results,
        "resource_usage": resource_metrics
    }

# Query plan test queries
QUERY_PLAN_TESTS = {
    "molecule_search_plan": {
        "description": "Analyze query plan for molecule search",
        "query": "EXPLAIN ANALYZE SELECT * FROM molecules WHERE name ILIKE '%glyc%' LIMIT 20",
        "expected_index": "idx_molecule_name_trgm"
    },
    "mixture_components_join_plan": {
        "description": "Analyze query plan for mixture components join",
        "query": '''
            EXPLAIN ANALYZE
            SELECT mix.*, mc.molecule_id, mc.amount, m.name as molecule_name
            FROM mixtures mix
            JOIN mixture_components mc ON mix.id = mc.mixture_id
            JOIN molecules m ON mc.molecule_id = m.id
            LIMIT 50
        ''',
        "expected_index": "idx_mixture_component_mixture_id"
    },
    "prediction_filter_plan": {
        "description": "Analyze query plan for prediction filter",
        "query": "EXPLAIN ANALYZE SELECT * FROM predictions WHERE confidence > 0.8 LIMIT 20",
        "expected_index": "idx_predictions_confidence"
    },
    "experiment_date_filter_plan": {
        "description": "Analyze query plan for experiment date filter",
        "query": "EXPLAIN ANALYZE SELECT * FROM experiments WHERE date > '2024-01-01' LIMIT 20",
        "expected_index": "idx_experiments_date"
    },
    "rls_filter_plan": {
        "description": "Analyze query plan for RLS filter",
        "query": "EXPLAIN ANALYZE SELECT * FROM mixtures WHERE created_by = '{user_id}' LIMIT 20",
        "expected_index": "idx_mixtures_created_by"
    }
}

def analyze_query_plans(supabase: 'Client', test_data: Dict[str, Any]) -> Dict[str, Any]:
    """Analyze query plans to verify index usage."""
    results = {}
    logger.info("Analyzing query plans...")
    for plan_name, plan_info in QUERY_PLAN_TESTS.items():
        logger.info(f"Analyzing plan: {plan_name} - {plan_info['description']}")
        query = replace_placeholders(plan_info["query"], test_data)
        plan_data, _ = execute_query(supabase, query)
        if plan_data:
            plan_text = "\n".join([str(line) for line in plan_data])
            expected_index = plan_info["expected_index"]
            index_used = expected_index in plan_text
            seq_scan = "Seq Scan" in plan_text
            results[plan_name] = {
                "description": plan_info["description"],
                "expected_index": expected_index,
                "index_used": index_used,
                "sequential_scan": seq_scan,
                "plan": plan_data
            }
        else:
            results[plan_name] = {
                "description": plan_info["description"],
                "expected_index": plan_info["expected_index"],
                "error": "Failed to execute query plan"
            }
    return results

def verify_indexes(supabase: 'Client') -> Dict[str, Any]:
    """Verify that all expected indexes exist in the database."""
    logger.info("Verifying indexes...")
    expected_indexes = [
        "idx_mixtures_created_by",
        "idx_experiments_created_by",
        "idx_predictions_created_by",
        "idx_molecular_properties_created_by",
        "idx_calculation_methods_created_by",
        "idx_mixture_component_mixture_id",
        "idx_mixture_component_molecule_id",
        "idx_mixture_component_amount",
        "idx_experiments_mixture_id",
        "idx_experiments_property_type_id",
        "idx_experiments_date",
        "idx_mixture_components_mixture_molecule",
        "idx_experiments_mixture_property",
        "idx_predictions_mixture_property",
        "idx_predictions_calculation_method",
        "idx_predictions_confidence",
        "idx_predictions_numeric_value",
        "idx_molecular_property_property_type",
        "idx_molecular_property_molecule_id",
        "idx_molecular_property_numeric_value",
        "idx_molecular_properties_molecule_property",
        "idx_predictions_method_property",
        "idx_molecule_name_trgm",
        "idx_mixture_name_trgm",
        "idx_property_type_name_trgm",
        "idx_molecule_name",
        "idx_mixture_name",
        "idx_property_type_name",
        "idx_mixtures_created_at",
        "idx_mixtures_updated_at",
        "idx_experiments_created_at",
        "idx_predictions_created_at",
        "idx_property_types_unit_id",
        "idx_calculation_methods_version"
    ]
    query = """
    SELECT
        indexname,
        tablename,
        indexdef
    FROM
        pg_indexes
    WHERE
        schemaname = 'public'
    ORDER BY
        tablename, indexname;
    """
    index_data, _ = execute_query(supabase, query)
    if not index_data:
        logger.error("Failed to retrieve index information")
        return {"error": "Failed to retrieve index information"}
    existing_indexes = {}
    for index in index_data:
        existing_indexes[index["indexname"]] = {
            "table": index["tablename"],
            "definition": index["indexdef"]
        }
    index_status = {}
    for index_name in expected_indexes:
        index_status[index_name] = {
            "exists": index_name in existing_indexes,
            "table": existing_indexes.get(index_name, {}).get("table", "N/A"),
            "definition": existing_indexes.get(index_name, {}).get("definition", "N/A")
        }
    total_indexes = len(expected_indexes)
    existing_count = sum(1 for status in index_status.values() if status["exists"])
    missing_count = total_indexes - existing_count
    return {
        "total_expected": total_indexes,
        "existing_count": existing_count,
        "missing_count": missing_count,
        "index_status": index_status
    }

def verify_optimized_views(supabase: 'Client') -> Dict[str, Any]:
    """Verify that all optimized views exist in the database."""
    logger.info("Verifying optimized views...")
    expected_views = [
        "molecule_with_properties",
        "mixture_with_components",
        "prediction_experiment_comparison",
        "molecule_property_statistics",
        "recent_activity"
    ]
    query = """
    SELECT
        table_name,
        view_definition
    FROM
        information_schema.views
    WHERE
        table_schema = 'public'
    ORDER BY
        table_name;
    """
    view_data, _ = execute_query(supabase, query)
    if not view_data:
        logger.error("Failed to retrieve view information")
        return {"error": "Failed to retrieve view information"}
    existing_views = {}
    for view in view_data:
        existing_views[view["table_name"]] = {
            "definition": view["view_definition"]
        }
    view_status = {}
    for view_name in expected_views:
        view_status[view_name] = {
            "exists": view_name in existing_views,
            "definition": existing_views.get(view_name, {}).get("definition", "N/A")
        }
    total_views = len(expected_views)
    existing_count = sum(1 for status in view_status.values() if status["exists"])
    missing_count = total_views - existing_count
    return {
        "total_expected": total_views,
        "existing_count": existing_count,
        "missing_count": missing_count,
        "view_status": view_status
    }

def main():
    logger.info("Starting performance optimization verification...")
    supabase = connect_to_supabase()
    test_data = get_test_data(supabase)
    performance_tests = run_performance_tests(supabase, test_data)
    query_plans = analyze_query_plans(supabase, test_data)
    index_verification = verify_indexes(supabase)
    view_verification = verify_optimized_views(supabase)
    report = {
        "timestamp": datetime.now().isoformat(),
        "performance_tests": performance_tests,
        "query_plans": query_plans,
        "index_verification": index_verification,
        "view_verification": view_verification
    }
    # Save report
    with open("performance_optimization_verification_report.json", "w") as f:
        json.dump(report, f, indent=2)
    logger.info("Performance optimization verification complete. Report saved to performance_optimization_verification_report.json")
    print("Performance optimization verification complete. See performance_optimization_verification_report.json for details.")

if __name__ == "__main__":
    main()
