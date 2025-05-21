#!/usr/bin/env python3
"""
CryoProtect v2 - Database Performance Optimization

This script identifies and creates necessary indexes for database performance optimization
using the direct PostgreSQL connection. It measures query performance before and after
applying the indexes and generates a detailed report.

Features:
1. Identifies tables and columns that would benefit from indexing
2. Creates appropriate indexes for commonly queried fields
3. Measures query performance before and after index creation
4. Provides transaction support for safe operations
5. Includes detailed logging and reporting
6. Implements error handling and rollback mechanisms
"""

import os
import sys
import json
import time
import logging
import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set
import statistics
from contextlib import contextmanager

# Import the direct PostgreSQL connection helper and SQL executor
from postgres_direct import PostgresDirectConnection
import sql_executor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("performance_optimization.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("performance_optimization")

# Common queries for performance testing
COMMON_QUERIES = [
    {
        "name": "get_molecules",
        "query": """
            SELECT * FROM molecules 
            ORDER BY created_at DESC 
            LIMIT 100
        """
    },
    {
        "name": "get_molecule_by_id",
        "query": """
            SELECT m.*, mp.property_type_id, mp.numeric_value, mp.text_value, mp.boolean_value
            FROM molecules m
            LEFT JOIN molecular_properties mp ON m.id = mp.molecule_id
            WHERE m.id = %s
        """,
        "params": ["random_molecule_id"]
    },
    {
        "name": "get_mixtures",
        "query": """
            SELECT * FROM mixtures 
            ORDER BY created_at DESC 
            LIMIT 100
        """
    },
    {
        "name": "get_mixture_by_id",
        "query": """
            SELECT mix.*, mc.molecule_id, mc.concentration, mc.concentration_unit
            FROM mixtures mix
            LEFT JOIN mixture_components mc ON mix.id = mc.mixture_id
            WHERE mix.id = %s
        """,
        "params": ["random_mixture_id"]
    },
    {
        "name": "get_predictions",
        "query": """
            SELECT p.*, pt.name as property_name, cm.name as method_name
            FROM predictions p
            JOIN property_types pt ON p.property_type_id = pt.id
            JOIN calculation_methods cm ON p.calculation_method_id = cm.id
            WHERE p.mixture_id = %s
            LIMIT 100
        """,
        "params": ["random_mixture_id"]
    },
    {
        "name": "get_experiments",
        "query": """
            SELECT e.*, pt.name as property_name
            FROM experiments e
            JOIN property_types pt ON e.property_type_id = pt.id
            WHERE e.mixture_id = %s
        """,
        "params": ["random_mixture_id"]
    },
    {
        "name": "search_molecules_by_name",
        "query": """
            SELECT * FROM molecules
            WHERE name ILIKE %s
            LIMIT 50
        """,
        "params": ["%cryo%"]
    },
    {
        "name": "filter_by_property",
        "query": """
            SELECT m.*, mp.numeric_value
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE pt.name = %s AND mp.numeric_value > %s
            LIMIT 100
        """,
        "params": ["LogP", "2.0"]
    }
]

class PerformanceMetrics:
    """Class for collecting and analyzing performance metrics."""
    
    def __init__(self):
        self.operation_times = {}
        self.start_time = None
        self.end_time = None
    
    def start_test(self):
        """Start the performance test."""
        self.start_time = datetime.now()
    
    def end_test(self):
        """End the performance test."""
        self.end_time = datetime.now()
    
    def record_operation_time(self, operation, time_ms):
        """Record the time taken for an operation."""
        if operation not in self.operation_times:
            self.operation_times[operation] = []
        self.operation_times[operation].append(time_ms)
    
    def get_summary(self):
        """Get a summary of the performance metrics."""
        summary = {
            "test_duration": (self.end_time - self.start_time).total_seconds() if self.end_time else None,
            "operations": {}
        }
        
        for operation, times in self.operation_times.items():
            if times:
                sorted_times = sorted(times)
                p95_index = int(len(times) * 0.95)
                p99_index = int(len(times) * 0.99)
                
                summary["operations"][operation] = {
                    "count": len(times),
                    "min_ms": min(times),
                    "max_ms": max(times),
                    "avg_ms": statistics.mean(times),
                    "median_ms": statistics.median(times),
                    "p95_ms": sorted_times[p95_index] if p95_index < len(times) else sorted_times[-1],
                    "p99_ms": sorted_times[p99_index] if p99_index < len(times) else sorted_times[-1]
                }
        
        return summary

def transaction():
    """Context manager for database transactions.
    
    This function uses the transaction context manager from sql_executor
    to ensure compatibility with all database adapter types.
    """
    return sql_executor.transaction()

def get_random_ids():
    """Get random IDs from the database for testing."""
    random_ids = {
        "random_molecule_id": None,
        "random_mixture_id": None
    }
    
    try:
        # Get a random molecule ID
        result = sql_executor.execute_query(
            "SELECT id FROM molecules ORDER BY RANDOM() LIMIT 1"
        )
        if result and len(result) > 0:
            random_ids["random_molecule_id"] = result[0]["id"]
        
        # Get a random mixture ID
        result = sql_executor.execute_query(
            "SELECT id FROM mixtures ORDER BY RANDOM() LIMIT 1"
        )
        if result and len(result) > 0:
            random_ids["random_mixture_id"] = result[0]["id"]
    
    except Exception as e:
        logger.warning(f"Could not get random IDs: {str(e)}")
    
    return random_ids

def run_performance_test():
    """Run performance tests on common queries."""
    metrics = PerformanceMetrics()
    metrics.start_test()
    
    # Get random IDs for parameterized queries
    random_ids = get_random_ids()
    
    # Run each query multiple times
    for query_info in COMMON_QUERIES:
        query_name = query_info["name"]
        query = query_info["query"]
        
        # Prepare parameters if needed
        params = None
        if "params" in query_info:
            params = []
            for param in query_info["params"]:
                if param in random_ids:
                    params.append(random_ids[param])
                else:
                    params.append(param)
        
        # Run the query multiple times
        for _ in range(5):  # Run each query 5 times
            try:
                start_time = time.time()
                sql_executor.execute_query(query, params)
                end_time = time.time()
                
                # Record time in milliseconds
                elapsed_ms = (end_time - start_time) * 1000
                metrics.record_operation_time(query_name, elapsed_ms)
                
            except Exception as e:
                logger.error(f"Error running query {query_name}: {str(e)}")
    
    metrics.end_test()
    return metrics

def analyze_metrics(metrics):
    """Analyze the performance metrics and identify bottlenecks."""
    analysis = {
        "bottlenecks": [],
        "slow_operations": []
    }
    
    # Check for slow operations (p95 > 500ms)
    for operation, stats in metrics["operations"].items():
        if stats["p95_ms"] > 500:
            analysis["slow_operations"].append({
                "operation": operation,
                "p95_ms": stats["p95_ms"],
                "avg_ms": stats["avg_ms"]
            })
    
    # Identify potential bottlenecks
    for operation, stats in metrics["operations"].items():
        # Operations with high variability (max/avg > 5)
        if stats["max_ms"] / stats["avg_ms"] > 5:
            analysis["bottlenecks"].append({
                "operation": operation,
                "type": "high_variability",
                "max_ms": stats["max_ms"],
                "avg_ms": stats["avg_ms"],
                "ratio": stats["max_ms"] / stats["avg_ms"]
            })
    
    return analysis

def generate_recommendations(metrics):
    """Generate recommendations based on the performance metrics."""
    recommendations = []
    
    # Check for slow read operations
    for operation, stats in metrics["operations"].items():
        if operation.startswith("get_") and stats["p95_ms"] > 500:
            if "molecules" in operation:
                recommendations.append({
                    "operation": operation,
                    "issue": "Slow read performance",
                    "recommendation": "Consider adding indexes on commonly queried fields in the molecule and molecular_property tables."
                })
            elif "mixtures" in operation:
                recommendations.append({
                    "operation": operation,
                    "issue": "Slow read performance",
                    "recommendation": "Consider adding indexes on commonly queried fields in the mixture and mixture_component tables."
                })
            elif "predictions" in operation:
                recommendations.append({
                    "operation": operation,
                    "issue": "Slow read performance",
                    "recommendation": "Consider adding indexes on mixture_id, property_type_id, and calculation_method_id in the predictions table."
                })
            elif "experiments" in operation:
                recommendations.append({
                    "operation": operation,
                    "issue": "Slow read performance",
                    "recommendation": "Consider adding indexes on mixture_id and property_type_id in the experiments table."
                })
    
    # General recommendations
    recommendations.append({
        "issue": "General optimization",
        "recommendation": "Consider implementing query caching for frequently accessed data."
    })
    
    return recommendations

def get_performance_indexes_sql(schema="public"):
    """
    Generate SQL statements for creating performance indexes.
    
    Args:
        schema (str): The database schema to modify
        
    Returns:
        dict: SQL statements for creating performance indexes
    """
    # Batch 1: RLS performance indexes
    rls_indexes = f"""
    -- Add index for RLS performance
    CREATE INDEX IF NOT EXISTS idx_mixtures_created_by ON {schema}.mixtures(created_by);
    CREATE INDEX IF NOT EXISTS idx_experiments_created_by ON {schema}.experiments(created_by);
    CREATE INDEX IF NOT EXISTS idx_predictions_created_by ON {schema}.predictions(created_by);
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_created_by ON {schema}.molecular_properties(created_by);
    CREATE INDEX IF NOT EXISTS idx_calculation_methods_created_by ON {schema}.calculation_methods(created_by);
    """
    
    # Batch 2: Mixture and component indexes
    mixture_indexes = f"""
    -- Add index for mixture components
    CREATE INDEX IF NOT EXISTS idx_mixture_components_mixture_id ON {schema}.mixture_components(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_mixture_components_molecule_id ON {schema}.mixture_components(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_mixture_components_concentration ON {schema}.mixture_components(concentration);
    
    -- Add index for experiments
    CREATE INDEX IF NOT EXISTS idx_experiments_mixture_id ON {schema}.experiments(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_experiments_property_type_id ON {schema}.experiments(property_type_id);
    CREATE INDEX IF NOT EXISTS idx_experiments_date ON {schema}.experiments(date_performed);
    
    -- Add composite indexes for common query patterns
    CREATE INDEX IF NOT EXISTS idx_mixture_components_mixture_molecule ON {schema}.mixture_components(mixture_id, molecule_id);
    CREATE INDEX IF NOT EXISTS idx_experiments_mixture_property ON {schema}.experiments(mixture_id, property_type_id);
    """
    
    # Batch 3: Prediction and property indexes
    prediction_indexes = f"""
    -- Add composite index for predictions
    CREATE INDEX IF NOT EXISTS idx_predictions_mixture_property ON {schema}.predictions(mixture_id, property_type_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_calculation_method ON {schema}.predictions(calculation_method_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_confidence ON {schema}.predictions(confidence);
    CREATE INDEX IF NOT EXISTS idx_predictions_numeric_value ON {schema}.predictions(numeric_value);
    
    -- Add index for property types
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_type ON {schema}.molecular_properties(property_type_id);
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON {schema}.molecular_properties(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_numeric_value ON {schema}.molecular_properties(numeric_value);
    
    -- Add composite indexes for common query patterns
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_property ON {schema}.molecular_properties(molecule_id, property_type_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_method_property ON {schema}.predictions(calculation_method_id, property_type_id);
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
    """
    
    # Batch 5: Additional performance indexes
    additional_indexes = f"""
    -- Add indexes for timestamp columns for time-based queries
    CREATE INDEX IF NOT EXISTS idx_mixtures_created_at ON {schema}.mixtures(created_at);
    CREATE INDEX IF NOT EXISTS idx_mixtures_updated_at ON {schema}.mixtures(updated_at);
    CREATE INDEX IF NOT EXISTS idx_experiments_created_at ON {schema}.experiments(created_at);
    CREATE INDEX IF NOT EXISTS idx_predictions_created_at ON {schema}.predictions(created_at);
    
    -- Add indexes for foreign keys that might not be automatically indexed
    CREATE INDEX IF NOT EXISTS idx_property_types_units ON {schema}.property_types(units);
    CREATE INDEX IF NOT EXISTS idx_calculation_methods_version ON {schema}.calculation_methods(version);
    """
    
    # Batch 6: Performance optimization indexes for verification queries
    performance_indexes = f"""
    -- Add indexes for verification queries
    CREATE INDEX IF NOT EXISTS idx_molecules_created_at ON {schema}.molecules(created_at);
    CREATE INDEX IF NOT EXISTS idx_molecules_lower_name ON {schema}.molecules(LOWER(name));
    
    -- Add composite indexes for join queries
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_property_value ON {schema}.molecular_properties(molecule_id, property_type_id, numeric_value);
    CREATE INDEX IF NOT EXISTS idx_property_types_name_lower ON {schema}.property_types(LOWER(name));
    
    -- Add covering indexes for common query patterns
    CREATE INDEX IF NOT EXISTS idx_molecules_all_fields ON {schema}.molecules(id, name, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid, chembl_id);
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_all_fields ON {schema}.molecular_properties(id, molecule_id, property_type_id, numeric_value, text_value, boolean_value);
    """
    
    return {
        "rls_indexes": rls_indexes,
        "mixture_indexes": mixture_indexes,
        "prediction_indexes": prediction_indexes,
        "text_search_indexes": text_search_indexes,
        "additional_indexes": additional_indexes,
        "performance_indexes": performance_indexes
    }

def verify_index_creation(schema, index_name):
    """
    Verify that an index was created successfully.
    
    Args:
        schema (str): The database schema
        index_name (str): The name of the index to verify
        
    Returns:
        bool: True if the index exists, False otherwise
    """
    try:
        query = """
        SELECT 1 FROM pg_indexes 
        WHERE schemaname = %s AND indexname = %s
        """
        result = sql_executor.execute_query(query, [schema, index_name])
        return result and len(result) > 0
    except Exception as e:
        logger.error(f"Error verifying index {index_name}: {str(e)}")
        return False

def apply_performance_indexes(schema="public", test_only=False):
    """
    Apply performance indexes to the database.
    
    Args:
        schema (str): The database schema to modify
        test_only (bool): If True, only test the indexes without creating them
        
    Returns:
        dict: Results of the operation
    """
    results = {
        "before_metrics": None,
        "after_metrics": None,
        "analysis_before": None,
        "analysis_after": None,
        "recommendations_before": None,
        "created_indexes": [],
        "failed_indexes": [],
        "errors": []
    }
    
    try:
        logger.info("Starting performance optimization")
        
        # Run performance tests before applying indexes
        logger.info("Running performance tests before applying indexes")
        before_metrics = run_performance_test()
        results["before_metrics"] = before_metrics.get_summary()
        results["analysis_before"] = analyze_metrics(results["before_metrics"])
        results["recommendations_before"] = generate_recommendations(results["before_metrics"])
        
        if test_only:
            logger.info("Test only mode, skipping index creation")
            return results
        
        # Get SQL statements for creating indexes
        index_batches = get_performance_indexes_sql(schema)
        
        # Log the number of index batches
        logger.info(f"Applying {len(index_batches)} index batches")
        
        # Apply each batch of indexes
        for batch_name, batch_sql in index_batches.items():
            logger.info(f"Applying {batch_name}")
            
            try:
                # Execute the batch in a transaction
                with transaction():
                    sql_executor.execute_query(batch_sql)
                
                # Extract index names from the SQL for verification
                index_names = []
                for line in batch_sql.split("\n"):
                    if "CREATE INDEX" in line and "idx_" in line:
                        parts = line.split("idx_")
                        if len(parts) > 1:
                            index_name = "idx_" + parts[1].split(" ")[0].strip()
                            if index_name.endswith(";"):
                                index_name = index_name[:-1]
                            index_names.append(index_name)
                
                # Verify each index was created
                for index_name in index_names:
                    if verify_index_creation(schema, index_name):
                        results["created_indexes"].append(index_name)
                        logger.info(f"Successfully created index {index_name}")
                    else:
                        results["failed_indexes"].append(index_name)
                        logger.warning(f"Failed to verify index {index_name}")
            
            except Exception as e:
                error_msg = f"Error applying {batch_name}: {str(e)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)
        
        # Run performance tests after applying indexes
        logger.info("Running performance tests after applying indexes")
        after_metrics = run_performance_test()
        results["after_metrics"] = after_metrics.get_summary()
        results["analysis_after"] = analyze_metrics(results["after_metrics"])
        
        # Generate performance improvement report
        generate_performance_report(results)
        
        logger.info("Performance optimization completed")
        
    except Exception as e:
        error_msg = f"Error during performance optimization: {str(e)}"
        logger.error(error_msg)
        results["errors"].append(error_msg)
    
    return results

def generate_performance_report(results):
    """
    Generate a detailed performance report.
    
    Args:
        results (dict): Results of the performance optimization
    """
    report = {
        "timestamp": datetime.now().isoformat(),
        "schema": "public",
        "created_indexes": results["created_indexes"],
        "failed_indexes": results["failed_indexes"],
        "errors": results["errors"],
        "performance": {
            "before": results["before_metrics"],
            "after": results["after_metrics"]
        },
        "analysis": {
            "before": results["analysis_before"],
            "after": results["analysis_after"]
        },
        "recommendations_before": results["recommendations_before"],
        "improvements": calculate_improvements(results)
    }
    
    # Save the report to a file
    with open("performance_optimization_report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    # Generate a text report
    text_report = generate_text_report(report)
    with open("performance_optimization_report.txt", "w") as f:
        f.write(text_report)
    
    logger.info("Performance report generated")

def calculate_improvements(results):
    """
    Calculate performance improvements between before and after metrics.
    
    Args:
        results (dict): Results of the performance optimization
        
    Returns:
        dict: Performance improvements
    """
    improvements = {}
    
    if not results["before_metrics"] or not results["after_metrics"]:
        return improvements
    
    before_ops = results["before_metrics"]["operations"]
    after_ops = results["after_metrics"]["operations"]
    
    for op_name in before_ops:
        if op_name in after_ops:
            before = before_ops[op_name]
            after = after_ops[op_name]
            
            # Calculate percentage improvements
            avg_improvement = ((before["avg_ms"] - after["avg_ms"]) / before["avg_ms"]) * 100 if before["avg_ms"] > 0 else 0
            p95_improvement = ((before["p95_ms"] - after["p95_ms"]) / before["p95_ms"]) * 100 if before["p95_ms"] > 0 else 0
            
            improvements[op_name] = {
                "avg_ms_before": before["avg_ms"],
                "avg_ms_after": after["avg_ms"],
                "avg_improvement_pct": avg_improvement,
                "p95_ms_before": before["p95_ms"],
                "p95_ms_after": after["p95_ms"],
                "p95_improvement_pct": p95_improvement
            }
    
    return improvements

def generate_text_report(report):
    """
    Generate a text report from the JSON report.
    
    Args:
        report (dict): The JSON report
        
    Returns:
        str: The text report
    """
    text_report = []
    
    text_report.append("=" * 80)
    text_report.append("CryoProtect v2 Database Performance Optimization Report")
    text_report.append("=" * 80)
    text_report.append(f"Timestamp: {report['timestamp']}")
    text_report.append("")
    
    # Created indexes
    text_report.append("Created Indexes:")
    if report["created_indexes"]:
        for index in report["created_indexes"]:
            text_report.append(f"  - {index}")
    else:
        text_report.append("  None")
    text_report.append("")
    
    # Failed indexes
    if report["failed_indexes"]:
        text_report.append("Failed Indexes:")
        for index in report["failed_indexes"]:
            text_report.append(f"  - {index}")
        text_report.append("")
    
    # Errors
    if report["errors"]:
        text_report.append("Errors:")
        for error in report["errors"]:
            text_report.append(f"  - {error}")
        text_report.append("")
    
    # Performance improvements
    text_report.append("Performance Improvements:")
    for op_name, improvement in report["improvements"].items():
        text_report.append(f"  {op_name}:")
        text_report.append(f"    Avg: {improvement['avg_ms_before']:.2f}ms -> {improvement['avg_ms_after']:.2f}ms ({improvement['avg_improvement_pct']:.2f}% improvement)")
        text_report.append(f"    P95: {improvement['p95_ms_before']:.2f}ms -> {improvement['p95_ms_after']:.2f}ms ({improvement['p95_improvement_pct']:.2f}% improvement)")
    text_report.append("")
    
    # Before/After operation performance
    text_report.append("Operation Performance (Before):")
    for operation, stats in report["performance"]["before"]["operations"].items():
        text_report.append(f"  {operation}:")
        text_report.append(f"    Avg: {stats['avg_ms']:.2f}ms")
        text_report.append(f"    P95: {stats['p95_ms']:.2f}ms")
    text_report.append("")
    
    text_report.append("Operation Performance (After):")
    for operation, stats in report["performance"]["after"]["operations"].items():
        text_report.append(f"  {operation}:")
        text_report.append(f"    Avg: {stats['avg_ms']:.2f}ms")
        text_report.append(f"    P95: {stats['p95_ms']:.2f}ms")
    text_report.append("")
    
    # Bottlenecks before
    if report["analysis"]["before"]["bottlenecks"]:
        text_report.append("Bottlenecks (Before):")
        for bottleneck in report["analysis"]["before"]["bottlenecks"]:
            text_report.append(f"  {bottleneck['operation']} - {bottleneck['type']}")
            text_report.append(f"    Max: {bottleneck['max_ms']:.2f}ms, Avg: {bottleneck['avg_ms']:.2f}ms, Ratio: {bottleneck['ratio']:.2f}")
        text_report.append("")
    
    # Bottlenecks after
    if report["analysis"]["after"]["bottlenecks"]:
        text_report.append("Bottlenecks (After):")
        for bottleneck in report["analysis"]["after"]["bottlenecks"]:
            text_report.append(f"  {bottleneck['operation']} - {bottleneck['type']}")
            text_report.append(f"    Max: {bottleneck['max_ms']:.2f}ms, Avg: {bottleneck['avg_ms']:.2f}ms, Ratio: {bottleneck['ratio']:.2f}")
        text_report.append("")
    
    # Recommendations
    text_report.append("Recommendations:")
    for rec in report["recommendations_before"]:
        text_report.append(f"  Issue: {rec['issue']}")
        text_report.append(f"  Recommendation: {rec['recommendation']}")
        text_report.append("")
    
    text_report.append("=" * 80)
    
    return "\n".join(text_report)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Apply database performance optimizations using direct PostgreSQL connection")
    parser.add_argument("--schema", default="public",
                        help="Database schema to modify (default: public)")
    parser.add_argument("--test-only", action="store_true",
                        help="Only test the performance without applying indexes")
    return parser.parse_args()

def main():
    """Main function."""
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Apply performance indexes
        results = apply_performance_indexes(
            schema=args.schema,
            test_only=args.test_only
        )
        
        # Print summary
        print("\n" + "=" * 80)
        print("CryoProtect v2 Database Performance Optimization")
        print("=" * 80)
        
        if results["errors"]:
            print("Errors occurred during optimization:")
            for error in results["errors"]:
                print(f"  - {error}")
        
        print(f"Created {len(results['created_indexes'])} indexes")
        print(f"Failed to create {len(results['failed_indexes'])} indexes")
        
        if results["before_metrics"] and results["after_metrics"]:
            print("\nPerformance Improvements:")
            improvements = calculate_improvements(results)
            for op_name, improvement in improvements.items():
                print(f"  {op_name}:")
                print(f"    Avg: {improvement['avg_ms_before']:.2f}ms -> {improvement['avg_ms_after']:.2f}ms ({improvement['avg_improvement_pct']:.2f}% improvement)")
        
        print("\nDetailed reports saved to:")
        print("  - performance_optimization_report.json")
        print("  - performance_optimization_report.txt")
        print("=" * 80)
        
        return 0
    
    except Exception as e:
        logger.error(f"Error in main function: {str(e)}", exc_info=True)
        return 1
    finally:
        # Close database connections
        sql_executor.close_connections()

if __name__ == "__main__":
    sys.exit(main())