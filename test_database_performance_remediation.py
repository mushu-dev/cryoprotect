#!/usr/bin/env python3
"""
CryoProtect v2 - Database Performance Validation for Remediation

This script validates the performance of the remediated database, focusing on:
1. Query performance on tables with RLS enabled
2. Query performance on tables with foreign key relationships
3. Query performance on junction tables
"""

import os
import time
import json
import uuid
import random
import statistics
import concurrent.futures
import psutil
from datetime import datetime, timedelta
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("database_performance_validation.log"),
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

# Test configuration
TEST_CONFIG = {
    "concurrent_users": [1, 5, 10],  # Number of concurrent users to simulate
    "iterations_per_user": 5,        # Number of iterations per user
    "rls_operations": [
        "get_user_mixtures",
        "get_user_experiments",
        "get_user_predictions",
        "create_user_mixture"
    ],
    "foreign_key_operations": [
        "get_molecules_with_properties",
        "get_mixtures_with_components",
        "get_predictions_with_methods",
        "get_experiments_with_properties"
    ],
    "junction_table_operations": [
        "get_molecule_proteins",
        "get_molecule_experiments",
        "get_mixture_components"
    ]
}

# Performance metrics
class PerformanceMetrics:
    def __init__(self):
        self.operation_times = {}
        self.cpu_usage = []
        self.memory_usage = []
        self.start_time = None
        self.end_time = None
    
    def start_test(self):
        self.start_time = datetime.now()
    
    def end_test(self):
        self.end_time = datetime.now()
    
    def record_operation_time(self, operation, time_ms):
        if operation not in self.operation_times:
            self.operation_times[operation] = []
        self.operation_times[operation].append(time_ms)
    
    def record_resource_usage(self):
        self.cpu_usage.append(psutil.cpu_percent())
        self.memory_usage.append(psutil.virtual_memory().percent)
    
    def get_summary(self):
        summary = {
            "test_duration": (self.end_time - self.start_time).total_seconds(),
            "operations": {},
            "resource_usage": {
                "cpu": {
                    "min": min(self.cpu_usage) if self.cpu_usage else 0,
                    "max": max(self.cpu_usage) if self.cpu_usage else 0,
                    "avg": statistics.mean(self.cpu_usage) if self.cpu_usage else 0
                },
                "memory": {
                    "min": min(self.memory_usage) if self.memory_usage else 0,
                    "max": max(self.memory_usage) if self.memory_usage else 0,
                    "avg": statistics.mean(self.memory_usage) if self.memory_usage else 0
                }
            }
        }
        
        for operation, times in self.operation_times.items():
            if not times:
                continue
                
            summary["operations"][operation] = {
                "count": len(times),
                "min_ms": min(times),
                "max_ms": max(times),
                "avg_ms": statistics.mean(times),
                "median_ms": statistics.median(times),
                "p95_ms": sorted(times)[int(len(times) * 0.95)] if len(times) >= 20 else max(times),
                "p99_ms": sorted(times)[int(len(times) * 0.99)] if len(times) >= 100 else max(times)
            }
        
        return summary

# Database operations
class DatabaseOperations:
    def __init__(self, supabase):
        self.supabase = supabase
        self.metrics = PerformanceMetrics()
        self.test_data = {
            "molecule_ids": [],
            "mixture_ids": [],
            "property_type_ids": [],
            "calculation_method_ids": [],
            "protein_ids": [],
            "experiment_ids": [],
            "user_id": None
        }
    
    def load_test_data(self):
        """Load test data IDs from the database."""
        logger.info("Loading test data IDs...")
        
        # Get molecule IDs
        response = self.supabase.table("molecules").select("id").limit(50).execute()
        if hasattr(response, 'data') and response.data:
            self.test_data["molecule_ids"] = [mol["id"] for mol in response.data]
        
        # Get mixture IDs
        response = self.supabase.table("mixtures").select("id").limit(50).execute()
        if hasattr(response, 'data') and response.data:
            self.test_data["mixture_ids"] = [mix["id"] for mix in response.data]
        
        # Get property type IDs
        response = self.supabase.table("property_types").select("id").limit(50).execute()
        if hasattr(response, 'data') and response.data:
            self.test_data["property_type_ids"] = [prop["id"] for prop in response.data]
        
        # Get calculation method IDs
        response = self.supabase.table("calculation_methods").select("id").limit(50).execute()
        if hasattr(response, 'data') and response.data:
            self.test_data["calculation_method_ids"] = [method["id"] for method in response.data]
        
        # Get protein IDs if the table exists
        try:
            response = self.supabase.table("proteins").select("id").limit(50).execute()
            if hasattr(response, 'data') and response.data:
                self.test_data["protein_ids"] = [protein["id"] for protein in response.data]
        except Exception as e:
            logger.warning(f"Could not load protein IDs: {str(e)}")
        
        # Get experiment IDs
        response = self.supabase.table("experiments").select("id").limit(50).execute()
        if hasattr(response, 'data') and response.data:
            self.test_data["experiment_ids"] = [exp["id"] for exp in response.data]
        
        # Create a test user ID
        self.test_data["user_id"] = str(uuid.uuid4())
        
        logger.info(f"Loaded {len(self.test_data['molecule_ids'])} molecules, {len(self.test_data['mixture_ids'])} mixtures, "
                   f"{len(self.test_data['property_type_ids'])} property types, {len(self.test_data['calculation_method_ids'])} calculation methods, "
                   f"{len(self.test_data['protein_ids'])} proteins, and {len(self.test_data['experiment_ids'])} experiments")
    
    def time_operation(self, operation_name, func, *args, **kwargs):
        """Time an operation and record the metrics."""
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        
        # Record time in milliseconds
        elapsed_ms = (end_time - start_time) * 1000
        self.metrics.record_operation_time(operation_name, elapsed_ms)
        
        # Record resource usage
        self.metrics.record_resource_usage()
        
        return result
    
    # RLS Operations
    def get_user_mixtures(self, user_id=None):
        """Get mixtures created by a specific user (tests RLS)."""
        user_id = user_id or self.test_data["user_id"]
        
        return self.time_operation(
            "get_user_mixtures",
            self.supabase.table("mixtures").select("*").eq("created_by", user_id).execute
        )
    
    def get_user_experiments(self, user_id=None):
        """Get experiments created by a specific user (tests RLS)."""
        user_id = user_id or self.test_data["user_id"]
        
        return self.time_operation(
            "get_user_experiments",
            self.supabase.table("experiments").select("*").eq("created_by", user_id).execute
        )
    
    def get_user_predictions(self, user_id=None):
        """Get predictions created by a specific user (tests RLS)."""
        user_id = user_id or self.test_data["user_id"]
        
        return self.time_operation(
            "get_user_predictions",
            self.supabase.table("predictions").select("*").eq("created_by", user_id).execute
        )
    
    def create_user_mixture(self, user_id=None):
        """Create a mixture for a specific user (tests RLS)."""
        user_id = user_id or self.test_data["user_id"]
        
        if not self.test_data["molecule_ids"]:
            logger.warning("No molecule IDs available for create_user_mixture")
            return None
        
        # Create a unique name for the test mixture
        name = f"Test Mixture {int(time.time())}"
        description = f"Test mixture created for performance testing"
        
        # Create mixture
        mixture_data = {
            "name": name,
            "description": description,
            "created_by": user_id
        }
        
        def create_mixture_with_components():
            # Insert mixture
            response = self.supabase.table("mixtures").insert(mixture_data).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error creating mixture: {response.error}")
                return None
            
            mixture_id = response.data[0]["id"]
            
            # Insert components
            component_inserts = []
            for i in range(min(3, len(self.test_data["molecule_ids"]))):
                molecule_id = self.test_data["molecule_ids"][i]
                component_inserts.append({
                    "mixture_id": mixture_id,
                    "molecule_id": molecule_id,
                    "amount": random.uniform(1.0, 50.0),
                    "amount_unit": random.choice(["%w/v", "%v/v", "M"]),
                    "role": "cryoprotectant",
                    "created_by": user_id
                })
            
            response = self.supabase.table("mixture_components").insert(component_inserts).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error adding components: {response.error}")
                # Rollback mixture creation
                self.supabase.table("mixtures").delete().eq("id", mixture_id).execute()
                return None
            
            return mixture_id
        
        return self.time_operation("create_user_mixture", create_mixture_with_components)
    
    # Foreign Key Operations
    def get_molecules_with_properties(self, limit=50):
        """Get molecules with their properties (tests foreign key relationships)."""
        return self.time_operation(
            "get_molecules_with_properties",
            self.supabase.from_("molecules").select("""
                *,
                molecular_properties(*)
            """).limit(limit).execute
        )
    
    def get_mixtures_with_components(self, limit=50):
        """Get mixtures with their components (tests foreign key relationships)."""
        return self.time_operation(
            "get_mixtures_with_components",
            self.supabase.from_("mixtures").select("""
                *,
                mixture_components(*)
            """).limit(limit).execute
        )
    
    def get_predictions_with_methods(self, limit=50):
        """Get predictions with their calculation methods (tests foreign key relationships)."""
        return self.time_operation(
            "get_predictions_with_methods",
            self.supabase.from_("predictions").select("""
                *,
                calculation_methods(*),
                property_types(*)
            """).limit(limit).execute
        )
    
    def get_experiments_with_properties(self, limit=50):
        """Get experiments with their properties (tests foreign key relationships)."""
        return self.time_operation(
            "get_experiments_with_properties",
            self.supabase.from_("experiments").select("""
                *,
                property_types(*)
            """).limit(limit).execute
        )
    
    # Junction Table Operations
    def get_molecule_proteins(self, limit=50):
        """Get molecule-protein relationships (tests junction tables)."""
        try:
            return self.time_operation(
                "get_molecule_proteins",
                self.supabase.from_("molecule_proteins").select("""
                    *,
                    molecules(*),
                    proteins(*)
                """).limit(limit).execute
            )
        except Exception as e:
            logger.warning(f"Error in get_molecule_proteins: {str(e)}")
            return None
    
    def get_molecule_experiments(self, limit=50):
        """Get molecule-experiment relationships (tests junction tables)."""
        try:
            return self.time_operation(
                "get_molecule_experiments",
                self.supabase.from_("molecule_experiments").select("""
                    *,
                    molecules(*),
                    experiments(*)
                """).limit(limit).execute
            )
        except Exception as e:
            logger.warning(f"Error in get_molecule_experiments: {str(e)}")
            return None
    
    def get_mixture_components(self, limit=50):
        """Get mixture components (tests junction tables)."""
        return self.time_operation(
            "get_mixture_components",
            self.supabase.from_("mixture_components").select("""
                *,
                mixtures(*),
                molecules(*)
            """).limit(limit).execute
        )

# Test execution
def run_rls_test(db_ops, operation, concurrent_users, iterations_per_user):
    """Run a test on RLS operations with the specified number of concurrent users."""
    logger.info(f"Running RLS test for {operation} with {concurrent_users} concurrent users...")
    
    def user_task():
        for _ in range(iterations_per_user):
            if operation == "get_user_mixtures":
                db_ops.get_user_mixtures()
            elif operation == "get_user_experiments":
                db_ops.get_user_experiments()
            elif operation == "get_user_predictions":
                db_ops.get_user_predictions()
            elif operation == "create_user_mixture":
                db_ops.create_user_mixture()
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrent_users) as executor:
        futures = [executor.submit(user_task) for _ in range(concurrent_users)]
        concurrent.futures.wait(futures)

def run_foreign_key_test(db_ops, operation, concurrent_users, iterations_per_user):
    """Run a test on foreign key operations with the specified number of concurrent users."""
    logger.info(f"Running foreign key test for {operation} with {concurrent_users} concurrent users...")
    
    def user_task():
        for _ in range(iterations_per_user):
            if operation == "get_molecules_with_properties":
                db_ops.get_molecules_with_properties()
            elif operation == "get_mixtures_with_components":
                db_ops.get_mixtures_with_components()
            elif operation == "get_predictions_with_methods":
                db_ops.get_predictions_with_methods()
            elif operation == "get_experiments_with_properties":
                db_ops.get_experiments_with_properties()
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrent_users) as executor:
        futures = [executor.submit(user_task) for _ in range(concurrent_users)]
        concurrent.futures.wait(futures)

def run_junction_table_test(db_ops, operation, concurrent_users, iterations_per_user):
    """Run a test on junction table operations with the specified number of concurrent users."""
    logger.info(f"Running junction table test for {operation} with {concurrent_users} concurrent users...")
    
    def user_task():
        for _ in range(iterations_per_user):
            if operation == "get_molecule_proteins":
                db_ops.get_molecule_proteins()
            elif operation == "get_molecule_experiments":
                db_ops.get_molecule_experiments()
            elif operation == "get_mixture_components":
                db_ops.get_mixture_components()
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrent_users) as executor:
        futures = [executor.submit(user_task) for _ in range(concurrent_users)]
        concurrent.futures.wait(futures)

def run_performance_tests():
    """Run all performance tests."""
    logger.info("Starting database performance validation tests...")
    
    # Connect to Supabase
    supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
    logger.info("Connected to Supabase")
    
    # Initialize database operations
    db_ops = DatabaseOperations(supabase)
    db_ops.load_test_data()
    
    # Start metrics collection
    db_ops.metrics.start_test()
    
    # Run tests with increasing concurrent users
    for concurrent_users in TEST_CONFIG["concurrent_users"]:
        logger.info(f"Running tests with {concurrent_users} concurrent users...")
        
        # Run RLS tests
        for operation in TEST_CONFIG["rls_operations"]:
            run_rls_test(db_ops, operation, concurrent_users, TEST_CONFIG["iterations_per_user"])
        
        # Run foreign key tests
        for operation in TEST_CONFIG["foreign_key_operations"]:
            run_foreign_key_test(db_ops, operation, concurrent_users, TEST_CONFIG["iterations_per_user"])
        
        # Run junction table tests
        for operation in TEST_CONFIG["junction_table_operations"]:
            run_junction_table_test(db_ops, operation, concurrent_users, TEST_CONFIG["iterations_per_user"])
    
    # End metrics collection
    db_ops.metrics.end_test()
    
    # Generate report
    metrics_summary = db_ops.metrics.get_summary()
    
    # Save report to file
    report = {
        "timestamp": datetime.now().isoformat(),
        "test_config": TEST_CONFIG,
        "metrics": metrics_summary,
        "analysis": analyze_metrics(metrics_summary),
        "recommendations": generate_recommendations(metrics_summary)
    }
    
    with open("database_performance_validation_report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info("Performance validation tests completed. Report saved to database_performance_validation_report.json")
    
    return report

def analyze_metrics(metrics):
    """Analyze the performance metrics and identify bottlenecks."""
    analysis = {
        "bottlenecks": [],
        "slow_operations": [],
        "resource_issues": [],
        "rls_performance": {
            "operations": [],
            "avg_response_time": 0,
            "status": "GOOD"
        },
        "foreign_key_performance": {
            "operations": [],
            "avg_response_time": 0,
            "status": "GOOD"
        },
        "junction_table_performance": {
            "operations": [],
            "avg_response_time": 0,
            "status": "GOOD"
        }
    }
    
    # Check for slow operations (p95 > 500ms)
    for operation, stats in metrics["operations"].items():
        if stats["p95_ms"] > 500:
            analysis["slow_operations"].append({
                "operation": operation,
                "p95_ms": stats["p95_ms"],
                "avg_ms": stats["avg_ms"]
            })
    
    # Check for high CPU usage (avg > 80%)
    if metrics["resource_usage"]["cpu"]["avg"] > 80:
        analysis["resource_issues"].append({
            "type": "high_cpu_usage",
            "value": metrics["resource_usage"]["cpu"]["avg"]
        })
    
    # Check for high memory usage (avg > 80%)
    if metrics["resource_usage"]["memory"]["avg"] > 80:
        analysis["resource_issues"].append({
            "type": "high_memory_usage",
            "value": metrics["resource_usage"]["memory"]["avg"]
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
    
    # Analyze RLS performance
    rls_times = []
    for operation, stats in metrics["operations"].items():
        if operation.startswith("get_user_") or operation == "create_user_mixture":
            analysis["rls_performance"]["operations"].append({
                "operation": operation,
                "avg_ms": stats["avg_ms"],
                "p95_ms": stats["p95_ms"]
            })
            rls_times.append(stats["avg_ms"])
    
    if rls_times:
        avg_rls_time = statistics.mean(rls_times)
        analysis["rls_performance"]["avg_response_time"] = avg_rls_time
        if avg_rls_time > 500:
            analysis["rls_performance"]["status"] = "POOR"
        elif avg_rls_time > 200:
            analysis["rls_performance"]["status"] = "FAIR"
    
    # Analyze foreign key performance
    fk_times = []
    for operation, stats in metrics["operations"].items():
        if operation.startswith("get_") and "_with_" in operation:
            analysis["foreign_key_performance"]["operations"].append({
                "operation": operation,
                "avg_ms": stats["avg_ms"],
                "p95_ms": stats["p95_ms"]
            })
            fk_times.append(stats["avg_ms"])
    
    if fk_times:
        avg_fk_time = statistics.mean(fk_times)
        analysis["foreign_key_performance"]["avg_response_time"] = avg_fk_time
        if avg_fk_time > 500:
            analysis["foreign_key_performance"]["status"] = "POOR"
        elif avg_fk_time > 200:
            analysis["foreign_key_performance"]["status"] = "FAIR"
    
    # Analyze junction table performance
    jt_times = []
    for operation, stats in metrics["operations"].items():
        if operation in ["get_molecule_proteins", "get_molecule_experiments", "get_mixture_components"]:
            analysis["junction_table_performance"]["operations"].append({
                "operation": operation,
                "avg_ms": stats["avg_ms"],
                "p95_ms": stats["p95_ms"]
            })
            jt_times.append(stats["avg_ms"])
    
    if jt_times:
        avg_jt_time = statistics.mean(jt_times)
        analysis["junction_table_performance"]["avg_response_time"] = avg_jt_time
        if avg_jt_time > 500:
            analysis["junction_table_performance"]["status"] = "POOR"
        elif avg_jt_time > 200:
            analysis["junction_table_performance"]["status"] = "FAIR"
    
    return analysis

def generate_recommendations(metrics):
    """Generate recommendations based on the performance metrics."""
    recommendations = []
    
    # Check for slow RLS operations
    for operation, stats in metrics["operations"].items():
        if operation.startswith("get_user_") and stats["p95_ms"] > 500:
            recommendations.append({
                "operation": operation,
                "issue": "Slow RLS performance",
                "recommendation": f"Consider adding an index on the created_by column in the table used by {operation}."
            })
    
    # Check for slow foreign key operations
    for operation, stats in metrics["operations"].items():
        if operation.startswith("get_") and "_with_" in operation and stats["p95_ms"] > 500:
            table_name = operation.replace("get_", "").replace("_with_", "_")
            recommendations.append({
                "operation": operation,
                "issue": "Slow foreign key join performance",
                "recommendation": f"Consider adding indexes on the foreign key columns in the {table_name} table."
            })
    
    # Check for slow junction table operations
    for operation, stats in metrics["operations"].items():
        if operation in ["get_molecule_proteins", "get_molecule_experiments", "get_mixture_components"] and stats["p95_ms"] > 500:
            table_name = operation.replace("get_", "")
            recommendations.append({
                "operation": operation,
                "issue": "Slow junction table performance",
                "recommendation": f"Consider adding composite indexes on the foreign key columns in the {table_name} table."
            })
    
    # Check for high resource usage
    if any(issue["type"] == "high_cpu_usage" for issue in analyze_metrics(metrics)["resource_issues"]):
        recommendations.append({
            "issue": "High CPU usage",
            "recommendation": "Consider optimizing complex queries or adding more indexes to reduce CPU load."
        })
    
    if any(issue["type"] == "high_memory_usage" for issue in analyze_metrics(metrics)["resource_issues"]):
        recommendations.append({
            "issue": "High memory usage",
            "recommendation": "Consider optimizing queries that return large result sets or implementing pagination."
        })
    
    # General recommendations
    recommendations.append({
        "issue": "General optimization",
        "recommendation": "Consider implementing connection pooling to improve connection reuse and reduce overhead."
    })
    
    recommendations.append({
        "issue": "General optimization",
        "recommendation": "Consider implementing query caching for frequently accessed data."
    })
    
    return recommendations

def generate_text_report(report):
    """Generate a text report from the JSON report."""
    text_report = []
    
    text_report.append("=" * 80)
    text_report.append("CryoProtect v2 Database Performance Validation Report")
    text_report.append("=" * 80)
    text_report.append(f"Timestamp: {report['timestamp']}")
    text_report.append(f"Test Duration: {report['metrics']['test_duration']:.2f} seconds")
    text_report.append("")
    
    # Overall status
    overall_status = "SUCCESS"
    if (report['analysis']['rls_performance']['status'] == "POOR" or
        report['analysis']['foreign_key_performance']['status'] == "POOR" or
        report['analysis']['junction_table_performance']['status'] == "POOR"):
        overall_status = "ERROR"
    elif (report['analysis']['rls_performance']['status'] == "FAIR" or
          report['analysis']['foreign_key_performance']['status'] == "FAIR" or
          report['analysis']['junction_table_performance']['status'] == "FAIR"):
        overall_status = "COMPLETED_WITH_WARNINGS"
    
    text_report.append(f"Status: {overall_status}")
    text_report.append("")
    
    # Resource usage
    text_report.append("Resource Usage:")
    text_report.append(f"  CPU: Min={report['metrics']['resource_usage']['cpu']['min']:.2f}%, "
                      f"Max={report['metrics']['resource_usage']['cpu']['max']:.2f}%, "
                      f"Avg={report['metrics']['resource_usage']['cpu']['avg']:.2f}%")
    text_report.append(f"  Memory: Min={report['metrics']['resource_usage']['memory']['min']:.2f}%, "
                      f"Max={report['metrics']['resource_usage']['memory']['max']:.2f}%, "
                      f"Avg={report['metrics']['resource_usage']['memory']['avg']:.2f}%")
    text_report.append("")
    
    # RLS Performance
    text_report.append("RLS Performance:")
    text_report.append(f"  Status: {report['analysis']['rls_performance']['status']}")
    text_report.append(f"  Average Response Time: {report['analysis']['rls_performance']['avg_response_time']:.2f} ms")
    text_report.append("  Operations:")
    for op in report['analysis']['rls_performance']['operations']:
        text_report.append(f"    - {op['operation']}: Avg={op['avg_ms']:.2f} ms, P95={op['p95_ms']:.2f} ms")
    text_report.append("")
    
    # Foreign Key Performance
    text_report.append("Foreign Key Performance:")
    text_report.append(f"  Status: {report['analysis']['foreign_key_performance']['status']}")
    text_report.append(f"  Average Response Time: {report['analysis']['foreign_key_performance']['avg_response_time']:.2f} ms")
    text_report.append("  Operations:")
    for op in report['analysis']['foreign_key_performance']['operations']:
        text_report.append(f"    - {op['operation']}: Avg={op['avg_ms']:.2f} ms, P95={op['p95_ms']:.2f} ms")
    text_report.append("")
    
    # Junction Table Performance
    text_report.append("Junction Table Performance:")
    text_report.append(f"  Status: {report['analysis']['junction_table_performance']['status']}")
    text_report.append(f"  Average Response Time: {report['analysis']['junction_table_performance']['avg_response_time']:.2f} ms")
    text_report.append("  Operations:")
    for op in report['analysis']['junction_table_performance']['operations']:
        text_report.append(f"    - {op['operation']}: Avg={op['avg_ms']:.2f} ms, P95={op['p95_ms']:.2f} ms")
    text_report.append("")
    
    # Bottlenecks
    if report['analysis']['bottlenecks']:
        text_report.append("Bottlenecks:")
        for bottleneck in report['analysis']['bottlenecks']:
            text_report.append(f"  {bottleneck['operation']} - {bottleneck['type']}")
            text_report.append(f"    Max: {bottleneck['max_ms']:.2f} ms, Avg: {bottleneck['avg_ms']:.2f} ms, Ratio: {bottleneck['ratio']:.2f}")
        text_report.append("")
    
    # Slow Operations
    if report['analysis']['slow_operations']:
        text_report.append("Slow Operations:")
        for op in report['analysis']['slow_operations']:
            text_report.append(f"  {op['operation']} - P95: {op['p95_ms']:.2f} ms, Avg: {op['avg_ms']:.2f} ms")
        text_report.append("")
    
    # Resource Issues
    if report['analysis']['resource_issues']:
        text_report.append("Resource Issues:")
        for issue in report['analysis']['resource_issues']:
            text_report.append(f"  {issue['type']} - {issue['value']:.2f}%")
        text_report.append("")
    
    # Recommendations
    text_report.append("Recommendations:")
    for rec in report['recommendations']:
        if "operation" in rec:
            text_report.append(f"  {rec['operation']} - {rec['issue']}")
        else:
            text_report.append(f"  {rec['issue']}")
        text_report.append(f"    {rec['recommendation']}")
    text_report.append("")
    
    text_report.append("=" * 80)
    
    return "\n".join(text_report)

def main():
    """Main function."""
    try:
        # Run performance tests
        report = run_performance_tests()
        
        # Generate text report
        text_report = generate_text_report(report)
        
        # Save text report to file
        with open("database_performance_validation_report.txt", "w") as f:
            f.write(text_report)
        
        print("\n" + "=" * 80)
        print("CryoProtect v2 Database Performance Validation Completed")
        print("=" * 80)
        print(f"Report saved to database_performance_validation_report.txt and database_performance_validation_report.json")
        print("=" * 80)
        print(text_report)
        
        return 0
    except Exception as e:
        logger.error(f"Error running performance validation: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    exit(main())
