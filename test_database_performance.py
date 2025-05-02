#!/usr/bin/env python3
"""
CryoProtect v2 - Database Performance Testing

This script tests the database performance under expected production load.
It measures query performance, response times, and resource utilization.
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
        logging.FileHandler("database_performance_test.log"),
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
    "concurrent_users": [1, 5, 10, 20],  # Number of concurrent users to simulate
    "iterations_per_user": 5,            # Number of iterations per user
    "read_operations": [
        "get_molecules",
        "get_molecule_by_id",
        "get_mixtures",
        "get_mixture_by_id",
        "get_predictions",
        "get_experiments",
        "compare_predictions_experiments"
    ],
    "write_operations": [
        "create_mixture",
        "add_prediction",
        "record_experiment"
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
                    "min": min(self.cpu_usage),
                    "max": max(self.cpu_usage),
                    "avg": statistics.mean(self.cpu_usage)
                },
                "memory": {
                    "min": min(self.memory_usage),
                    "max": max(self.memory_usage),
                    "avg": statistics.mean(self.memory_usage)
                }
            }
        }
        
        for operation, times in self.operation_times.items():
            summary["operations"][operation] = {
                "count": len(times),
                "min_ms": min(times),
                "max_ms": max(times),
                "avg_ms": statistics.mean(times),
                "median_ms": statistics.median(times),
                "p95_ms": sorted(times)[int(len(times) * 0.95)] if times else 0,
                "p99_ms": sorted(times)[int(len(times) * 0.99)] if times else 0
            }
        
        return summary

# Database operations
class DatabaseOperations:
    def __init__(self, supabase):
        self.supabase = supabase
        self.metrics = PerformanceMetrics()
        self.molecule_ids = []
        self.mixture_ids = []
        self.property_type_ids = []
        self.calculation_method_ids = []
    
    def load_test_data(self):
        """Load test data IDs from the database."""
        logger.info("Loading test data IDs...")
        
        # Get molecule IDs
        response = self.supabase.table("molecule").select("id").execute()
        if hasattr(response, 'data') and response.data:
            self.molecule_ids = [mol["id"] for mol in response.data]
        
        # Get mixture IDs
        response = self.supabase.table("mixture").select("id").execute()
        if hasattr(response, 'data') and response.data:
            self.mixture_ids = [mix["id"] for mix in response.data]
        
        # Get property type IDs
        response = self.supabase.table("property_types").select("id").execute()
        if hasattr(response, 'data') and response.data:
            self.property_type_ids = [prop["id"] for prop in response.data]
        
        # Get calculation method IDs
        response = self.supabase.table("calculation_methods").select("id").execute()
        if hasattr(response, 'data') and response.data:
            self.calculation_method_ids = [method["id"] for method in response.data]
        
        logger.info(f"Loaded {len(self.molecule_ids)} molecules, {len(self.mixture_ids)} mixtures, "
                   f"{len(self.property_type_ids)} property types, and {len(self.calculation_method_ids)} calculation methods")
    
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
    
    # Read operations
    def get_molecules(self, limit=100, offset=0):
        """Get a list of molecules with their properties."""
        return self.time_operation(
            "get_molecules",
            self.supabase.table("molecule_with_properties").select("*").range(offset, offset + limit - 1).execute
        )
    
    def get_molecule_by_id(self, molecule_id):
        """Get a molecule with its properties by ID."""
        return self.time_operation(
            "get_molecule_by_id",
            self.supabase.table("molecule_with_properties").select("*").eq("id", molecule_id).execute
        )
    
    def get_mixtures(self, limit=100, offset=0):
        """Get a list of mixtures with their components."""
        return self.time_operation(
            "get_mixtures",
            self.supabase.table("mixture_with_components").select("*").range(offset, offset + limit - 1).execute
        )
    
    def get_mixture_by_id(self, mixture_id):
        """Get a mixture with its components by ID."""
        return self.time_operation(
            "get_mixture_by_id",
            self.supabase.table("mixture_with_components").select("*").eq("id", mixture_id).execute
        )
    
    def get_predictions(self, mixture_id, limit=100, offset=0):
        """Get predictions for a mixture."""
        return self.time_operation(
            "get_predictions",
            self.supabase.from_("predictions").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value,
                confidence, created_at,
                property_types(name), calculation_methods(name)
            """).eq("mixture_id", mixture_id).range(offset, offset + limit - 1).execute
        )
    
    def get_experiments(self, mixture_id):
        """Get experiments for a mixture."""
        return self.time_operation(
            "get_experiments",
            self.supabase.from_("experiments").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value, 
                experimental_conditions, date_performed, created_at, 
                property_types(name)
            """).eq("mixture_id", mixture_id).execute
        )
    
    def compare_predictions_experiments(self, mixture_id, property_type_id=None):
        """Compare predictions with experiments for a mixture."""
        params = {"p_mixture_id": mixture_id}
        if property_type_id:
            params["p_property_type_id"] = property_type_id
        
        return self.time_operation(
            "compare_predictions_experiments",
            self.supabase.rpc("compare_predictions_with_experiments", params).execute
        )
    
    # Write operations
    def create_mixture(self, name, description, components):
        """Create a new mixture with components."""
        # Create a unique user ID for testing
        user_id = str(uuid.uuid4())
        
        # Create mixture
        mixture_data = {
            "name": name,
            "description": description,
            "created_by": user_id
        }
        
        def create_mixture_with_components():
            # Insert mixture
            response = self.supabase.table("mixture").insert(mixture_data).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error creating mixture: {response.error}")
                return None
            
            mixture_id = response.data[0]["id"]
            
            # Insert components
            component_inserts = []
            for comp in components:
                component_inserts.append({
                    "mixture_id": mixture_id,
                    "molecule_id": comp["molecule_id"],
                    "amount": comp["amount"],
                    "amount_unit": comp["amount_unit"],
                    "role": "cryoprotectant",
                    "created_by": user_id
                })
            
            response = self.supabase.table("mixture_component").insert(component_inserts).execute()
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error adding components: {response.error}")
                # Rollback mixture creation
                self.supabase.table("mixture").delete().eq("id", mixture_id).execute()
                return None
            
            return mixture_id
        
        return self.time_operation("create_mixture", create_mixture_with_components)
    
    def add_prediction(self, mixture_id, property_type_id, calculation_method_id, value, confidence=0.9):
        """Add a prediction for a mixture."""
        # Create a unique user ID for testing
        user_id = str(uuid.uuid4())
        
        # Get property type data type
        response = self.supabase.table("property_types").select("data_type").eq("id", property_type_id).execute()
        if not response.data:
            logger.error(f"Property type {property_type_id} not found")
            return None
        
        data_type = response.data[0]["data_type"]
        
        # Prepare prediction data
        prediction_data = {
            "mixture_id": mixture_id,
            "property_type_id": property_type_id,
            "calculation_method_id": calculation_method_id,
            "confidence": confidence,
            "created_by": user_id
        }
        
        # Set the appropriate value field based on data type
        if data_type == "numeric":
            prediction_data["numeric_value"] = float(value)
        elif data_type == "text":
            prediction_data["text_value"] = str(value)
        elif data_type == "boolean":
            prediction_data["boolean_value"] = bool(value)
        
        return self.time_operation(
            "add_prediction",
            self.supabase.table("predictions").upsert(
                prediction_data,
                on_conflict="mixture_id,property_type_id,calculation_method_id"
            ).execute
        )
    
    def record_experiment(self, mixture_id, property_type_id, value, conditions="Standard conditions"):
        """Record an experiment for a mixture."""
        # Create a unique user ID for testing
        user_id = str(uuid.uuid4())
        
        # Get property type data type
        response = self.supabase.table("property_types").select("data_type").eq("id", property_type_id).execute()
        if not response.data:
            logger.error(f"Property type {property_type_id} not found")
            return None
        
        data_type = response.data[0]["data_type"]
        
        # Prepare experiment data
        experiment_data = {
            "mixture_id": mixture_id,
            "property_type_id": property_type_id,
            "experimental_conditions": conditions,
            "date_performed": datetime.now().isoformat(),
            "created_by": user_id
        }
        
        # Set the appropriate value field based on data type
        if data_type == "numeric":
            experiment_data["numeric_value"] = float(value)
        elif data_type == "text":
            experiment_data["text_value"] = str(value)
        elif data_type == "boolean":
            experiment_data["boolean_value"] = bool(value)
        
        return self.time_operation(
            "record_experiment",
            self.supabase.table("experiments").insert(experiment_data).execute
        )

# Test execution
def run_read_test(db_ops, operation, concurrent_users, iterations_per_user):
    """Run a read test with the specified number of concurrent users."""
    logger.info(f"Running read test for {operation} with {concurrent_users} concurrent users...")
    
    def user_task():
        for _ in range(iterations_per_user):
            if operation == "get_molecules":
                db_ops.get_molecules(limit=random.randint(10, 100), offset=random.randint(0, 5))
            elif operation == "get_molecule_by_id":
                if db_ops.molecule_ids:
                    molecule_id = random.choice(db_ops.molecule_ids)
                    db_ops.get_molecule_by_id(molecule_id)
            elif operation == "get_mixtures":
                db_ops.get_mixtures(limit=random.randint(10, 100), offset=random.randint(0, 5))
            elif operation == "get_mixture_by_id":
                if db_ops.mixture_ids:
                    mixture_id = random.choice(db_ops.mixture_ids)
                    db_ops.get_mixture_by_id(mixture_id)
            elif operation == "get_predictions":
                if db_ops.mixture_ids:
                    mixture_id = random.choice(db_ops.mixture_ids)
                    db_ops.get_predictions(mixture_id, limit=random.randint(10, 100), offset=random.randint(0, 5))
            elif operation == "get_experiments":
                if db_ops.mixture_ids:
                    mixture_id = random.choice(db_ops.mixture_ids)
                    db_ops.get_experiments(mixture_id)
            elif operation == "compare_predictions_experiments":
                if db_ops.mixture_ids and db_ops.property_type_ids:
                    mixture_id = random.choice(db_ops.mixture_ids)
                    property_type_id = random.choice(db_ops.property_type_ids) if random.random() > 0.5 else None
                    db_ops.compare_predictions_experiments(mixture_id, property_type_id)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrent_users) as executor:
        futures = [executor.submit(user_task) for _ in range(concurrent_users)]
        concurrent.futures.wait(futures)

def run_write_test(db_ops, operation, concurrent_users, iterations_per_user):
    """Run a write test with the specified number of concurrent users."""
    logger.info(f"Running write test for {operation} with {concurrent_users} concurrent users...")
    
    def user_task():
        for i in range(iterations_per_user):
            if operation == "create_mixture":
                if db_ops.molecule_ids:
                    # Create a mixture with 2-3 random components
                    name = f"Test Mixture {datetime.now().strftime('%Y%m%d%H%M%S')}-{i}"
                    description = f"Test mixture created for performance testing"
                    components = []
                    for _ in range(random.randint(2, 3)):
                        molecule_id = random.choice(db_ops.molecule_ids)
                        components.append({
                            "molecule_id": molecule_id,
                            "amount": random.uniform(1.0, 50.0),
                            "amount_unit": random.choice(["%w/v", "%v/v", "M"])
                        })
                    mixture_id = db_ops.create_mixture(name, description, components)
                    if mixture_id:
                        db_ops.mixture_ids.append(mixture_id)
            elif operation == "add_prediction":
                if db_ops.mixture_ids and db_ops.property_type_ids and db_ops.calculation_method_ids:
                    mixture_id = random.choice(db_ops.mixture_ids)
                    property_type_id = random.choice(db_ops.property_type_ids)
                    calculation_method_id = random.choice(db_ops.calculation_method_ids)
                    value = random.uniform(-100.0, 100.0)  # Assuming numeric for simplicity
                    confidence = random.uniform(0.5, 1.0)
                    db_ops.add_prediction(mixture_id, property_type_id, calculation_method_id, value, confidence)
            elif operation == "record_experiment":
                if db_ops.mixture_ids and db_ops.property_type_ids:
                    mixture_id = random.choice(db_ops.mixture_ids)
                    property_type_id = random.choice(db_ops.property_type_ids)
                    value = random.uniform(-100.0, 100.0)  # Assuming numeric for simplicity
                    conditions = f"Test conditions at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
                    db_ops.record_experiment(mixture_id, property_type_id, value, conditions)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrent_users) as executor:
        futures = [executor.submit(user_task) for _ in range(concurrent_users)]
        concurrent.futures.wait(futures)

def run_performance_tests():
    """Run all performance tests."""
    logger.info("Starting database performance tests...")
    
    # Connect to Supabase
    supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
    logger.info("Connected to Supabase")
    
    # Initialize database operations
    db_ops = DatabaseOperations(supabase)
    db_ops.load_test_data()
    
    # Start metrics collection
    db_ops.metrics.start_test()
    
    # Run baseline tests (single user)
    logger.info("Running baseline tests (single user)...")
    for operation in TEST_CONFIG["read_operations"]:
        run_read_test(db_ops, operation, 1, 5)
    
    for operation in TEST_CONFIG["write_operations"]:
        run_write_test(db_ops, operation, 1, 2)
    
    # Run load tests with increasing concurrent users
    for concurrent_users in TEST_CONFIG["concurrent_users"][1:]:  # Skip the first one (1 user) as we already did baseline
        logger.info(f"Running load tests with {concurrent_users} concurrent users...")
        
        # Run read tests
        for operation in TEST_CONFIG["read_operations"]:
            run_read_test(db_ops, operation, concurrent_users, TEST_CONFIG["iterations_per_user"])
        
        # Run write tests
        for operation in TEST_CONFIG["write_operations"]:
            run_write_test(db_ops, operation, concurrent_users, TEST_CONFIG["iterations_per_user"] // 2)  # Fewer write operations
    
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
    
    with open("database_performance_report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info("Performance tests completed. Report saved to database_performance_report.json")
    
    return report

def analyze_metrics(metrics):
    """Analyze the performance metrics and identify bottlenecks."""
    analysis = {
        "bottlenecks": [],
        "slow_operations": [],
        "resource_issues": []
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
    
    return analysis

def generate_recommendations(metrics):
    """Generate recommendations based on the performance metrics."""
    recommendations = []
    
    # Check for slow read operations
    for operation, stats in metrics["operations"].items():
        if operation.startswith("get_") and stats["p95_ms"] > 500:
            if "molecule_with_properties" in operation:
                recommendations.append({
                    "operation": operation,
                    "issue": "Slow read performance",
                    "recommendation": "Consider adding indexes on commonly queried fields in the molecule and molecular_property tables."
                })
            elif "mixture_with_components" in operation:
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
    
    # Check for slow write operations
    for operation, stats in metrics["operations"].items():
        if operation.startswith("create_") or operation.startswith("add_") or operation.startswith("record_") and stats["p95_ms"] > 1000:
            recommendations.append({
                "operation": operation,
                "issue": "Slow write performance",
                "recommendation": "Consider optimizing the database schema or using batch operations for better write performance."
            })
    
    # Check for high resource usage
    if metrics["resource_usage"]["cpu"]["avg"] > 80:
        recommendations.append({
            "issue": "High CPU usage",
            "recommendation": "Consider scaling up the database server or optimizing queries to reduce CPU load."
        })
    
    if metrics["resource_usage"]["memory"]["avg"] > 80:
        recommendations.append({
            "issue": "High memory usage",
            "recommendation": "Consider scaling up the database server memory or optimizing queries to reduce memory usage."
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
    text_report.append("CryoProtect v2 Database Performance Test Report")
    text_report.append("=" * 80)
    text_report.append(f"Timestamp: {report['timestamp']}")
    text_report.append(f"Test Duration: {report['metrics']['test_duration']:.2f} seconds")
    text_report.append("")
    
    text_report.append("Resource Usage:")
    text_report.append(f"  CPU: Min={report['metrics']['resource_usage']['cpu']['min']:.2f}%, "
                      f"Max={report['metrics']['resource_usage']['cpu']['max']:.2f}%, "
                      f"Avg={report['metrics']['resource_usage']['cpu']['avg']:.2f}%")
    text_report.append(f"  Memory: Min={report['metrics']['resource_usage']['memory']['min']:.2f}%, "
                      f"Max={report['metrics']['resource_usage']['memory']['max']:.2f}%, "
                      f"Avg={report['metrics']['resource_usage']['memory']['avg']:.2f}%")
    text_report.append("")
    
    text_report.append("Operation Performance:")
    for operation, stats in report['metrics']['operations'].items():
        text_report.append(f"  {operation}:")
        text_report.append(f"    Count: {stats['count']}")
        text_report.append(f"    Min: {stats['min_ms']:.2f} ms")
        text_report.append(f"    Max: {stats['max_ms']:.2f} ms")
        text_report.append(f"    Avg: {stats['avg_ms']:.2f} ms")
        text_report.append(f"    Median: {stats['median_ms']:.2f} ms")
        text_report.append(f"    P95: {stats['p95_ms']:.2f} ms")
        text_report.append(f"    P99: {stats['p99_ms']:.2f} ms")
        text_report.append("")
    
    if report['analysis']['bottlenecks']:
        text_report.append("Bottlenecks:")
        for bottleneck in report['analysis']['bottlenecks']:
            text_report.append(f"  {bottleneck['operation']} - {bottleneck['type']}")
            text_report.append(f"    Max: {bottleneck['max_ms']:.2f} ms, Avg: {bottleneck['avg_ms']:.2f} ms, Ratio: {bottleneck['ratio']:.2f}")
        text_report.append("")
    
    if report['analysis']['slow_operations']:
        text_report.append("Slow Operations:")
        for op in report['analysis']['slow_operations']:
            text_report.append(f"  {op['operation']} - P95: {op['p95_ms']:.2f} ms, Avg: {op['avg_ms']:.2f} ms")
        text_report.append("")
    
    if report['analysis']['resource_issues']:
        text_report.append("Resource Issues:")
        for issue in report['analysis']['resource_issues']:
            text_report.append(f"  {issue['type']} - {issue['value']:.2f}%")
        text_report.append("")
    
    text_report.append("Recommendations:")
    for rec in report['recommendations']:
        text_report.append(f"  Issue: {rec['issue']}")
        text_report.append(f"  Recommendation: {rec['recommendation']}")
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
        with open("database_performance_report.txt", "w") as f:
            f.write(text_report)
        
        print("\n" + "=" * 80)
        print("CryoProtect v2 Database Performance Testing Completed")
        print("=" * 80)
        print(f"Report saved to database_performance_report.txt and database_performance_report.json")
        print("=" * 80)
        
        return 0
    except Exception as e:
        logger.error(f"Error running performance tests: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    exit(main())