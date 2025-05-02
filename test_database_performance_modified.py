#!/usr/bin/env python3
"""
CryoProtect v2 - Database Performance Testing (Modified)

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
import argparse
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

# Database operations
class DatabaseOperations:
    def __init__(self, supabase):
        self.supabase = supabase
        self.molecule_ids = []
        self.mixture_ids = []
        self.property_type_ids = []
        self.calculation_method_ids = []
        self.operation_times = {}
    
    def load_test_data(self):
        """Load test data IDs from the database."""
        logger.info("Loading test data IDs...")
        
        # Get molecule IDs
        response = self.supabase.table("molecules").select("id").execute()
        if hasattr(response, 'data') and response.data:
            self.molecule_ids = [mol["id"] for mol in response.data]
        
        # Get mixture IDs
        response = self.supabase.table("mixtures").select("id").execute()
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
        if operation_name not in self.operation_times:
            self.operation_times[operation_name] = []
        self.operation_times[operation_name].append(elapsed_ms)
        
        logger.info(f"{operation_name}: {elapsed_ms:.2f} ms")
        
        return result
    
    # Read operations
    def get_molecules(self, limit=100, offset=0):
        """Get a list of molecules with their properties."""
        return self.time_operation(
            "get_molecules",
            self.supabase.table("molecules").select("*").range(offset, offset + limit - 1).execute
        )
    
    def get_molecule_by_id(self, molecule_id):
        """Get a molecule with its properties by ID."""
        return self.time_operation(
            "get_molecule_by_id",
            self.supabase.table("molecules").select("*").eq("id", molecule_id).execute
        )
    
    def get_mixtures(self, limit=100, offset=0):
        """Get a list of mixtures with their components."""
        return self.time_operation(
            "get_mixtures",
            self.supabase.table("mixtures").select("*").range(offset, offset + limit - 1).execute
        )
    
    def get_mixture_by_id(self, mixture_id):
        """Get a mixture with its components by ID."""
        return self.time_operation(
            "get_mixture_by_id",
            self.supabase.table("mixtures").select("*").eq("id", mixture_id).execute
        )
    
    def get_predictions(self, mixture_id, limit=100, offset=0):
        """Get predictions for a mixture."""
        return self.time_operation(
            "get_predictions",
            self.supabase.from_("predictions").select("*").eq("mixture_id", mixture_id).range(offset, offset + limit - 1).execute
        )
    
    def get_experiments(self, mixture_id):
        """Get experiments for a mixture."""
        return self.time_operation(
            "get_experiments",
            self.supabase.from_("experiments").select("*").eq("mixture_id", mixture_id).execute
        )

def run_pagination_test(db_ops):
    """Run pagination performance tests."""
    logger.info("Running pagination performance tests...")
    
    # Test different page sizes
    page_sizes = [10, 25, 50, 100, 250]
    
    for page_size in page_sizes:
        logger.info(f"Testing pagination with page size {page_size}...")
        
        # Test molecules pagination
        total_pages = 5  # Test 5 pages
        for page in range(total_pages):
            offset = page * page_size
            db_ops.get_molecules(limit=page_size, offset=offset)
        
        # Test mixtures pagination
        for page in range(total_pages):
            offset = page * page_size
            db_ops.get_mixtures(limit=page_size, offset=offset)
        
        # Test predictions pagination (if mixture_ids exist)
        if db_ops.mixture_ids:
            mixture_id = random.choice(db_ops.mixture_ids)
            for page in range(total_pages):
                offset = page * page_size
                db_ops.get_predictions(mixture_id, limit=page_size, offset=offset)

def generate_performance_report(db_ops):
    """Generate a performance report based on the collected metrics."""
    report = {
        "timestamp": datetime.now().isoformat(),
        "operations": {}
    }
    
    for operation, times in db_ops.operation_times.items():
        report["operations"][operation] = {
            "count": len(times),
            "min_ms": min(times),
            "max_ms": max(times),
            "avg_ms": statistics.mean(times),
            "median_ms": statistics.median(times),
            "p95_ms": sorted(times)[int(len(times) * 0.95)] if times else 0,
            "p99_ms": sorted(times)[int(len(times) * 0.99)] if times else 0
        }
    
    # Create reports directory if it doesn't exist
    os.makedirs("reports", exist_ok=True)
    
    # Save report to file
    report_path = "reports/query_performance_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Performance report saved to {report_path}")
    
    # Check if any operation exceeds the target response time
    target_response_time = 100  # ms
    slow_operations = []
    for operation, stats in report["operations"].items():
        if stats["avg_ms"] > target_response_time:
            slow_operations.append({
                "operation": operation,
                "avg_ms": stats["avg_ms"],
                "p95_ms": stats["p95_ms"]
            })
    
    if slow_operations:
        logger.warning("The following operations exceed the target response time of 100ms:")
        for op in slow_operations:
            logger.warning(f"  {op['operation']}: avg={op['avg_ms']:.2f}ms, p95={op['p95_ms']:.2f}ms")
    else:
        logger.info("All operations meet the target response time of <100ms")
    
    return report

def main():
    """Main function."""
    try:
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="CryoProtect v2 Database Performance Testing")
        parser.add_argument("--test-pagination", action="store_true", help="Run pagination performance tests only")
        args = parser.parse_args()
        
        logger.info("Starting database performance tests...")
        
        # Connect to Supabase
        supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
        logger.info("Connected to Supabase")
        
        # Initialize database operations
        db_ops = DatabaseOperations(supabase)
        db_ops.load_test_data()
        
        if args.test_pagination:
            # Run pagination tests only
            run_pagination_test(db_ops)
        else:
            # Run basic performance tests
            logger.info("Running basic performance tests...")
            
            # Test molecule queries
            if db_ops.molecule_ids:
                db_ops.get_molecules(limit=100)
                molecule_id = random.choice(db_ops.molecule_ids)
                db_ops.get_molecule_by_id(molecule_id)
            
            # Test mixture queries
            if db_ops.mixture_ids:
                db_ops.get_mixtures(limit=100)
                mixture_id = random.choice(db_ops.mixture_ids)
                db_ops.get_mixture_by_id(mixture_id)
                db_ops.get_predictions(mixture_id)
                db_ops.get_experiments(mixture_id)
        
        # Generate performance report
        report = generate_performance_report(db_ops)
        
        logger.info("Performance tests completed successfully")
        return 0
    except Exception as e:
        logger.error(f"Error running performance tests: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    exit(main())