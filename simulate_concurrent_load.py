#!/usr/bin/env python3
"""
CryoProtect v2 - Database Concurrent Load Simulation

This script simulates concurrent user load on the database to test how it performs under stress.
It simulates multiple users performing various operations simultaneously.
"""

import os
import time
import json
import uuid
import random
import statistics
import concurrent.futures
import threading
from datetime import datetime, timedelta
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("concurrent_load_test.log"),
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
    "duration_seconds": 60,           # Duration of the test in seconds
    "ramp_up_seconds": 10,            # Time to ramp up to full user load
    "concurrent_users": 20,           # Maximum number of concurrent users
    "user_think_time_ms": [500, 3000],  # Random think time between operations (min, max)
    "operation_weights": {            # Relative frequency of operations
        "browse_molecules": 30,
        "view_molecule_details": 15,
        "browse_mixtures": 20,
        "view_mixture_details": 15,
        "search_molecules": 10,
        "create_mixture": 5,
        "add_prediction": 3,
        "record_experiment": 2
    }
}

# Global test data
test_data = {
    "molecule_ids": [],
    "mixture_ids": [],
    "property_type_ids": [],
    "calculation_method_ids": []
}

# Global metrics
class ConcurrentLoadMetrics:
    def __init__(self):
        self.lock = threading.Lock()
        self.operation_times = {}
        self.operation_counts = {}
        self.errors = {}
        self.start_time = None
        self.end_time = None
    
    def start_test(self):
        self.start_time = datetime.now()
    
    def end_test(self):
        self.end_time = datetime.now()
    
    def record_operation(self, operation, time_ms, success=True):
        with self.lock:
            if operation not in self.operation_times:
                self.operation_times[operation] = []
                self.operation_counts[operation] = {"success": 0, "error": 0}
            
            self.operation_times[operation].append(time_ms)
            if success:
                self.operation_counts[operation]["success"] += 1
            else:
                self.operation_counts[operation]["error"] += 1
    
    def record_error(self, operation, error_message):
        with self.lock:
            if operation not in self.errors:
                self.errors[operation] = []
            self.errors[operation].append(error_message)
    
    def get_summary(self):
        summary = {
            "test_duration": (self.end_time - self.start_time).total_seconds(),
            "operations": {},
            "errors": self.errors
        }
        
        total_operations = 0
        total_time_ms = 0
        
        for operation, times in self.operation_times.items():
            if not times:
                continue
                
            op_summary = {
                "count": len(times),
                "success_count": self.operation_counts[operation]["success"],
                "error_count": self.operation_counts[operation]["error"],
                "min_ms": min(times),
                "max_ms": max(times),
                "avg_ms": statistics.mean(times),
                "median_ms": statistics.median(times),
                "p95_ms": sorted(times)[int(len(times) * 0.95)] if len(times) >= 20 else max(times),
                "throughput_per_second": len(times) / summary["test_duration"]
            }
            
            summary["operations"][operation] = op_summary
            total_operations += len(times)
            total_time_ms += sum(times)
        
        summary["total_operations"] = total_operations
        summary["avg_response_time_ms"] = total_time_ms / total_operations if total_operations > 0 else 0
        summary["throughput_per_second"] = total_operations / summary["test_duration"]
        
        return summary

# Global metrics instance
metrics = ConcurrentLoadMetrics()

# User simulation
class SimulatedUser:
    def __init__(self, user_id):
        self.user_id = user_id
        self.supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
        self.session_id = str(uuid.uuid4())
        logger.debug(f"User {user_id} initialized with session {self.session_id}")
    
    def run(self, duration_seconds):
        """Run the simulated user for the specified duration."""
        end_time = time.time() + duration_seconds
        
        while time.time() < end_time:
            # Select an operation based on weights
            operation = self._select_weighted_operation()
            
            # Execute the operation
            try:
                if operation == "browse_molecules":
                    self._browse_molecules()
                elif operation == "view_molecule_details":
                    self._view_molecule_details()
                elif operation == "browse_mixtures":
                    self._browse_mixtures()
                elif operation == "view_mixture_details":
                    self._view_mixture_details()
                elif operation == "search_molecules":
                    self._search_molecules()
                elif operation == "create_mixture":
                    self._create_mixture()
                elif operation == "add_prediction":
                    self._add_prediction()
                elif operation == "record_experiment":
                    self._record_experiment()
            except Exception as e:
                logger.error(f"User {self.user_id} error executing {operation}: {str(e)}")
                metrics.record_error(operation, str(e))
            
            # Simulate user think time
            think_time = random.randint(
                TEST_CONFIG["user_think_time_ms"][0],
                TEST_CONFIG["user_think_time_ms"][1]
            ) / 1000.0  # Convert to seconds
            time.sleep(think_time)
    
    def _select_weighted_operation(self):
        """Select an operation based on the configured weights."""
        operations = list(TEST_CONFIG["operation_weights"].keys())
        weights = list(TEST_CONFIG["operation_weights"].values())
        return random.choices(operations, weights=weights, k=1)[0]
    
    def _time_operation(self, operation, func, *args, **kwargs):
        """Time an operation and record metrics."""
        start_time = time.time()
        success = True
        result = None
        
        try:
            result = func(*args, **kwargs)
        except Exception as e:
            success = False
            logger.error(f"User {self.user_id} error in {operation}: {str(e)}")
            metrics.record_error(operation, str(e))
            raise
        finally:
            end_time = time.time()
            elapsed_ms = (end_time - start_time) * 1000
            metrics.record_operation(operation, elapsed_ms, success)
        
        return result
    
    # User operations
    def _browse_molecules(self):
        """Browse molecules with pagination."""
        limit = random.randint(10, 50)
        offset = random.randint(0, 10) * limit
        
        return self._time_operation(
            "browse_molecules",
            self.supabase.table("molecule_with_properties").select("*").range(offset, offset + limit - 1).execute
        )
    
    def _view_molecule_details(self):
        """View details of a specific molecule."""
        if not test_data["molecule_ids"]:
            logger.warning("No molecule IDs available for view_molecule_details")
            return None
        
        molecule_id = random.choice(test_data["molecule_ids"])
        
        return self._time_operation(
            "view_molecule_details",
            self.supabase.table("molecule_with_properties").select("*").eq("id", molecule_id).execute
        )
    
    def _browse_mixtures(self):
        """Browse mixtures with pagination."""
        limit = random.randint(10, 50)
        offset = random.randint(0, 10) * limit
        
        return self._time_operation(
            "browse_mixtures",
            self.supabase.table("mixture_with_components").select("*").range(offset, offset + limit - 1).execute
        )
    
    def _view_mixture_details(self):
        """View details of a specific mixture."""
        if not test_data["mixture_ids"]:
            logger.warning("No mixture IDs available for view_mixture_details")
            return None
        
        mixture_id = random.choice(test_data["mixture_ids"])
        
        return self._time_operation(
            "view_mixture_details",
            self.supabase.table("mixture_with_components").select("*").eq("id", mixture_id).execute
        )
    
    def _search_molecules(self):
        """Search for molecules by name."""
        search_terms = ["glyc", "eth", "meth", "suc", "tre", "prop"]
        search_term = random.choice(search_terms)
        
        return self._time_operation(
            "search_molecules",
            self.supabase.table("molecule").select("*").ilike("name", f"%{search_term}%").execute
        )
    
    def _create_mixture(self):
        """Create a new mixture with components."""
        if not test_data["molecule_ids"] or len(test_data["molecule_ids"]) < 2:
            logger.warning("Not enough molecule IDs available for create_mixture")
            return None
        
        # Create a unique name for the test mixture
        name = f"Test Mixture {self.user_id}-{int(time.time())}"
        description = f"Test mixture created by simulated user {self.user_id}"
        
        # Select 2-3 random molecules for components
        num_components = random.randint(2, min(3, len(test_data["molecule_ids"])))
        selected_molecules = random.sample(test_data["molecule_ids"], num_components)
        
        # Create components with random concentrations
        components = []
        for molecule_id in selected_molecules:
            components.append({
                "molecule_id": molecule_id,
                "amount": round(random.uniform(1.0, 50.0), 2),
                "amount_unit": random.choice(["%w/v", "%v/v", "M"]),
                "role": "cryoprotectant"
            })
        
        def create_mixture_operation():
            # Insert mixture
            response = self.supabase.table("mixture").insert({
                "name": name,
                "description": description,
                "created_by": self.user_id
            }).execute()
            
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
                    "role": comp["role"],
                    "created_by": self.user_id
                })
            
            response = self.supabase.table("mixture_component").insert(component_inserts).execute()
            
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error adding components: {response.error}")
                # Rollback mixture creation
                self.supabase.table("mixture").delete().eq("id", mixture_id).execute()
                return None
            
            # Add the new mixture ID to the test data
            with threading.Lock():
                test_data["mixture_ids"].append(mixture_id)
            
            return mixture_id
        
        return self._time_operation("create_mixture", create_mixture_operation)
    
    def _add_prediction(self):
        """Add a prediction for a mixture."""
        if (not test_data["mixture_ids"] or not test_data["property_type_ids"] or 
            not test_data["calculation_method_ids"]):
            logger.warning("Missing IDs for add_prediction")
            return None
        
        mixture_id = random.choice(test_data["mixture_ids"])
        property_type_id = random.choice(test_data["property_type_ids"])
        calculation_method_id = random.choice(test_data["calculation_method_ids"])
        
        # Generate a random numeric value for the prediction
        value = round(random.uniform(-100.0, 100.0), 2)
        confidence = round(random.uniform(0.5, 1.0), 2)
        
        # Prepare prediction data
        prediction_data = {
            "mixture_id": mixture_id,
            "property_type_id": property_type_id,
            "calculation_method_id": calculation_method_id,
            "numeric_value": value,  # Assuming numeric for simplicity
            "confidence": confidence,
            "created_by": self.user_id
        }
        
        return self._time_operation(
            "add_prediction",
            self.supabase.table("predictions").upsert(
                prediction_data,
                on_conflict="mixture_id,property_type_id,calculation_method_id"
            ).execute
        )
    
    def _record_experiment(self):
        """Record an experiment for a mixture."""
        if not test_data["mixture_ids"] or not test_data["property_type_ids"]:
            logger.warning("Missing IDs for record_experiment")
            return None
        
        mixture_id = random.choice(test_data["mixture_ids"])
        property_type_id = random.choice(test_data["property_type_ids"])
        
        # Generate a random numeric value for the experiment
        value = round(random.uniform(-100.0, 100.0), 2)
        conditions = f"Test conditions at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        
        # Prepare experiment data
        experiment_data = {
            "mixture_id": mixture_id,
            "property_type_id": property_type_id,
            "numeric_value": value,  # Assuming numeric for simplicity
            "experimental_conditions": conditions,
            "date_performed": datetime.now().isoformat(),
            "created_by": self.user_id
        }
        
        return self._time_operation(
            "record_experiment",
            self.supabase.table("experiments").insert(experiment_data).execute
        )

def load_test_data():
    """Load test data from the database."""
    logger.info("Loading test data...")
    
    supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
    
    # Get molecule IDs
    response = supabase.table("molecule").select("id").execute()
    if hasattr(response, 'data') and response.data:
        test_data["molecule_ids"] = [mol["id"] for mol in response.data]
    
    # Get mixture IDs
    response = supabase.table("mixture").select("id").execute()
    if hasattr(response, 'data') and response.data:
        test_data["mixture_ids"] = [mix["id"] for mix in response.data]
    
    # Get property type IDs
    response = supabase.table("property_types").select("id").execute()
    if hasattr(response, 'data') and response.data:
        test_data["property_type_ids"] = [prop["id"] for prop in response.data]
    
    # Get calculation method IDs
    response = supabase.table("calculation_methods").select("id").execute()
    if hasattr(response, 'data') and response.data:
        test_data["calculation_method_ids"] = [method["id"] for method in response.data]
    
    logger.info(f"Loaded {len(test_data['molecule_ids'])} molecules, {len(test_data['mixture_ids'])} mixtures, "
               f"{len(test_data['property_type_ids'])} property types, and {len(test_data['calculation_method_ids'])} calculation methods")

def run_concurrent_load_test():
    """Run the concurrent load test."""
    logger.info("Starting concurrent load test...")
    
    # Load test data
    load_test_data()
    
    # Start metrics collection
    metrics.start_test()
    
    # Calculate how many users to add during each ramp-up interval
    ramp_up_interval = TEST_CONFIG["ramp_up_seconds"] / TEST_CONFIG["concurrent_users"]
    
    # Create and start user threads
    with concurrent.futures.ThreadPoolExecutor(max_workers=TEST_CONFIG["concurrent_users"]) as executor:
        futures = []
        
        # Ramp up users
        for i in range(TEST_CONFIG["concurrent_users"]):
            user = SimulatedUser(f"user-{i+1}")
            
            # Calculate how long this user should run
            user_duration = TEST_CONFIG["duration_seconds"] - (i * ramp_up_interval)
            
            # Submit the user task
            future = executor.submit(user.run, user_duration)
            futures.append(future)
            
            # Wait for the ramp-up interval before adding the next user
            if i < TEST_CONFIG["concurrent_users"] - 1:  # Don't wait after the last user
                time.sleep(ramp_up_interval)
        
        # Wait for all users to complete
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logger.error(f"User thread error: {str(e)}")
    
    # End metrics collection
    metrics.end_test()
    
    # Generate report
    summary = metrics.get_summary()
    
    # Save report to file
    report = {
        "timestamp": datetime.now().isoformat(),
        "config": TEST_CONFIG,
        "summary": summary
    }
    
    with open("concurrent_load_test_report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    # Generate text report
    text_report = []
    
    text_report.append("=" * 80)
    text_report.append("CryoProtect v2 Database Concurrent Load Test Report")
    text_report.append("=" * 80)
    text_report.append(f"Timestamp: {report['timestamp']}")
    text_report.append(f"Test Duration: {summary['test_duration']:.2f} seconds")
    text_report.append(f"Concurrent Users: {TEST_CONFIG['concurrent_users']}")
    text_report.append(f"Ramp-up Time: {TEST_CONFIG['ramp_up_seconds']} seconds")
    text_report.append("")
    
    text_report.append("Overall Performance:")
    text_report.append(f"  Total Operations: {summary['total_operations']}")
    text_report.append(f"  Average Response Time: {summary['avg_response_time_ms']:.2f} ms")
    text_report.append(f"  Overall Throughput: {summary['throughput_per_second']:.2f} operations/second")
    text_report.append("")
    
    text_report.append("Operation Performance:")
    for operation, stats in summary['operations'].items():
        text_report.append(f"  {operation}:")
        text_report.append(f"    Count: {stats['count']} (Success: {stats['success_count']}, Error: {stats['error_count']})")
        text_report.append(f"    Min: {stats['min_ms']:.2f} ms")
        text_report.append(f"    Max: {stats['max_ms']:.2f} ms")
        text_report.append(f"    Avg: {stats['avg_ms']:.2f} ms")
        text_report.append(f"    Median: {stats['median_ms']:.2f} ms")
        text_report.append(f"    P95: {stats['p95_ms']:.2f} ms")
        text_report.append(f"    Throughput: {stats['throughput_per_second']:.2f} operations/second")
        text_report.append("")
    
    if summary['errors']:
        text_report.append("Errors:")
        for operation, errors in summary['errors'].items():
            text_report.append(f"  {operation}:")
            for error in errors[:5]:  # Show only the first 5 errors
                text_report.append(f"    - {error}")
            if len(errors) > 5:
                text_report.append(f"    ... and {len(errors) - 5} more errors")
            text_report.append("")
    
    text_report.append("=" * 80)
    
    # Save text report to file
    with open("concurrent_load_test_report.txt", "w") as f:
        f.write("\n".join(text_report))
    
    logger.info("Concurrent load test completed. Reports saved to concurrent_load_test_report.json and concurrent_load_test_report.txt")
    
    return summary

def main():
    """Main function."""
    try:
        # Run concurrent load test
        summary = run_concurrent_load_test()
        
        print("\n" + "=" * 80)
        print("CryoProtect v2 Database Concurrent Load Testing Completed")
        print("=" * 80)
        print(f"Total Operations: {summary['total_operations']}")
        print(f"Average Response Time: {summary['avg_response_time_ms']:.2f} ms")
        print(f"Overall Throughput: {summary['throughput_per_second']:.2f} operations/second")
        print("=" * 80)
        print(f"Reports saved to concurrent_load_test_report.txt and concurrent_load_test_report.json")
        print("=" * 80)
        
        return 0
    except Exception as e:
        logger.error(f"Error running concurrent load test: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    exit(main())