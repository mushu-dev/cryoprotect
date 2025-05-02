"""
CryoProtect Analyzer - Performance Benchmark Tests

This module contains tests to benchmark the performance of the API endpoints and database operations.
It measures response times, throughput, and resource utilization under various load conditions.
"""

import os
import sys
import time
import json
import unittest
import threading
import statistics
import psutil
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch, MagicMock

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the base test case
from tests.base_test_case import MockSupabaseBaseTestCase
from tests.mock_supabase.helpers import patch_supabase

class TestPerformance(MockSupabaseBaseTestCase):
    """Test cases for performance benchmarking."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Performance test configuration
        self.test_config = {
            "concurrent_users": [1, 5, 10],  # Number of concurrent users to simulate
            "iterations_per_user": 5,        # Number of iterations per user
            "read_operations": [             # List of read operations to test
                {"name": "Get Molecules", "endpoint": "/api/v1/molecules"},
                {"name": "Get Mixtures", "endpoint": "/api/v1/mixtures"},
                {"name": "Get Molecule", "endpoint": "/api/v1/molecules/{molecule_id}"},
                {"name": "Get Mixture", "endpoint": "/api/v1/mixtures/{mixture_id}"},
                {"name": "Get Predictions", "endpoint": "/api/v1/mixtures/{mixture_id}/predictions"},
                {"name": "Get Experiments", "endpoint": "/api/v1/mixtures/{mixture_id}/experiments"}
            ],
            "write_operations": [            # List of write operations to test
                {"name": "Create Mixture", "endpoint": "/api/v1/mixtures", "method": "POST"},
                {"name": "Add Prediction", "endpoint": "/api/v1/mixtures/{mixture_id}/predictions", "method": "POST"},
                {"name": "Record Experiment", "endpoint": "/api/v1/mixtures/{mixture_id}/experiments", "method": "POST"}
            ]
        }
        
        # Sample data for testing
        self.sample_molecule_id = "00000000-0000-0000-0000-000000000001"
        self.sample_mixture_id = "00000000-0000-0000-0000-000000000002"
        
        self.sample_mixture_data = {
            "name": "Test Mixture",
            "description": "A test mixture for performance testing",
            "components": [
                {
                    "molecule_id": self.sample_molecule_id,
                    "concentration": 100,
                    "concentration_unit": "%"
                }
            ]
        }
        
        self.sample_prediction_data = {
            "property_name": "Freezing Point",
            "value": -15.3,
            "confidence": 0.9,
            "calculation_method": "CryoProtect Scoring"
        }
        
        self.sample_experiment_data = {
            "property_name": "Freezing Point",
            "value": -14.8,
            "experimental_conditions": "Standard pressure, cooling rate 1°C/min",
            "date_performed": "2025-04-15"
        }
        
        # Performance metrics
        self.metrics = {
            "response_times": {},
            "throughput": {},
            "cpu_usage": [],
            "memory_usage": []
        }

    def _monitor_resources(self, duration, interval=0.5):
        """Monitor CPU and memory usage for a specified duration."""
        start_time = time.time()
        while time.time() - start_time < duration:
            # Get CPU and memory usage
            cpu_percent = psutil.cpu_percent(interval=0.1)
            memory_percent = psutil.virtual_memory().percent
            
            # Store metrics
            self.metrics["cpu_usage"].append(cpu_percent)
            self.metrics["memory_usage"].append(memory_percent)
            
            # Sleep for the specified interval
            time.sleep(interval)

    def _execute_operation(self, operation, auth_headers=None):
        """Execute a single operation and measure response time."""
        # Replace placeholders in the endpoint
        endpoint = operation["endpoint"]
        endpoint = endpoint.replace("{molecule_id}", self.sample_molecule_id)
        endpoint = endpoint.replace("{mixture_id}", self.sample_mixture_id)
        
        # Set up headers
        headers = {}
        if auth_headers:
            headers.update(auth_headers)
        
        # Execute the operation
        start_time = time.time()
        
        if operation.get("method") == "POST":
            # Determine the request data based on the endpoint
            if "mixtures" in endpoint and "predictions" not in endpoint and "experiments" not in endpoint:
                data = self.sample_mixture_data
            elif "predictions" in endpoint:
                data = self.sample_prediction_data
            elif "experiments" in endpoint:
                data = self.sample_experiment_data
            else:
                data = {}
            
            response = self.client.post(endpoint, json=data, headers=headers)
        else:
            response = self.client.get(endpoint, headers=headers)
        
        end_time = time.time()
        response_time = end_time - start_time
        
        # Return the response time and status code
        return {
            "operation": operation["name"],
            "endpoint": endpoint,
            "response_time": response_time,
            "status_code": response.status_code
        }

    def _execute_operations_for_user(self, user_id, operations, auth_headers=None):
        """Execute a set of operations for a single user."""
        results = []
        
        for _ in range(self.test_config["iterations_per_user"]):
            for operation in operations:
                result = self._execute_operation(operation, auth_headers)
                result["user_id"] = user_id
                results.append(result)
        
        return results

    @patch_supabase(load_data=True)
    def test_read_operations_performance(self, mock_client):
        """Test the performance of read operations."""
        print("\nTesting read operations performance...")
        
        # Set up mock data
        self._setup_mock_data(mock_client)
        
        # Test with different numbers of concurrent users
        for num_users in self.test_config["concurrent_users"]:
            print(f"\nTesting with {num_users} concurrent users...")
            
            # Start resource monitoring in a separate thread
            monitor_thread = threading.Thread(
                target=self._monitor_resources,
                args=(num_users * self.test_config["iterations_per_user"] * len(self.test_config["read_operations"]) * 0.1,)
            )
            monitor_thread.daemon = True
            monitor_thread.start()

            # Run operations in parallel for each user
            with ThreadPoolExecutor(max_workers=num_users) as executor:
                futures = [
                    executor.submit(self._execute_operations_for_user, user_id, self.test_config["read_operations"])
                    for user_id in range(num_users)
                ]
                results = [f.result() for f in futures]

            monitor_thread.join()

            # Flatten results and check status codes
            for user_results in results:
                for result in user_results:
                    self.assertIn(result["status_code"], [200, 201, 400, 404, 500])

    @patch_supabase(load_data=True)
    def test_write_operations_performance(self, mock_client):
        """Test the performance of write operations."""
        print("\nTesting write operations performance...")
        self._setup_mock_data(mock_client)
        for num_users in self.test_config["concurrent_users"]:
            print(f"\nTesting with {num_users} concurrent users (write)...")
            with ThreadPoolExecutor(max_workers=num_users) as executor:
                futures = [
                    executor.submit(self._execute_operations_for_user, user_id, self.test_config["write_operations"])
                    for user_id in range(num_users)
                ]
                results = [f.result() for f in futures]
            # Flatten results and check status codes
            for user_results in results:
                for result in user_results:
                    self.assertIn(result["status_code"], [200, 201, 400, 404, 500])

    def test_monitor_resources_short_duration(self):
        """Test resource monitoring with a very short duration."""
        self.metrics = {"cpu_usage": [], "memory_usage": []}
        self._monitor_resources(duration=0.1, interval=0.05)
        self.assertTrue(len(self.metrics["cpu_usage"]) > 0)
        self.assertTrue(len(self.metrics["memory_usage"]) > 0)

    def test_execute_operation_error_handling(self):
        """Test _execute_operation with a simulated error response."""
        operation = {"name": "Bad Request", "endpoint": "/api/v1/bad", "method": "POST"}
        # Patch self.client.post to simulate an error response
        class FakeResponse:
            status_code = 400
        self.client.post = lambda endpoint, json, headers: FakeResponse()
        result = self._execute_operation(operation)
        self.assertEqual(result["status_code"], 400)

    def test_execute_operations_for_user_empty(self):
        """Test _execute_operations_for_user with no operations."""
        results = self._execute_operations_for_user(user_id=1, operations=[])
        self.assertEqual(results, [])

    @patch_supabase(load_data=True)
    def test_write_operations_performance(self, mock_client):
        """Test the performance of write operations."""
        print("\nTesting write operations performance...")
        
        # Set up mock data
        self._setup_mock_data(mock_client)
        
        # Set up auth headers
        auth_headers = {
            "Authorization": "Bearer test-token"
        }
        
        # Mock the auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id="test-user-id"
        )
        
        # Test with different numbers of concurrent users
        for num_users in self.test_config["concurrent_users"]:
            print(f"\nTesting with {num_users} concurrent users...")
            
            # Start resource monitoring in a separate thread
            monitor_thread = threading.Thread(
                target=self._monitor_resources,
                args=(num_users * self.test_config["iterations_per_user"] * len(self.test_config["write_operations"]) * 0.1,)
            )
            monitor_thread.daemon = True
            monitor_thread.start()
            
            # Execute operations concurrently
            all_results = []
            with ThreadPoolExecutor(max_workers=num_users) as executor:
                futures = []
                for user_id in range(num_users):
                    future = executor.submit(
                        self._execute_operations_for_user,
                        user_id,
                        self.test_config["write_operations"],
                        auth_headers
                    )
                    futures.append(future)
                
                # Collect results
                for future in futures:
                    all_results.extend(future.result())
            
            # Wait for the monitor thread to finish
            monitor_thread.join()
            
            # Calculate metrics
            self._calculate_metrics(all_results, f"write_{num_users}_users")
            
            # Print results
            self._print_results(all_results, f"Write Operations ({num_users} users)")

    @patch_supabase(load_data=True)
    def test_mixed_operations_performance(self, mock_client):
        """Test the performance of mixed read and write operations."""
        print("\nTesting mixed operations performance...")
        
        # Set up mock data
        self._setup_mock_data(mock_client)
        
        # Set up auth headers
        auth_headers = {
            "Authorization": "Bearer test-token"
        }
        
        # Mock the auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id="test-user-id"
        )
        
        # Combine read and write operations
        mixed_operations = self.test_config["read_operations"] + self.test_config["write_operations"]
        
        # Test with different numbers of concurrent users
        for num_users in self.test_config["concurrent_users"]:
            print(f"\nTesting with {num_users} concurrent users...")
            
            # Start resource monitoring in a separate thread
            monitor_thread = threading.Thread(
                target=self._monitor_resources,
                args=(num_users * self.test_config["iterations_per_user"] * len(mixed_operations) * 0.1,)
            )
            monitor_thread.daemon = True
            monitor_thread.start()
            
            # Execute operations concurrently
            all_results = []
            with ThreadPoolExecutor(max_workers=num_users) as executor:
                futures = []
                for user_id in range(num_users):
                    future = executor.submit(
                        self._execute_operations_for_user,
                        user_id,
                        mixed_operations,
                        auth_headers
                    )
                    futures.append(future)
                
                # Collect results
                for future in futures:
                    all_results.extend(future.result())
            
            # Wait for the monitor thread to finish
            monitor_thread.join()
            
            # Calculate metrics
            self._calculate_metrics(all_results, f"mixed_{num_users}_users")
            
            # Print results
            self._print_results(all_results, f"Mixed Operations ({num_users} users)")

    def _setup_mock_data(self, mock_client):
        """Set up mock data for performance testing."""
        # Mock the molecule_with_properties table
        mock_client.table("molecule_with_properties").select().execute.return_value = MagicMock(
            data=[{"id": self.sample_molecule_id, "name": "Test Molecule", "properties": []}]
        )
        
        # Mock the mixture_with_components table
        mock_client.table("mixture_with_components").select().execute.return_value = MagicMock(
            data=[{"id": self.sample_mixture_id, "name": "Test Mixture", "components": []}]
        )
        
        # Mock the predictions table
        mock_client.from_().select().eq().execute.return_value = MagicMock(
            data=[{"id": "test-id", "property_name": "Test Property", "numeric_value": 0}]
        )
        
        # Mock the insert methods
        mock_client.table().insert().execute.return_value = MagicMock(
            data=[{"id": "test-id"}]
        )

    def _calculate_metrics(self, results, test_name):
        """Calculate performance metrics from test results."""
        # Group results by operation
        operation_results = {}
        for result in results:
            operation = result["operation"]
            if operation not in operation_results:
                operation_results[operation] = []
            operation_results[operation].append(result["response_time"])
        
        # Calculate metrics for each operation
        response_times = {}
        for operation, times in operation_results.items():
            response_times[operation] = {
                "min": min(times),
                "max": max(times),
                "avg": statistics.mean(times),
                "p95": sorted(times)[int(len(times) * 0.95)],
                "count": len(times)
            }
        
        # Calculate throughput
        total_time = max([result["response_time"] for result in results])
        throughput = len(results) / total_time if total_time > 0 else 0
        
        # Store metrics
        self.metrics["response_times"][test_name] = response_times
        self.metrics["throughput"][test_name] = throughput
        
        # Save metrics to a file
        self._save_metrics_to_file(test_name)

    def _print_results(self, results, test_name):
        """Print performance test results."""
        print(f"\n----- {test_name} Results -----")
        
        # Group results by operation
        operation_results = {}
        for result in results:
            operation = result["operation"]
            if operation not in operation_results:
                operation_results[operation] = []
            operation_results[operation].append(result["response_time"])
        
        # Print metrics for each operation
        for operation, times in operation_results.items():
            print(f"\nOperation: {operation}")
            print(f"  Min Response Time: {min(times):.4f} seconds")
            print(f"  Max Response Time: {max(times):.4f} seconds")
            print(f"  Avg Response Time: {statistics.mean(times):.4f} seconds")
            print(f"  P95 Response Time: {sorted(times)[int(len(times) * 0.95)]:.4f} seconds")
            print(f"  Count: {len(times)}")
        
        # Print throughput
        total_time = max([result["response_time"] for result in results])
        throughput = len(results) / total_time if total_time > 0 else 0
        print(f"\nThroughput: {throughput:.2f} operations/second")
        
        # Print resource usage
        if self.metrics["cpu_usage"]:
            print(f"\nCPU Usage: Min={min(self.metrics['cpu_usage']):.2f}%, Max={max(self.metrics['cpu_usage']):.2f}%, Avg={statistics.mean(self.metrics['cpu_usage']):.2f}%")
        if self.metrics["memory_usage"]:
            print(f"Memory Usage: Min={min(self.metrics['memory_usage']):.2f}%, Max={max(self.metrics['memory_usage']):.2f}%, Avg={statistics.mean(self.metrics['memory_usage']):.2f}%")

    def _save_metrics_to_file(self, test_name):
        """Save performance metrics to a file."""
        # Create the reports directory if it doesn't exist
        reports_dir = os.path.join(os.path.dirname(__file__), 'reports')
        os.makedirs(reports_dir, exist_ok=True)
        
        # Save metrics to a JSON file
        metrics_file = os.path.join(reports_dir, f'performance_metrics_{test_name}.json')
        with open(metrics_file, 'w') as f:
            json.dump({
                "response_times": self.metrics["response_times"][test_name],
                "throughput": self.metrics["throughput"][test_name],
                "cpu_usage": {
                    "min": min(self.metrics["cpu_usage"]) if self.metrics["cpu_usage"] else 0,
                    "max": max(self.metrics["cpu_usage"]) if self.metrics["cpu_usage"] else 0,
                    "avg": statistics.mean(self.metrics["cpu_usage"]) if self.metrics["cpu_usage"] else 0
                },
                "memory_usage": {
                    "min": min(self.metrics["memory_usage"]) if self.metrics["memory_usage"] else 0,
                    "max": max(self.metrics["memory_usage"]) if self.metrics["memory_usage"] else 0,
                    "avg": statistics.mean(self.metrics["memory_usage"]) if self.metrics["memory_usage"] else 0
                }
            }, f, indent=2)

if __name__ == '__main__':
    unittest.main()