#!/usr/bin/env python3
"""
Example of using the Unified Monitoring Module for CryoProtect v2

This script demonstrates how to use the unified_monitoring module in various scenarios:
1. Basic monitoring setup
2. Database operation tracking
3. API endpoint monitoring
4. Progress tracking for batch operations
5. Custom metrics collection
6. Error and alert handling
7. Monitoring dashboard

Usage:
    python monitoring_example.py

Dependencies:
    - Flask (for web dashboard)
    - psutil (for system monitoring)
"""

import time
import random
import threading
from flask import Flask, jsonify

from unified_monitoring import (
    MonitoringService, 
    track_operation, 
    track_database_operation, 
    track_api_operation, 
    start_monitoring
)

# Create a Flask application
app = Flask(__name__)

# Initialize the monitoring service with the Flask app
monitor = start_monitoring(app, dashboard_port=5001)

# Example routes to monitor
@app.route('/api/health')
def health_check():
    """Simple health check endpoint."""
    return jsonify({"status": "healthy"})

@app.route('/api/examples')
@track_api_operation(name="list_examples")
def get_examples():
    """Example API endpoint with performance tracking."""
    # Simulate some processing
    time.sleep(0.05)
    return jsonify({"examples": ["example1", "example2", "example3"]})

@app.route('/api/simulate-error')
def simulate_error():
    """Endpoint that simulates an error for testing error tracking."""
    try:
        # Deliberately cause an error
        result = 1 / 0
    except Exception as e:
        # Record the error with the monitoring service
        monitor.record_error(
            "division_by_zero", 
            str(e), 
            None,  # No stack trace needed for this example
            {"endpoint": "/api/simulate-error"}
        )
        return jsonify({"error": "Internal server error"}), 500
    
    return jsonify({"result": result})

# Example database operations
@track_database_operation(items_count=1)
def get_database_record(record_id):
    """Example database query function with performance tracking."""
    # Simulate database query
    time.sleep(0.02 + (random.random() * 0.05))
    return {"id": record_id, "name": f"Record {record_id}", "value": random.randint(1, 100)}

@track_database_operation(items_count=10)
def batch_database_operation():
    """Example batch database operation with performance tracking."""
    # Simulate batch database operation
    time.sleep(0.1 + (random.random() * 0.1))
    return [{"id": i, "processed": True} for i in range(10)]

# Example long-running process with progress tracking
def long_running_process():
    """Example of a long-running process with progress tracking."""
    # Create a progress tracker
    total_items = 100
    progress = monitor.track_progress("data_import", total_items)
    
    for i in range(0, total_items, 5):
        # Simulate processing a batch of items
        time.sleep(0.2)
        
        # Randomly simulate some failures
        successful = 5
        failed = 0
        if random.random() < 0.1:  # 10% chance of failure
            successful = 4
            failed = 1
        
        # Update progress
        progress.update(5, successful=successful, failed=failed, batch_size=5)
        
        # Occasionally add custom metrics
        if i % 20 == 0:
            monitor.performance_metrics["database"].record_custom_metric(
                f"batch_{i}_special_metric", 
                random.random() * 100
            )
    
    return progress.generate_report()

# Example of using context manager for operation tracking
def manual_operation_tracking():
    """Example of manually tracking an operation with context manager."""
    with monitor.track_operation("manual_complex_operation", items_count=3):
        # First step
        time.sleep(0.1)
        
        # Second step
        get_database_record(1)
        
        # Third step
        time.sleep(0.05)
        
        # Log a custom metric
        monitor.performance_metrics["api"].record_custom_metric(
            "special_timing", time.time()
        )

# Function to simulate database health changes
def simulate_database_health_changes():
    """Simulates changes in database health for testing."""
    while True:
        # Simulate an occasional database issue
        if random.random() < 0.05:  # 5% chance of database issue
            monitor._trigger_alert(
                "database_connection_issue",
                "Simulated database connection issue for testing",
                "database",
                "warning"
            )
            # Override health status for demo purposes
            monitor.health_status["database"]["status"] = "Warning"
            time.sleep(20)  # Keep in warning state for 20 seconds
            monitor.health_status["database"]["status"] = "OK"
        
        time.sleep(10)

# Example alert handler
def example_alert_handler(alert):
    """Example handler function for alerts."""
    print(f"ALERT RECEIVED: {alert['severity']} - {alert['message']}")
    
    # In a real application, you might send an email, SMS, or webhook notification
    # This is just a simple example

# Main example function
def run_example():
    """Run the full monitoring example."""
    print("Starting Unified Monitoring Example...")
    print(f"Dashboard available at: http://localhost:5001/monitoring")
    
    # Register an alert handler
    monitor.add_alert_handler(example_alert_handler)
    
    # Start a thread to simulate database health changes
    db_simulator = threading.Thread(
        target=simulate_database_health_changes, 
        daemon=True
    )
    db_simulator.start()
    
    # Run examples in a loop
    try:
        print("Running examples (press Ctrl+C to stop)...")
        iteration = 0
        
        while True:
            iteration += 1
            print(f"Iteration {iteration}")
            
            # Run database operations
            for i in range(5):
                get_database_record(random.randint(1, 100))
            
            # Run batch operation
            batch_database_operation()
            
            # Manual tracking example
            manual_operation_tracking()
            
            # Every 5 iterations, run the long process
            if iteration % 5 == 0:
                print("Starting long-running process...")
                progress_report = long_running_process()
                print("Long-running process complete!")
                print(progress_report)
            
            # Every 10 iterations, trigger a test error
            if iteration % 10 == 0:
                try:
                    # Simulate an application error
                    raise ValueError("Example error for testing error recording")
                except Exception as e:
                    monitor.record_error(
                        "example_error", 
                        str(e), 
                        None,
                        {"iteration": iteration}
                    )
                    print(f"Simulated error recorded in iteration {iteration}")
            
            # Sleep between iterations
            time.sleep(2)
            
    except KeyboardInterrupt:
        print("\nStopping example...")
    finally:
        # Stop the monitoring service
        monitor.stop()
        print("Monitoring service stopped.")

# Run the Flask app in a separate thread
def run_flask_app():
    app.run(host='0.0.0.0', port=5000, debug=False)

if __name__ == "__main__":
    # Start Flask in a separate thread
    flask_thread = threading.Thread(target=run_flask_app)
    flask_thread.daemon = True
    flask_thread.start()
    
    # Run the main example
    run_example()