# Unified Monitoring System for CryoProtect v2

This document provides an overview of the unified monitoring system implemented in the CryoProtect v2 application. The monitoring system combines various monitoring capabilities into a single, easy-to-use interface.

## Features

- **Database Connection Monitoring**: Track database connection health and performance
- **API Endpoint Health Checking**: Monitor API endpoint status and response times
- **System Resource Monitoring**: Track CPU, memory, and disk usage
- **Performance Metrics Collection**: Gather detailed metrics on operations
- **Progress Tracking**: Monitor and estimate completion time for batch operations
- **Observability and Tracing**: Implement distributed tracing for requests
- **Error Aggregation**: Collect and analyze errors
- **Alerting System**: Generate alerts based on configurable thresholds
- **Web Dashboard**: Visualize monitoring data in real-time
- **Metrics Persistence**: Save monitoring data to disk for historical analysis

## Quick Start

### Basic Setup

```python
from unified_monitoring import start_monitoring

# Initialize with Flask app (optional)
monitor = start_monitoring(app, dashboard_port=5001)

# Start monitoring
monitor.start()
```

### Track Database Operations

```python
from unified_monitoring import track_database_operation

@track_database_operation(items_count=10)
def fetch_data_from_database():
    # Database operations here
    return results
```

### Track API Endpoints

```python
from unified_monitoring import track_api_operation

@app.route('/api/resource')
@track_api_operation(name="get_resource")
def get_resource():
    # API endpoint logic here
    return jsonify(resource)
```

### Track Progress for Batch Operations

```python
# Create a progress tracker
progress = monitor.track_progress("data_import", total_items=1000)

# Update as processing occurs
for batch in batches:
    process_batch(batch)
    progress.update(
        items_processed=len(batch),
        successful=success_count,
        failed=fail_count
    )

# Get a progress report
report = progress.generate_report()
```

### Manual Operation Tracking

```python
# Use the context manager for custom operations
with monitor.track_operation("complex_task", items_count=10):
    # Your complex task code here
    pass
```

### Error Recording

```python
try:
    # Risky operation
    result = perform_risky_operation()
except Exception as e:
    # Record the error
    monitor.record_error(
        "operation_failure", 
        str(e),
        traceback.format_exc(),
        {"context": "additional context data"}
    )
```

### Handle Alerts

```python
def alert_handler(alert):
    # The alert object contains details about the alert
    severity = alert["severity"]
    message = alert["message"]
    source = alert["source"]
    
    # Send email, SMS, or other notification
    send_notification(f"{severity} alert from {source}: {message}")

# Register the alert handler
monitor.add_alert_handler(alert_handler)
```

## Monitoring Dashboard

The monitoring system includes a web dashboard that provides real-time visualization of monitoring data. To access the dashboard:

1. Start the monitoring service with dashboard enabled
2. Navigate to `http://localhost:5001/monitoring`

The dashboard displays:
- Overall system health status
- Database and API endpoint status
- System resource usage
- Recent errors and alerts
- Performance metrics
- Progress status for batch operations

## API Reference

### Monitoring Service

The `MonitoringService` class is the main entry point for the monitoring system:

```python
from unified_monitoring import MonitoringService

# Get the singleton instance
monitor = MonitoringService.get_instance()

# Or create and initialize with an app
monitor = MonitoringService(app)
```

Important methods:
- `start()`: Start all monitoring services
- `stop()`: Stop all monitoring services
- `track_operation(name, items_count=0, context=None)`: Context manager for operation tracking
- `track_progress(name, total_items, checkpoint_file=None)`: Create a progress tracker
- `record_error(error_type, message, stack_trace=None, context=None)`: Record an error
- `add_alert_handler(handler)`: Add an alert handler
- `get_monitoring_data()`: Get all monitoring data as a dictionary

### Decorators

The module provides decorators for easy integration:

- `@track_operation(operation_type, items_count=0, context=None)`: General operation tracking
- `@track_database_operation(items_count=0)`: Database operation tracking
- `@track_api_operation(name=None)`: API endpoint tracking

### Progress Tracker

The `ProgressTracker` class tracks progress of batch operations:

Important methods:
- `update(items_processed, successful=None, failed=None, batch_size=None)`: Update progress
- `get_status()`: Get current progress status
- `generate_report()`: Generate a human-readable progress report
- `reset()`: Reset progress to initial state

### Performance Metrics

The `PerformanceMetrics` class collects performance metrics:

Important methods:
- `record_operation(operation_type, success, execution_time, items_processed=0)`: Record an operation
- `record_custom_metric(metric_name, value)`: Record a custom metric
- `get_metrics()`: Get current metrics
- `generate_report()`: Generate a human-readable report
- `reset()`: Reset metrics to initial state

## Example Usage

See `monitoring_example.py` for a complete example of using the unified monitoring system. This example demonstrates:

1. Setting up the monitoring service
2. Tracking database operations
3. Monitoring API endpoints
4. Tracking progress for batch operations
5. Recording errors and handling alerts
6. Using the monitoring dashboard

## Configuration

The monitoring system can be configured through the `config` dictionary on the `MonitoringService` instance:

```python
monitor = MonitoringService.get_instance()
monitor.config.update({
    "monitoring_dir": "monitoring",
    "metrics_retention_days": 14,
    "health_check_intervals": {
        "database": 30,  # seconds
        "api": 60,       # seconds
        "system": 15     # seconds
    },
    "alert_thresholds": {
        "database_failures": 5,
        "api_errors": 3,
        "memory_usage": 85,  # percentage
        "cpu_usage": 85,     # percentage
        "disk_usage": 90     # percentage
    }
})
```

## Dependencies

The unified monitoring system has the following dependencies:

- **Flask**: For the web dashboard (optional)
- **psutil**: For system resource monitoring
- **threading**: For concurrent monitoring
- **json**: For metrics persistence
- **datetime**: For time-based operations

## Architecture

The unified monitoring system follows a modular architecture:

1. **MonitoringService**: Central coordinator for all monitoring activities
2. **PerformanceMetrics**: Collects and stores performance metrics
3. **ProgressTracker**: Tracks progress of batch operations
4. **Health Checkers**: Monitor database, API, and system health
5. **Dashboard**: Web interface for visualization of monitoring data
6. **Alerting System**: Generates and dispatches alerts

Monitoring data flows from various parts of the application through the `MonitoringService`, which aggregates and persists the data.

## Best Practices

1. **Start Early**: Initialize monitoring at application startup
2. **Track Appropriately**: Use the right tracking method for each operation
3. **Handle Alerts**: Register alert handlers for important alerts
4. **Tune Alert Thresholds**: Adjust thresholds to avoid alert fatigue
5. **Review Dashboards**: Regularly check the monitoring dashboard
6. **Persist Metrics**: Configure metrics persistence for historical analysis

## Troubleshooting

### Common Issues

1. **Dashboard not showing**: Ensure Flask is installed and the application is running
2. **Missing system metrics**: Ensure psutil is installed
3. **High memory usage**: Adjust history size for errors, alerts, and metrics
4. **Alert fatigue**: Tune alert thresholds and use silencing

For more detailed troubleshooting, check the logs generated by the monitoring system.