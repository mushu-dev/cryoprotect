# CryoProtect Enhanced Monitoring System

This document provides comprehensive information about the enhanced monitoring system implemented for the CryoProtect application. The monitoring system combines various monitoring capabilities into a unified solution, providing real-time visibility into the application's health and performance.

## System Components

The enhanced monitoring system consists of the following components:

1. **Unified Monitoring Module** (`unified_monitoring.py`): Core monitoring functionality including database connection health, API endpoint health, system resource monitoring, performance metrics collection, progress tracking, error and alert management, and web dashboard.

2. **Monitoring Integrations** (`monitoring_integrations.py`): Integrations with external notification services including email, Slack, PagerDuty, SMS (via Twilio), and webhooks.

3. **Configuration Loader** (`monitoring_config_loader.py`): Configuration management from multiple sources including default values, JSON configuration file, environment variables, and command-line arguments.

4. **Command-Line Interface** (`monitoring_cli.py`): CLI tool for interacting with the monitoring system, viewing status and metrics, managing alerts, and generating reports.

5. **Monitoring Launcher** (`monitoring_launcher.py`): Main entry point for starting the monitoring system, setting up logging, initializing services, and providing a graceful shutdown mechanism.

6. **Dashboard Template** (`monitoring/dashboard_template.html`): HTML template for the web dashboard with responsive design, real-time updates, and visualizations.

## Getting Started

### Installation

The monitoring system is integrated with the CryoProtect application and doesn't require separate installation. However, it has the following dependencies:

- Python 3.6+
- Flask (for web dashboard)
- psutil (for system resource monitoring)
- requests (for HTTP requests to notification services)
- Additional dependencies for specific notification integrations:
  - twilio (for SMS notifications)

### Basic Usage

1. **Start the monitoring system**:

```bash
# Start with default configuration
./monitoring_launcher.py

# Start with custom configuration file
./monitoring_launcher.py --config custom_config.json

# Start with dashboard on specific port
./monitoring_launcher.py --dashboard-port 5002

# Start and open dashboard in browser
./monitoring_launcher.py --open-browser
```

2. **Access the dashboard**:

Open your web browser and navigate to:
```
http://localhost:5001/monitoring
```

3. **Check status using CLI**:

```bash
# Show current system health status
./monitoring_cli.py status

# Show detailed status
./monitoring_cli.py status --verbose
```

## Configuration

The monitoring system can be configured through the following methods:

1. **Configuration file** (`monitoring_config.json`): JSON file containing configuration settings.

2. **Environment variables**: Environment variables in the format `CRYOPROTECT_SECTION_KEY=value`.

3. **Command-line arguments**: Arguments in the format `--section-key=value`.

### Configuration Structure

The configuration is organized into the following sections:

- **monitoring**: Core monitoring settings (enabled, dashboard port, intervals, thresholds)
- **notifications**: Settings for notification integrations (email, Slack, PagerDuty, SMS, webhooks)
- **logging**: Logging configuration (level, file, format)
- **database_monitoring**: Database monitoring settings
- **api_monitoring**: API monitoring settings (endpoints, thresholds)
- **system_monitoring**: System monitoring settings (CPU, memory, disk thresholds)
- **progress_trackers**: Progress tracking settings
- **observability**: Distributed tracing settings

### Example Configuration

See `monitoring_config.json` for a complete example configuration.

## Monitoring Features

### Database Monitoring

- Connection health monitoring
- Pool size monitoring
- Query performance tracking
- Automatic fallback to alternative connection methods
- Circuit breaker pattern to prevent cascading failures

### API Monitoring

- Endpoint health checking
- Response time tracking
- Error rate monitoring
- Customizable thresholds for warnings and errors

### System Resource Monitoring

- CPU usage monitoring
- Memory usage monitoring
- Disk usage monitoring
- Network usage monitoring
- Process-specific metrics
- Customizable thresholds for warnings and errors

### Performance Metrics

- Operation timing statistics (min, max, average)
- Throughput tracking (operations per second)
- Memory usage tracking
- Custom metrics support

### Progress Tracking

- Batch operation progress tracking
- Completion time estimation
- Success/failure tracking
- Checkpoint support for resumable operations

### Observability

- Request tracing across components
- Correlation ID propagation
- Timing information for requests and operations
- Context propagation for debugging

### Error and Alert Management

- Error aggregation and classification
- Alert generation based on thresholds
- Alert silencing functionality
- Integration with external notification services

## Command-Line Interface

The CLI provides the following commands:

- **status**: Show current system health status
- **metrics**: Show performance metrics
- **progress**: Show progress of ongoing operations
- **alerts**: Show recent alerts
- **errors**: Show recent errors
- **silence**: Silence alerts by type
- **unsilence**: Unsilence alerts by type
- **report**: Generate monitoring report
- **dashboard**: Open the monitoring dashboard in a browser

### Examples

```bash
# Show current status
./monitoring_cli.py status

# Show database metrics
./monitoring_cli.py metrics --type database

# Show progress of ongoing operations
./monitoring_cli.py progress

# Show recent alerts
./monitoring_cli.py alerts

# Show recent errors
./monitoring_cli.py errors

# Silence an alert for 1 hour
./monitoring_cli.py silence --alert-type database_connection_failure --duration 3600

# Unsilence an alert
./monitoring_cli.py unsilence --alert-type database_connection_failure

# Generate a report
./monitoring_cli.py report --output report.txt

# Generate a JSON report
./monitoring_cli.py report --format json --output report.json

# Open dashboard in browser
./monitoring_cli.py dashboard
```

## Web Dashboard

The web dashboard provides a visual representation of the monitoring data, including:

- System health overview
- Resource usage charts
- API endpoint status
- Performance metrics
- Recent alerts and errors
- Progress of ongoing operations

The dashboard automatically refreshes every 30 seconds, but can also be manually refreshed.

## Notification Integrations

The monitoring system can send alerts to the following services:

- **Email**: Send alerts via SMTP
- **Slack**: Send alerts to Slack channels via webhooks
- **PagerDuty**: Create incidents in PagerDuty
- **SMS**: Send text messages via Twilio
- **Webhooks**: Send alerts to custom HTTP endpoints

Each integration can be configured with:
- Minimum severity level for alerts
- Throttling settings to prevent alert storms
- Custom formatting options

### Configuring Notifications

To enable a notification integration, update the corresponding section in the configuration file:

```json
"notifications": {
  "slack": {
    "enabled": true,
    "webhook_url": "https://hooks.slack.com/services/XXX/YYY/ZZZ",
    "channel": "#monitoring",
    "min_severity": "warning"
  }
}
```

## API Integration

The monitoring system can be integrated with your Flask application:

```python
from flask import Flask
from unified_monitoring import start_monitoring

# Create Flask app
app = Flask(__name__)

# Initialize monitoring
monitor = start_monitoring(app)

# Define routes
@app.route('/')
def index():
    return "Hello, World!"

if __name__ == '__main__':
    app.run()
```

## Code Integration

### Track Database Operations

```python
from unified_monitoring import track_database_operation

@track_database_operation(items_count=10)
def fetch_data_from_database():
    # Database operations
    return results
```

### Track API Operations

```python
from unified_monitoring import track_api_operation

@app.route('/api/resource')
@track_api_operation(name="get_resource")
def get_resource():
    # API endpoint logic
    return jsonify(resource)
```

### Track Progress

```python
from unified_monitoring import MonitoringService

# Get monitoring service
monitor = MonitoringService.get_instance()

# Create progress tracker
progress = monitor.track_progress("data_import", total_items=1000)

# Update progress
for batch in batches:
    # Process batch
    process_batch(batch)
    
    # Update progress
    progress.update(
        items_processed=len(batch),
        successful=success_count,
        failed=fail_count,
        batch_size=len(batch)
    )
```

### Record Errors

```python
from unified_monitoring import MonitoringService

# Get monitoring service
monitor = MonitoringService.get_instance()

try:
    # Risky operation
    result = perform_risky_operation()
except Exception as e:
    # Record error
    monitor.record_error(
        "operation_failure",
        str(e),
        traceback.format_exc(),
        {"context": "additional context"}
    )
```

### Manual Operation Tracking

```python
from unified_monitoring import MonitoringService

# Get monitoring service
monitor = MonitoringService.get_instance()

# Track operation
with monitor.track_operation("complex_operation", items_count=10):
    # Operation steps
    step1()
    step2()
    step3()
```

## Best Practices

1. **Start Early**: Initialize monitoring at application startup to capture all activity.

2. **Use the Right Tools**: Use decorators for automatic tracking when possible, and context managers for manual tracking when needed.

3. **Track Progress**: Use progress trackers for long-running operations to provide visibility.

4. **Add Context**: Include relevant context in errors and alerts to aid debugging.

5. **Set Appropriate Thresholds**: Configure appropriate thresholds for alerts to avoid false positives or alert fatigue.

6. **Review Dashboard Regularly**: Check the monitoring dashboard regularly to identify trends and potential issues.

7. **Use Notifications Judiciously**: Configure notifications only for important alerts to avoid overwhelming notification channels.

8. **Customize for Your Application**: Adjust monitoring parameters to suit your application's specific requirements and characteristics.

## Troubleshooting

### Dashboard Not Available

- Check if monitoring service is running
- Verify dashboard port is correct
- Ensure there are no other services using the same port
- Check for exceptions in the monitoring service logs

### Missing Metrics

- Ensure the corresponding monitoring component is enabled
- Check if the component is correctly integrated with your code
- Verify that operations are being tracked using the appropriate methods

### Notification Issues

- Check notification service credentials
- Verify network connectivity to the notification service
- Check notification service logs for errors
- Verify that the alert meets the minimum severity threshold

### Performance Impact

- Adjust health check intervals for less critical components
- Reduce the detail level of tracking for performance-sensitive operations
- Limit the number of endpoints monitored
- Use sampling for high-volume operations

## Contributing

Contributions to the monitoring system are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Write tests for your changes
5. Submit a pull request

## License

The monitoring system is licensed under the same license as the CryoProtect application.