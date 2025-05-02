# CryoProtect v2 Operational Runbooks

This document contains detailed step-by-step runbooks for common operational tasks related to monitoring, logging, backup, and observability components of CryoProtect v2.

## Table of Contents

1. [Monitoring Runbooks](#monitoring-runbooks)
2. [Logging Runbooks](#logging-runbooks)
3. [Backup Runbooks](#backup-runbooks)
4. [Observability Runbooks](#observability-runbooks)
5. [Maintenance Runbooks](#maintenance-runbooks)

## Monitoring Runbooks

### Setting Up New Grafana Dashboards

**Purpose**: Create new dashboards in Grafana to visualize system metrics.

**Steps**:

1. Log in to Grafana (http://localhost:3000)
2. Click "+" icon in the left sidebar and select "Dashboard"
3. Click "Add new panel"
4. Select Prometheus as the data source
5. Enter a PromQL query (e.g., `rate(http_requests_total{job="cryoprotect"}[5m])`)
6. Set visualization type (Graph, Gauge, etc.)
7. Configure axes, legends, and thresholds
8. Add more panels as needed
9. Click "Save dashboard" and enter a name and description

### Configuring Alert Rules

**Purpose**: Set up alerting rules in Prometheus to notify when metrics exceed thresholds.

**Steps**:

1. Create or edit alert rules file in `monitoring/prometheus/rules/`
2. Define alert rules with expressions, labels, and annotations
3. Validate the rules file with `promtool check rules`
4. Update Prometheus configuration to include the rules file
5. Reload Prometheus configuration
6. Configure Alertmanager for notification routing
7. Verify alerts in Prometheus UI

### Responding to Alerts

**Purpose**: Effectively respond to and resolve alerts from the monitoring system.

**Steps**:

1. Acknowledge the alert in Alertmanager
2. Investigate using Grafana dashboards
3. Check logs in Kibana for related errors
4. Determine root cause
5. Implement resolution (restart service, scale, fix application issue)
6. Verify resolution in Prometheus UI
7. Document incident and resolution
8. Implement preventive measures

## Logging Runbooks

### Configuring Log Rotation

**Purpose**: Configure log rotation to manage disk space and maintain log history.

**Steps**:

1. Review current log configuration in docker-compose.yml
2. Modify log rotation settings (size and backup count)
3. Apply changes by restarting the service
4. Verify configuration by checking log directory structure
5. Optionally configure system logrotate for additional management

### Setting Up Log Filters in Kibana

**Purpose**: Create saved searches and filters in Kibana for efficient log analysis.

**Steps**:

1. Access Kibana and create index pattern for logs
2. Navigate to Discover view
3. Create basic filters for error logs, API requests, and performance issues
4. Create advanced filters for specific scenarios
5. Save searches for future use
6. Create dashboards combining saved searches
7. Set up automated reports if needed

### Troubleshooting Missing Logs

**Purpose**: Diagnose and resolve issues with missing logs in the logging system.

**Steps**:

1. Verify application is generating logs
2. Check Filebeat status and configuration
3. Verify Logstash status and configuration
4. Check Elasticsearch status and indices
5. Restart components if needed
6. Generate test logs and verify flow
7. Check for index mapping issues
8. Rebuild indices if necessary

## Backup Runbooks

### Configuring Backup Schedule

**Purpose**: Configure automated backup schedule for database and application data.

**Steps**:

1. Review current backup configuration in `backup/config.json`
2. Modify backup schedule (daily, weekly, monthly)
3. Update retention policy
4. Configure cross-region backup if needed
5. Apply configuration by restarting backup scheduler
6. Verify configuration in logs
7. Test backup creation

### Performing Database Restore

**Purpose**: Restore database from a backup in case of data loss or corruption.

**Steps**:

1. List available backups
2. Stop the application
3. Create a backup of the current state if possible
4. Perform the restore operation
5. Verify restore success
6. Start the application
7. Verify application functionality
8. Document the restore operation

### Verifying Backup Integrity

**Purpose**: Verify the integrity of backups to ensure they can be used for restoration.

**Steps**:

1. List available backups
2. Verify specific backups or all backups
3. Check verification results
4. Examine verification details in logs
5. Address any verification failures
6. Perform test restore for critical backups
7. Document verification results

## Observability Runbooks

### Tracing Requests Across Services

**Purpose**: Trace requests as they flow through different services to diagnose issues.

**Steps**:

1. Identify the correlation ID from logs or response headers
2. Search logs by correlation ID in Kibana
3. Analyze the request flow from start to completion
4. Identify performance bottlenecks
5. Create visualizations if needed
6. Check for errors or exceptions
7. Document findings

### Analyzing Performance Bottlenecks

**Purpose**: Identify and resolve performance bottlenecks in the application.

**Steps**:

1. Identify slow endpoints in Grafana dashboards
2. Check system resource usage (CPU, memory, disk, network)
3. Analyze database performance metrics
4. Look for slow functions in logs
5. Profile the application if needed
6. Identify optimization opportunities
7. Implement improvements
8. Measure impact of changes
9. Document optimizations

### Setting Up Custom Metrics

**Purpose**: Create custom metrics to monitor application-specific behavior.

**Steps**:

1. Identify metrics to track (business, application-specific, performance)
2. Add metric definitions to the code
3. Instrument the code to record metrics
4. Expose metrics endpoint
5. Update Prometheus configuration to scrape new metrics
6. Create Grafana dashboards for new metrics
7. Set up alerts for new metrics if needed
8. Document new metrics and their purpose

## Maintenance Runbooks

### Performing System Updates

**Purpose**: Safely update system components without disrupting service.

**Steps**:

1. Create a backup before updates
2. Update Docker images in docker-compose.yml
3. Pull new images
4. Stop and restart services one at a time
5. Verify each service after update
6. Check logs for errors
7. Verify system functionality
8. Document updates and any issues

### Scaling System Components

**Purpose**: Scale system components to handle increased load.

**Steps**:

1. Identify components that need scaling
2. Update resource limits in docker-compose.yml
3. Scale horizontally by increasing replica count
4. Update load balancer configuration if needed
5. Monitor performance during scaling
6. Adjust scaling based on metrics
7. Document scaling changes and impact

### Database Maintenance

**Purpose**: Perform regular database maintenance to ensure optimal performance.

**Steps**:

1. Schedule maintenance window
2. Create a backup before maintenance
3. Run database vacuum and analyze
4. Update statistics
5. Check for and fix index fragmentation
6. Verify query performance
7. Document maintenance activities and results
