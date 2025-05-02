# CryoProtect v2 Operations Guide

This document provides comprehensive operational documentation for the CryoProtect v2 system, focusing on monitoring, logging, backup, and observability components.

## Table of Contents

1. [Architecture and Design Decisions](#architecture-and-design-decisions)
2. [Configuration Options](#configuration-options)
3. [Operational Procedures](#operational-procedures)
4. [Integration with Other Components](#integration-with-other-components)
5. [Runbooks for Common Tasks](#runbooks-for-common-tasks)

## Architecture and Design Decisions

### Monitoring and Observability

CryoProtect v2 implements a comprehensive monitoring and observability stack based on industry best practices:

- **Prometheus**: Used for metrics collection and storage, providing a time-series database for all system metrics.
- **Grafana**: Provides visualization and dashboarding capabilities for metrics stored in Prometheus.
- **Alertmanager**: Handles alert routing, grouping, and notifications based on rules defined in Prometheus.
- **Node Exporter**: Collects system-level metrics from the host machine, including CPU, memory, disk, and network usage.
- **Custom Metrics**: Application-specific metrics are exposed via a Prometheus endpoint in the API.

The observability middleware in the application provides:
- Request tracing with correlation IDs
- Performance timing for API requests
- Detailed error reporting with full context
- Integration with the structured logging system

### Logging System

The logging system is designed to provide comprehensive visibility into application behavior:

- **JSON-structured logging**: All logs are formatted as JSON for easier parsing and analysis.
- **Correlation IDs**: Each request is assigned a unique correlation ID that is propagated through all logs.
- **ELK Stack**: Elasticsearch, Logstash, and Kibana are used for log storage, processing, and visualization.
- **Filebeat**: Collects log files and forwards them to Logstash.
- **Log Rotation**: Logs are automatically rotated based on size and age to prevent disk space issues.
- **Contextual Information**: Logs include detailed contextual information such as request details, user information, and system state.

### Backup System

The backup system ensures data durability and disaster recovery capabilities:

- **Scheduled Backups**: Automated backups are performed on daily, weekly, monthly, and yearly schedules.
- **Retention Policies**: Old backups are automatically pruned based on configurable retention policies.
- **Verification**: Backups are verified for integrity after creation.
- **Compression**: Backups are compressed to reduce storage requirements.
- **Encryption**: Optional encryption for sensitive data.
- **Cross-Region Storage**: Backups can be stored in multiple regions for disaster recovery.

## Configuration Options

### Monitoring Configuration

The monitoring system can be configured through environment variables in the Docker Compose file:

```yaml
environment:
  - PROMETHEUS_METRICS=1  # Enable Prometheus metrics (0 to disable)
  - PROMETHEUS_ENDPOINT=/metrics  # Endpoint for Prometheus metrics
```

Prometheus configuration is stored in `monitoring/prometheus/prometheus.yml`:

```yaml
global:
  scrape_interval: 15s
  evaluation_interval: 15s

alerting:
  alertmanagers:
    - static_configs:
        - targets:
            - alertmanager:9093

rule_files:
  - "rules/api_alerts.yml"
  - "rules/system_alerts.yml"

scrape_configs:
  - job_name: 'prometheus'
    static_configs:
      - targets: ['localhost:9090']

  - job_name: 'cryoprotect'
    static_configs:
      - targets: ['cryoprotect:5000']

  - job_name: 'node-exporter'
    static_configs:
      - targets: ['node-exporter:9100']
```

### Logging Configuration

Logging can be configured through environment variables:

```yaml
environment:
  - LOG_LEVEL=INFO  # Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
  - LOG_TO_FILE=1  # Enable file logging (0 to disable)
  - LOG_TO_CONSOLE=1  # Enable console logging (0 to disable)
  - LOG_TO_ELK=1  # Enable ELK logging (0 to disable)
  - LOG_JSON_FORMAT=1  # Enable JSON formatting (0 to disable)
  - ELASTICSEARCH_HOST=elasticsearch:9200  # Elasticsearch host
  - ELASTICSEARCH_INDEX=cryoprotect-logs  # Elasticsearch index
```

### Backup Configuration

Backup configuration is stored in `backup/config.json`:

```json
{
  "backup_dir": "backup/data",
  "retention": {
    "daily": 7,
    "weekly": 4,
    "monthly": 6,
    "yearly": 2
  },
  "schedule": {
    "daily": "02:00",
    "weekly": "sunday",
    "monthly": 1
  },
  "cross_region": {
    "enabled": false,
    "method": "s3",
    "config": {
      "bucket": "cryoprotect-backups",
      "region": "us-west-2",
      "access_key": "${AWS_ACCESS_KEY}",
      "secret_key": "${AWS_SECRET_KEY}"
    }
  },
  "verification": {
    "enabled": true,
    "methods": ["checksum", "restore_test"]
  },
  "compression": {
    "enabled": true,
    "method": "zip"
  },
  "encryption": {
    "enabled": false,
    "method": "aes256",
    "key_file": ""
  }
}
```

## Operational Procedures

### Setup and Deployment

1. **Initial Setup**:
   ```bash
   # Clone the repository
   git clone https://github.com/your-org/cryoprotect-v2.git
   cd cryoprotect-v2

   # Create .env file with required environment variables
   cp .env.template .env
   # Edit .env file with your configuration

   # Start the application with all components
   docker-compose --profile all up -d
   ```

2. **Production Deployment**:
   ```bash
   # Start only production components
   docker-compose --profile prod up -d
   ```

3. **Development Setup**:
   ```bash
   # Start only development components
   docker-compose --profile dev up -d
   ```

### Maintenance

#### Regular Maintenance Tasks

1. **Log Rotation**:
   - Logs are automatically rotated based on size and age.
   - Old logs are compressed and stored in the `logs` directory.
   - Configure rotation settings in the Docker Compose file.

2. **Backup Verification**:
   ```bash
   # Manually verify backups
   docker-compose exec backup-scheduler python -m backup.backup_manager --verify
   ```

3. **Database Maintenance**:
   ```bash
   # Run database maintenance tasks
   docker-compose exec cryoprotect python -m database.maintenance
   ```

4. **Monitoring System Updates**:
   ```bash
   # Update Prometheus configuration
   docker-compose exec prometheus promtool check config /etc/prometheus/prometheus.yml
   # Reload Prometheus configuration
   curl -X POST http://localhost:9090/-/reload
   ```

### Troubleshooting

#### Common Issues and Solutions

1. **API Service Not Starting**:
   - Check logs: `docker-compose logs cryoprotect`
   - Verify environment variables in `.env` file
   - Check database connection: `docker-compose exec cryoprotect python -m check_supabase_connection`

2. **Monitoring System Issues**:
   - Check Prometheus logs: `docker-compose logs prometheus`
   - Verify Prometheus targets: `http://localhost:9090/targets`
   - Check Grafana dashboards: `http://localhost:3000`

3. **Logging System Issues**:
   - Check Elasticsearch status: `curl -X GET "localhost:9200/_cluster/health"`
   - Check Logstash logs: `docker-compose logs logstash`
   - Verify Kibana connection: `http://localhost:5601`

4. **Backup System Issues**:
   - Check backup logs: `docker-compose logs backup-scheduler`
   - Verify backup directory permissions
   - Test manual backup: `docker-compose exec backup-scheduler python -m backup.backup_manager --backup`

## Integration with Other Components

### Integration with Supabase

CryoProtect v2 integrates with Supabase for database storage and authentication:

1. **Database Connection**:
   - Connection is managed through the Supabase client in `api/utils.py`
   - Connection pooling is implemented for performance
   - Monitoring includes database connection metrics

2. **Authentication**:
   - Supabase Auth is used for user authentication
   - JWT tokens are validated and used for authorization
   - Session management is integrated with the monitoring system

### Integration with External Systems

1. **Email Notifications**:
   - Alertmanager can send email notifications for alerts
   - Configure email settings in `monitoring/alertmanager/alertmanager.yml`

2. **Slack Notifications**:
   - Alertmanager can send Slack notifications for alerts
   - Configure Slack webhook in `monitoring/alertmanager/alertmanager.yml`

3. **Cloud Storage for Backups**:
   - Backups can be stored in S3, Azure Blob Storage, or Google Cloud Storage
   - Configure in `backup/config.json`

## Runbooks for Common Tasks

### Runbook: Monitoring System Management

#### Accessing Monitoring Dashboards

1. **Prometheus UI**:
   - URL: `http://localhost:9090`
   - Use for querying metrics and checking alert status

2. **Grafana Dashboards**:
   - URL: `http://localhost:3000`
   - Default credentials: admin/cryoprotect
   - Main dashboards:
     - API Performance: Overall API performance metrics
     - Database Performance: Database query performance
     - System Resources: Host system resource usage

#### Managing Alerts

1. **Viewing Active Alerts**:
   - Prometheus UI: `http://localhost:9090/alerts`
   - Alertmanager UI: `http://localhost:9093`

2. **Silencing Alerts**:
   - Alertmanager UI: `http://localhost:9093/#/silences`
   - Create a new silence with appropriate matcher

3. **Modifying Alert Rules**:
   - Edit files in `monitoring/prometheus/rules/`
   - Reload Prometheus configuration: `curl -X POST http://localhost:9090/-/reload`

### Runbook: Logging System Management

#### Accessing Logs

1. **Kibana UI**:
   - URL: `http://localhost:5601`
   - Navigate to "Discover" to search logs
   - Use saved searches for common queries

2. **Log Files**:
   - Location: `logs/` directory
   - Format: JSON or plain text depending on configuration

#### Creating Visualizations

1. **Kibana Visualizations**:
   - Navigate to "Visualize" in Kibana
   - Create new visualization based on log data
   - Add to dashboard for regular monitoring

2. **Log Analysis**:
   - Use Kibana's query language to filter logs
   - Example: `event_type:error AND user_id:123`
   - Save searches for common patterns

### Runbook: Backup and Restore

#### Creating Manual Backups

```bash
# Create a manual backup
docker-compose exec backup-scheduler python -m backup.backup_manager --backup

# Create a specific type of backup
docker-compose exec backup-scheduler python -m backup.backup_manager --backup-type daily
```

#### Restoring from Backup

```bash
# List available backups
docker-compose exec backup-scheduler python -m backup.backup_manager --list

# Restore from a specific backup
docker-compose exec backup-scheduler python -m backup.restore_backup --backup-dir backup/data/daily_backup_20250422_020000
```

#### Verifying Backup Integrity

```bash
# Verify all backups
docker-compose exec backup-scheduler python -m backup.backup_manager --verify-all

# Verify a specific backup
docker-compose exec backup-scheduler python -m backup.backup_manager --verify --backup-dir backup/data/daily_backup_20250422_020000
```

### Runbook: System Scaling

#### Scaling the Application

1. **Horizontal Scaling**:
   ```bash
   # Scale the API service
   docker-compose up -d --scale cryoprotect=3
   ```

2. **Vertical Scaling**:
   - Modify resource limits in Docker Compose file
   - Update and restart the service

#### Scaling the Monitoring System

1. **Prometheus Storage**:
   - Modify storage retention in `monitoring/prometheus/prometheus.yml`
   - Increase volume size for `prometheus_data`

2. **Elasticsearch Resources**:
   - Modify JVM heap size in Docker Compose file
   - Increase volume size for `elasticsearch-data`

### Runbook: Disaster Recovery

#### Complete System Failure

1. **Restore from Backups**:
   ```bash
   # Restore database from latest backup
   docker-compose exec backup-scheduler python -m backup.restore_backup --latest
   ```

2. **Rebuild Infrastructure**:
   ```bash
   # Rebuild all containers
   docker-compose down
   docker-compose --profile prod up -d
   ```

3. **Verify System Integrity**:
   ```bash
   # Run system verification
   docker-compose exec cryoprotect python -m system_verification
   ```

#### Data Corruption

1. **Identify Corrupted Data**:
   - Use monitoring alerts to identify issues
   - Check database consistency

2. **Restore Specific Data**:
   ```bash
   # Restore specific tables
   docker-compose exec backup-scheduler python -m backup.restore_backup --tables table1,table2
   ```

3. **Verify Data Integrity**:
   ```bash
   # Run data verification
   docker-compose exec cryoprotect python -m data_verification