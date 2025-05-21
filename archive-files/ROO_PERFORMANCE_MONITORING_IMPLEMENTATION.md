# Performance Monitoring Implementation for CryoProtect v2

## Overview

As part of Phase 3.2 of the CryoProtect v2 project, we need to implement a comprehensive performance monitoring system. While we've already made progress on the centralized logging system, we now need to focus on monitoring application performance, system resources, database performance, and setting up alerting.

## Current Status

- Enhanced logging system implemented with JSON structured logs and ELK stack integration (README_Enhanced_Logging.md)
- Blue/green deployment completed (Phase 3.1, Stage 2)
- Database foundation solidified (Phase 3.1, Stage 1)

## Implementation Goals

1. Application Performance Monitoring
2. System Resource Monitoring
3. Database Performance Tracking
4. Real-time Health Dashboard
5. Performance Alerting System

## Detailed Tasks

### 1. Application Performance Monitoring (APM)

Implement a comprehensive APM solution to track application performance metrics:

- **Response time tracking**: Measure response times of all API endpoints
- **Endpoint usage metrics**: Track the most frequently used endpoints
- **Error rate monitoring**: Monitor application errors and exceptions
- **Transaction tracing**: Trace requests across different components
- **Bottleneck identification**: Identify performance bottlenecks

Implementation approach:
- Use OpenTelemetry for instrumentation
- Integrate with Prometheus for metrics collection
- Set up Grafana dashboards for visualization

Files to create/modify:
- `monitoring/apm.py`: APM configuration and middleware
- `monitoring/middleware.py`: Request timing middleware for Flask
- `monitoring/exporters.py`: Exporters for metrics
- `app.py`: Add APM initialization
- `docker-compose.monitoring.yml`: Add monitoring stack containers

### 2. System Resource Monitoring

Implement system-level monitoring for all infrastructure components:

- **CPU usage**: Monitor CPU utilization across all services
- **Memory usage**: Track memory consumption
- **Disk usage and I/O**: Monitor disk space and I/O operations
- **Network traffic**: Monitor network bandwidth utilization
- **Container health**: Monitor Docker container status and health

Implementation approach:
- Use Prometheus Node Exporter for system metrics
- Configure cAdvisor for container metrics
- Set up alerting thresholds for resource usage

Files to create/modify:
- `monitoring/system_metrics.py`: System metrics collection
- `docker-compose.monitoring.yml`: Add Node Exporter and cAdvisor
- `prometheus/prometheus.yml`: Configure Prometheus for system metrics
- `grafana/dashboards/system_resources.json`: System resources dashboard

### 3. Database Performance Tracking

Implement database performance monitoring:

- **Query performance**: Monitor slow queries and query execution times
- **Connection pool metrics**: Track database connection usage
- **Index usage**: Monitor index utilization
- **Table statistics**: Track table growth and usage patterns
- **Database load**: Monitor overall database load

Implementation approach:
- Use pg_stat_statements and custom queries for Postgres monitoring
- Set up Postgres Exporter for Prometheus integration
- Create database-specific dashboards in Grafana

Files to create/modify:
- `monitoring/db_metrics.py`: Database metrics collection
- `docker-compose.monitoring.yml`: Add Postgres Exporter
- `prometheus/prometheus.yml`: Configure Prometheus for database metrics
- `grafana/dashboards/database_performance.json`: Database dashboard

### 4. Real-time Health Dashboard

Create comprehensive dashboards for monitoring system health:

- **Overview dashboard**: High-level system health metrics
- **Application dashboard**: Detailed application performance metrics
- **Database dashboard**: Database performance metrics
- **Infrastructure dashboard**: System resource utilization metrics
- **Alert dashboard**: Active alerts and historical alert data

Implementation approach:
- Use Grafana for dashboards
- Set up dashboard provisioning for automated setup
- Create custom panels for key metrics

Files to create/modify:
- `grafana/dashboards/`: Dashboard JSON files
- `grafana/provisioning/`: Dashboard provisioning configuration
- `monitoring/dashboard.py`: Dashboard setup and configuration

### 5. Performance Alerting System

Implement a comprehensive alerting system:

- **Alert rules**: Define thresholds for various metrics
- **Alert channels**: Configure multiple notification channels (email, Slack, PagerDuty)
- **Alert grouping**: Group related alerts to reduce noise
- **Alert history**: Track historical alerts for analysis
- **Escalation policies**: Define escalation procedures for critical alerts

Implementation approach:
- Use Alertmanager for alert management
- Set up Prometheus alert rules
- Configure notification channels and templates

Files to create/modify:
- `prometheus/alert_rules.yml`: Alert rule definitions
- `alertmanager/alertmanager.yml`: Alertmanager configuration
- `monitoring/alerts.py`: Alert configuration helpers
- `templates/alerts/`: Alert notification templates

## Implementation Details

### 1. Prometheus Configuration

```yaml
# prometheus/prometheus.yml
global:
  scrape_interval: 15s
  evaluation_interval: 15s

rule_files:
  - "alert_rules.yml"

alerting:
  alertmanagers:
    - static_configs:
        - targets: ['alertmanager:9093']

scrape_configs:
  - job_name: 'prometheus'
    static_configs:
      - targets: ['localhost:9090']

  - job_name: 'cryoprotect'
    metrics_path: '/metrics'
    static_configs:
      - targets: ['cryoprotect:5000']

  - job_name: 'node-exporter'
    static_configs:
      - targets: ['node-exporter:9100']

  - job_name: 'cadvisor'
    static_configs:
      - targets: ['cadvisor:8080']

  - job_name: 'postgres-exporter'
    static_configs:
      - targets: ['postgres-exporter:9187']
```

### 2. Alert Rules

```yaml
# prometheus/alert_rules.yml
groups:
  - name: cryoprotect_alerts
    rules:
      - alert: HighCPUUsage
        expr: 100 - (avg by(instance) (irate(node_cpu_seconds_total{mode="idle"}[5m])) * 100) > 80
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "High CPU usage detected"
          description: "CPU usage is above 80% for more than 5 minutes on {{ $labels.instance }}"

      - alert: HighMemoryUsage
        expr: (node_memory_MemTotal_bytes - node_memory_MemAvailable_bytes) / node_memory_MemTotal_bytes * 100 > 85
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "High memory usage detected"
          description: "Memory usage is above 85% for more than 5 minutes on {{ $labels.instance }}"

      - alert: DiskSpaceRunningOut
        expr: (node_filesystem_size_bytes - node_filesystem_free_bytes) / node_filesystem_size_bytes * 100 > 85
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "Disk space running out"
          description: "Disk usage is above 85% for more than 5 minutes on {{ $labels.instance }}"

      - alert: HighAPILatency
        expr: http_request_duration_seconds{quantile="0.9"} > 1
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "High API latency detected"
          description: "90th percentile API response time is above 1 second for more than 5 minutes"

      - alert: HighErrorRate
        expr: sum(rate(http_requests_total{status=~"5.."}[5m])) / sum(rate(http_requests_total[5m])) * 100 > 5
        for: 5m
        labels:
          severity: critical
        annotations:
          summary: "High error rate detected"
          description: "Error rate is above 5% for more than 5 minutes"

      - alert: DatabaseHighConnections
        expr: pg_stat_activity_count{datname="cryoprotect"} > 20
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "High database connections"
          description: "Database has more than 20 connections for more than 5 minutes"

      - alert: SlowDatabaseQueries
        expr: pg_stat_activity_max_tx_duration{datname="cryoprotect"} > 30
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "Slow database queries detected"
          description: "Database has queries running for more than 30 seconds"
```

### 3. Alertmanager Configuration

```yaml
# alertmanager/alertmanager.yml
global:
  resolve_timeout: 5m
  slack_api_url: 'https://hooks.slack.com/services/TXXXXXXXX/BXXXXXXXX/XXXXXXXXXX'

route:
  group_by: ['alertname', 'job']
  group_wait: 30s
  group_interval: 5m
  repeat_interval: 4h
  receiver: 'slack-notifications'
  routes:
  - match:
      severity: critical
    receiver: 'pagerduty-critical'
    continue: true

receivers:
- name: 'slack-notifications'
  slack_configs:
  - channel: '#cryoprotect-alerts'
    send_resolved: true
    title: '{{ .GroupLabels.alertname }}'
    text: >
      {{ range .Alerts }}
        *Alert:* {{ .Annotations.summary }}
        *Description:* {{ .Annotations.description }}
        *Severity:* {{ .Labels.severity }}
        *Time:* {{ .StartsAt }}
      {{ end }}

- name: 'email-notifications'
  email_configs:
  - to: 'alerts@example.com'
    from: 'cryoprotect-alerts@example.com'
    smarthost: 'smtp.example.com:587'
    auth_username: 'alerts@example.com'
    auth_password: 'password'
    send_resolved: true

- name: 'pagerduty-critical'
  pagerduty_configs:
  - service_key: 'your_pagerduty_service_key'
    send_resolved: true
```

### 4. Docker Compose Configuration

```yaml
# docker-compose.monitoring.yml
version: '3.8'

services:
  prometheus:
    image: prom/prometheus:v2.37.0
    container_name: prometheus
    volumes:
      - ./prometheus:/etc/prometheus
      - prometheus_data:/prometheus
    command:
      - '--config.file=/etc/prometheus/prometheus.yml'
      - '--storage.tsdb.path=/prometheus'
      - '--web.console.libraries=/etc/prometheus/console_libraries'
      - '--web.console.templates=/etc/prometheus/consoles'
      - '--web.enable-lifecycle'
    ports:
      - "9090:9090"
    networks:
      - cryoprotect-network
    restart: unless-stopped

  node-exporter:
    image: prom/node-exporter:v1.3.1
    container_name: node-exporter
    volumes:
      - /proc:/host/proc:ro
      - /sys:/host/sys:ro
      - /:/rootfs:ro
    command:
      - '--path.procfs=/host/proc'
      - '--path.sysfs=/host/sys'
      - '--collector.filesystem.mount-points-exclude=^/(sys|proc|dev|host|etc)($$|/)'
    ports:
      - "9100:9100"
    networks:
      - cryoprotect-network
    restart: unless-stopped

  cadvisor:
    image: gcr.io/cadvisor/cadvisor:v0.45.0
    container_name: cadvisor
    volumes:
      - /:/rootfs:ro
      - /var/run:/var/run:ro
      - /sys:/sys:ro
      - /var/lib/docker/:/var/lib/docker:ro
      - /dev/disk/:/dev/disk:ro
    ports:
      - "8080:8080"
    networks:
      - cryoprotect-network
    restart: unless-stopped

  postgres-exporter:
    image: prometheuscommunity/postgres-exporter:v0.10.1
    container_name: postgres-exporter
    environment:
      DATA_SOURCE_NAME: "postgresql://postgres:password@postgres:5432/cryoprotect?sslmode=disable"
    ports:
      - "9187:9187"
    networks:
      - cryoprotect-network
    restart: unless-stopped

  grafana:
    image: grafana/grafana:9.1.0
    container_name: grafana
    volumes:
      - ./grafana/provisioning:/etc/grafana/provisioning
      - ./grafana/dashboards:/var/lib/grafana/dashboards
      - grafana_data:/var/lib/grafana
    environment:
      - GF_SECURITY_ADMIN_USER=admin
      - GF_SECURITY_ADMIN_PASSWORD=password
      - GF_USERS_ALLOW_SIGN_UP=false
    ports:
      - "3000:3000"
    networks:
      - cryoprotect-network
    restart: unless-stopped

  alertmanager:
    image: prom/alertmanager:v0.24.0
    container_name: alertmanager
    volumes:
      - ./alertmanager:/etc/alertmanager
      - alertmanager_data:/alertmanager
    command:
      - '--config.file=/etc/alertmanager/alertmanager.yml'
      - '--storage.path=/alertmanager'
    ports:
      - "9093:9093"
    networks:
      - cryoprotect-network
    restart: unless-stopped

volumes:
  prometheus_data:
  grafana_data:
  alertmanager_data:

networks:
  cryoprotect-network:
    external: true
```

### 5. Application Metrics Middleware

```python
# monitoring/middleware.py
import time
from functools import wraps
from flask import request, g
from prometheus_client import Counter, Histogram, Gauge

# Request count metric
REQUEST_COUNT = Counter(
    'http_requests_total',
    'Total number of HTTP requests',
    ['method', 'endpoint', 'status']
)

# Request latency metric
REQUEST_LATENCY = Histogram(
    'http_request_duration_seconds',
    'HTTP request latency in seconds',
    ['method', 'endpoint']
)

# Active requests gauge
ACTIVE_REQUESTS = Gauge(
    'http_requests_active',
    'Number of active HTTP requests',
    ['method', 'endpoint']
)

def timing_middleware():
    """
    Middleware to time requests and record metrics.
    """
    def decorator(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            # Record start time
            request_start_time = time.time()
            
            # Increment active requests
            endpoint = request.endpoint or 'unknown'
            ACTIVE_REQUESTS.labels(
                method=request.method,
                endpoint=endpoint
            ).inc()
            
            try:
                # Execute request handler
                response = f(*args, **kwargs)
                
                # Record metrics
                request_end_time = time.time()
                request_latency = request_end_time - request_start_time
                
                # Update Prometheus metrics
                REQUEST_COUNT.labels(
                    method=request.method,
                    endpoint=endpoint,
                    status=response.status_code
                ).inc()
                
                REQUEST_LATENCY.labels(
                    method=request.method,
                    endpoint=endpoint
                ).observe(request_latency)
                
                # Add timing header for debugging
                response.headers['X-Response-Time'] = f"{request_latency:.6f}s"
                
                return response
            finally:
                # Decrement active requests
                ACTIVE_REQUESTS.labels(
                    method=request.method,
                    endpoint=endpoint
                ).dec()
        
        return wrapped
    
    return decorator
```

### 6. Database Metrics Collector

```python
# monitoring/db_metrics.py
import time
import logging
import threading
from typing import Dict, Any, List
import psycopg2
from prometheus_client import Gauge, Counter

# Set up logging
logger = logging.getLogger(__name__)

# Database metrics
DB_CONNECTION_COUNT = Gauge(
    'database_connections',
    'Number of active database connections',
    ['database']
)

DB_TRANSACTION_COUNT = Counter(
    'database_transactions_total',
    'Total number of database transactions',
    ['database', 'status']
)

DB_QUERY_LATENCY = Gauge(
    'database_query_latency_seconds',
    'Database query latency in seconds',
    ['database', 'query_type']
)

DB_TABLE_SIZE = Gauge(
    'database_table_size_bytes',
    'Size of database tables in bytes',
    ['database', 'table']
)

DB_INDEX_SIZE = Gauge(
    'database_index_size_bytes',
    'Size of database indexes in bytes',
    ['database', 'index']
)

DB_TABLE_ROWS = Gauge(
    'database_table_rows',
    'Number of rows in database tables',
    ['database', 'table']
)

class DatabaseMetricsCollector:
    """Collects metrics from the database."""

    def __init__(self, connection_string: str, database: str = "cryoprotect", collection_interval: int = 60):
        """
        Initialize the database metrics collector.
        
        Args:
            connection_string: PostgreSQL connection string
            database: Database name
            collection_interval: How often to collect metrics (in seconds)
        """
        self.connection_string = connection_string
        self.database = database
        self.collection_interval = collection_interval
        self.running = False
        self.collection_thread = None
    
    def start(self):
        """Start collecting metrics."""
        if self.running:
            return
        
        self.running = True
        self.collection_thread = threading.Thread(target=self._collect_metrics_loop)
        self.collection_thread.daemon = True
        self.collection_thread.start()
        logger.info(f"Started database metrics collection for {self.database}")
    
    def stop(self):
        """Stop collecting metrics."""
        self.running = False
        if self.collection_thread:
            self.collection_thread.join(timeout=5)
        logger.info(f"Stopped database metrics collection for {self.database}")
    
    def _collect_metrics_loop(self):
        """Loop to continuously collect metrics."""
        while self.running:
            try:
                self.collect_metrics()
            except Exception as e:
                logger.error(f"Error collecting database metrics: {str(e)}")
            
            # Sleep until next collection
            time.sleep(self.collection_interval)
    
    def collect_metrics(self):
        """Collect all database metrics."""
        try:
            conn = psycopg2.connect(self.connection_string)
            
            # Collect connection metrics
            self._collect_connection_metrics(conn)
            
            # Collect table metrics
            self._collect_table_metrics(conn)
            
            # Collect index metrics
            self._collect_index_metrics(conn)
            
            # Collect query metrics
            self._collect_query_metrics(conn)
            
            conn.close()
        except Exception as e:
            logger.error(f"Error collecting metrics: {str(e)}")
            raise
    
    def _collect_connection_metrics(self, conn):
        """Collect connection metrics."""
        with conn.cursor() as cur:
            cur.execute("""
                SELECT count(*) FROM pg_stat_activity 
                WHERE datname = %s
            """, (self.database,))
            count = cur.fetchone()[0]
            DB_CONNECTION_COUNT.labels(database=self.database).set(count)
    
    def _collect_table_metrics(self, conn):
        """Collect table metrics."""
        with conn.cursor() as cur:
            # Table sizes
            cur.execute("""
                SELECT relname, pg_total_relation_size(c.oid) as total_size
                FROM pg_class c
                JOIN pg_namespace n ON n.oid = c.relnamespace
                WHERE relkind = 'r' 
                AND n.nspname = 'public'
            """)
            for table, size in cur.fetchall():
                DB_TABLE_SIZE.labels(database=self.database, table=table).set(size)
            
            # Table row counts
            cur.execute("""
                SELECT relname, reltuples::bigint as rows
                FROM pg_class c
                JOIN pg_namespace n ON n.oid = c.relnamespace
                WHERE relkind = 'r' 
                AND n.nspname = 'public'
            """)
            for table, rows in cur.fetchall():
                DB_TABLE_ROWS.labels(database=self.database, table=table).set(rows)
    
    def _collect_index_metrics(self, conn):
        """Collect index metrics."""
        with conn.cursor() as cur:
            # Index sizes
            cur.execute("""
                SELECT relname, pg_relation_size(c.oid) as index_size
                FROM pg_class c
                JOIN pg_namespace n ON n.oid = c.relnamespace
                WHERE relkind = 'i' 
                AND n.nspname = 'public'
            """)
            for index, size in cur.fetchall():
                DB_INDEX_SIZE.labels(database=self.database, index=index).set(size)
    
    def _collect_query_metrics(self, conn):
        """Collect query metrics."""
        with conn.cursor() as cur:
            # Check if pg_stat_statements is installed
            cur.execute("""
                SELECT EXISTS (
                    SELECT 1 FROM pg_extension WHERE extname = 'pg_stat_statements'
                ) as has_extension
            """)
            has_extension = cur.fetchone()[0]
            
            if has_extension:
                # Query latency from pg_stat_statements
                cur.execute("""
                    SELECT 
                        CASE 
                            WHEN query ~ '^SELECT' THEN 'SELECT'
                            WHEN query ~ '^INSERT' THEN 'INSERT'
                            WHEN query ~ '^UPDATE' THEN 'UPDATE'
                            WHEN query ~ '^DELETE' THEN 'DELETE'
                            ELSE 'OTHER'
                        END as query_type,
                        avg(mean_exec_time) / 1000 as avg_time_seconds
                    FROM pg_stat_statements
                    GROUP BY 1
                """)
                for query_type, avg_time in cur.fetchall():
                    DB_QUERY_LATENCY.labels(
                        database=self.database, 
                        query_type=query_type
                    ).set(avg_time)
```

## Testing Strategy

1. **Unit Testing**: Write unit tests for all monitoring components using pytest.
2. **Integration Testing**: Test the entire monitoring system in a staging environment.
3. **Load Testing**: Generate synthetic load to test alerting thresholds.
4. **End-to-end Testing**: Verify that alerts are properly sent to all configured channels.

## Documentation

Comprehensive documentation should be created for the monitoring system:

1. **System Architecture**: Document the overall monitoring architecture and components
2. **Dashboard Guide**: Guide to using and interpreting the Grafana dashboards
3. **Alert Reference**: Document all alert rules and their meaning
4. **Troubleshooting Guide**: Guide for diagnosing and fixing common issues
5. **Runbook**: Step-by-step procedures for common operational tasks

## Implementation Plan

### Week 1: Infrastructure Setup

1. Set up Prometheus, Grafana, and Alertmanager
2. Configure basic system monitoring
3. Set up Docker Compose for the monitoring stack
4. Create foundation of monitoring dashboard

### Week 2: Application Instrumentation

1. Implement application metrics middleware
2. Set up database metrics collection
3. Configure application-level metrics in Prometheus
4. Create detailed application dashboards

### Week 3: Alerting and Documentation

1. Configure comprehensive alert rules
2. Set up notification channels
3. Create alerting dashboards
4. Write monitoring system documentation

### Week 4: Testing and Refinement

1. Test monitoring system under load
2. Refine alert thresholds
3. Optimize metric collection
4. Complete documentation and runbooks

## Success Criteria

The performance monitoring implementation will be considered successful when:

1. All key metrics are collected and displayed on dashboards
2. Alerts are triggered when thresholds are exceeded
3. Notification channels deliver alerts correctly
4. The monitoring system has minimal performance impact
5. Documentation is complete and clear
6. Operations team can effectively use the monitoring tools

## Resources Required

1. Server resources for monitoring stack
2. Prometheus, Grafana, Alertmanager configuration
3. Developer time for implementation and testing
4. Operations team time for review and feedback
5. External services for notifications (email, Slack, PagerDuty)

## Next Steps

After implementing performance monitoring, the project should move on to implementing:

1. Scheduled backups
2. Maintenance runbooks
3. Security monitoring
4. Disaster recovery planning