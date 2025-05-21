# Monitoring Module Reference

## Overview

The monitoring module provides comprehensive observability features for the CryoProtect application, integrating metrics collection, logging, alerting, and performance tracking capabilities.

## Key Components

### Middleware
- **Request Tracking**: `middleware.py` for HTTP request monitoring
- **Performance Metrics**: Response time and throughput collection
- **Error Tracking**: Exception monitoring and categorization
- **Resource Utilization**: Memory and CPU usage tracking

### Prometheus Integration
- **Metrics Export**: `prometheus_metrics.py` for exporting metrics
- **Custom Counters**: Track business-specific operations
- **Histograms**: Measure distribution of response times
- **Gauges**: Monitor current state values (e.g., connection pool size)
- **Endpoint**: `/metrics` for Prometheus scraping

### Configuration
- **Metric Naming**: Standardized naming convention
- **Label Structure**: Consistent labeling approach
- **Aggregation Rules**: How metrics are combined
- **Retention Policies**: How long metrics are stored
- **Scrape Intervals**: How frequently metrics are collected

## Monitoring Categories

### Application Health
- **API Health**: Endpoint availability and response times
- **Database Health**: Connection pool status and query times
- **External Services**: Integration point health (ChEMBL, PubChem)
- **Authentication**: JWT verification success rates

### Business Metrics
- **User Activity**: Active users and sessions
- **Operation Counts**: Molecular analysis operations
- **Data Volume**: Size of datasets being processed
- **Calculation Distribution**: Types of properties being calculated

### System Resources
- **Memory Usage**: Heap and non-heap memory consumption
- **CPU Utilization**: Process and system CPU usage
- **Disk Activity**: Storage read/write patterns
- **Network Traffic**: Inbound and outbound data rates

### Error Tracking
- **Exception Rates**: Count of exceptions by type
- **Status Codes**: Distribution of HTTP status codes
- **Query Failures**: Database operation errors
- **Integration Failures**: External API communication issues

## Alert Configuration

Alerts are configured for various thresholds:

- **High Latency**: Response times exceed acceptable limits
- **Error Rate Spike**: Sudden increase in application errors
- **Resource Exhaustion**: Memory/CPU approaching limits
- **Integration Failure**: External API dependency issues
- **Security Events**: Authentication or authorization failures

## Dashboards

Predefined Grafana dashboards are available for:

- **API Overview**: General API performance and health
- **Database Performance**: Query performance and connection stats
- **User Activity**: User engagement and operation patterns
- **System Resources**: Infrastructure utilization
- **Integration Status**: External service health and performance

## Best Practices

1. **Use Standardized Metrics**: Follow the naming conventions
2. **Add Context Labels**: Include relevant dimensions for filtering
3. **Focus on Actionable Metrics**: Measure what matters for operations
4. **Set Appropriate Thresholds**: Avoid alert fatigue
5. **Document Dashboard Usage**: Help others understand metrics
6. **Monitor Business Outcomes**: Track metrics that reflect user success

## Integration Points

- **Prometheus**: For metrics collection and storage
- **Grafana**: For visualization and dashboarding
- **AlertManager**: For alert routing and notification
- **ELK Stack**: For log aggregation and analysis
- **Jaeger/OpenTelemetry**: For distributed tracing

## Usage Examples

```python
# Instrumenting code with metrics
from monitoring.prometheus_metrics import molecule_counter, query_latency

# Count molecule processing
molecule_counter.inc()

# Measure query time
with query_latency.time():
    results = db.execute_query("SELECT * FROM molecules")

# Custom metrics with labels
from monitoring.prometheus_metrics import calculation_counter

calculation_counter.labels(
    property_type="logP",
    calculation_method="rdkit"
).inc()

# Middleware application
from monitoring.middleware import MetricsMiddleware

app.wsgi_app = MetricsMiddleware(app.wsgi_app)
```