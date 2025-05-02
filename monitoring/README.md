# CryoProtect v2 - Monitoring System

This directory contains the monitoring and metrics collection system for the CryoProtect API. It provides comprehensive monitoring capabilities using Prometheus and Grafana.

## Features

- Prometheus metrics collection with custom application metrics
- Exposed `/metrics` endpoint for Prometheus scraping
- Grafana dashboards for:
  - API latency and throughput
  - Error rates and status codes
  - System resources (CPU, memory, disk)
  - Database performance
- Alerting rules for critical conditions (e.g., high error rates, resource exhaustion)
- Unit tests for all monitoring components

## Directory Structure

```
monitoring/
├── __init__.py                      # Package initialization
├── prometheus_metrics.py            # Prometheus metrics definitions
├── middleware.py                    # Flask middleware for metrics collection
├── docker-compose.yml               # Docker Compose configuration for the monitoring stack
├── prometheus/                      # Prometheus configuration
│   ├── prometheus.yml               # Main Prometheus configuration
│   └── rules/                       # Alerting rules
│       ├── api_alerts.yml           # API-related alerting rules
│       └── system_alerts.yml        # System-related alerting rules
├── grafana/                         # Grafana configuration
│   ├── dashboards/                  # Grafana dashboards
│   │   ├── api_performance.json     # API performance dashboard
│   │   ├── system_resources.json    # System resources dashboard
│   │   ├── database_performance.json # Database performance dashboard
│   │   └── provisioning/            # Dashboard provisioning configuration
│   │       └── dashboards.yml       # Dashboard provisioning configuration
│   └── datasources/                 # Grafana datasources
│       └── prometheus.yml           # Prometheus datasource configuration
├── alertmanager/                    # Alertmanager configuration
│   └── alertmanager.yml             # Alertmanager configuration
└── tests/                           # Unit tests
    ├── __init__.py                  # Test package initialization
    ├── test_prometheus_metrics.py   # Tests for prometheus_metrics.py
    └── test_middleware.py           # Tests for middleware.py
```

## Setup and Usage

### Prerequisites

- Docker and Docker Compose
- Python 3.8 or higher
- Flask application (CryoProtect API)

### Installation

1. Install the required Python packages:

```bash
pip install prometheus-client flask psutil
```

2. Integrate the monitoring system with your Flask application:

```python
from flask import Flask
from monitoring import init_metrics, PrometheusMiddleware

app = Flask(__name__)

# Initialize Prometheus metrics
init_metrics(app)

# Add Prometheus middleware
prometheus_middleware = PrometheusMiddleware(app)
```

### Starting the Monitoring Stack

1. Start the monitoring stack using Docker Compose:

```bash
cd monitoring
docker-compose up -d
```

2. Access the Grafana dashboard at http://localhost:3000 (default credentials: admin/cryoprotect)

3. Access the Prometheus dashboard at http://localhost:9090

4. Access the Alertmanager dashboard at http://localhost:9093

### Running Tests

Run the unit tests for the monitoring system:

```bash
python -m unittest discover -s monitoring/tests
```

## Metrics Collected

### API Metrics

- Request count by endpoint, method, and status code
- Request latency by endpoint and method
- Error count by type and endpoint
- In-progress requests by endpoint and method

### Database Metrics

- Query count by operation and table
- Query latency by operation and table
- Connection pool size by state (active, idle, max)

### System Metrics

- CPU usage
- Memory usage (total, available, used)
- Disk usage (total, free, used)

### Business Metrics

- Molecule count
- Mixture count
- Prediction count
- Experiment count

## Alerting Rules

### API Alerts

- High error rate (>5% for 2 minutes)
- High latency (95th percentile >1s for 2 minutes)
- Endpoint down (>80% error rate for 2 minutes)
- High request rate (>100 req/s for 2 minutes)
- High database latency (95th percentile >0.5s for 2 minutes)
- Database errors (>10 in 5 minutes)
- Too many in-progress requests (>50 for 2 minutes)

### System Alerts

- High CPU usage (>80% for 5 minutes)
- Critical CPU usage (>95% for 2 minutes)
- High memory usage (>80% for 5 minutes)
- Critical memory usage (>95% for 2 minutes)
- High disk usage (>80% for 5 minutes)
- Critical disk usage (>95% for 2 minutes)
- Connection pool exhaustion (>90% for 2 minutes)
- Node down (for 1 minute)
- API service down (for 1 minute)

## Customization

### Adding Custom Metrics

To add custom metrics, modify the `prometheus_metrics.py` file:

```python
# Add a new counter
MY_COUNTER = Counter(
    'cryoprotect_my_counter_total',
    'Description of my counter',
    ['label1', 'label2']
)

# Add a new gauge
MY_GAUGE = Gauge(
    'cryoprotect_my_gauge',
    'Description of my gauge',
    ['label1', 'label2']
)

# Add a new histogram
MY_HISTOGRAM = Histogram(
    'cryoprotect_my_histogram_seconds',
    'Description of my histogram',
    ['label1', 'label2'],
    buckets=(0.01, 0.1, 1.0, 10.0)
)
```

### Adding Custom Alerting Rules

To add custom alerting rules, create a new file in the `prometheus/rules/` directory or modify the existing files.

### Adding Custom Dashboards

To add custom dashboards, create a new JSON file in the `grafana/dashboards/` directory and update the provisioning configuration if needed.