# Phase 2 Implementation Patterns and Examples

This document provides implementation patterns and code examples to guide the Roo PM agent in completing Phase 2 tasks. These patterns follow the established coding style and architecture of CryoProtect v2.

## Database Patterns

### Foreign Key Constraint Pattern

When adding foreign key constraints in SQL migrations, follow this pattern:

```sql
-- Check if constraint doesn't exist before adding
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.table_constraints 
        WHERE constraint_type = 'FOREIGN KEY' 
        AND table_name = 'table_name' 
        AND constraint_name = 'constraint_name'
    ) THEN
        ALTER TABLE table_name
        ADD CONSTRAINT constraint_name
        FOREIGN KEY (column_name) REFERENCES referenced_table(referenced_column) 
        ON DELETE CASCADE; -- or RESTRICT/SET NULL as appropriate
        
        RAISE NOTICE 'Added foreign key constraint: constraint_name';
    END IF;
END $$;
```

### RLS Policy Pattern

When implementing RLS policies, follow this pattern:

```sql
-- Enable RLS on table
ALTER TABLE table_name ENABLE ROW LEVEL SECURITY;

-- Select policy (who can view records)
CREATE POLICY table_name_select_policy ON table_name
    FOR SELECT
    USING (
        auth.uid() = created_by
        OR
        EXISTS (
            SELECT 1 FROM user_team_membership utm
            JOIN team_roles tr ON utm.role_id = tr.id
            WHERE utm.user_id = auth.uid() AND tr.can_view = true
        )
    );

-- Insert policy (who can create records)
CREATE POLICY table_name_insert_policy ON table_name
    FOR INSERT
    WITH CHECK (auth.uid() = created_by);

-- Update policy (who can modify records)
CREATE POLICY table_name_update_policy ON table_name
    FOR UPDATE
    USING (auth.uid() = created_by);

-- Delete policy (who can remove records)
CREATE POLICY table_name_delete_policy ON table_name
    FOR DELETE
    USING (auth.uid() = created_by);

-- Service role bypass policy (for admin operations)
CREATE POLICY table_name_service_role_policy ON table_name
    USING (auth.role() = 'service_role')
    WITH CHECK (auth.role() = 'service_role');
```

## Python Module Patterns

### Database Utility Module Pattern

For database utility modules like `health_check.py`, follow this pattern:

```python
"""
Module name and purpose.

Detailed description of the module's functionality.
"""

import logging
import json
from typing import Dict, List, Any, Optional
from datetime import datetime

from database.utils.connection import supabase_connection

# Set up logging
logger = logging.getLogger(__name__)

def main_function() -> Dict[str, Any]:
    """
    Main function description.
    
    Returns:
        Dict with results
    """
    results = {
        'status': 'passed',
        'details': {}
    }
    
    try:
        with supabase_connection() as conn:
            # Database operations
            pass
    except Exception as e:
        logger.error(f"Error in main_function: {str(e)}", exc_info=True)
        results['status'] = 'failed'
        results['details']['error'] = str(e)
    
    return results

# Additional functions as needed

if __name__ == "__main__":
    # Code to run when module is executed directly
    pass
```

### API Resource Module Pattern

For API resources like `predictive_models_resources.py`, follow this pattern:

```python
"""
API resource module name and purpose.

Detailed description of the API endpoints.
"""

import logging
from flask import request
from flask_restful import Resource, reqparse
from marshmallow import Schema, fields, validate, ValidationError

from api.utils import token_required, handle_error, marshal_with
from api.models import Model1, Model2

# Set up logging
logger = logging.getLogger(__name__)

# Define request/response schemas
class RequestSchema(Schema):
    field1 = fields.String(required=True)
    field2 = fields.Integer(required=False)

class ResponseSchema(Schema):
    id = fields.String()
    name = fields.String()
    created_at = fields.DateTime()

class ResourceName(Resource):
    """Resource for specific functionality."""
    
    @token_required
    def get(self, resource_id=None):
        """Get a resource or list of resources."""
        try:
            if resource_id:
                # Get specific resource
                pass
            else:
                # Get list of resources
                pass
        except Exception as e:
            return handle_error(e)
    
    @token_required
    def post(self):
        """Create a new resource."""
        try:
            # Validate request
            schema = RequestSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return {'error': err.messages}, 400
            
            # Create resource
            pass
        except Exception as e:
            return handle_error(e)
    
    # Additional methods as needed

def register_resources(api):
    """
    Register resources with the API.
    
    Args:
        api: Flask-RESTful API instance
    """
    api.add_resource(ResourceName, '/api/v1/resource-path', '/api/v1/resource-path/<string:resource_id>')
    # Register additional resources
```

## JavaScript Patterns

### Visualization Component Pattern

For visualization components like those in `predictive-models.js`, follow this pattern:

```javascript
/**
 * Create a specific chart type
 * 
 * @param {string} elementId - HTML element ID to render the chart
 * @param {Object} data - Data for the chart
 * @param {string} title - Chart title
 * @param {Object} options - Additional options
 * @returns {Object} Chart instance
 */
function createChart(elementId, data, title, options = {}) {
    const ctx = document.getElementById(elementId).getContext('2d');
    
    // Process and format data
    const processedData = processData(data);
    
    // Default options that can be overridden
    const defaultOptions = {
        responsive: true,
        maintainAspectRatio: false,
        title: {
            display: true,
            text: title || 'Chart Title'
        },
        // Additional default options
    };
    
    // Merge default options with provided options
    const chartOptions = Object.assign({}, defaultOptions, options);
    
    // Create the chart
    const chart = new Chart(ctx, {
        type: options.type || 'line',
        data: processedData,
        options: chartOptions
    });
    
    return chart;
}

/**
 * Process raw data for chart
 * 
 * @param {Object} rawData - Raw data from API
 * @returns {Object} Processed data for Chart.js
 */
function processData(rawData) {
    // Process data for chart
    return {
        labels: [],
        datasets: []
    };
}

// Export functions for external use
window.chartModule = {
    createChart,
    // Additional exported functions
};
```

### API Interaction Pattern

For API interactions in JavaScript, follow this pattern:

```javascript
/**
 * Function to interact with API
 * 
 * @param {string} resourceId - Optional resource ID
 * @param {Object} params - Optional query parameters
 * @returns {Promise} Promise that resolves with the data
 */
function fetchResourceFromApi(resourceId = null, params = {}) {
    // Build URL
    let url = '/api/v1/resource';
    if (resourceId) {
        url += `/${resourceId}`;
    }
    
    // Add query parameters
    if (Object.keys(params).length > 0) {
        const queryParams = new URLSearchParams();
        for (const [key, value] of Object.entries(params)) {
            queryParams.append(key, value);
        }
        url += `?${queryParams.toString()}`;
    }
    
    // Show loading indicator
    showLoading();
    
    // Make request
    return fetch(url)
        .then(response => {
            if (!response.ok) {
                throw new Error(`HTTP error! Status: ${response.status}`);
            }
            return response.json();
        })
        .then(data => {
            // Process data
            return data;
        })
        .catch(error => {
            console.error('Error fetching resource:', error);
            showError(`Error: ${error.message}`);
            throw error;
        })
        .finally(() => {
            // Hide loading indicator
            hideLoading();
        });
}
```

## CI/CD Patterns

### GitHub Actions Workflow Pattern

For GitHub Actions in `.github/workflows/ci-cd.yml`, follow this pattern:

```yaml
name: CI/CD Pipeline

on:
  push:
    branches: [ main, staging ]
    tags:
      - 'v*.*.*'
  pull_request:
    branches: [ main, staging ]
  workflow_dispatch:

jobs:
  test:
    name: Run Tests
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        test-group: [python, javascript, integration]
    
    steps:
    - name: Checkout code
      uses: actions/checkout@v3
    
    - name: Set up environment
      uses: actions/setup-python@v4
      with:
        python-version: '3.9'
    
    - name: Install dependencies
      run: |
        pip install -r requirements.txt
    
    - name: Run tests
      run: |
        python run_tests.py
    
    - name: Upload test results
      uses: actions/upload-artifact@v3
      with:
        name: test-results
        path: test-results/

  build:
    name: Build and Test Docker
    runs-on: ubuntu-latest
    needs: test
    
    steps:
    - name: Checkout code
      uses: actions/checkout@v3
    
    - name: Build Docker image
      run: |
        docker build -t app:test .
    
    - name: Test Docker image
      run: |
        docker run app:test python -m pytest
  
  deploy:
    name: Deploy to Environment
    runs-on: ubuntu-latest
    needs: build
    if: github.event_name == 'push' && (github.ref == 'refs/heads/main' || startsWith(github.ref, 'refs/tags/v'))
    
    steps:
    - name: Checkout code
      uses: actions/checkout@v3
    
    - name: Deploy
      run: |
        ./scripts/deploy.sh
```

## Monitoring Patterns

### Prometheus Metrics Pattern

For Prometheus metrics in `monitoring/prometheus_metrics.py`, follow this pattern:

```python
"""
Prometheus metrics for application monitoring.

This module defines Prometheus metrics for monitoring application performance.
"""

import time
from prometheus_client import Counter, Histogram, Gauge, Summary

# Request metrics
REQUEST_COUNT = Counter(
    'app_request_count_total',
    'Total request count',
    ['method', 'endpoint', 'status']
)

REQUEST_LATENCY = Histogram(
    'app_request_latency_seconds',
    'Request latency in seconds',
    ['method', 'endpoint'],
    buckets=[0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0]
)

# Database metrics
DB_QUERY_COUNT = Counter(
    'app_db_query_count_total',
    'Total database query count',
    ['query_type']
)

DB_QUERY_LATENCY = Histogram(
    'app_db_query_latency_seconds',
    'Database query latency in seconds',
    ['query_type'],
    buckets=[0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0]
)

# System metrics
MEMORY_USAGE = Gauge(
    'app_memory_usage_bytes',
    'Memory usage in bytes'
)

CPU_USAGE = Gauge(
    'app_cpu_usage_percent',
    'CPU usage in percent'
)

# Error metrics
ERROR_COUNT = Counter(
    'app_error_count_total',
    'Total error count',
    ['type', 'endpoint']
)

# Request tracking decorator
def track_request_latency(func):
    """Decorator to track request latency."""
    def wrapper(*args, **kwargs):
        method = request.method
        endpoint = request.endpoint or 'unknown'
        
        start_time = time.time()
        try:
            response = func(*args, **kwargs)
            status = response[1] if isinstance(response, tuple) and len(response) > 1 else 200
            REQUEST_COUNT.labels(method=method, endpoint=endpoint, status=status).inc()
            return response
        except Exception as e:
            status = 500
            REQUEST_COUNT.labels(method=method, endpoint=endpoint, status=status).inc()
            ERROR_COUNT.labels(type=type(e).__name__, endpoint=endpoint).inc()
            raise
        finally:
            duration = time.time() - start_time
            REQUEST_LATENCY.labels(method=method, endpoint=endpoint).observe(duration)
    
    return wrapper
```

### Alert Rules Pattern

For Prometheus alert rules, follow this pattern:

```yaml
groups:
- name: api_alerts
  rules:
  - alert: HighErrorRate
    expr: sum(rate(app_error_count_total[5m])) / sum(rate(app_request_count_total[5m])) > 0.05
    for: 2m
    labels:
      severity: critical
    annotations:
      summary: "High Error Rate"
      description: "Error rate is above 5% for the last 5 minutes"

  - alert: SlowEndpoint
    expr: avg(rate(app_request_latency_seconds_sum[5m]) / rate(app_request_latency_seconds_count[5m])) by (endpoint) > 1.0
    for: 5m
    labels:
      severity: warning
    annotations:
      summary: "Slow Endpoint: {{ $labels.endpoint }}"
      description: "Endpoint {{ $labels.endpoint }} has average latency > 1s for the last 5 minutes"

  - alert: ApiDown
    expr: up{job="api"} == 0
    for: 1m
    labels:
      severity: critical
    annotations:
      summary: "API Down"
      description: "API has been down for more than 1 minute"
```

## Advanced Implementation Techniques

### Batched Database Operations

For efficient database operations, use this pattern:

```python
def batch_insert(items: List[Dict[str, Any]], table_name: str, batch_size: int = 100) -> List[Dict[str, Any]]:
    """
    Insert items in batches for better performance.
    
    Args:
        items: List of items to insert
        table_name: Name of the table
        batch_size: Size of each batch
        
    Returns:
        List of inserted items with IDs
    """
    results = []
    
    for i in range(0, len(items), batch_size):
        batch = items[i:i+batch_size]
        
        with supabase_connection() as conn:
            response = conn.table(table_name).insert(batch).execute()
            
            if hasattr(response, 'data'):
                results.extend(response.data)
    
    return results
```

### Feature Flag Implementation

For feature flag implementation, use this pattern:

```python
class FeatureFlags:
    """Feature flag management."""
    
    @staticmethod
    def is_enabled(feature_name: str, user_id: str = None) -> bool:
        """
        Check if a feature is enabled.
        
        Args:
            feature_name: Name of the feature
            user_id: Optional user ID for user-specific features
            
        Returns:
            True if feature is enabled, False otherwise
        """
        try:
            with supabase_connection() as conn:
                query = conn.table('feature_flags').select('*').eq('name', feature_name)
                
                if user_id:
                    # Check for user-specific override
                    user_query = query.or_(f"user_ids.cs.{{'{user_id}'}}")
                    response = user_query.execute()
                else:
                    response = query.execute()
                
                if hasattr(response, 'data') and response.data:
                    return response.data[0].get('enabled', False)
                
                return False
        except Exception as e:
            logger.error(f"Error checking feature flag {feature_name}: {str(e)}")
            # Default to disabled if error
            return False
```

### Connection Pooling

For database connection pooling, use this pattern:

```python
from functools import wraps
import time

# Pooled connections cache
_connection_pool = {}
_max_pool_size = 10
_connection_timeout = 300  # 5 minutes

def get_pooled_connection():
    """
    Get a connection from the pool or create a new one.
    
    Returns:
        Database connection
    """
    current_time = time.time()
    
    # Clean up expired connections
    expired_keys = []
    for key, conn_data in _connection_pool.items():
        if current_time - conn_data['last_used'] > _connection_timeout:
            expired_keys.append(key)
    
    for key in expired_keys:
        del _connection_pool[key]
    
    # If pool has space, create new connection
    if len(_connection_pool) < _max_pool_size:
        conn = create_connection()
        conn_id = id(conn)
        _connection_pool[conn_id] = {
            'connection': conn,
            'last_used': current_time,
            'in_use': True
        }
        return conn
    
    # Find available connection
    for conn_id, conn_data in _connection_pool.items():
        if not conn_data['in_use']:
            conn_data['last_used'] = current_time
            conn_data['in_use'] = True
            return conn_data['connection']
    
    # All connections in use, wait or error
    raise ConnectionError("Connection pool exhausted")

def release_connection(conn):
    """
    Release a connection back to the pool.
    
    Args:
        conn: Connection to release
    """
    conn_id = id(conn)
    if conn_id in _connection_pool:
        _connection_pool[conn_id]['in_use'] = False
```

## Conclusion

These implementation patterns provide a comprehensive guide for the Roo PM agent to complete Phase 2 tasks. By following these patterns, the implementation will maintain consistency with the existing codebase while adding new functionality efficiently.