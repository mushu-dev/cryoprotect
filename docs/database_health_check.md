# Database Health Check System Documentation

## Overview

The CryoProtect v2 Database Health Check System provides comprehensive monitoring, validation, and diagnostics for the database layer. It enables administrators to proactively identify and resolve issues before they impact application performance or data integrity.

### Objectives

- **Ensure Schema Integrity**: Validate table structures, column definitions, constraints, and views
- **Verify Data Integrity**: Detect orphaned records, constraint violations, and data inconsistencies
- **Monitor Performance**: Identify slow queries, missing indexes, and optimization opportunities
- **Provide Actionable Insights**: Generate recommendations for resolving detected issues
- **Support Multiple Output Formats**: Deliver reports in JSON, Markdown, or HTML formats
- **Enable Automated Monitoring**: Schedule regular health checks with configurable intervals

### Architecture

The health check system follows a modular architecture with specialized validator components:

- **SchemaValidator**: Validates database schema integrity (tables, columns, constraints, views, RLS policies)
- **IntegrityChecker**: Verifies data integrity (orphaned records, unique constraints, required fields)
- **PerformanceAnalyzer**: Diagnoses performance issues (query times, indexes, connection pool)
- **DatabaseHealthCheck**: Orchestrates the validators and generates comprehensive reports

## Core Modules

### DatabaseHealthCheck Orchestrator

The `DatabaseHealthCheck` class serves as the main entry point for running health checks and generating reports. It integrates the specialized validator modules and provides a unified interface.

#### Basic Usage

```python
from database.utils.health_check import DatabaseHealthCheck
from database.utils.connection import supabase_connection

# Create a connection to the database
with supabase_connection() as conn:
    # Initialize the health check orchestrator
    health_check = DatabaseHealthCheck(conn)
    
    # Run a comprehensive health check
    results = health_check.run_health_check()
    
    # Generate a report in the desired format
    report = health_check.generate_report(results, format='markdown')
    
    # Print or save the report
    print(report)
```

#### Running Specific Checks

You can run checks for specific categories:

```python
# Run only schema validation checks
schema_results = health_check.run_health_check(categories=['schema'])

# Run integrity and performance checks
partial_results = health_check.run_health_check(categories=['integrity', 'performance'])
```

### SchemaValidator

The `SchemaValidator` verifies the database schema structure against expected definitions.

#### Checks Performed

- Table existence verification
- Column definition and type checking
- Constraint validation (primary keys, foreign keys, unique)
- View definition verification
- RLS policy verification

### IntegrityChecker

The `IntegrityChecker` validates data integrity across the database.

#### Checks Performed

- Orphaned records detection (foreign key constraint violations)
- Unique constraint validation
- Data consistency across related tables
- Required field validation
- Data corruption detection

### PerformanceAnalyzer

The `PerformanceAnalyzer` diagnoses database performance issues.

#### Checks Performed

- Query execution time measurement
- Index usage and effectiveness analysis
- Slow query identification
- Connection pool utilization monitoring
- Table statistics and growth pattern analysis

## API Endpoints

The health check system exposes several REST API endpoints for integration with monitoring tools and dashboards.

### `/health` Endpoint

A public endpoint that provides basic system health status without authentication.

#### Request

```
GET /health
```

#### Response

```json
{
  "status": "ok",
  "version": "1.0.0",
  "uptime": 86400,
  "uptime_formatted": "1 day, 0:00:00",
  "database_status": "connected",
  "timestamp": "2025-04-23T17:30:00.000Z"
}
```

#### Status Values

- `ok`: System is fully operational
- `degraded`: System is operational but with reduced functionality
- `error`: System is experiencing critical issues

### `/health/database` Endpoint

An authenticated endpoint that provides detailed database health information.

#### Request

```
GET /health/database
```

#### Query Parameters

- `verbosity`: Level of detail in the response (`minimal`, `normal`, `detailed`; default: `normal`)
- `format`: Response format (`json`, `markdown`, `html`; default: `json`)

#### Response (JSON format)

```json
{
  "status": "warning",
  "overall_health": "warning",
  "schema_status": "passed",
  "integrity_status": "warning",
  "performance_status": "not_checked",
  "issues_count": 3,
  "timestamp": "2025-04-23T17:30:00.000Z",
  "details": {
    "categories_checked": ["schema", "integrity"],
    "recommendations": [
      "Fix orphaned records in 'mixture_components' that reference non-existent 'molecules' records",
      "Create index on mixture_components(molecule_id) to improve join performance"
    ]
  }
}
```

#### Authentication

This endpoint requires:
- A valid authentication token
- Admin privileges

### `/health/performance` Endpoint

An authenticated endpoint that provides detailed performance metrics.

#### Request

```
GET /health/performance
```

#### Query Parameters

- `period`: Time period for metrics (`minute`, `hour`, `day`, `week`; default: `hour`)
- `metrics`: Comma-separated list of metric categories to include (`cpu`, `memory`, `database`, `api`; default: all)

#### Response

```json
{
  "status": "ok",
  "metrics": {
    "cpu": {
      "current": 45.2,
      "history": {
        "labels": ["0", "1", "2", "..."],
        "values": [42.5, 43.1, 44.8, "..."],
        "unit": "%",
        "interval": "minute"
      }
    },
    "memory": {
      "current": {
        "total": 16777216000,
        "available": 8388608000,
        "used": 8388608000,
        "percent": 50.0
      },
      "history": {
        "labels": ["0", "1", "2", "..."],
        "values": [48.5, 49.2, 50.1, "..."],
        "unit": "%",
        "interval": "minute"
      }
    }
  },
  "database_metrics": {
    "status": "passed",
    "query_performance": {
      "get_molecules_with_properties": {
        "execution_time_ms": 45.2,
        "threshold_ms": 100,
        "status": "ok"
      }
    },
    "connection_stats": {
      "total_connections": 10,
      "active_connections": 3,
      "idle_connections": 7
    }
  },
  "api_metrics": {
    "requests": {
      "labels": ["0", "1", "2", "..."],
      "values": [120, 125, 118, "..."],
      "unit": "count",
      "interval": "minute"
    },
    "response_time": {
      "labels": ["0", "1", "2", "..."],
      "values": [210.5, 215.2, 208.7, "..."],
      "unit": "ms",
      "interval": "minute"
    },
    "error_rate": {
      "labels": ["0", "1", "2", "..."],
      "values": [1.2, 1.5, 1.1, "..."],
      "unit": "%",
      "interval": "minute"
    }
  },
  "timestamp": "2025-04-23T17:30:00.000Z"
}
```

#### Authentication

This endpoint requires:
- A valid authentication token
- Admin privileges

## Configuration Options

### Health Check Configuration

The `DatabaseHealthCheck` constructor accepts an optional configuration dictionary:

```python
config = {
    'report_dir': '/custom/path/to/reports',
    'thresholds': {
        'query_execution': {
            'warning_ms': 100,
            'critical_ms': 500
        },
        'connection_pool': {
            'max_idle_transactions': 5,
            'max_connections_percent': 80
        }
    },
    'scheduled_checks': {
        'enabled': True,
        'interval_hours': 24,
        'alert_on_failure': True,
        'alert_email': 'admin@example.com'
    }
}

health_check = DatabaseHealthCheck(conn, config)
```

### Validator Configuration

Each validator module can be configured with custom validation rules:

#### SchemaValidator

```python
# Define custom expected tables
schema_validator = SchemaValidator()
schema_validator.expected_tables = {
    'custom_table': [
        {'name': 'id', 'type': 'uuid', 'nullable': False, 'primary_key': True},
        {'name': 'name', 'type': 'character varying', 'nullable': False}
    ]
}
```

#### IntegrityChecker

```python
# Define custom table relationships
integrity_checker = IntegrityChecker()
integrity_checker.table_relationships = {
    'parent_table': {
        'dependent_tables': [
            {'table': 'child_table', 'fk_column': 'parent_id'}
        ]
    }
}
```

#### PerformanceAnalyzer

```python
# Define custom test queries
performance_analyzer = PerformanceAnalyzer()
performance_analyzer.test_queries = [
    {
        'name': 'custom_query',
        'query': 'SELECT * FROM custom_table LIMIT 10;',
        'threshold_ms': 50
    }
]
```

## Integration Guidance

### Integrating with Monitoring Systems

The health check system can be integrated with external monitoring systems:

```python
# Run health check and send results to monitoring system
def send_health_metrics_to_monitoring():
    with supabase_connection() as conn:
        health_check = DatabaseHealthCheck(conn)
        results = health_check.run_health_check()
        
        # Extract metrics
        metrics = {
            'overall_status': results['overall_status'],
            'issues_count': results['total_issues_found'],
            'schema_status': results['results_by_category']['schema']['status'],
            'integrity_status': results['results_by_category']['integrity']['status'],
            'performance_status': results['results_by_category']['performance']['status']
        }
        
        # Send to monitoring system (example)
        monitoring_client.send_metrics('database_health', metrics)
```

### Scheduled Health Checks

The system supports scheduled health checks:

```python
from database.utils.health_check import run_scheduled_health_check

# Run health check every 12 hours
run_scheduled_health_check(interval_hours=12)
```

### Alert Integration

Configure alerts for critical issues:

```python
def health_check_with_alerts():
    with supabase_connection() as conn:
        health_check = DatabaseHealthCheck(conn)
        results = health_check.run_health_check()
        
        if results['overall_status'] == 'failed':
            # Generate alert message
            alert_message = f"Database health check failed with {results['total_issues_found']} issues"
            
            # Send alert (example)
            send_alert('DATABASE_HEALTH_CRITICAL', alert_message, results)
```

## Output Formats

The health check system supports multiple output formats for reports.

### JSON Format

The JSON format provides structured data suitable for programmatic processing:

```python
# Generate JSON report
json_report = health_check.generate_report(results, format='json')
```

Example JSON output:
```json
{
  "timestamp": "2025-04-23T17:30:00.000Z",
  "overall_status": "warning",
  "categories_checked": ["schema", "integrity", "performance"],
  "total_issues_found": 3,
  "results_by_category": {
    "schema": {
      "status": "passed",
      "issues_found": 0,
      "details": {
        "tables": {"status": "passed", "issues_found": 0},
        "columns": {"status": "passed", "issues_found": 0},
        "constraints": {"status": "passed", "issues_found": 0},
        "views": {"status": "passed", "issues_found": 0},
        "rls_policies": {"status": "passed", "issues_found": 0}
      }
    },
    "integrity": {
      "status": "warning",
      "issues_found": 2,
      "details": {
        "orphaned_records": {
          "status": "warning",
          "issues_found": 2,
          "orphaned_records": {
            "mixture_components_to_molecules": [
              {"id": "123e4567-e89b-12d3-a456-426614174000", "molecule_id": "123e4567-e89b-12d3-a456-426614174111"}
            ]
          }
        }
      },
      "recommendations": [
        "Fix orphaned records in 'mixture_components' that reference non-existent 'molecules' records"
      ]
    },
    "performance": {
      "status": "warning",
      "issues_found": 1,
      "details": {
        "missing_indexes": [
          {"table_name": "mixture_components", "column_name": "molecule_id"}
        ]
      },
      "recommendations": [
        "Create index on mixture_components(molecule_id) to improve join performance"
      ]
    }
  },
  "recommendations": [
    "Fix orphaned records in 'mixture_components' that reference non-existent 'molecules' records",
    "Create index on mixture_components(molecule_id) to improve join performance"
  ],
  "execution_time_seconds": 1.25
}
```

### Markdown Format

The Markdown format provides human-readable reports suitable for documentation:

```python
# Generate Markdown report
markdown_report = health_check.generate_report(results, format='markdown')
```

Example Markdown output:
```markdown
# Database Health Check Report

**Generated:** 2025-04-23T17:30:00.000Z
**Overall Status:** WARNING
**Total Issues Found:** 3
**Categories Checked:** schema, integrity, performance
**Execution Time:** 1.25 seconds

## Summary by Category
- **Schema**: PASSED (0 issues)
- **Integrity**: WARNING (2 issues)
- **Performance**: WARNING (1 issues)

## Recommendations
1. Fix orphaned records in 'mixture_components' that reference non-existent 'molecules' records
2. Create index on mixture_components(molecule_id) to improve join performance

## Detailed Results

### Schema Checks
**Status:** PASSED
**Issues Found:** 0

### Integrity Checks
**Status:** WARNING
**Issues Found:** 2

#### Details

##### Orphaned Records
Status: warning
- Orphaned Records: 
  - mixture_components_to_molecules: 1 record(s)

### Performance Checks
**Status:** WARNING
**Issues Found:** 1

#### Details

##### Missing Indexes
Status: warning
- Missing Indexes: 
  - mixture_components.molecule_id
```

### HTML Format

The HTML format provides visually formatted reports suitable for web display:

```python
# Generate HTML report
html_report = health_check.generate_report(results, format='html')
```

The HTML report includes the same information as the Markdown report but with additional styling for better readability in web browsers.

## Verbosity Levels

The health check system supports different verbosity levels for API responses:

### Minimal Verbosity

Only includes schema validation checks:

```
GET /health/database?verbosity=minimal
```

### Normal Verbosity (Default)

Includes schema and integrity checks:

```
GET /health/database?verbosity=normal
```

### Detailed Verbosity

Includes all checks (schema, integrity, performance) with full details:

```
GET /health/database?verbosity=detailed
```

## Example Health Reports

### Healthy System Report

```json
{
  "timestamp": "2025-04-23T17:30:00.000Z",
  "overall_status": "passed",
  "categories_checked": ["schema", "integrity", "performance"],
  "total_issues_found": 0,
  "results_by_category": {
    "schema": {"status": "passed", "issues_found": 0},
    "integrity": {"status": "passed", "issues_found": 0},
    "performance": {"status": "passed", "issues_found": 0}
  },
  "recommendations": [],
  "execution_time_seconds": 0.95
}
```

### Warning Report

```json
{
  "timestamp": "2025-04-23T17:30:00.000Z",
  "overall_status": "warning",
  "categories_checked": ["schema", "integrity", "performance"],
  "total_issues_found": 2,
  "results_by_category": {
    "schema": {"status": "passed", "issues_found": 0},
    "integrity": {"status": "passed", "issues_found": 0},
    "performance": {
      "status": "warning",
      "issues_found": 2,
      "details": {
        "query_execution_times": {
          "get_mixtures_with_components": {
            "execution_time_ms": 120.5,
            "threshold_ms": 100,
            "status": "slow"
          }
        },
        "missing_indexes": [
          {"table_name": "mixture_components", "column_name": "molecule_id"}
        ]
      }
    }
  },
  "recommendations": [
    "Optimize query 'get_mixtures_with_components' which took 120.50ms (threshold: 100ms)",
    "Create index on mixture_components(molecule_id) to improve join performance"
  ],
  "execution_time_seconds": 1.05
}
```

### Critical Report

```json
{
  "timestamp": "2025-04-23T17:30:00.000Z",
  "overall_status": "failed",
  "categories_checked": ["schema", "integrity", "performance"],
  "total_issues_found": 5,
  "results_by_category": {
    "schema": {
      "status": "failed",
      "issues_found": 2,
      "details": {
        "missing_tables": ["molecular_properties"],
        "column_issues": {
          "molecules": {
            "missing_columns": ["smiles"]
          }
        }
      }
    },
    "integrity": {
      "status": "failed",
      "issues_found": 3,
      "details": {
        "orphaned_records": {
          "mixture_components_to_molecules": [
            {"id": "123e4567-e89b-12d3-a456-426614174000", "molecule_id": "123e4567-e89b-12d3-a456-426614174111"}
          ]
        },
        "missing_values": {
          "molecules_name": [
            {"id": "123e4567-e89b-12d3-a456-426614174222"}
          ]
        }
      }
    },
    "performance": {"status": "passed", "issues_found": 0}
  },
  "recommendations": [
    "Create missing table 'molecular_properties'",
    "Add missing column 'smiles' to table 'molecules'",
    "Fix orphaned records in 'mixture_components' that reference non-existent 'molecules' records",
    "Update record with ID '123e4567-e89b-12d3-a456-426614174222' in table 'molecules' to provide a value for 'name'"
  ],
  "execution_time_seconds": 1.15
}
```

## Troubleshooting

### Common Issues and Solutions

#### Schema Validation Failures

**Issue**: Missing tables or columns
**Solution**: Run database migrations to create the missing schema objects:

```sql
-- Example: Create missing table
CREATE TABLE IF NOT EXISTS molecular_properties (
    id UUID PRIMARY KEY,
    molecule_id UUID NOT NULL REFERENCES molecules(id),
    property_name VARCHAR NOT NULL,
    numeric_value DOUBLE PRECISION,
    text_value TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Example: Add missing column
ALTER TABLE molecules ADD COLUMN IF NOT EXISTS smiles VARCHAR;
```

#### Data Integrity Issues

**Issue**: Orphaned records
**Solution**: Either delete the orphaned records or create the missing parent records:

```sql
-- Option 1: Delete orphaned records
DELETE FROM mixture_components
WHERE molecule_id NOT IN (SELECT id FROM molecules);

-- Option 2: Create missing parent records
INSERT INTO molecules (id, name)
SELECT DISTINCT mc.molecule_id, 'Unknown Molecule ' || mc.molecule_id
FROM mixture_components mc
LEFT JOIN molecules m ON mc.molecule_id = m.id
WHERE m.id IS NULL;
```

#### Performance Issues

**Issue**: Missing indexes
**Solution**: Create indexes on frequently queried columns:

```sql
-- Create index on foreign key column
CREATE INDEX IF NOT EXISTS idx_mixture_components_molecule_id
ON mixture_components(molecule_id);
```

**Issue**: Slow queries
**Solution**: Analyze and optimize query execution plans:

```sql
-- Analyze query plan
EXPLAIN ANALYZE
SELECT m.id, m.name, mc.molecule_id, mc.concentration
FROM mixtures m
JOIN mixture_components mc ON m.id = mc.mixture_id
LIMIT 10;

-- Optimize by adding indexes or rewriting the query
```

### Diagnostic Commands

Use these SQL commands to diagnose specific issues:

```sql
-- Check for orphaned records
SELECT mc.id, mc.mixture_id, mc.molecule_id
FROM mixture_components mc
LEFT JOIN molecules m ON mc.molecule_id = m.id
WHERE m.id IS NULL;

-- Check for missing required values
SELECT id, name
FROM molecules
WHERE name IS NULL OR name = '';

-- Check index usage
SELECT
    relname AS table_name,
    indexrelname AS index_name,
    idx_scan AS index_scans,
    idx_tup_read AS tuples_read,
    idx_tup_fetch AS tuples_fetched
FROM pg_stat_user_indexes
JOIN pg_statio_user_indexes USING (indexrelid)
ORDER BY idx_scan DESC;
```

## Best Practices

### Scheduled Health Checks

- Run comprehensive health checks during off-peak hours
- Schedule more frequent lightweight checks during business hours
- Configure alerts for critical issues

```python
# Example: Configure scheduled health checks
run_scheduled_health_check(interval_hours=24)  # Daily comprehensive check
```

### Proactive Monitoring

- Monitor trends in database size and performance metrics
- Set up alerts for sudden changes in error rates or response times
- Regularly review health check reports even when no critical issues are reported

### Regular Maintenance

- Apply recommended fixes promptly
- Run VACUUM and ANALYZE regularly to maintain statistics
- Implement index maintenance procedures

```sql
-- Regular maintenance commands
VACUUM ANALYZE;
REINDEX DATABASE cryoprotect;
```

### Integration with CI/CD

- Run health checks as part of deployment verification
- Include database schema validation in CI/CD pipelines
- Verify data integrity after major migrations

```python
# Example: Health check in deployment script
def verify_deployment():
    with supabase_connection() as conn:
        health_check = DatabaseHealthCheck(conn)
        results = health_check.run_health_check()
        
        if results['overall_status'] != 'passed':
            print("Deployment verification failed!")
            print(f"Issues found: {results['total_issues_found']}")
            for rec in results['recommendations']:
                print(f"- {rec}")
            return False
        
        return True
```

## Remediation Recommendations

### Schema Issues

- Always use migrations for schema changes
- Test schema changes in development/staging environments first
- Document schema changes and their impact

### Data Integrity Issues

- Implement application-level validation to prevent integrity issues
- Use database constraints (foreign keys, unique, not null) to enforce integrity
- Regularly audit data for consistency

### Performance Issues

- Create indexes for frequently queried columns and join conditions
- Optimize slow queries by rewriting or adding appropriate indexes
- Monitor query execution plans for changes after database updates

## Conclusion

The Database Health Check System provides comprehensive monitoring and diagnostics for the CryoProtect v2 database. By regularly running health checks and addressing identified issues, you can ensure optimal database performance, data integrity, and application reliability.

For additional assistance or to report issues with the health check system itself, please contact the database administration team.