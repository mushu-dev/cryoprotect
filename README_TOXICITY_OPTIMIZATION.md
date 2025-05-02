# Toxicity Query Optimization

This module provides tools to optimize database queries for toxicity data in the CryoProtect v2 system.

## Overview

The `optimize_toxicity_queries.py` script enhances database performance for toxicity-related queries by:

1. Adding performance indexes for commonly queried toxicity data
2. Creating materialized views for frequently accessed data
3. Verifying query performance with EXPLAIN ANALYZE

## Features

### Performance Indexes

The script creates the following indexes to improve query performance:

- `idx_toxicity_data_molecule_hit_call`: Optimizes queries filtering by molecule ID and hit call
- `idx_toxicity_data_assay_activity`: Optimizes queries filtering by assay ID and activity value
- `idx_toxicity_score_molecule_type`: Optimizes queries filtering by molecule ID and score type
- `idx_assay_endpoint_relevance`: Optimizes queries filtering by assay ID, endpoint ID, and relevance score
- `idx_toxicity_assay_endpoint`: Optimizes queries filtering by toxicological endpoint

### Materialized Views

The script creates the following materialized views to improve query performance:

- `toxicity_molecule_summary`: Provides a pre-computed summary of toxicity data for each molecule
- `toxicity_endpoint_summary`: Provides a pre-computed summary of toxicity data for each endpoint

### Query Performance Testing

The script tests the performance of common toxicity data queries using EXPLAIN ANALYZE to verify that the indexes and materialized views are working correctly.

## Usage

### Basic Usage

```bash
python optimize_toxicity_queries.py
```

This will:
1. Add all performance indexes
2. Create all materialized views
3. Test query performance

### Command Line Options

The script supports the following command line options:

- `--skip-indexes`: Skip creating indexes
- `--skip-views`: Skip creating materialized views
- `--skip-tests`: Skip query performance tests
- `--refresh-only`: Only refresh existing materialized views

### Examples

Refresh materialized views only:
```bash
python optimize_toxicity_queries.py --refresh-only
```

Add indexes only:
```bash
python optimize_toxicity_queries.py --skip-views --skip-tests
```

Create materialized views only:
```bash
python optimize_toxicity_queries.py --skip-indexes --skip-tests
```

## Testing

A test script is provided to verify the functionality of the optimization script:

```bash
python test_optimize_toxicity_queries.py
```

This runs unit tests to ensure that all functions in the optimization script work correctly.

## Integration with CryoProtect v2

This script is part of the toxicity data integration for CryoProtect v2. It works with the following components:

- Database schema in `migrations/012_toxicity_schema.sql`
- Toxicity data client in `chemical_data/toxicity/tox21_client.py`
- Toxicity scorer in `chemical_data/toxicity/toxicity_scorer.py`
- API resources in `api/toxicity_resources.py`

## Maintenance

The materialized views should be refreshed periodically to ensure they contain up-to-date data. This can be done by running:

```bash
python optimize_toxicity_queries.py --refresh-only
```

Consider scheduling this command to run automatically using a cron job or similar scheduling mechanism.