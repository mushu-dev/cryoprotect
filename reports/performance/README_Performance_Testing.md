# CryoProtect v2 Database Performance Testing

This document provides instructions for running database performance tests on the CryoProtect v2 application and interpreting the results.

## Overview

The performance testing suite consists of several scripts designed to evaluate different aspects of database performance:

1. **test_database_performance.py**: Tests query performance, response times, and resource utilization under various load conditions.
2. **analyze_query_plans.py**: Analyzes the execution plans of common database queries to identify optimization opportunities.
3. **simulate_concurrent_load.py**: Simulates concurrent user load to test how the database performs under stress.
4. **run_all_performance_tests.py**: Runs all the above tests and generates a comprehensive report.

## Prerequisites

Before running the performance tests, ensure you have:

1. Python 3.7+ installed
2. Required Python packages:
   ```
   pip install supabase psutil python-dotenv
   ```
3. A valid `.env` file with Supabase credentials:
   ```
   SUPABASE_URL=your_supabase_url
   SUPABASE_KEY=your_supabase_service_role_key
   ```
4. The database has been remediated and populated with test data

## Running the Tests

### Option 1: Run All Tests

To run all performance tests and generate a comprehensive report:

```bash
python run_all_performance_tests.py
```

This will execute all three test scripts sequentially and combine the results into a comprehensive report.

### Option 2: Run Individual Tests

You can also run each test script individually:

```bash
# Test database performance
python test_database_performance.py

# Analyze query plans
python analyze_query_plans.py

# Simulate concurrent load
python simulate_concurrent_load.py
```

## Test Configuration

Each test script has configuration parameters that can be adjusted:

### test_database_performance.py

```python
TEST_CONFIG = {
    "concurrent_users": [1, 5, 10, 20],  # Number of concurrent users to simulate
    "iterations_per_user": 5,            # Number of iterations per user
    "read_operations": [...],            # List of read operations to test
    "write_operations": [...]            # List of write operations to test
}
```

### analyze_query_plans.py

```python
QUERIES = [
    {
        "name": "Get all molecules with properties",
        "query": "SELECT * FROM molecule_with_properties LIMIT 100"
    },
    # Additional queries...
]
```

### simulate_concurrent_load.py

```python
TEST_CONFIG = {
    "duration_seconds": 60,           # Duration of the test in seconds
    "ramp_up_seconds": 10,            # Time to ramp up to full user load
    "concurrent_users": 20,           # Maximum number of concurrent users
    "user_think_time_ms": [500, 3000],  # Random think time between operations
    "operation_weights": {            # Relative frequency of operations
        "browse_molecules": 30,
        # Additional operations...
    }
}
```

## Understanding the Results

### Performance Test Report

The performance test report (`database_performance_report.txt`) includes:

- **Resource Usage**: CPU and memory utilization during the test
- **Operation Performance**: Response times for each database operation
- **Bottlenecks**: Identified performance bottlenecks
- **Recommendations**: Suggested optimizations

Key metrics to look for:
- **P95 Response Time**: 95% of requests complete within this time
- **Max Response Time**: Worst-case response time
- **Throughput**: Operations per second

### Query Plan Analysis Report

The query plan analysis report (`query_plan_analysis_report.txt`) includes:

- **Execution Plans**: How the database executes each query
- **Issues**: Identified issues in query execution
- **Optimization Suggestions**: Recommended indexes and query optimizations

Look for:
- **Sequential Scans**: Indicates missing indexes
- **Hash Joins**: May indicate missing indexes on join columns
- **High Cost Operations**: Operations with high computational cost

### Concurrent Load Test Report

The concurrent load test report (`concurrent_load_test_report.txt`) includes:

- **Overall Performance**: Total operations, average response time, throughput
- **Operation Performance**: Performance metrics for each operation type
- **Errors**: Any errors encountered during the test

Key metrics:
- **Average Response Time**: Average time to complete operations
- **Throughput**: Operations per second under load
- **Error Rate**: Percentage of operations that failed

### Comprehensive Report

The comprehensive report (`comprehensive_performance_report.txt`) combines results from all tests and provides:

- **Overall Assessment**: Whether the database is ready for production
- **Critical Issues**: Issues that must be addressed before production
- **Key Recommendations**: Most important optimization recommendations

## Interpreting the Assessment

The comprehensive report will indicate whether the database is ready for production:

- **READY**: The database can handle the expected production load
- **NOT READY**: Critical issues need to be addressed

Even if the database is deemed ready, consider implementing the recommended optimizations to improve performance further.

## Common Performance Issues and Solutions

### Slow Queries

- **Issue**: Queries with high response times (>500ms)
- **Solution**: Add indexes on commonly queried fields, optimize query structure

### High Resource Usage

- **Issue**: CPU or memory usage consistently above 80%
- **Solution**: Scale up database resources, optimize queries, implement caching

### Concurrency Issues

- **Issue**: Performance degradation with increasing concurrent users
- **Solution**: Implement connection pooling, optimize transaction handling, add indexes

### High Error Rates

- **Issue**: Operations frequently failing under load
- **Solution**: Implement retry logic, optimize database configuration, check for deadlocks

## Next Steps

After running the performance tests:

1. Review the comprehensive report and address any critical issues
2. Implement the recommended optimizations
3. Re-run the tests to verify improvements
4. Adjust database configuration based on the results
5. Document the performance characteristics for operational planning

## Troubleshooting

If you encounter issues running the tests:

- **Connection Errors**: Verify Supabase credentials in the `.env` file
- **Permission Errors**: Ensure you're using a service role key with appropriate permissions
- **Memory Errors**: Reduce the number of concurrent users or iterations
- **Timeout Errors**: Increase timeout settings or optimize the database

For additional help, consult the CryoProtect v2 documentation or contact the development team.