# CryoProtect v2 API Performance Benchmarking

This document provides instructions for running performance benchmarks on the CryoProtect v2 API endpoints and interpreting the results.

## Overview

The API performance benchmarking suite consists of scripts designed to evaluate different aspects of API performance:

1. **benchmark_api_endpoints.py**: Tests API endpoint performance, measuring response times, throughput, and error rates under various load conditions.
2. **run_api_benchmark.bat**: Windows batch script to run the benchmark.
3. **run_api_benchmark.sh**: Linux/macOS shell script to run the benchmark.

## Prerequisites

Before running the performance benchmarks, ensure you have:

1. Python 3.7+ installed
2. Required Python packages:
   ```
   pip install requests psutil
   ```
3. The CryoProtect v2 API server running (default: http://localhost:5000)

## Running the Benchmarks

### On Windows

1. Start the API server:
   ```
   .\run_app.bat
   ```

2. In a separate terminal, run the benchmark:
   ```
   .\run_api_benchmark.bat
   ```

### On Linux/macOS

1. Start the API server:
   ```
   ./run_app.sh
   ```

2. In a separate terminal, run the benchmark:
   ```
   chmod +x run_api_benchmark.sh
   ./run_api_benchmark.sh
   ```

### Command Line Options

The benchmark script supports several command line options:

```
usage: benchmark_api_endpoints.py [-h] [--base-url BASE_URL] [--num-requests NUM_REQUESTS]
                                 [--concurrency CONCURRENCY] [--timeout TIMEOUT]
                                 [--warmup-requests WARMUP_REQUESTS] [--auth-token AUTH_TOKEN]
                                 [--endpoints-file ENDPOINTS_FILE] [--output-file OUTPUT_FILE]

options:
  -h, --help            show this help message and exit
  --base-url BASE_URL   Base URL of the API (default: http://localhost:5000)
  --num-requests NUM_REQUESTS
                        Number of requests per endpoint (default: 50)
  --concurrency CONCURRENCY
                        Number of concurrent requests (default: 10)
  --timeout TIMEOUT     Request timeout in seconds (default: 10)
  --warmup-requests WARMUP_REQUESTS
                        Number of warmup requests (default: 5)
  --auth-token AUTH_TOKEN
                        Authentication token (if required)
  --endpoints-file ENDPOINTS_FILE
                        Path to endpoints file (default: memory-bank/api_endpoints.json)
  --output-file OUTPUT_FILE
                        Output report file (default: API_Performance_Benchmark_Report.md)
```

For example, to run the benchmark with 100 requests per endpoint and 20 concurrent requests:

```
.\run_api_benchmark.bat --num-requests 100 --concurrency 20
```

## Understanding the Results

The benchmark generates a comprehensive Markdown report (`API_Performance_Benchmark_Report.md`) with detailed performance metrics for each API endpoint.

### Performance Ratings

Endpoints are rated based on their average response time:

- **Excellent**: < 100ms
- **Good**: 100-300ms
- **Acceptable**: 300-1000ms
- **Poor**: > 1000ms

### Key Metrics

For each endpoint, the report includes:

1. **Response Time Metrics**:
   - Minimum response time
   - Maximum response time
   - Average response time
   - Median response time
   - 95th percentile response time

2. **Throughput and Error Metrics**:
   - Throughput (requests per second)
   - Success count
   - Error count
   - Error rate

3. **Status Code Distribution**:
   - Count and percentage of each status code returned

### System Resource Usage

The report also includes system resource usage metrics:

- CPU usage (minimum, maximum, average)
- Memory usage (minimum, maximum, average)

## Interpreting the Assessment

The report provides an overall assessment of API performance and identifies:

1. **Endpoints Requiring Performance Optimization**: Endpoints with poor performance (average response time > 1000ms)
2. **Endpoints with High Error Rates**: Endpoints with error rates > 10%
3. **General Recommendations**: Suggestions for improving overall API performance

## Common Performance Issues and Solutions

### Slow Endpoints

- **Issue**: Endpoints with high response times (>500ms)
- **Solution**: 
  - Add indexes on commonly queried fields
  - Optimize query structure
  - Implement caching for frequently accessed data
  - Check for N+1 query problems

### High Error Rates

- **Issue**: Endpoints with high error rates (>10%)
- **Solution**:
  - Review error handling
  - Improve input validation
  - Verify authentication and authorization logic
  - Add more comprehensive request validation

### Concurrency Issues

- **Issue**: Performance degradation with increasing concurrent users
- **Solution**:
  - Implement connection pooling
  - Optimize transaction handling
  - Add indexes
  - Consider asynchronous processing for long-running operations

## Next Steps

After running the performance benchmarks:

1. Review the comprehensive report and address any critical issues
2. Implement the recommended optimizations
3. Re-run the benchmarks to verify improvements
4. Adjust API configuration based on the results
5. Document the performance characteristics for operational planning

## Troubleshooting

If you encounter issues running the benchmarks:

- **Connection Errors**: Verify that the API server is running and accessible
- **Authentication Errors**: Check that the authentication token is valid
- **Timeout Errors**: Increase the timeout setting using the `--timeout` option
- **Memory Errors**: Reduce the number of concurrent requests using the `--concurrency` option

For additional help, consult the CryoProtect v2 documentation or contact the development team.