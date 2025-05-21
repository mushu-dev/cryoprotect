# Connection Pool Tuning Guide

This guide provides comprehensive instructions for optimizing the database connection pool based on stress test results. It outlines key parameters, monitoring strategies, and recommended configurations for various workloads.

## Connection Pool Parameters

| Parameter | Description | Default | Recommended Range |
|-----------|-------------|---------|------------------|
| `min_connections` | Minimum number of connections to maintain | 1-2 | 2-5 |
| `max_connections` | Maximum number of connections allowed | 10 | 10-30 |
| `connection_timeout` | Connection acquisition timeout (seconds) | 30 | 5-30 |
| `connection_lifetime` | Maximum connection lifetime (seconds) | 3600 | 1800-3600 |
| `idle_timeout` | Maximum idle time before recycling (seconds) | 300 | 180-600 |
| `health_check_interval` | Interval between health checks (seconds) | 60 | 30-120 |
| `retry_attempts` | Number of retry attempts | 3 | 3-5 |
| `initial_retry_delay` | Initial delay before retry (seconds) | 0.2 | 0.1-0.5 |
| `max_retry_delay` | Maximum delay between retries (seconds) | 10 | 5-15 |
| `retry_jitter_factor` | Random jitter factor to avoid thundering herd | 0.1 | 0.05-0.2 |
| `validation_query` | Query used to validate connections | "SELECT 1" | Simple query |
| `circuit_breaker_threshold` | Failures to trigger circuit breaker | 5 | 3-10 |
| `circuit_breaker_timeout` | Time circuit stays open (seconds) | 30 | 10-60 |
| `circuit_breaker_reset` | Successes needed to close circuit | 2 | 1-3 |

## Workload-Based Configurations

### Light Workload (< 10 concurrent users)
```python
connection_pool_config = {
    'min_connections': 2,
    'max_connections': 10,
    'connection_timeout': 20,
    'connection_lifetime': 3600,  # 1 hour
    'idle_timeout': 300,  # 5 minutes
    'health_check_interval': 60,  # 1 minute
    'retry_attempts': 3,
}
```

### Medium Workload (10-50 concurrent users)
```python
connection_pool_config = {
    'min_connections': 3,
    'max_connections': 20,
    'connection_timeout': 15,
    'connection_lifetime': 1800,  # 30 minutes
    'idle_timeout': 240,  # 4 minutes
    'health_check_interval': 45,  # 45 seconds
    'retry_attempts': 3,
    'initial_retry_delay': 0.2,
    'max_retry_delay': 5,
    'retry_jitter_factor': 0.1,
}
```

### Heavy Workload (50+ concurrent users)
```python
connection_pool_config = {
    'min_connections': 5,
    'max_connections': 30,
    'connection_timeout': 10,
    'connection_lifetime': 1800,  # 30 minutes
    'idle_timeout': 180,  # 3 minutes
    'health_check_interval': 30,  # 30 seconds
    'retry_attempts': 5,
    'initial_retry_delay': 0.1,
    'max_retry_delay': 5,
    'retry_jitter_factor': 0.15,
    'validation_query': 'SELECT 1',
    'validation_timeout': 5,
    'circuit_breaker_threshold': 5,
    'circuit_breaker_timeout': 30,
    'circuit_breaker_reset': 2,
}
```

## Stress Testing

The stress test script provides comprehensive testing of the connection pool's performance and reliability. Run it with different configurations to determine the optimal settings for your workload.

```bash
# Basic test
./stress_test_connection_pool.py

# Comprehensive test suite
./stress_test_connection_pool.py --full-suite

# Custom test configuration
./stress_test_connection_pool.py --workers 25 --duration 60 --min-connections 3 --max-connections 25
```

### Key Metrics to Monitor

1. **Queries Per Second (QPS)**: Overall throughput of the system
2. **Success Rate**: Percentage of successful queries (target: >99%)
3. **Query Time**: Average and 95th percentile query execution time
4. **Connection Acquisition Time**: Time to acquire a connection from the pool
5. **Pool Utilization**: Ratio of active to total connections
6. **Error Rates**: Connection errors, timeouts, and other failures

### Common Performance Issues

| Issue | Symptoms | Solution |
|-------|----------|----------|
| Pool too small | High connection acquisition times, timeouts | Increase `max_connections` |
| Pool too large | High memory usage, database connection errors | Decrease `max_connections` |
| Connection leaks | Steadily increasing active connections | Reduce `connection_lifetime`, ensure proper connection release |
| Frequent reconnects | High CPU usage, connection spikes | Increase `connection_lifetime` and `idle_timeout` |
| Cascading failures | Multiple failures in short time | Enable circuit breaker, adjust retry parameters |

## Implementation Steps

1. **Baseline Test**: Run the stress test with default settings to establish a baseline
2. **Parameter Tuning**: Adjust parameters based on baseline results
3. **Verification Test**: Run the stress test with new parameters to verify improvements
4. **Monitoring Setup**: Configure monitoring to track connection pool metrics
5. **Production Deployment**: Gradually roll out changes to production

## Example Optimized Configuration for CryoProtect

Based on stress test results, the following configuration is recommended for the CryoProtect production environment:

```python
# In config.py
SUPABASE_MIN_CONNECTIONS = 3
SUPABASE_MAX_CONNECTIONS = 20
SUPABASE_CONNECTION_TIMEOUT = 20
SUPABASE_CONNECTION_LIFETIME = 1800
SUPABASE_IDLE_TIMEOUT = 180
SUPABASE_HEALTH_CHECK_INTERVAL = 30
SUPABASE_RETRY_ATTEMPTS = 4
SUPABASE_INITIAL_RETRY_DELAY = 0.2
SUPABASE_MAX_RETRY_DELAY = 5
SUPABASE_RETRY_JITTER_FACTOR = 0.1
SUPABASE_VALIDATION_QUERY = 'SELECT 1'
SUPABASE_VALIDATION_TIMEOUT = 5
SUPABASE_CIRCUIT_BREAKER_THRESHOLD = 5
SUPABASE_CIRCUIT_BREAKER_TIMEOUT = 30
SUPABASE_CIRCUIT_BREAKER_RESET = 2
```

## Monitoring

Set up monitoring to continuously track the following metrics:

1. **Active Connections**: Number of connections currently in use
2. **Total Connections**: Total connections in the pool
3. **Pool Utilization**: Percentage of connections in use
4. **Connection Acquisition Time**: Time to get a connection from the pool
5. **Query Execution Time**: Time to execute queries
6. **Connection Errors**: Count of connection-related errors
7. **Circuit Breaker Status**: Whether the circuit breaker is open, closed, or half-open

### Monitoring Implementation

For comprehensive monitoring, consider:

1. **Logging**: Configure detailed logging of connection pool events
2. **Metrics Collection**: Use Prometheus or similar tools to collect metrics
3. **Dashboards**: Create dashboards to visualize metrics
4. **Alerts**: Set up alerts for unusual conditions

## Conclusion

Proper connection pool tuning is critical for database performance and reliability. Regular stress testing and monitoring will help maintain optimal configuration as your application evolves. Use the provided stress test script to validate changes before deploying to production.