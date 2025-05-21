# Enhanced Connection Pool

The Enhanced Connection Pool provides advanced connection management with features for reliability, performance optimization, and monitoring.

## Features

- **Health Checking**: Automatic validation of connections to detect stale or broken connections
- **Dynamic Pool Sizing**: Automatically adjusts pool size based on load and utilization
- **Circuit Breaker Pattern**: Prevents cascading failures when database connections fail
- **Detailed Metrics**: Comprehensive statistics on connection usage and performance
- **Async Support**: Full asynchronous support for high-performance applications
- **Connection Validation**: Verifies connections before returning them to clients

## Configuration

### Basic Configuration

```python
from unified_importer.core.enhanced_connection_pool import EnhancedConnectionPool

# Create a connection pool with default settings
pool = EnhancedConnectionPool(
    min_size=5,                # Minimum number of connections to maintain
    max_size=20,               # Maximum number of connections allowed
    target_utilization=0.7,    # Target pool utilization (70%)
    connection_type="direct",  # "direct" or "supabase"
    connection_params={
        'host': 'localhost',
        'port': 5432,
        'database': 'mydb',
        'user': 'myuser',
        'password': 'mypassword'
    }
)
```

### Advanced Configuration

```python
# Create a pool with advanced settings
pool = EnhancedConnectionPool(
    min_size=5,
    max_size=20,
    target_utilization=0.7,
    timeout=30.0,               # Connection acquisition timeout
    max_lifetime=3600.0,        # Maximum connection lifetime (1 hour)
    validation_interval=300.0,  # Validate connections every 5 minutes
    health_check_interval=60.0, # Run health check every minute
    connection_type="direct",
    connection_params={
        'host': 'localhost',
        'port': 5432,
        'database': 'mydb',
        'user': 'myuser',
        'password': 'mypassword'
    },
    logger=custom_logger        # Optional custom logger
)
```

## Usage

### Synchronous Usage

```python
# Get a connection from the pool
conn, conn_id, is_new = pool.get_connection()

try:
    # Use the connection
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM molecules LIMIT 10")
    results = cursor.fetchall()
finally:
    # Release the connection back to the pool
    pool.release_connection(conn, conn_id)
```

### Asynchronous Usage

```python
# Get a connection asynchronously
conn, conn_id, is_new = await pool.get_connection_async()

try:
    # Use the connection
    # Note: For direct PostgreSQL connections, you'll need to use asyncio.to_thread
    # to run blocking operations in a thread
    results = await asyncio.to_thread(
        lambda: conn.cursor().execute("SELECT * FROM molecules LIMIT 10").fetchall()
    )
finally:
    # Release the connection back to the pool
    await pool.release_connection_async(conn, conn_id)
```

## Monitoring and Metrics

The pool provides detailed metrics to monitor its health and performance:

```python
# Get pool statistics
stats = pool.get_stats()
print(f"Pool size: {stats['pool_size']}")
print(f"Idle connections: {stats['idle_connections']}")
print(f"Active connections: {stats['active_connections']}")
print(f"Connection acquisition count: {stats['acquired']}")
print(f"Connection release count: {stats['released']}")
print(f"Error count: {stats['errors']}")

# For async operations
stats = await pool.get_stats_async()
```

## Health Checks

The pool automatically runs health checks on a configurable interval, but you can also manually trigger them:

```python
# Run a manual health check
pool._run_health_check()

# Or asynchronously
health_info = await pool.run_health_check_async()
print(f"Pool health: {health_info}")
```

## Circuit Breaker

The pool includes a circuit breaker pattern to prevent cascading failures when the database is experiencing issues:

- **Closed State**: Normal operation, connections are created as needed
- **Open State**: After multiple failures, connection creation is blocked to prevent hammering the database
- **Half-Open State**: After a timeout, a single connection attempt is allowed to test if the database is back

The circuit breaker is automatically managed, but you can check its state:

```python
circuit_state = pool.conn_circuit_breaker.get_state()
print(f"Circuit breaker state: {circuit_state}")
```

## Dynamic Pool Sizing

The pool automatically adjusts its size based on utilization:

- If utilization is high (above target_utilization), the pool will grow up to max_size
- If utilization is low (significantly below target_utilization), the pool will shrink down toward min_size

This ensures optimal resource usage based on application demand.

## Connection Validation

Connections are validated before being returned to clients:

- Periodic validation based on validation_interval
- Validation after connection has been idle for a long time
- Validation when connections are returned to the pool

This ensures that clients never receive stale or broken connections.

## Error Handling

The pool includes robust error handling:

- Failed connection attempts are tracked and retried with exponential backoff
- Connections with repeated errors are discarded and replaced
- Health checks automatically repair the pool by replacing bad connections

## Best Practices

1. **Sizing Guidelines**:
   - Set min_size to cover your typical load
   - Set max_size to handle peak load without overwhelming the database
   - A good rule of thumb is max_size = 2 * number_of_cores

2. **Connection Lifetime**:
   - Use a reasonable max_lifetime (1-4 hours) to prevent resource leaks
   - Don't set it too low or you'll create unnecessary overhead

3. **Health Check Interval**:
   - For production, 30-60 seconds is a good balance
   - For development, longer intervals (5+ minutes) reduce overhead

4. **Target Utilization**:
   - 0.7 (70%) is a good default for most applications
   - Lower values create more idle connections but improve response time
   - Higher values save resources but might cause waits during usage spikes

## Thread Safety

The enhanced connection pool is fully thread-safe:

- Synchronized access to shared state with locks
- Safe to use from multiple threads concurrently
- Asynchronous operations are also thread-safe with async locks

## Cleanup

Always clean up the pool when shutting down your application:

```python
# Synchronously
pool.close_all()

# Asynchronously
await pool.close_all_async()
```