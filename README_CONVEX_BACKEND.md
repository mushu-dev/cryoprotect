# Convex Backend Integration Guide

This README provides an overview of the Convex backend integration for CryoProtect, explaining how to set up, configure, and use the enhanced database adapter.

## Overview

The Convex backend integration provides a robust way to use Convex as a database backend for the CryoProtect API. This integration supports:

- Seamless database adapter pattern to work with both Supabase and Convex
- Bidirectional data synchronization between databases
- Enhanced authentication with JWT tokens
- Transaction support for data consistency
- Circuit breaker pattern for resilience
- Comprehensive error handling

## Quick Start

### 1. Configure Environment Variables

Create a `.env` file with the following variables:

```env
# Convex Connection
CONVEX_URL=https://your-convex-deployment.convex.cloud
CONVEX_DEPLOYMENT_KEY=your-deployment-key
CONVEX_DB_ENABLED=true

# Connection Settings
CONVEX_TIMEOUT=30
CONVEX_RETRY_COUNT=3
CONVEX_CIRCUIT_BREAKER_THRESHOLD=5
CONVEX_CIRCUIT_BREAKER_TIMEOUT=60

# JWT Authentication
JWT_SECRET=your-jwt-secret-key
JWT_ALGORITHM=HS256
JWT_EXPIRY_SECONDS=86400

# Sync Settings (optional)
CONVEX_SYNC_INTERVAL=60
CONVEX_SYNC_TABLES=molecules,experiments,protocols
CONVEX_CONFLICT_RESOLUTION=last_modified
```

### 2. Initialize in Flask App

```python
from flask import Flask
from database.flask_convex_integration import init_flask_convex

app = Flask(__name__)

# Initialize Flask-Convex integration
flask_convex = init_flask_convex(app)

# Get Convex adapter
convex = flask_convex.get_convex_adapter()

@app.route('/api/data')
def get_data():
    # Use Convex adapter to execute queries
    data = convex.execute_query('api.data.list', {
        'limit': 10
    })
    
    return jsonify(data)
```

### 3. Use in Connection Manager

The ConvexAdapter is already integrated with the ConnectionManager. Configure the adapter order to use Convex:

```env
# Use Convex first, then fall back to other adapters
DB_ADAPTER_ORDER=convex,supabase,local,mcp
```

Then use the ConnectionManager as usual:

```python
from database.connection_manager import ConnectionManager

# Get connection manager instance
manager = ConnectionManager.get_instance()

# Connect to database (will try Convex first)
if manager.connect():
    # Get active adapter
    adapter = manager.get_active_adapter()
    
    # Execute query
    result = adapter.execute_query("SELECT * FROM molecules LIMIT 10")
```

## Modules

### Enhanced Convex Adapter

The `enhanced_convex_adapter.py` module provides a comprehensive implementation of the DatabaseAdapter interface:

- **Features**:
  - Real-time data subscriptions
  - Full transaction support
  - Circuit breaker pattern for resilience
  - Comprehensive error handling
  - Supabase-like query interface

- **Usage Example**:
  ```python
  from database.enhanced_convex_adapter import ConvexAdapter
  
  # Create adapter
  convex = ConvexAdapter({
      'url': 'https://your-deployment.convex.cloud',
      'key': 'your-deployment-key'
  })
  
  # Connect to Convex
  if convex.connect():
      # Execute a query
      data = convex.execute_query('api.data.list', {
          'limit': 10
      })
  ```

### Flask Convex Integration

The `flask_convex_integration.py` module provides a convenient wrapper for integrating Convex with a Flask application:

- **Features**:
  - Flask middleware for JWT processing
  - Integration with Flask app context
  - Transaction management
  - Simplified API for database operations

- **Usage Example**:
  ```python
  from flask import Flask
  from database.flask_convex_integration import init_flask_convex
  
  app = Flask(__name__)
  
  # Initialize Flask-Convex integration
  flask_convex = init_flask_convex(app)
  
  # Use in routes
  @app.route('/api/data')
  def get_data():
      with flask_convex.transaction():
          # Operations within a transaction
          pass
  ```

### Convex Sync Manager

The `convex_sync.py` module provides bidirectional data synchronization between Supabase and Convex:

- **Features**:
  - Change detection
  - Conflict resolution
  - Bidirectional sync
  - Background sync worker

- **Usage Example**:
  ```python
  from database.enhanced_convex_adapter import ConvexAdapter
  from database.supabase_adapter import SupabaseDirectAdapter
  from database.convex_sync import create_sync_manager
  
  # Create adapters
  supabase = SupabaseDirectAdapter(...)
  convex = ConvexAdapter(...)
  
  # Create sync manager
  sync_manager = create_sync_manager(supabase, convex)
  
  # Start sync in background
  sync_manager.start_sync(background=True)
  ```

## Authentication

The `auth_bridge.py` module provides authentication integration between Flask and Convex:

- **Features**:
  - JWT generation and validation
  - Role-based access control
  - Compatible with both Supabase and Convex

- **Usage Example**:
  ```python
  from database.auth_bridge import AuthBridge
  
  # Create AuthBridge
  auth_bridge = AuthBridge(
      secret_key='your-jwt-secret',
      algorithm='HS256'
  )
  
  # Generate token
  token = auth_bridge.generate_token('user_123', {
      'email': 'user@example.com',
      'role': 'admin'
  })
  
  # Create Convex identity for frontend
  convex_identity = auth_bridge.create_convex_identity('user_123', {
      'email': 'user@example.com',
      'role': 'admin'
  })
  ```

## Example Application

A complete example application is available in `examples/enhanced_convex_example.py`. This example demonstrates:

- Flask integration with Convex
- Authentication with JWT tokens
- CRUD operations with Convex
- Transaction management
- Error handling

To run the example:

```bash
python examples/enhanced_convex_example.py
```

## Bidirectional Sync

The bidirectional sync feature keeps data consistent between Supabase and Convex:

- **Configuration**:
  ```env
  CONVEX_SYNC_INTERVAL=60
  CONVEX_SYNC_TABLES=molecules,experiments,protocols
  CONVEX_CONFLICT_RESOLUTION=last_modified
  ```

- **Conflict Resolution Strategies**:
  - `last_modified`: Use record with later modification time
  - `supabase_wins`: Always use Supabase record
  - `convex_wins`: Always use Convex record
  - Custom callback: Provide your own resolution function

- **Start Sync**:
  ```python
  # Start sync in background
  sync_manager.start_sync(background=True)
  
  # Or run sync immediately
  sync_manager.sync_all()
  ```

## Transactions

Transactions ensure data consistency across related operations:

- **Using Flask Integration**:
  ```python
  with flask_convex.transaction():
      # Operations within transaction
      convex.execute_query('api.data.create', {
          'data': {'name': 'Example'}
      })
  ```

- **Using Adapter Directly**:
  ```python
  # Begin transaction
  tx = convex.begin_transaction()
  
  try:
      # Execute operations
      convex_adapter._execute_query_in_transaction(
          tx, "INSERT INTO data (name) VALUES ('Example')"
      )
      
      # Commit transaction
      convex.commit_transaction(tx)
      
  except Exception:
      # Rollback on error
      convex.rollback_transaction(tx)
      raise
  ```

## Resilience Features

The Convex integration includes several resilience features:

- **Circuit Breaker Pattern**:
  ```python
  convex = ConvexAdapter({
      'circuit_breaker_threshold': 5,  # Failures before opening
      'circuit_breaker_timeout': 60    # Seconds before half-open
  })
  ```

- **Retry Logic**:
  ```python
  convex = ConvexAdapter({
      'retry_count': 3  # Number of retries
  })
  ```

- **Timeout Handling**:
  ```python
  convex = ConvexAdapter({
      'timeout': 30  # Seconds
  })
  ```

## Monitoring

Monitor the Convex integration:

```python
# Get connection information
connection_info = convex.get_connection_info()
print(f"Convex status: {connection_info['connected']}")
print(f"Circuit breaker state: {connection_info['circuit_state']}")

# Test connection
success, message = convex.test_connection()
print(f"Connection test: {message}")
```

## Advanced Configuration

### Adapter Order

Configure the order in which adapters are tried:

```env
# Try Convex first, then Supabase, then local
DB_ADAPTER_ORDER=convex,supabase,local,mcp
```

### Custom Auth Bridge

```python
from database.auth_bridge import init_auth_bridge

auth_bridge = init_auth_bridge(app, 
    secret_key='custom-secret',
    algorithm='HS384',
    expiry=12 * 60 * 60  # 12 hours
)
```

### Custom Sync Callbacks

```python
def on_conflict(table, record_id, supabase_data, convex_data):
    # Custom conflict resolution logic
    if table == 'important_table':
        return supabase_data  # Supabase wins for important table
    else:
        # Merge data for other tables
        merged = convex_data.copy()
        for key, value in supabase_data.items():
            if key not in convex_data:
                merged[key] = value
        return merged

# Set callbacks
sync_manager.set_callbacks(
    on_conflict=on_conflict,
    on_sync_error=lambda e: logger.error(f"Sync error: {str(e)}"),
    on_sync_start=lambda: logger.info("Sync started"),
    on_sync_complete=lambda: logger.info("Sync completed")
)
```

## Troubleshooting

### Connection Issues

- Check Convex URL and deployment key
- Ensure network connectivity to Convex
- Check if circuit breaker is open
- Try increasing timeout value

### Authentication Problems

- Verify JWT secret is consistent
- Check token expiry and roles
- Ensure Convex claims are properly structured
- Test token with Convex directly

### Sync Errors

- Check for schema differences between Supabase and Convex
- Verify record IDs are consistent across both systems
- Look for data type mismatches
- Check conflict resolution strategy

### Transaction Failures

- Ensure transaction IDs are tracked
- Check for duplicate operations
- Look for concurrent transaction conflicts
- Verify transaction state before commit/rollback

## Further Documentation

For more detailed information, see:

- [CONVEX_BACKEND_INTEGRATION.md](CONVEX_BACKEND_INTEGRATION.md): Comprehensive backend integration guide
- [CONVEX_FRONTEND_INTEGRATION.md](CONVEX_FRONTEND_INTEGRATION.md): Frontend integration guide
- [API Documentation](https://your-convex-deployment.convex.cloud/api)

## License

This project is licensed under the same terms as the main CryoProtect project.