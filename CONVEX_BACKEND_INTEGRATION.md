# Convex Backend Integration

This document describes the backend integration of Convex with CryoProtect, focusing on how the Flask API seamlessly works with both Supabase and Convex databases.

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Components](#components)
4. [Configuration](#configuration)
5. [Usage Examples](#usage-examples)
6. [Bidirectional Sync](#bidirectional-sync)
7. [Resilience Features](#resilience-features)
8. [Authentication](#authentication)
9. [Transactions](#transactions)
10. [Monitoring](#monitoring)
11. [Best Practices](#best-practices)

## Overview

The Convex backend integration provides a robust and resilient way to use Convex as a database backend for the CryoProtect API. This integration supports the following key features:

- Seamless database adapter pattern to work with both Supabase and Convex
- Bidirectional data synchronization between Supabase and Convex
- Enhanced authentication with JWT tokens that work across systems
- Transaction support for data consistency
- Circuit breaker pattern for resilience
- Comprehensive error handling
- Monitoring and observability hooks

## Architecture

The backend integration is built around a modular, adapter-based architecture:

```
┌────────────────────────────────────────────────────┐
│ Flask API                                          │
│                                                    │
│  ┌────────────────────┐    ┌────────────────────┐  │
│  │ Auth Bridge        │    │ Connection Manager │  │
│  │                    │    │                    │  │
│  │  • JWT Generation  │    │  • Adapter Factory │  │
│  │  • Token Validation│    │  • Connection Pool │  │
│  │  • Role-Based Auth │    │  • Failover Logic  │  │
│  └────────────────────┘    └────────────────────┘  │
│                                                    │
│  ┌────────────────────┐    ┌────────────────────┐  │
│  │ Database Adapters  │    │ Sync Manager       │  │
│  │                    │    │                    │  │
│  │  • Supabase Adapter│    │  • Data Sync       │  │
│  │  • Convex Adapter  │    │  • Conflict Resolver│  │
│  │  • Local Adapter   │    │  • Change Detection │  │
│  └────────────────────┘    └────────────────────┘  │
│                                                    │
│  ┌────────────────────┐    ┌────────────────────┐  │
│  │ Request/Response   │    │ Error Handling     │  │
│  │ Middleware         │    │                    │  │
│  │  • JWT Processing  │    │  • Circuit Breaker │  │
│  │  • Convex Token    │    │  • Retry Logic     │  │
│  │  • Transaction     │    │  • Error Logging   │  │
│  └────────────────────┘    └────────────────────┘  │
│                                                    │
└────────────────────────────────────────────────────┘
```

## Components

### 1. Enhanced Convex Adapter

The `ConvexAdapter` (in `enhanced_convex_adapter.py`) provides a comprehensive implementation of the `DatabaseAdapter` interface, enabling seamless integration with the existing database connection infrastructure. It includes:

- Real-time data subscriptions
- Full transaction support
- Comprehensive error handling
- Circuit breaker pattern for resilience
- Supabase-like query interface for compatibility

### 2. Auth Bridge

The `AuthBridge` (in `auth_bridge.py`) handles authentication between Flask and Convex, including:

- JWT generation and validation compatible with both systems
- Role-based access control
- Token refresh and management
- Middleware for JWT processing

### 3. Flask-Convex Integration

The `FlaskConvexIntegration` (in `flask_convex_integration.py`) provides a convenient wrapper for integrating Convex with a Flask application, including:

- Flask middleware for JWT processing
- Integration with the Flask app context
- Transaction management
- Simplified API for database operations
- Context managers for transactions

### 4. Convex Sync Manager

The `ConvexSyncManager` (in `convex_sync.py`) handles bidirectional data synchronization between Supabase and Convex, including:

- Change detection
- Conflict resolution
- Bidirectional sync
- Background sync worker
- Sync progress tracking

### 5. Connection Manager Integration

The existing `ConnectionManager` (in `connection_manager.py`) has been updated to include the `ConvexAdapter`, allowing it to be used as a database backend alongside Supabase and other adapters.

## Configuration

Configure the Convex integration using environment variables or direct configuration:

```ini
# Convex Connection
CONVEX_URL=https://dynamic-mink-63.convex.cloud
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

# Sync Settings
CONVEX_SYNC_INTERVAL=60
CONVEX_SYNC_TABLES=molecules,experiments,protocols
CONVEX_CONFLICT_RESOLUTION=last_modified
```

## Usage Examples

### Basic Connection and Query

```python
from database.enhanced_convex_adapter import ConvexAdapter

# Create adapter
convex = ConvexAdapter({
    'url': 'https://dynamic-mink-63.convex.cloud',
    'key': 'your-deployment-key'
})

# Connect to Convex
if convex.connect():
    # Execute a query
    molecules = convex.execute_query('api.molecules.list', {
        'limit': 10
    })
    
    print(f"Found {len(molecules)} molecules")
```

### Flask Integration

```python
from flask import Flask, jsonify
from database.flask_convex_integration import init_flask_convex

app = Flask(__name__)

# Initialize Flask-Convex integration
flask_convex = init_flask_convex(app, {
    'convex_url': 'https://dynamic-mink-63.convex.cloud',
    'convex_key': 'your-deployment-key'
})

# Get Convex adapter
convex = flask_convex.get_convex_adapter()

@app.route('/api/molecules')
def get_molecules():
    molecules = convex.execute_query('api.molecules.list', {
        'limit': 10
    })
    
    return jsonify(molecules)

if __name__ == '__main__':
    app.run(debug=True)
```

### Transactions

```python
from database.flask_convex_integration import init_flask_convex

# Initialize Flask-Convex integration
flask_convex = init_flask_convex(app)

# Use transaction context manager
with flask_convex.transaction() as tx:
    # Execute queries within transaction
    convex.execute_query('api.molecules.create', {
        'data': {'name': 'Water', 'formula': 'H2O'}
    })
    
    convex.execute_query('api.experiments.create', {
        'data': {'name': 'Water Experiment'}
    })
    
    # Transaction is automatically committed on success
    # or rolled back on exception
```

### Connection Manager

```python
from database.connection_manager import ConnectionManager

# Get connection manager instance
manager = ConnectionManager.get_instance()

# Connect to database (will try adapters in configured order)
if manager.connect():
    # Get active adapter
    adapter = manager.get_active_adapter()
    
    # Execute query
    result = adapter.execute_query("SELECT * FROM molecules LIMIT 10")
    
    print(f"Found {len(result)} molecules")
```

## Bidirectional Sync

The bidirectional sync provides data consistency between Supabase and Convex:

```python
from database.enhanced_convex_adapter import ConvexAdapter
from database.supabase_adapter import SupabaseDirectAdapter
from database.convex_sync import create_sync_manager

# Create adapters
supabase = SupabaseDirectAdapter({
    'host': 'your-supabase-host',
    'user': 'your-supabase-user',
    'password': 'your-supabase-password'
})

convex = ConvexAdapter({
    'url': 'https://dynamic-mink-63.convex.cloud',
    'key': 'your-deployment-key'
})

# Connect adapters
supabase.connect()
convex.connect()

# Create sync manager
sync_manager = create_sync_manager(supabase, convex, {
    'sync_interval': 60,
    'sync_tables': ['molecules', 'experiments'],
    'conflict_resolution': 'last_modified'
})

# Start sync in background
sync_manager.start_sync(background=True)

# Force sync a specific record
sync_manager.force_sync_record('molecules', 'molecule-id-123')

# Stop sync when done
sync_manager.stop_sync_thread()
```

## Resilience Features

The Convex integration includes several resilience features:

### Circuit Breaker Pattern

The ConvexAdapter implements the circuit breaker pattern to prevent cascading failures:

- **Closed State**: Normal operation, requests flow through
- **Open State**: After repeated failures, requests are blocked
- **Half-Open State**: Testing if the service is back

### Retry Logic

Automatic retry with exponential backoff for transient errors:

```python
# Configure retry count
convex = ConvexAdapter({
    'url': 'https://dynamic-mink-63.convex.cloud',
    'key': 'your-deployment-key',
    'retry_count': 3  # Number of retries
})
```

### Timeout Handling

Configurable timeouts to prevent hanging requests:

```python
# Configure timeout
convex = ConvexAdapter({
    'url': 'https://dynamic-mink-63.convex.cloud',
    'key': 'your-deployment-key',
    'timeout': 30  # Seconds
})
```

## Authentication

The authentication system uses JWT tokens that work with both Supabase and Convex:

```python
from database.auth_bridge import AuthBridge

# Create AuthBridge
auth_bridge = AuthBridge(
    secret_key='your-jwt-secret',
    algorithm='HS256',
    expiry=86400  # 24 hours
)

# Generate token with user data
user_id = 'user_123'
user_data = {
    'email': 'user@example.com',
    'role': 'admin',
    'name': 'Example User'
}

token = auth_bridge.generate_token(user_id, user_data)

# Create Convex-compatible identity for frontend
convex_identity = auth_bridge.create_convex_identity(user_id, user_data)
```

## Transactions

Transactions ensure data consistency across related operations:

```python
# Begin transaction
tx = convex.begin_transaction()

try:
    # Execute operations within transaction
    convex_adapter._execute_query_in_transaction(
        tx, "INSERT INTO molecules (name) VALUES ('Water')"
    )
    
    convex_adapter._execute_query_in_transaction(
        tx, "INSERT INTO experiments (name) VALUES ('Water Experiment')"
    )
    
    # Commit transaction
    convex.commit_transaction(tx)
    
except Exception as e:
    # Rollback transaction on error
    convex.rollback_transaction(tx)
    raise
```

## Monitoring

The integration includes monitoring and observability features:

```python
# Get connection information
connection_info = convex.get_connection_info()
print(f"Convex connection status: {connection_info['connected']}")
print(f"Circuit breaker state: {connection_info['circuit_state']}")
print(f"Active transactions: {connection_info['active_transactions']}")

# Test connection
success, message = convex.test_connection()
print(f"Connection test: {message}")
```

## Best Practices

1. **Use the Adapter Pattern**: Always use the database adapter pattern to ensure compatibility with both Supabase and Convex.

2. **Always Use Transactions**: Wrap related operations in transactions to ensure data consistency.

3. **Handle Authentication Tokens**: Ensure JWT tokens are properly managed and include Convex-specific claims.

4. **Error Handling**: Implement proper error handling and use the circuit breaker pattern for resilience.

5. **Data Synchronization**: Use the bidirectional sync manager to keep data consistent between Supabase and Convex.

6. **Configuration**: Use environment variables for configuration to make deployment easier.

7. **Monitoring**: Monitor connection status, circuit breaker state, and transaction activity.

8. **Testing**: Test with both Supabase and Convex databases to ensure compatibility.