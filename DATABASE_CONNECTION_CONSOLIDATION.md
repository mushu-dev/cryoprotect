# Database Connection Consolidation Plan

This document outlines the approach for consolidating the multiple Supabase connection mechanisms in the CryoProtect project to reduce technical debt and streamline the database interaction.

## Current State Analysis

The codebase currently has multiple overlapping approaches to connect to Supabase:

1. **Direct REST API Connection** (simplified_app.py)
   - Uses requests library
   - Simple HTTP calls to Supabase REST API
   - Minimal dependencies
   - Currently working in Fedora environment

2. **Adapter Pattern** (database/adapter.py, database/supabase_adapter.py)
   - Abstract DatabaseAdapter interface
   - SupabaseDirectAdapter implementation with connection pooling
   - More robust error handling and reconnection logic
   - Complex but comprehensive

3. **Connection Factory** (database/connection.py)
   - More sophisticated with multiple adapter types
   - Fallback logic between connection methods
   - Health checks and monitoring
   - Possibly overengineered for current needs

4. **MCP Adapter** (database/mcp_adapter.py)
   - Uses Supabase MCP for direct SQL execution
   - Simpler approach with stateless connections
   - Works well with Cursor IDE environment

5. **Legacy/Miscellaneous Approaches**
   - Various direct connection scripts scattered throughout the codebase
   - Multiple connection utilities and helpers
   - Redundant and potentially conflicting implementations

## Consolidation Strategy

### 1. Primary Connection Approach

For the Fedora migration, we'll standardize on a hybrid approach combining:

- **Direct REST API** for simple operations (already working in simplified_app.py)
- **MCP Adapter** for complex SQL operations (leveraging Cursor IDE capabilities)

### 2. Implementation Steps

#### Phase 1: Establish Core Connection Module

1. Create a new `database/core.py` module that provides:
   - Simple REST API connection functions
   - MCP SQL execution functions
   - Configuration loading from environment
   - Basic connection validation functions

2. Create convenience functions for common operations:
   - `execute_query(query, params=None)`
   - `get_table_data(table_name, filters=None, limit=100)`
   - `insert_data(table_name, data)`
   - `update_data(table_name, id, data)`
   - `delete_data(table_name, id)`

#### Phase 2: Update Simplified App

1. Refactor simplified_app.py to use the new core module
2. Add more endpoints that demonstrate the consolidated approach:
   - `/api/tables/{table_name}` - List records from specified table
   - `/api/execute-sql` - Execute SQL query through MCP
   - `/api/connection-test` - Run comprehensive connection test

#### Phase 3: Create Transition Plan for Existing Code

1. Add backward compatibility helpers to support existing code
2. Create utility functions to assist migration:
   - `get_legacy_adapter()` - Returns adapter compatible with old code
   - `get_connection_pool()` - Returns a connection pool for batch operations
   - `get_transaction_manager()` - Provides transaction capabilities

#### Phase 4: Documentation and Testing

1. Create comprehensive documentation:
   - Usage examples for different scenarios
   - Migration guide for existing code
   - Best practices for different operation types

2. Develop test suite for connection module:
   - Connection validation tests
   - Performance benchmarks
   - Error handling tests
   - Concurrent operation tests

## Implementation Details

### Core Connection Module Structure

```python
# database/core.py

import os
import logging
import requests
import json
from typing import Any, Dict, List, Optional, Union, Tuple

# Configure logging
logger = logging.getLogger(__name__)

# Load configuration from environment
def load_config():
    """Load database configuration from environment variables."""
    return {
        'supabase_url': os.environ.get('SUPABASE_URL'),
        'supabase_key': os.environ.get('SUPABASE_KEY'),
        'supabase_service_key': os.environ.get('SUPABASE_SERVICE_KEY'),
        'project_id': os.environ.get('MCP_PROJECT_ID')
    }

# REST API functions
def rest_request(method, path, data=None, params=None, headers=None):
    """Make a request to Supabase REST API."""
    config = load_config()
    url = f"{config['supabase_url']}/rest/v1/{path}"
    
    if headers is None:
        headers = {}
    
    headers.update({
        'apikey': config['supabase_key'],
        'Authorization': f"Bearer {config['supabase_key']}"
    })
    
    return requests.request(
        method=method,
        url=url,
        json=data,
        params=params,
        headers=headers
    )

# MCP SQL execution
def execute_sql(query, params=None, project_id=None):
    """Execute SQL query through MCP."""
    config = load_config()
    proj_id = project_id or config['project_id']
    
    if not proj_id:
        raise ValueError("Project ID not provided")
    
    # Handle query parameters by formatting the query
    if params:
        # Implementation similar to MCPAdapter.execute_query
        pass
    
    try:
        # Import mcp_tools
        from supabase_mcp_tools import execute_sql_through_mcp
        return execute_sql_through_mcp(query, proj_id)
    except ImportError:
        logger.error("Failed to import MCP tools. MCP SQL execution will not work.")
        raise

# Convenience functions
def get_table_data(table_name, filters=None, limit=100):
    """Get data from a table with optional filters."""
    path = table_name
    params = {'limit': limit}
    
    if filters:
        params.update(filters)
    
    response = rest_request('GET', path, params=params)
    response.raise_for_status()
    return response.json()

def insert_data(table_name, data):
    """Insert data into a table."""
    path = table_name
    response = rest_request('POST', path, data=data)
    response.raise_for_status()
    return response.json()

def update_data(table_name, id, data):
    """Update data in a table."""
    path = f"{table_name}?id=eq.{id}"
    headers = {'Prefer': 'return=representation'}
    response = rest_request('PATCH', path, data=data, headers=headers)
    response.raise_for_status()
    return response.json()

def delete_data(table_name, id):
    """Delete data from a table."""
    path = f"{table_name}?id=eq.{id}"
    response = rest_request('DELETE', path)
    response.raise_for_status()
    return True

# Legacy adapter compatibility
def get_legacy_adapter():
    """Return an adapter compatible with legacy code."""
    from .adapter import DatabaseAdapter
    
    # Implementation that wraps the core functions in adapter interface
    pass
```

### Migration Guidelines

1. **Direct REST API Operations**
   - Use for simple CRUD operations
   - Good for small data volumes
   - Example: `get_table_data('molecules', {'name': 'eq.Glycerol'})`

2. **MCP SQL Execution**
   - Use for complex queries or bulk operations
   - Better for large data volumes or joins
   - Example: `execute_sql("SELECT * FROM molecules WHERE name LIKE '%glycol%'")`

3. **Transactions**
   - For operations that require ACID properties
   - Example: Using the transaction manager

4. **Legacy Code Compatibility**
   - Use the adapter compatibility functions
   - Gradually migrate code to the new approach
   - Example: `adapter = get_legacy_adapter()`

## Advantages of this Approach

1. **Simplicity**: The direct REST API is simple to understand and maintain
2. **Power**: MCP provides SQL execution for complex operations
3. **Minimal Dependencies**: Reduced dependency on third-party libraries
4. **Fedora Compatibility**: Tested and working in Fedora environment
5. **Backward Compatibility**: Supports existing code through adapter pattern
6. **Flexibility**: Can be extended for different operation types
7. **Performance**: Uses the right tool for each job
8. **Maintainability**: Centralized connection logic

## Timeline

1. **Week 1**: Create core module and test in Fedora environment
2. **Week 2**: Update simplified app and implement REST endpoints
3. **Week 3**: Develop compatibility layer and migration utilities
4. **Week 4**: Documentation and testing

## Conclusion

This consolidation plan provides a pragmatic approach to managing database connections in the CryoProtect project. By focusing on a hybrid approach that combines simple REST API operations with powerful MCP SQL execution, we can reduce technical debt while ensuring compatibility with the Fedora environment.