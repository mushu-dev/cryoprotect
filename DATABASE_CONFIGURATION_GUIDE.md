# CryoProtect Database Configuration Guide

This document explains the new unified database configuration system implemented in CryoProtect v2. The system provides a centralized, type-safe, and flexible approach to configuring and managing database connections.

## Table of Contents

1. [Overview](#overview)
2. [Key Features](#key-features)
3. [Configuration Structure](#configuration-structure)
4. [Connection Modes](#connection-modes)
5. [Adapter Types](#adapter-types)
6. [Environment Variables](#environment-variables)
7. [JSON Configuration](#json-configuration)
8. [Configuration Initialization](#configuration-initialization)
9. [Integration with Existing Code](#integration-with-existing-code)
10. [Testing and Verification](#testing-and-verification)
11. [Best Practices](#best-practices)
12. [Troubleshooting](#troubleshooting)

## Overview

The database configuration system provides a centralized way to manage all database connection settings across the application. It supports multiple connection methods (local PostgreSQL, Supabase) and provides a consistent interface for all database operations.

The system is designed to be:
- **Centralized**: All configuration is managed in one place
- **Type-safe**: All configuration values have defined types and validation
- **Flexible**: Supports multiple connection methods with automatic fallback
- **Environment-aware**: Automatically detects and adapts to different environments
- **Backward compatible**: Works with existing code with minimal changes

## Key Features

- **Unified Configuration**: Single source of truth for all database connection parameters
- **Multiple Adapters**: Support for local PostgreSQL and Supabase connections
- **Connection Order Preference**: Configure the order in which connection methods are tried
- **Type Validation**: Built-in validation for configuration values
- **Auto-fallback**: Automatically tries alternative connection methods if preferred one fails
- **Connection Pooling**: Integrated connection pool configuration
- **Service Role Support**: Special handling for service role connections to bypass RLS
- **Public API Access**: Support for public authenticated access with RLS

## Configuration Structure

The configuration is structured as a hierarchical JSON object with the following structure:

```json
{
  "database": {
    "connection_mode": "auto",
    "adapter_order": ["local", "supabase"],
    "connection_timeout": 30,
    "connection_lifetime": 3600,
    "idle_timeout": 300,
    "application_name": "CryoProtect",
    "local": {
      "enabled": true,
      "host": "localhost",
      "port": 5432,
      "database": "postgres",
      "user": "postgres",
      "password": "",
      "min_connections": 1,
      "max_connections": 5,
      "use_ssl": false
    },
    "supabase": {
      "enabled": true,
      "url": "",
      "key": "",
      "service_key": "",
      "project_id": "",
      "host": "",
      "port": 5432,
      "database": "postgres",
      "user": "",
      "password": "",
      "min_connections": 1,
      "max_connections": 10
    }
  }
}
```

## Connection Modes

The `connection_mode` setting determines which connection method to use:

- **auto**: Automatically try connection methods in the order specified by `adapter_order`
- **local**: Use only local PostgreSQL connection
- **supabase**: Use only Supabase connection

## Adapter Types

The system supports the following adapter types:

### Local PostgreSQL

Direct connection to a PostgreSQL database, typically running locally.

```json
"local": {
  "enabled": true,
  "host": "localhost",
  "port": 5432,
  "database": "postgres",
  "user": "postgres",
  "password": "",
  "min_connections": 1,
  "max_connections": 5,
  "use_ssl": false
}
```

### Supabase

Connection to a Supabase project's PostgreSQL database, supporting both the API and direct connections.

```json
"supabase": {
  "enabled": true,
  "url": "https://your-project.supabase.co",
  "key": "your-anon-key",
  "service_key": "your-service-key",
  "project_id": "your-project-id",
  "host": "db.your-project.supabase.co",
  "port": 5432,
  "database": "postgres",
  "user": "postgres",
  "password": "your-password",
  "min_connections": 1,
  "max_connections": 10
}
```

## Environment Variables

The configuration system can be configured using environment variables. The following environment variables are supported:

### General Configuration

- `DB_CONNECTION_MODE`: Connection mode (local, supabase, auto)
- `DB_ADAPTER_ORDER`: Comma-separated list of adapters in order of preference
- `DB_CONNECTION_TIMEOUT`: Connection timeout in seconds
- `DB_CONNECTION_LIFETIME`: Connection lifetime in seconds
- `DB_IDLE_TIMEOUT`: Idle timeout in seconds
- `DB_APPLICATION_NAME`: Application name for database connections

### Local PostgreSQL Configuration

- `LOCAL_DB_ENABLED`: Enable local PostgreSQL connection (true, false)
- `LOCAL_DB_HOST`: Local PostgreSQL host
- `LOCAL_DB_PORT`: Local PostgreSQL port
- `LOCAL_DB_NAME`: Local PostgreSQL database name
- `LOCAL_DB_USER`: Local PostgreSQL user
- `LOCAL_DB_PASSWORD`: Local PostgreSQL password
- `LOCAL_DB_MIN_CONNECTIONS`: Minimum connections in pool
- `LOCAL_DB_MAX_CONNECTIONS`: Maximum connections in pool
- `LOCAL_DB_USE_SSL`: Use SSL for local connection (true, false)

### Supabase Configuration

- `SUPABASE_ENABLED`: Enable Supabase connection (true, false)
- `SUPABASE_URL`: Supabase URL
- `SUPABASE_KEY`: Supabase anon key
- `SUPABASE_SERVICE_KEY`: Supabase service key
- `SUPABASE_PROJECT_ID`: Supabase project ID
- `SUPABASE_DB_HOST`: Supabase database host
- `SUPABASE_DB_PORT`: Supabase database port
- `SUPABASE_DB_NAME`: Supabase database name
- `SUPABASE_DB_USER`: Supabase database user
- `SUPABASE_DB_PASSWORD`: Supabase database password
- `SUPABASE_DB_MIN_CONNECTIONS`: Minimum connections in pool
- `SUPABASE_DB_MAX_CONNECTIONS`: Maximum connections in pool

## JSON Configuration

Besides environment variables, the configuration can also be loaded from a JSON file. The default location is `./config/database.json`, but this can be configured using the `CONFIG_FILE` environment variable.

## Configuration Initialization

The configuration can be initialized using the provided `initialize_db_config.py` script:

```bash
python initialize_db_config.py [--config-file CONFIG_FILE] [--env-file ENV_FILE] [--validate-only] [--override]
```

This script will:
1. Read configuration from environment variables
2. Validate the configuration
3. Write the configuration to a JSON file
4. Generate an .env file with all environment variables set

### Options

- `--config-file`: Path to configuration file to generate (default: `./config/database.json`)
- `--env-file`: Path to .env file to generate (default: `./.env.database`)
- `--validate-only`: Only validate the configuration, don't write to file
- `--override`: Override existing files

## Integration with Existing Code

The configuration system is designed to work with existing code with minimal changes. The following modules have been updated to use the new configuration system:

- Database connection factory (`database/connection.py`)
- Database adapters (local, pooler, supabase)
- Database operation modules (db.py, db_service_role.py, db_public.py)
- Database population scripts
- Database migration scripts
- Supabase client for API access

### Example Usage

```python
from database.connection_config import (
    validate_config,
    get_connection_config,
    test_adapter_configuration
)

# Validate configuration
validate_config()

# Get configuration for specific adapter
config = get_connection_config('local')

# Test adapter configuration
is_valid, message = test_adapter_configuration('local')
```

### Backward Compatibility

For backward compatibility, the system will fall back to environment variables if the configuration file is not found, or the `connection_config` module is not available. This ensures that existing code will continue to work without changes.

## Testing and Verification

A test script is provided to verify the configuration and connection functionality:

```bash
python test_new_configuration.py
```

This script will:
1. Validate the configuration
2. Test getting connection configurations for all adapters
3. Test the connection factory
4. Test each database module (db.py, db_service_role.py, db_public.py)
5. Test the Supabase client for API access
6. Generate a report with the test results

## Best Practices

1. **Use the configuration initialization script**: Always use the provided script to initialize the configuration, as it performs validation and ensures consistency.

2. **Prefer central configuration over environment variables**: While environment variables are supported for backward compatibility, prefer using the central configuration file for new code.

3. **Validate configuration before use**: Always call `validate_config()` before using any configuration values to ensure they are valid.

4. **Use typed configuration values**: When accessing configuration values, use the appropriate type (e.g., `int`, `bool`) to ensure type safety.

5. **Handle connection failures gracefully**: Always handle the case where a connection cannot be established, and provide appropriate fallback behavior.

6. **Test all connection methods**: When testing your code, ensure it works with all supported connection methods.

7. **Use connection pooling for performance**: Configure appropriate connection pool parameters for your application's workload.

## Troubleshooting

### Configuration Validation Errors

If you encounter validation errors when using the configuration system, check the following:

1. Make sure all required environment variables are set correctly.
2. Ensure the configuration file exists and has the correct format.
3. Verify that the adapter you're trying to use is enabled.

### Connection Failures

If you encounter connection failures, check the following:

1. Verify that the database server is running and accessible.
2. Check that the connection parameters (host, port, user, password) are correct.
3. Ensure that network connectivity to the database server is available.
4. For Supabase connections, verify that the project is active and the database is not paused.

### Module Import Errors

If you encounter module import errors, check the following:

1. Ensure you're importing the correct modules from the right locations.
2. Verify that the directory structure is correct.
3. Check that all required dependencies are installed.