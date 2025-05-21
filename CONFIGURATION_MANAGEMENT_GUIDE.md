# CryoProtect Configuration Management Guide

This guide provides detailed instructions for managing configuration in the CryoProtect application across different environments, from development to production deployment.

## Table of Contents

1. [Overview](#overview)
2. [Configuration Architecture](#configuration-architecture)
3. [Getting Started](#getting-started)
4. [Working with the Schema](#working-with-the-schema)
5. [Environment-Specific Configuration](#environment-specific-configuration)
6. [Secret Management](#secret-management)
7. [Configuration Overrides](#configuration-overrides)
8. [Validation](#validation)
9. [Deployment Workflow](#deployment-workflow)
10. [Troubleshooting](#troubleshooting)

## Overview

The CryoProtect configuration system is designed to provide:

- A **single source of truth** for all configuration parameters
- **Type safety** across both backend (Python) and frontend (TypeScript)
- **Environment-specific** configuration for development, testing, staging, and production
- **Secure handling of secrets** and sensitive data
- **Validation** to catch misconfigurations early
- **Override mechanisms** for local development and testing

## Configuration Architecture

The configuration system consists of these key components:

1. **Schema Definition**: A `schema.json` file that defines all configuration parameters, their types, defaults, and environment overrides.

2. **Code Generators**: Tools that generate language-specific configuration files from the schema:
   - Python configuration for the backend
   - TypeScript configuration for the frontend

3. **Validation Tools**: Scripts to validate configuration against the schema and check for inconsistencies.

4. **Secret Management**: Tools for securely handling sensitive configuration values.

5. **Override System**: Mechanisms to override configuration values without modifying the main configuration files.

## Getting Started

### Installation

The configuration system is located in the `/config` directory. To get started:

```bash
# Navigate to the config directory
cd /path/to/CryoProtect/config

# Make sure scripts are executable
chmod +x *.py
```

### Basic Usage

1. **Generate configuration files**:

```bash
./generate_configs.py
```

This will:
- Read the schema from `schema.json`
- Generate `../config.py.new` for the backend
- Generate `../frontend/src/config/config.ts` for the frontend

2. **Validate the configuration**:

```bash
./validate_config.py
```

3. **Manage secrets**:

```bash
./secret_manager.py generate
```

4. **Create an override file**:

```bash
./config_override.py generate --output my_overrides.json
```

## Working with the Schema

The schema (`schema.json`) is the central definition of all configuration parameters. It follows this structure:

```json
{
  "properties": {
    "section_name": {
      "type": "object",
      "description": "Section description",
      "properties": {
        "property_name": {
          "type": "string",
          "description": "Property description",
          "default": "default value"
        }
      }
    }
  },
  "environmentMapping": {
    "development": {
      "section_name.property_name": "development value"
    }
  }
}
```

### Adding New Configuration Parameters

To add a new configuration parameter:

1. Add it to the appropriate section in `schema.json`
2. Specify its type, description, and default value (if any)
3. Add environment-specific overrides if needed
4. Regenerate the configuration files

Example:

```json
"app": {
  "properties": {
    "new_feature_enabled": {
      "type": "boolean",
      "description": "Enable the new feature",
      "default": false
    }
  }
}
```

Then in `environmentMapping`:

```json
"development": {
  "app.new_feature_enabled": true
}
```

### Removing Configuration Parameters

To remove a parameter:

1. Remove it from `schema.json`
2. Regenerate the configuration files
3. Update any code that might be using the parameter

## Environment-Specific Configuration

The configuration system supports different environments through the `environmentMapping` section in the schema:

- `development`: Local development environment
- `testing`: Test/CI environment
- `staging`: Pre-production environment
- `production`: Production environment

### Setting the Current Environment

For the **backend**:

- Set the `FLASK_ENV` or `APP_ENV` environment variable to one of: `development`, `testing`, `staging`, `production`
- If not set, it defaults to `development`

For the **frontend**:

- Set the `NODE_ENV` environment variable
- In the React application, this is typically set automatically by the build system

## Secret Management

The configuration system includes tools for managing secrets securely:

### Listing Secrets

To see all secrets defined in the schema:

```bash
./secret_manager.py list
```

### Generating Secrets

To generate and store secrets for all properties marked as sensitive:

```bash
./secret_manager.py generate
```

This will:
- Identify all secret properties from the schema
- Generate random values or prompt for sensitive ones
- Store them in the `.env` file

### Using Docker Secrets

For production deployments with Docker:

```bash
./secret_manager.py generate --docker-secrets /path/to/secrets
```

This will create secret files for Docker Compose or Kubernetes.

### Rotating Secrets

To rotate a specific secret:

```bash
./secret_manager.py rotate APP_SECRET_KEY
```

## Configuration Overrides

The configuration system provides several ways to override configuration values:

### Local Override Files

Create a JSON file with overrides:

```bash
./config_override.py generate --output local_overrides.json --env development
```

Then apply it:

```bash
./config_override.py apply --override local_overrides.json
```

### Environment Variable Overrides

Set environment variables with the `OVERRIDE_` prefix:

```bash
OVERRIDE_APP_DEBUG=false python app.py
```

### Command-Line Overrides

When running the application:

```bash
python app.py --app-debug false
```

## Validation

To validate your configuration:

```bash
./validate_config.py
```

This will:
1. Validate the schema for correctness
2. Check the Python configuration against the schema
3. Check the TypeScript configuration against the schema
4. Look for inconsistencies between the two

### Continuous Integration

Add this to your CI pipeline:

```yaml
- name: Validate Configuration
  run: ./config/validate_config.py --check-only
```

## Deployment Workflow

### Development Environment

1. Use the default configuration
2. Override with a local `.env` file for secrets
3. Use local override files for customization

### Testing/CI Environment

1. Set `FLASK_ENV=testing` / `NODE_ENV=test`
2. Inject secrets from CI secrets store
3. Validate configuration before tests

### Staging Environment

1. Set `FLASK_ENV=staging` / `NODE_ENV=staging`
2. Use Docker secrets or environment variables for secrets
3. Validate configuration on startup

### Production Environment

1. Set `FLASK_ENV=production` / `NODE_ENV=production`
2. Use Docker secrets or a secure secrets manager
3. Validate configuration on startup
4. No debugging or experimental features

## Troubleshooting

### Common Issues

#### "Configuration validation failed"

- Check the specific error messages
- Ensure the schema is valid
- Make sure all required parameters have values
- Check for type mismatches

#### "Secret generation failed"

- Check file permissions on the .env file
- Ensure required secrets have values
- Check Docker secret directory permissions

#### "Schema validation failed"

- Check the JSON syntax in schema.json
- Ensure all sections and properties have the required attributes
- Verify environment mappings reference existing properties

### Getting Help

If you encounter issues with the configuration system:

1. Run validation with verbose output: `./validate_config.py --verbose`
2. Check the documentation in the `config/README.md` file
3. Reach out to the CryoProtect development team

## Further Reading

- [Configuration System Technical Documentation](./config/README.md)
- [Python Configuration API Documentation](../docs/API_DOCUMENTATION.md)
- [TypeScript Configuration Interface Documentation](../frontend/docs/CONFIG.md)