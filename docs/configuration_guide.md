# Configuration Guide

This document describes how to configure CryoProtect v2 for different environments.

## Configuration Overview

CryoProtect v2 uses a multi-tiered configuration approach:

1. **Environment Variables**: Primary configuration method
2. **Configuration Files**: Default values and structure
3. **Database Configuration**: Stored settings

## Environment Variables

Environment variables can be set in multiple ways:
- `.env` file in the project root
- System environment variables
- Docker environment variables

### Critical Environment Variables

| Variable | Description | Default | Required |
|----------|-------------|---------|----------|
| `FLASK_ENV` | Environment name | `development` | Yes |
| `SECRET_KEY` | Flask secret key | None | Yes |
| `JWT_SECRET_KEY` | JWT authentication key | None | Yes |
| `DB_CONNECTION_MODE` | Database connection strategy | `auto` | Yes |

### Database Connection Modes

CryoProtect v2 supports multiple database connection modes:

- `auto`: Tries Supabase, then local, then MCP
- `local`: Only uses local PostgreSQL
- `supabase`: Only uses Supabase
- `mcp`: Only uses MCP

### Local Database Configuration

| Variable | Description | Default | Required |
|----------|-------------|---------|----------|
| `LOCAL_DB_HOST` | Database host | `localhost` | Yes |
| `LOCAL_DB_PORT` | Database port | `5432` | No |
| `LOCAL_DB_NAME` | Database name | `cryoprotect` | Yes |
| `LOCAL_DB_USER` | Database user | `postgres` | Yes |
| `LOCAL_DB_PASSWORD` | Database password | None | Yes |
| `LOCAL_DB_MIN_CONNECTIONS` | Min connection pool size | `1` | No |
| `LOCAL_DB_MAX_CONNECTIONS` | Max connection pool size | `5` | No |

### Supabase Configuration

| Variable | Description | Default | Required |
|----------|-------------|---------|----------|
| `SUPABASE_URL` | Supabase URL | None | Yes for Supabase |
| `SUPABASE_KEY` | Anon key | None | Yes for Supabase |
| `SUPABASE_SERVICE_KEY` | Service role key | None | Yes for Supabase |
| `SUPABASE_PROJECT_ID` | Project ID | None | Yes for Supabase |
| `SUPABASE_DB_HOST` | Database host | None | Yes for Supabase |
| `SUPABASE_DB_PORT` | Database port | `5432` | No |
| `SUPABASE_DB_NAME` | Database name | `postgres` | No |
| `SUPABASE_DB_USER` | Database user | `postgres` | Yes for Supabase |
| `SUPABASE_DB_PASSWORD` | Database password | None | Yes for Supabase |
| `SUPABASE_DB_IP_ADDRESS` | IP address fallback | None | No |

### API Configuration

| Variable | Description | Default | Required |
|----------|-------------|---------|----------|
| `RATE_LIMIT` | API rate limit | `100/hour` | No |

## Environment-Specific Configuration Files

- `.env.template`: Template for all environments
- `.env.production`: Production-specific settings
- `.env.staging`: Staging-specific settings

## Configuration Hierarchy

Configuration values are loaded in the following order (later overrides earlier):

1. Default values in `config.py`
2. Environment-specific configuration class
3. Environment variables

## Creating Custom Environments

To create a custom environment:

1. Create a new configuration class in `config.py`
2. Create a corresponding `.env.[environment]` file
3. Set `FLASK_ENV=[environment]`

## Secrets Management

Secrets should never be committed to version control:

1. Store secrets in environment variables
2. Use a secrets manager for production
3. For development, use `.env` files (excluded from git)

## Configuration Validation

The application validates its configuration at startup:

1. Required values are checked
2. Database connection is tested
3. Validation errors are logged

## Examples

### Development Configuration

```env
FLASK_ENV=development
SECRET_KEY=dev-secret-key
DB_CONNECTION_MODE=local
LOCAL_DB_HOST=localhost
LOCAL_DB_USER=postgres
LOCAL_DB_PASSWORD=postgres
```

### Production Configuration

```env
FLASK_ENV=production
SECRET_KEY=production-secret-key
DB_CONNECTION_MODE=auto
SUPABASE_URL=https://xxx.supabase.co
SUPABASE_KEY=xxx
SUPABASE_SERVICE_KEY=xxx
SUPABASE_DB_HOST=db.xxx.supabase.co
SUPABASE_DB_USER=postgres
SUPABASE_DB_PASSWORD=xxx