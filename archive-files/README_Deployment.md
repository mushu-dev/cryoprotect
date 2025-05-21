# CryoProtect Analyzer - Deployment Guide

This document provides information on the deployment process for the CryoProtect Analyzer project, including the CI/CD pipeline and environment configurations.

## CI/CD Pipeline

The CryoProtect Analyzer project uses GitHub Actions for continuous integration and deployment. The pipeline is defined in `.github/workflows/ci-cd.yml` and consists of the following stages:

### 1. Test Stage

- Runs on every push to main/master and on pull requests
- Sets up Python and Conda environment
- Installs dependencies from environment.yml and requirements_updated.txt
- Runs the test suite using tests/run_tests.py

### 2. Docker Stage

- Builds the Docker image using the Dockerfile
- Tests the Docker image by:
  - Running the container
  - Executing the RDKit verification script
  - Checking the Flask application health endpoint
  - Stopping the container

### 3. Staging Deployment

- Triggered only on pushes to main/master
- Builds and pushes the Docker image to GitHub Container Registry
- Tags the image with 'staging' and the commit SHA
- Deploys the image to the staging environment

### 4. Production Deployment

- Triggered only on pushes to main/master
- Requires manual approval (using GitHub Environments)
- Deploys the staging image to the production environment

## Environment Configurations

The project includes environment-specific configuration files:

- `config.py` - Base configuration
- `config_staging.py` - Staging environment configuration
- `config_production.py` - Production environment configuration

To use a specific configuration, set the `FLASK_ENV` environment variable:

```bash
# For development
export FLASK_ENV=development

# For staging
export FLASK_ENV=staging

# For production
export FLASK_ENV=production
```

## Docker Configuration

The Docker container includes:

- Python 3.9
- RDKit 2023.9.1
- IPython for enhanced visualization
- Curl for health checks
- All required Python packages

## Docker Compose Profiles

The project now uses Docker Compose profiles to manage different deployment environments:

- `dev` - Development environment with volume mounting for live code changes
- `staging` - Staging environment with production-like configuration but potentially different resource limits
- `prod` - Production environment with resource limits and logging configuration
- `all` - Includes all services for comprehensive testing

To use a specific profile, use the `--profile` flag with docker-compose:

```bash
# Start the development environment
docker-compose --profile dev up -d

# Start the staging environment
docker-compose --profile staging up -d

# Start the production environment
docker-compose --profile prod up -d

# Start all services
docker-compose --profile all up -d
```

## Resource Limits

Resource limits are configured for both development and production environments:

### Production/Staging Resource Limits

- Memory Limit: 512MB (default, can be overridden)
- CPU Limit: 0.5 cores (default, can be overridden)
- Memory Reservation: 256MB (default, can be overridden)
- CPU Reservation: 0.25 cores (default, can be overridden)

### Development Resource Limits

- Memory Limit: 1GB (default, can be overridden)
- CPU Limit: 1.0 core (default, can be overridden)

You can override these limits using environment variables:

```bash
# Override production resource limits
export MEM_LIMIT=1G
export CPU_LIMIT=1.0
export MEM_RESERVATION=512M
export CPU_RESERVATION=0.5

# Override development resource limits
export DEV_MEM_LIMIT=2G
export DEV_CPU_LIMIT=2.0

# Start with custom resource limits
docker-compose --profile prod up -d
```

## Logging Configuration

Logging is configured with the following defaults:

### Production/Staging Logging

- Driver: json-file
- Max Size: 10MB
- Max Files: 3
- Compression: enabled

### Development Logging

- Driver: json-file
- Max Size: 20MB
- Max Files: 5

You can override these settings using environment variables:

```bash
# Override production logging settings
export LOG_DRIVER=local
export LOG_MAX_SIZE=20m
export LOG_MAX_FILE=5
export LOG_COMPRESS=false

# Override development logging settings
export DEV_LOG_DRIVER=local
export DEV_LOG_MAX_SIZE=50m
export DEV_LOG_MAX_FILE=10

# Start with custom logging settings
docker-compose --profile prod up -d
```

## Manual Deployment

If you need to deploy manually, follow these steps:

### Staging Deployment

```bash
# Build the Docker image
docker build -t cryoprotect:staging .

# Tag the image for the registry
docker tag cryoprotect:staging ghcr.io/yourusername/cryoprotect:staging

# Push the image to the registry
docker push ghcr.io/yourusername/cryoprotect:staging

# Create Docker secrets on the staging server
ssh user@staging-server "
  # Create Docker secrets for sensitive values
  echo 'your-staging-supabase-url' | docker secret create cryoprotect_supabase_url - || docker secret rm cryoprotect_supabase_url && echo 'your-staging-supabase-url' | docker secret create cryoprotect_supabase_url - &&
  echo 'your-staging-supabase-key' | docker secret create cryoprotect_supabase_key - || docker secret rm cryoprotect_supabase_key && echo 'your-staging-supabase-key' | docker secret create cryoprotect_supabase_key - &&
  echo 'your-staging-secret-key' | docker secret create cryoprotect_secret_key - || docker secret rm cryoprotect_secret_key && echo 'your-staging-secret-key' | docker secret create cryoprotect_secret_key -
"

# Deploy to staging server with resource limits
ssh user@staging-server "
  docker pull ghcr.io/yourusername/cryoprotect:staging &&
  export FLASK_ENV=staging &&
  export USE_EXTERNAL_SECRETS=true &&
  # Set custom resource limits for staging if needed
  export MEM_LIMIT=768M &&
  export CPU_LIMIT=0.75 &&
  # Set custom logging options if needed
  export LOG_MAX_SIZE=15m &&
  export LOG_MAX_FILE=4 &&
  # Deploy using the staging profile
  docker-compose --profile staging up -d
"
```

### Production Deployment

```bash
# Tag the staging image for production
docker tag ghcr.io/yourusername/cryoprotect:staging ghcr.io/yourusername/cryoprotect:production

# Push the image to the registry
docker push ghcr.io/yourusername/cryoprotect:production

# Create Docker secrets on the production server
ssh user@production-server "
  # Create Docker secrets for sensitive values
  echo 'your-production-supabase-url' | docker secret create cryoprotect_supabase_url - || docker secret rm cryoprotect_supabase_url && echo 'your-production-supabase-url' | docker secret create cryoprotect_supabase_url - &&
  echo 'your-production-supabase-key' | docker secret create cryoprotect_supabase_key - || docker secret rm cryoprotect_supabase_key && echo 'your-production-supabase-key' | docker secret create cryoprotect_supabase_key - &&
  echo 'your-production-secret-key' | docker secret create cryoprotect_secret_key - || docker secret rm cryoprotect_secret_key && echo 'your-production-secret-key' | docker secret create cryoprotect_secret_key - &&
  echo 'your-redis-url' | docker secret create cryoprotect_redis_url - || docker secret rm cryoprotect_redis_url && echo 'your-redis-url' | docker secret create cryoprotect_redis_url -
"

# Deploy to production server with resource limits
ssh user@production-server "
  docker pull ghcr.io/yourusername/cryoprotect:production &&
  export FLASK_ENV=production &&
  export USE_EXTERNAL_SECRETS=true &&
  # Set production resource limits
  export MEM_LIMIT=1G &&
  export CPU_LIMIT=1.0 &&
  export MEM_RESERVATION=512M &&
  export CPU_RESERVATION=0.5 &&
  # Set production logging options
  export LOG_DRIVER=json-file &&
  export LOG_MAX_SIZE=10m &&
  export LOG_MAX_FILE=5 &&
  export LOG_COMPRESS=true &&
  # Deploy using the production profile
  docker-compose --profile prod up -d
"
```

> **Note**: Replace the placeholder values with your actual secret values. Never store these values in your codebase or commit them to version control.

### Development Deployment

For local development with live code changes:

```bash
# Set development resource limits (optional)
export DEV_MEM_LIMIT=2G
export DEV_CPU_LIMIT=2.0

# Set development logging options (optional)
export DEV_LOG_MAX_SIZE=50m
export DEV_LOG_MAX_FILE=10

# Start the development environment with the dev profile
docker-compose --profile dev up
```

## Secrets Management

For secure deployment, we use Docker Secrets to manage sensitive values instead of environment variables. See [README_Secrets_Management.md](./README_Secrets_Management.md) for detailed information.

### Required Secrets

#### Staging Environment

- `cryoprotect_supabase_url` - Supabase URL for staging
- `cryoprotect_supabase_key` - Supabase key for staging
- `cryoprotect_secret_key` - Secret key for session encryption and JWT signing

#### Production Environment

- `cryoprotect_supabase_url` - Supabase URL for production
- `cryoprotect_supabase_key` - Supabase key for production
- `cryoprotect_secret_key` - Secret key for session encryption and JWT signing
- `cryoprotect_redis_url` - Redis URL for caching (production only)

### Local Development

For local development, you can use the development profile:

```bash
# Use the development profile
docker-compose --profile dev up
```

This uses environment variables from your `.env` file instead of Docker Secrets and mounts the project directory as a volume for live code changes.

You can also customize the development environment:

```bash
# Set custom development resource limits
export DEV_MEM_LIMIT=4G
export DEV_CPU_LIMIT=2.0

# Set custom development logging options
export DEV_LOG_DRIVER=json-file
export DEV_LOG_MAX_SIZE=50m
export DEV_LOG_MAX_FILE=10

# Start with custom settings
docker-compose --profile dev up
```

## Monitoring and Logging

### Docker Logging

The Docker Compose configuration now includes explicit logging configuration:

- Application logs are sent to stdout/stderr and collected by Docker's logging driver
- The default logging driver is `json-file` with size and file rotation limits
- Logs can be viewed using `docker-compose logs cryoprotect` or `docker-compose logs cryoprotect-dev`
- Log settings can be customized using environment variables (see [Logging Configuration](#logging-configuration))

### Application Monitoring

- Health checks are available at the `/health` endpoint
- For production monitoring, consider setting up Prometheus and Grafana
- Resource usage can be monitored with `docker stats`

### Log Aggregation

For production environments, consider setting up a centralized logging solution:

```bash
# Example: Use fluentd logging driver
export LOG_DRIVER=fluentd
export LOG_OPTIONS="fluentd-address=localhost:24224"
docker-compose --profile prod up -d
```

## Rollback Procedure

If a deployment fails or causes issues:

1. Identify the last stable version (tag or commit SHA)
2. Pull that version's Docker image
3. Ensure Docker secrets are properly configured
4. Deploy the stable version to the affected environment
5. Investigate the issue in the staging environment

```bash
# Example rollback command
ssh user@production-server "
  # Verify Docker secrets exist
  docker secret ls | grep cryoprotect
  
  # Pull the last stable image
  docker pull ghcr.io/yourusername/cryoprotect:previous-stable-tag
  
  # Deploy with Docker secrets
  export FLASK_ENV=production
  export USE_EXTERNAL_SECRETS=true
  # Set appropriate resource limits and logging options
  export MEM_LIMIT=1G
  export CPU_LIMIT=1.0
  # Deploy using the production profile
  docker-compose --profile prod up -d
"

## Blue/Green Deployment

CryoProtect v2 now supports blue/green deployment for zero-downtime updates. This deployment strategy involves running two identical production environments (blue and green) with only one environment active at a time.

### Key Features

- Zero-downtime deployments
- Automated health checking
- Easy rollbacks
- Safe database migrations
- Traffic switching between environments

### Components

- NGINX load balancer for traffic routing
- Blue and green environments running identical application code
- Deployment scripts for automation
- Health check system

### Usage

For detailed instructions on using the blue/green deployment system, see [README_Blue_Green_Deployment.md](./README_Blue_Green_Deployment.md).

#### Quick Start

```bash
# Initialize blue/green deployment
./scripts/init-blue-green.sh --version v1.0.0

# Deploy to green environment
./scripts/deploy-green.sh --version v1.1.0

# Switch traffic to green environment
./scripts/switch-to-green.sh

# Check health of all environments
./scripts/check-health.sh
```

### Integration with CI/CD

The blue/green deployment system integrates with the existing CI/CD pipeline. When a new version is ready for deployment:

1. The CI/CD pipeline builds and tests the application
2. The Docker image is pushed to the registry
3. The deployment script deploys to the inactive environment
4. Health checks verify the deployment
5. Traffic is switched to the new environment

This ensures a smooth, zero-downtime transition between versions.