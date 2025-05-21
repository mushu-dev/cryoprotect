# Phase 3.1 Delivery Plan: CI/CD, Docker & Environment Configuration

## Project Status

Based on the current project state, we have successfully completed:

1. Database Connection Remediation (Adapter Pattern)
2. Local Database Setup
3. Database Utility Functions
4. Script Refactoring for PubChem, ChEMBL and Verification

We now need to implement the remaining tasks in the Deployment Infrastructure Enhancement phase:
- Task 3.1: CI/CD Pipeline Implementation
- Task 3.2: Docker Configuration Optimization
- Task 3.3: Environment Configuration Standardization
- Task 4.1: Implementation Documentation
- Task 4.2: Comprehensive Testing

## Implementation Plan

### Task 3.1: CI/CD Pipeline Implementation

**Description:** Implement a comprehensive CI/CD pipeline using GitHub Actions that automates testing, building, and deployment processes.

**Files to Modify:**
- `.github/workflows/deploy.yml`
- `.github/workflows/ci-cd.yml`

**Implementation Details:**

For `deploy.yml`:
```yaml
name: Deploy CryoProtect

on:
  push:
    branches: [master]
  workflow_dispatch:

jobs:
  test:
    name: Run Tests
    runs-on: ubuntu-latest
    services:
      postgres:
        image: postgres:13
        env:
          POSTGRES_PASSWORD: postgres
          POSTGRES_USER: postgres
          POSTGRES_DB: cryoprotect_test
        ports:
          - 5432:5432
        options: >-
          --health-cmd pg_isready
          --health-interval 10s
          --health-timeout 5s
          --health-retries 5

    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.9'
          
      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install -r requirements.txt
          
      - name: Initialize local database
        run: python database/init_local_db.py
        env:
          LOCAL_DB_HOST: localhost
          LOCAL_DB_USER: postgres
          LOCAL_DB_PASSWORD: postgres
          LOCAL_DB_NAME: cryoprotect_test
          
      - name: Run tests
        run: python -m pytest tests/
        env:
          DB_CONNECTION_MODE: local
          LOCAL_DB_HOST: localhost
          LOCAL_DB_USER: postgres
          LOCAL_DB_PASSWORD: postgres
          LOCAL_DB_NAME: cryoprotect_test

  build:
    name: Build and Push Docker Image
    needs: test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v2
        
      - name: Login to GitHub Container Registry
        uses: docker/login-action@v2
        with:
          registry: ghcr.io
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}
          
      - name: Extract metadata for Docker
        id: meta
        uses: docker/metadata-action@v4
        with:
          images: ghcr.io/${{ github.repository }}
          tags: |
            type=sha,format=long
            type=ref,event=branch
            latest
            
      - name: Build and push
        uses: docker/build-push-action@v4
        with:
          context: .
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=gha
          cache-to: type=gha,mode=max

  deploy-staging:
    name: Deploy to Staging
    needs: build
    runs-on: ubuntu-latest
    environment: staging
    steps:
      - uses: actions/checkout@v3
      
      - name: Install SSH key
        uses: shimataro/ssh-key-action@v2
        with:
          key: ${{ secrets.SSH_PRIVATE_KEY }}
          known_hosts: ${{ secrets.KNOWN_HOSTS }}
          
      - name: Deploy to staging server
        run: |
          ssh ${{ secrets.SSH_USER }}@${{ secrets.SSH_HOST }} "cd ${{ secrets.DEPLOY_PATH }} && \
          docker pull ghcr.io/${{ github.repository }}:latest && \
          docker-compose down && \
          docker-compose up -d"
```

For `ci-cd.yml`:
```yaml
name: CI Checks

on:
  pull_request:
    branches: [master]
  workflow_dispatch:

jobs:
  lint:
    name: Code Quality
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.9'
          
      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install flake8 pylint black mypy bandit safety
          pip install -r requirements.txt
          
      - name: Run flake8
        run: flake8 .
        
      - name: Run black
        run: black --check .
        
      - name: Run pylint
        run: pylint --disable=C0111,C0103,C0303,C0330,C0326 api database chembl pubchem
        
      - name: Run mypy
        run: mypy api database chembl pubchem
        
      - name: Run bandit
        run: bandit -r api database chembl pubchem -f json -o bandit-scan.json
        
      - name: Run safety
        run: safety check -r requirements.txt -o json -f safety-scan.json

  test:
    name: Run Tests
    runs-on: ubuntu-latest
    services:
      postgres:
        image: postgres:13
        env:
          POSTGRES_PASSWORD: postgres
          POSTGRES_USER: postgres
          POSTGRES_DB: cryoprotect_test
        ports:
          - 5432:5432
        options: >-
          --health-cmd pg_isready
          --health-interval 10s
          --health-timeout 5s
          --health-retries 5

    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.9'
          
      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install pytest pytest-cov
          pip install -r requirements.txt
          
      - name: Initialize local database
        run: python database/init_local_db.py
        env:
          LOCAL_DB_HOST: localhost
          LOCAL_DB_USER: postgres
          LOCAL_DB_PASSWORD: postgres
          LOCAL_DB_NAME: cryoprotect_test
          
      - name: Run tests with coverage
        run: python -m pytest tests/ --cov=. --cov-report=xml
        env:
          DB_CONNECTION_MODE: local
          LOCAL_DB_HOST: localhost
          LOCAL_DB_USER: postgres
          LOCAL_DB_PASSWORD: postgres
          LOCAL_DB_NAME: cryoprotect_test
          
      - name: Upload coverage to Codecov
        uses: codecov/codecov-action@v3
        with:
          file: ./coverage.xml
          fail_ci_if_error: true
```

**Success Criteria:**
- Pipeline automatically triggers on push to main branch
- Tests run successfully before deployment
- Docker image is built and pushed to registry
- Deployment to staging environment is automated
- Pull requests trigger CI checks before merging

### Task 3.2: Docker Configuration Optimization

**Description:** Optimize Docker configuration for production deployment, focusing on security, performance, and resource utilization.

**Files to Modify:**
- `Dockerfile`
- `docker-compose.yml`

**Implementation Details:**

For `Dockerfile`:
```dockerfile
# Build stage
FROM python:3.9-slim AS builder

WORKDIR /app

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements file
COPY requirements.txt .

# Install dependencies
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /app/wheels -r requirements.txt

# Final stage
FROM python:3.9-slim

WORKDIR /app

# Create non-root user
RUN groupadd -r cryoprotect && \
    useradd -r -g cryoprotect cryoprotect && \
    mkdir -p /app/data /app/logs /app/cache && \
    chown -R cryoprotect:cryoprotect /app

# Copy wheels from builder stage and install
COPY --from=builder /app/wheels /wheels
RUN pip install --no-cache /wheels/*

# Copy application code
COPY --chown=cryoprotect:cryoprotect . .

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app

# Expose port
EXPOSE 5000

# Switch to non-root user
USER cryoprotect

# Health check
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
    CMD curl -f http://localhost:5000/api/health || exit 1

# Set entrypoint and command
ENTRYPOINT ["./docker-entrypoint.sh"]
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "--workers", "4", "--threads", "2", "app:app"]
```

For `docker-compose.yml`:
```yaml
version: '3.8'

services:
  app:
    image: ghcr.io/organization/cryoprotect:latest
    restart: unless-stopped
    depends_on:
      - postgres
    environment:
      - DB_CONNECTION_MODE=local
      - LOCAL_DB_HOST=postgres
      - LOCAL_DB_PORT=5432
      - LOCAL_DB_NAME=cryoprotect
      - LOCAL_DB_USER=postgres
      - LOCAL_DB_PASSWORD=${POSTGRES_PASSWORD}
    ports:
      - "5000:5000"
    networks:
      - cryoprotect-network
    volumes:
      - app-data:/app/data
      - app-logs:/app/logs
      - app-cache:/app/cache
    deploy:
      resources:
        limits:
          cpus: '1'
          memory: 2G
        reservations:
          cpus: '0.5'
          memory: 1G
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:5000/api/health"]
      interval: 30s
      timeout: 5s
      retries: 3
      start_period: 30s

  postgres:
    image: postgres:13
    restart: unless-stopped
    environment:
      - POSTGRES_PASSWORD=${POSTGRES_PASSWORD}
      - POSTGRES_USER=postgres
      - POSTGRES_DB=cryoprotect
    volumes:
      - postgres-data:/var/lib/postgresql/data
    networks:
      - cryoprotect-network
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U postgres"]
      interval: 10s
      timeout: 5s
      retries: 5
    deploy:
      resources:
        limits:
          cpus: '1'
          memory: 1G
        reservations:
          cpus: '0.25'
          memory: 512M

networks:
  cryoprotect-network:
    driver: bridge

volumes:
  app-data:
  app-logs:
  app-cache:
  postgres-data:
```

**Success Criteria:**
- Docker image size reduced by at least 30%
- No critical or high security vulnerabilities in final image
- Application starts within 5 seconds in container
- Resource limits properly configured for container

### Task 3.3: Environment Configuration Standardization

**Description:** Create standardized environment configurations for development, staging, and production environments.

**Files to Create/Modify:**
- `config.py` (update)
- `.env.template` (update)
- `.env.production` (create)
- `.env.staging` (create)

**Implementation Details:**

For `config.py`:
```python
#!/usr/bin/env python3
"""
Configuration management for CryoProtect v2.

This module handles environment-specific configuration settings and
provides a unified interface for accessing configuration values.
"""

import os
import logging
from typing import Any, Dict, Optional
from dotenv import load_dotenv

logger = logging.getLogger(__name__)

# Load environment variables from .env file
load_dotenv()

class Config:
    """Base configuration class with common settings."""
    
    # Application settings
    APP_NAME = 'CryoProtect v2'
    DEBUG = False
    TESTING = False
    LOG_LEVEL = logging.INFO
    
    # Security settings
    SECRET_KEY = os.getenv('SECRET_KEY', 'development-key-not-secure')
    JWT_SECRET_KEY = os.getenv('JWT_SECRET_KEY', 'development-jwt-key-not-secure')
    JWT_ACCESS_TOKEN_EXPIRES = 3600  # 1 hour
    JWT_REFRESH_TOKEN_EXPIRES = 2592000  # 30 days
    
    # Database settings
    DB_CONNECTION_MODE = os.getenv('DB_CONNECTION_MODE', 'auto')
    
    # API settings
    RATE_LIMIT = os.getenv('RATE_LIMIT', '100/hour')
    
    @classmethod
    def get_database_config(cls) -> Dict[str, Any]:
        """Get database configuration based on connection mode."""
        load_environment_variables()
        
        connection_mode = os.getenv('DB_CONNECTION_MODE', 'auto').lower()
        
        if connection_mode == 'local' or connection_mode == 'auto':
            return {
                'local': {
                    'host': os.getenv('LOCAL_DB_HOST', 'localhost'),
                    'port': os.getenv('LOCAL_DB_PORT', '5432'),
                    'database': os.getenv('LOCAL_DB_NAME', 'cryoprotect'),
                    'user': os.getenv('LOCAL_DB_USER', 'postgres'),
                    'password': os.getenv('LOCAL_DB_PASSWORD', ''),
                    'min_connections': int(os.getenv('LOCAL_DB_MIN_CONNECTIONS', '1')),
                    'max_connections': int(os.getenv('LOCAL_DB_MAX_CONNECTIONS', '5'))
                }
            }
        
        if connection_mode == 'supabase' or connection_mode == 'auto':
            # Ensure we have the IP address resolved if available
            ip_address = os.getenv('SUPABASE_DB_IP_ADDRESS')
            
            return {
                'supabase': {
                    'host': os.getenv('SUPABASE_DB_HOST'),
                    'port': os.getenv('SUPABASE_DB_PORT', '5432'),
                    'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
                    'user': os.getenv('SUPABASE_DB_USER'),
                    'password': os.getenv('SUPABASE_DB_PASSWORD'),
                    'ip_address': ip_address,
                    'min_connections': int(os.getenv('SUPABASE_DB_MIN_CONNECTIONS', '1')),
                    'max_connections': int(os.getenv('SUPABASE_DB_MAX_CONNECTIONS', '10'))
                }
            }
        
        if connection_mode == 'mcp' or connection_mode == 'auto':
            return {
                'mcp': {
                    'project_id': os.getenv('SUPABASE_PROJECT_ID')
                }
            }
        
        # Default to empty config if no valid connection mode
        return {}

class DevelopmentConfig(Config):
    """Development environment configuration."""
    
    DEBUG = True
    LOG_LEVEL = logging.DEBUG

class TestingConfig(Config):
    """Testing environment configuration."""
    
    TESTING = True
    DEBUG = True
    LOG_LEVEL = logging.DEBUG
    
    # Use local database for testing
    DB_CONNECTION_MODE = 'local'

class StagingConfig(Config):
    """Staging environment configuration."""
    
    # Default to Supabase with local fallback
    DB_CONNECTION_MODE = 'auto'

class ProductionConfig(Config):
    """Production environment configuration."""
    
    # Stricter security settings
    JWT_ACCESS_TOKEN_EXPIRES = 1800  # 30 minutes
    
    # Higher rate limits for production
    RATE_LIMIT = os.getenv('RATE_LIMIT', '1000/hour')
    
    # Default to Supabase with MCP fallback
    DB_CONNECTION_MODE = 'auto'

def load_environment_variables() -> None:
    """
    Load and normalize environment variables.
    Ensures variables are consistently available regardless of naming convention.
    """
    # DB_* and SUPABASE_DB_* variables normalization
    db_vars = {
        'HOST': os.getenv('DB_HOST') or os.getenv('SUPABASE_DB_HOST'),
        'PORT': os.getenv('DB_PORT') or os.getenv('SUPABASE_DB_PORT') or '5432',
        'NAME': os.getenv('DB_NAME') or os.getenv('SUPABASE_DB_NAME') or 'postgres',
        'USER': os.getenv('DB_USER') or os.getenv('SUPABASE_DB_USER'),
        'PASSWORD': os.getenv('DB_PASSWORD') or os.getenv('SUPABASE_DB_PASSWORD'),
        'IP_ADDRESS': os.getenv('DB_IP_ADDRESS') or os.getenv('SUPABASE_DB_IP_ADDRESS'),
        'MIN_CONNECTIONS': os.getenv('DB_MIN_CONNECTIONS') or os.getenv('SUPABASE_DB_MIN_CONNECTIONS') or '1',
        'MAX_CONNECTIONS': os.getenv('DB_MAX_CONNECTIONS') or os.getenv('SUPABASE_DB_MAX_CONNECTIONS') or '10'
    }
    
    # Set both DB_* and SUPABASE_DB_* variables
    for key, value in db_vars.items():
        if value:
            if not os.getenv(f'DB_{key}'):
                os.environ[f'DB_{key}'] = value
            if not os.getenv(f'SUPABASE_DB_{key}'):
                os.environ[f'SUPABASE_DB_{key}'] = value

def get_config() -> Config:
    """
    Get the appropriate configuration based on environment.
    
    Returns:
        Config: Configuration object for the current environment
    """
    env = os.getenv('FLASK_ENV', 'development').lower()
    
    if env == 'production':
        return ProductionConfig()
    elif env == 'testing':
        return TestingConfig()
    elif env == 'staging':
        return StagingConfig()
    else:
        return DevelopmentConfig()

# Load environment variables when module is imported
load_environment_variables()
```

For `.env.template`:
```
# Application Settings
FLASK_ENV=development  # development, testing, staging, production
SECRET_KEY=your-secret-key-here
JWT_SECRET_KEY=your-jwt-secret-key-here

# Database Connection Mode
# Options: auto, local, supabase, mcp
DB_CONNECTION_MODE=auto

# Local Database Settings
LOCAL_DB_HOST=localhost
LOCAL_DB_PORT=5432
LOCAL_DB_NAME=cryoprotect
LOCAL_DB_USER=postgres
LOCAL_DB_PASSWORD=postgres
LOCAL_DB_MIN_CONNECTIONS=1
LOCAL_DB_MAX_CONNECTIONS=5

# Supabase Database Settings
SUPABASE_URL=https://xxxxxxxxxxxxxxxxxxxx.supabase.co
SUPABASE_KEY=your-anon-key
SUPABASE_SERVICE_KEY=your-service-role-key
SUPABASE_PROJECT_ID=your-project-id
SUPABASE_DB_HOST=db.xxxxxxxxxxxxxxxxxxxx.supabase.co
SUPABASE_DB_PORT=5432
SUPABASE_DB_NAME=postgres
SUPABASE_DB_USER=postgres
SUPABASE_DB_PASSWORD=your-db-password
SUPABASE_DB_IP_ADDRESS=  # Optional, for DNS fallback

# API Settings
RATE_LIMIT=100/hour
```

For `.env.production`:
```
# Production Environment Settings
FLASK_ENV=production
SECRET_KEY=generate-a-secure-random-key-here
JWT_SECRET_KEY=generate-another-secure-random-key-here

# Database Connection Mode
# Using auto for Supabase with MCP fallback
DB_CONNECTION_MODE=auto

# Supabase Database Settings
SUPABASE_URL=https://xxxxxxxxxxxxxxxxxxxx.supabase.co
SUPABASE_KEY=production-anon-key
SUPABASE_SERVICE_KEY=production-service-role-key
SUPABASE_PROJECT_ID=production-project-id
SUPABASE_DB_HOST=db.xxxxxxxxxxxxxxxxxxxx.supabase.co
SUPABASE_DB_PORT=5432
SUPABASE_DB_NAME=postgres
SUPABASE_DB_USER=postgres
SUPABASE_DB_PASSWORD=production-db-password
SUPABASE_DB_MIN_CONNECTIONS=5
SUPABASE_DB_MAX_CONNECTIONS=20

# API Settings
RATE_LIMIT=1000/hour

# Logging Settings
LOG_LEVEL=INFO
```

For `.env.staging`:
```
# Staging Environment Settings
FLASK_ENV=staging
SECRET_KEY=staging-secure-key
JWT_SECRET_KEY=staging-jwt-secure-key

# Database Connection Mode
# Using auto for Supabase with local fallback
DB_CONNECTION_MODE=auto

# Local Database Settings (Fallback)
LOCAL_DB_HOST=localhost
LOCAL_DB_PORT=5432
LOCAL_DB_NAME=cryoprotect_staging
LOCAL_DB_USER=postgres
LOCAL_DB_PASSWORD=staging-postgres-password
LOCAL_DB_MIN_CONNECTIONS=1
LOCAL_DB_MAX_CONNECTIONS=10

# Supabase Database Settings
SUPABASE_URL=https://xxxxxxxxxxxxxxxxxxxx.supabase.co
SUPABASE_KEY=staging-anon-key
SUPABASE_SERVICE_KEY=staging-service-role-key
SUPABASE_PROJECT_ID=staging-project-id
SUPABASE_DB_HOST=db.xxxxxxxxxxxxxxxxxxxx.supabase.co
SUPABASE_DB_PORT=5432
SUPABASE_DB_NAME=postgres
SUPABASE_DB_USER=postgres
SUPABASE_DB_PASSWORD=staging-db-password
SUPABASE_DB_MIN_CONNECTIONS=2
SUPABASE_DB_MAX_CONNECTIONS=15

# API Settings
RATE_LIMIT=500/hour

# Logging Settings
LOG_LEVEL=DEBUG
```

**Success Criteria:**
- All environment-specific values are externalized to configuration
- Secrets are properly managed and not exposed in code
- Configuration validation prevents startup with invalid settings
- Easy switching between environments during development

### Task 4.1: Implementation Documentation

**Description:** Create comprehensive documentation for the implementation.

**Files to Create:**
- `docs/deployment_guide.md`
- `docs/configuration_guide.md`
- `docs/database_guide.md`

**Implementation Details:**

For `docs/deployment_guide.md`:
```markdown
# Deployment Guide

This document provides comprehensive instructions for deploying CryoProtect v2 to various environments.

## Prerequisites

- Docker and Docker Compose installed
- Access to GitHub Container Registry (ghcr.io)
- PostgreSQL database (local deployment) or Supabase credentials (cloud deployment)
- SSH access to deployment server (for manual deployment)

## Environment Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/organization/cryoprotect.git
   cd cryoprotect
   ```

2. Create environment configuration:
   ```bash
   # For development
   cp .env.template .env
   
   # For staging
   cp .env.staging .env
   
   # For production
   cp .env.production .env
   ```

3. Edit the `.env` file with appropriate values for your environment.

## Local Deployment

1. Initialize the local database:
   ```bash
   python database/init_local_db.py
   ```

2. Run the application using Docker Compose:
   ```bash
   docker-compose up -d
   ```

3. Verify the deployment:
   ```bash
   curl http://localhost:5000/api/health
   ```

## Cloud Deployment

### Automated Deployment via CI/CD

1. Push changes to the `master` branch to trigger the deployment pipeline.
2. The CI/CD pipeline will:
   - Run tests
   - Build and push a Docker image
   - Deploy to the staging environment
   - (Manual approval required for production deployment)

### Manual Deployment

1. Build the Docker image:
   ```bash
   docker build -t cryoprotect:latest .
   ```

2. Push the image to a registry:
   ```bash
   docker tag cryoprotect:latest ghcr.io/organization/cryoprotect:latest
   docker push ghcr.io/organization/cryoprotect:latest
   ```

3. SSH into the deployment server:
   ```bash
   ssh user@server
   ```

4. Deploy using Docker Compose:
   ```bash
   cd /path/to/deployment
   docker-compose pull
   docker-compose down
   docker-compose up -d
   ```

## Blue-Green Deployment

For zero-downtime production deployments, we use a blue-green deployment strategy:

1. Set up the environment:
   ```bash
   ./scripts/init-blue-green.sh
   ```

2. Deploy to the inactive environment:
   ```bash
   # If blue is active, deploy to green
   ./scripts/deploy-green.sh
   
   # If green is active, deploy to blue
   ./scripts/deploy-blue.sh
   ```

3. Switch traffic after verifying the deployment:
   ```bash
   # Switch to blue
   ./scripts/switch-to-blue.sh
   
   # Switch to green
   ./scripts/switch-to-green.sh
   ```

4. In case of issues, rollback:
   ```bash
   ./scripts/rollback.sh
   ```

## Database Migration

1. Apply migrations:
   ```bash
   ./scripts/migrate-database.sh
   ```

## Health Checks

1. Check application health:
   ```bash
   ./scripts/check-health.sh
   ```

## Troubleshooting

### Common Issues

1. **Database Connection Failure**
   - Check the database credentials in the `.env` file
   - Verify network connectivity to the database
   - Check if the database service is running

2. **Image Pull Failure**
   - Verify Docker registry credentials
   - Check network connectivity to the registry

3. **Application Not Starting**
   - Check application logs: `docker-compose logs app`
   - Verify environment variables
   - Check resource availability on the host

For additional help, refer to the [Troubleshooting Guide](./troubleshooting_guide.md).
```

For `docs/configuration_guide.md`:
```markdown
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
```
```

For `docs/database_guide.md`:
```markdown
# Database Guide

This document describes the database architecture and operations for CryoProtect v2.

## Database Architecture

CryoProtect v2 uses a flexible database approach that supports:

1. **Local PostgreSQL**: For development and testing
2. **Supabase PostgreSQL**: For staging and production
3. **MCP API**: As a fallback mechanism

### Schema Overview

The core tables in the database are:

- `molecules`: Basic molecular information
- `molecular_properties`: Properties of molecules in a normalized structure
- `property_types`: Types of molecular properties
- `mixtures`: Cryoprotectant mixtures
- `mixture_components`: Components of mixtures with concentrations

## Database Adapters

The application uses a database adapter pattern to abstract connection details:

1. **DatabaseAdapter**: Abstract interface for all adapters
2. **LocalPostgreSQLAdapter**: Connects to local PostgreSQL
3. **SupabaseDirectAdapter**: Connects directly to Supabase PostgreSQL
4. **MCPAdapter**: Connects via MCP API

## Connection Management

The `ConnectionManager` handles database connections with features:

1. **Connection Pooling**: Efficiently manages database connections
2. **Fallback Mechanism**: Automatically tries alternative connection methods
3. **Transaction Support**: Manages database transactions
4. **Retry Logic**: Handles transient errors

## Local Database Setup

To set up the local database:

1. Install PostgreSQL
2. Configure environment variables:
   ```
   LOCAL_DB_HOST=localhost
   LOCAL_DB_PORT=5432
   LOCAL_DB_NAME=cryoprotect
   LOCAL_DB_USER=postgres
   LOCAL_DB_PASSWORD=your-password
   ```
3. Run the initialization script:
   ```bash
   python database/init_local_db.py
   ```

## Database Migrations

Database schema changes are managed through migrations:

1. Migration files are stored in `migrations/`
2. The format is `NNN_description.sql`
3. Migrations are applied in numerical order

To apply migrations:

```bash
python database/migrations/runner.py
```

## Data Population

To populate the database with scientific data:

1. Configure database connection
2. Run the population scripts:
   ```bash
   # For reference compounds
   python import_reference_compounds.py
   
   # For ChEMBL data
   python import_full_chembl.py
   
   # For PubChem data
   python PubChem_CryoProtectants_Supabase.py
   
   # For property enhancement
   python enhance_pubchem_properties.py
   ```

## Verification

To verify the database population:

```bash
python verify_imported_data.py
```

This will check:
- Molecule counts (>=500 expected)
- Property completeness
- Reference compound presence
- Cross-referencing
- Query performance

## Backup and Restore

### Creating a Backup

```bash
# Full backup
python database/utils/backup.py --full

# Partial backup (molecules only)
python database/utils/backup.py --tables molecules,molecular_properties
```

### Restoring from Backup

```bash
python database/utils/restore.py --file backup_20250430_120000.sql
```

## Performance Optimization

1. **Indexes**: Added to frequently queried columns
2. **Connection Pooling**: Reduces connection overhead
3. **Query Optimization**: Complex queries are optimized

## Troubleshooting

### Common Issues

1. **Connection Failures**
   - Check credentials in `.env`
   - Verify network connectivity
   - Check if PostgreSQL is running

2. **DNS Resolution Issues**
   - Use IP address fallback: Set `SUPABASE_DB_IP_ADDRESS`

3. **Performance Issues**
   - Check query plans: `EXPLAIN ANALYZE SELECT...`
   - Verify indexes

## Database Utilities

The `database/utils.py` module provides helper functions:

- `execute_query`: Execute SQL query
- `execute_batch`: Execute multiple queries
- `get_molecule_by_id`: Get molecule by ID
- `get_molecule_properties`: Get properties for a molecule
- `insert_molecule`: Insert new molecule
- `set_molecule_property`: Set molecule property
- `test_database_connection`: Test connection

## Best Practices

1. Use connection pool for multiple operations
2. Use transactions for related operations
3. Handle connection errors gracefully
4. Use parameterized queries to prevent SQL injection
5. Close connections when done
```

**Success Criteria:**
- Architecture clearly documented
- Setup procedures well explained
- Troubleshooting steps are actionable
- Documentation is up-to-date with code

### Task 4.2: Comprehensive Testing

**Description:** Create comprehensive tests for all components of the system.

**Files to Create/Modify:**
- `tests/test_database_adapters.py`
- `tests/test_connection_manager.py`
- `tests/test_database_utils.py`
- `tests/test_environment_config.py`
- `tests/test_integration.py`

**Implementation Details:**

For `tests/test_database_adapters.py`:
```python
#!/usr/bin/env python3
"""
Tests for database adapters.
"""

import os
import unittest
from unittest.mock import patch, MagicMock
import psycopg2

from database.adapter import DatabaseAdapter
from database.local_adapter import LocalPostgreSQLAdapter
from database.supabase_adapter import SupabaseDirectAdapter
from database.mcp_adapter import MCPAdapter

class TestLocalPostgreSQLAdapter(unittest.TestCase):
    """Test the LocalPostgreSQLAdapter."""
    
    def setUp(self):
        """Set up the test environment."""
        self.config = {
            'host': 'localhost',
            'port': 5432,
            'database': 'test_db',
            'user': 'postgres',
            'password': 'postgres',
            'min_connections': 1,
            'max_connections': 5
        }
        
        # Create adapter
        self.adapter = LocalPostgreSQLAdapter(self.config)
        
        # Mock the connection pool
        self.adapter.connection_pool = MagicMock()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_connect(self, mock_pool):
        """Test the connect method."""
        # Set up mock
        mock_pool.return_value = MagicMock()
        
        # Test connect
        result = self.adapter.connect()
        
        # Verify
        self.assertTrue(result)
        mock_pool.assert_called_once_with(
            minconn=1,
            maxconn=5,
            host='localhost',
            port=5432,
            dbname='test_db',
            user='postgres',
            password='postgres'
        )
    
    def test_disconnect(self):
        """Test the disconnect method."""
        # Test disconnect
        result = self.adapter.disconnect()
        
        # Verify
        self.assertTrue(result)
        self.adapter.connection_pool.closeall.assert_called_once()
    
    def test_execute_query(self):
        """Test the execute_query method."""
        # Set up mocks
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.fetchall.return_value = [{'id': 1, 'name': 'Test'}]
        mock_conn.cursor.return_value = mock_cursor
        self.adapter.connection_pool.getconn.return_value = mock_conn
        
        # Test execute_query
        result = self.adapter.execute_query("SELECT * FROM test")
        
        # Verify
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])
        mock_cursor.execute.assert_called_once_with("SELECT * FROM test", None)
        self.adapter.connection_pool.getconn.assert_called_once()
        self.adapter.connection_pool.putconn.assert_called_once_with(mock_conn)
    
    def test_execute_batch(self):
        """Test the execute_batch method."""
        # Set up mocks
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.fetchall.return_value = [{'id': 1, 'name': 'Test'}]
        mock_conn.cursor.return_value = mock_cursor
        self.adapter.connection_pool.getconn.return_value = mock_conn
        
        # Test execute_batch
        result = self.adapter.execute_batch(["SELECT * FROM test", "SELECT * FROM test2"])
        
        # Verify
        self.assertEqual(len(result), 2)
        self.assertEqual(mock_cursor.execute.call_count, 2)
        self.adapter.connection_pool.getconn.assert_called_once()
        self.adapter.connection_pool.putconn.assert_called_once_with(mock_conn)

class TestSupabaseDirectAdapter(unittest.TestCase):
    """Test the SupabaseDirectAdapter."""
    
    def setUp(self):
        """Set up the test environment."""
        self.config = {
            'host': 'db.example.supabase.co',
            'port': 5432,
            'database': 'postgres',
            'user': 'postgres',
            'password': 'password',
            'min_connections': 1,
            'max_connections': 5
        }
        
        # Create adapter
        self.adapter = SupabaseDirectAdapter(self.config)
        
        # Mock the connection pool
        self.adapter.connection_pool = MagicMock()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    @patch('database.supabase_adapter.SupabaseDirectAdapter._resolve_hostname')
    def test_connect_direct(self, mock_resolve, mock_pool):
        """Test the connect method with direct hostname connection."""
        # Set up mocks
        mock_pool.return_value = MagicMock()
        
        # Test connect
        result = self.adapter.connect()
        
        # Verify
        self.assertTrue(result)
        mock_pool.assert_called_once_with(
            minconn=1,
            maxconn=5,
            host='db.example.supabase.co',
            port=5432,
            dbname='postgres',
            user='postgres',
            password='password'
        )
        mock_resolve.assert_not_called()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    @patch('database.supabase_adapter.SupabaseDirectAdapter._resolve_hostname')
    def test_connect_fallback(self, mock_resolve, mock_pool):
        """Test the connect method with IP fallback."""
        # Set up mocks
        mock_pool.side_effect = [psycopg2.OperationalError("Test error"), MagicMock()]
        mock_resolve.return_value = "192.168.1.1"
        
        # Test connect
        result = self.adapter.connect()
        
        # Verify
        self.assertTrue(result)
        self.assertEqual(mock_pool.call_count, 2)
        self.assertEqual(self.adapter.ip_address, "192.168.1.1")
        mock_resolve.assert_called_once_with('db.example.supabase.co')

class TestMCPAdapter(unittest.TestCase):
    """Test the MCPAdapter."""
    
    def setUp(self):
        """Set up the test environment."""
        self.config = {
            'project_id': 'test-project'
        }
        
        # Create adapter
        with patch('database.mcp_adapter.MCPAdapter.execute_sql_through_mcp', MagicMock()):
            with patch('database.mcp_adapter.MCPAdapter.get_project_id_for_mcp', MagicMock(return_value='test-project')):
                self.adapter = MCPAdapter(self.config)
    
    @patch('database.mcp_adapter.MCPAdapter.execute_sql_through_mcp')
    def test_connect(self, mock_execute):
        """Test the connect method."""
        # Set up mock
        mock_execute.return_value = [{'test': 1}]
        
        # Test connect
        result = self.adapter.connect()
        
        # Verify
        self.assertTrue(result)
        mock_execute.assert_called_once_with("SELECT 1 as test", 'test-project')
        self.assertTrue(self.adapter.connected)
    
    def test_disconnect(self):
        """Test the disconnect method."""
        # Test disconnect
        self.adapter.connected = True
        result = self.adapter.disconnect()
        
        # Verify
        self.assertTrue(result)
        self.assertFalse(self.adapter.connected)
    
    @patch('database.mcp_adapter.MCPAdapter.execute_sql_through_mcp')
    def test_execute_query(self, mock_execute):
        """Test the execute_query method."""
        # Set up mock
        mock_execute.return_value = [{'id': 1, 'name': 'Test'}]
        
        # Test execute_query
        result = self.adapter.execute_query("SELECT * FROM test")
        
        # Verify
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])
        mock_execute.assert_called_once_with("SELECT * FROM test", 'test-project')

if __name__ == '__main__':
    unittest.main()
```

For `tests/test_connection_manager.py`:
```python
#!/usr/bin/env python3
"""
Tests for connection manager.
"""

import os
import unittest
from unittest.mock import patch, MagicMock

from database.connection_manager import ConnectionManager
from database.adapter import DatabaseAdapter
from database.local_adapter import LocalPostgreSQLAdapter
from database.supabase_adapter import SupabaseDirectAdapter
from database.mcp_adapter import MCPAdapter

class TestConnectionManager(unittest.TestCase):
    """Test the ConnectionManager."""
    
    def setUp(self):
        """Set up the test environment."""
        # Save original environment
        self.original_env = dict(os.environ)
        
        # Set test environment variables
        os.environ['DB_CONNECTION_MODE'] = 'auto'
        os.environ['LOCAL_DB_HOST'] = 'localhost'
        os.environ['LOCAL_DB_USER'] = 'postgres'
        os.environ['LOCAL_DB_PASSWORD'] = 'postgres'
        os.environ['SUPABASE_DB_HOST'] = 'db.example.supabase.co'
        os.environ['SUPABASE_DB_USER'] = 'postgres'
        os.environ['SUPABASE_DB_PASSWORD'] = 'password'
        os.environ['SUPABASE_PROJECT_ID'] = 'test-project'
        
        # Create connection manager
        with patch('database.connection_manager.LocalPostgreSQLAdapter', MagicMock()):
            with patch('database.connection_manager.SupabaseDirectAdapter', MagicMock()):
                with patch('database.connection_manager.MCPAdapter', MagicMock()):
                    self.manager = ConnectionManager()
                    
                    # Set up mock adapters
                    self.mock_local = MagicMock(spec=LocalPostgreSQLAdapter)
                    self.mock_supabase = MagicMock(spec=SupabaseDirectAdapter)
                    self.mock_mcp = MagicMock(spec=MCPAdapter)
                    
                    self.manager.adapters = {
                        'local': self.mock_local,
                        'supabase': self.mock_supabase,
                        'mcp': self.mock_mcp
                    }
    
    def tearDown(self):
        """Tear down the test environment."""
        # Restore original environment
        os.environ.clear()
        os.environ.update(self.original_env)
    
    def test_connect_primary_success(self):
        """Test connect when primary adapter succeeds."""
        # Set up mocks
        self.mock_supabase.connect.return_value = True
        
        # Test connect
        result = self.manager.connect()
        
        # Verify
        self.assertTrue(result)
        self.assertEqual(self.manager.active_adapter, 'supabase')
        self.mock_supabase.connect.assert_called_once()
        self.mock_local.connect.assert_not_called()
        self.mock_mcp.connect.assert_not_called()
    
    def test_connect_primary_fail_fallback_success(self):
        """Test connect when primary adapter fails but fallback succeeds."""
        # Set up mocks
        self.mock_supabase.connect.return_value = False
        self.mock_local.connect.return_value = True
        
        # Test connect
        result = self.manager.connect()
        
        # Verify
        self.assertTrue(result)
        self.assertEqual(self.manager.active_adapter, 'local')
        self.mock_supabase.connect.assert_called_once()
        self.mock_local.connect.assert_called_once()
        self.mock_mcp.connect.assert_not_called()
    
    def test_connect_all_fail(self):
        """Test connect when all adapters fail."""
        # Set up mocks
        self.mock_supabase.connect.return_value = False
        self.mock_local.connect.return_value = False
        self.mock_mcp.connect.return_value = False
        
        # Test connect
        result = self.manager.connect()
        
        # Verify
        self.assertFalse(result)
        self.assertIsNone(self.manager.active_adapter)
        self.mock_supabase.connect.assert_called_once()
        self.mock_local.connect.assert_called_once()
        self.mock_mcp.connect.assert_called_once()
    
    def test_disconnect(self):
        """Test disconnect."""
        # Set up state
        self.manager.active_adapter = 'supabase'
        
        # Test disconnect
        result = self.manager.disconnect()
        
        # Verify
        self.assertTrue(result)
        self.assertIsNone(self.manager.active_adapter)
        self.mock_supabase.disconnect.assert_called_once()
        self.mock_local.disconnect.assert_called_once()
        self.mock_mcp.disconnect.assert_called_once()
    
    def test_execute_query(self):
        """Test execute_query."""
        # Set up mocks
        self.manager.active_adapter = 'supabase'
        self.mock_supabase.execute_query.return_value = [{'id': 1, 'name': 'Test'}]
        
        # Test execute_query
        result = self.manager.execute_query("SELECT * FROM test")
        
        # Verify
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])
        self.mock_supabase.execute_query.assert_called_once_with("SELECT * FROM test", None)
    
    def test_execute_query_not_connected(self):
        """Test execute_query when not connected."""
        # Set up mocks
        self.manager.active_adapter = None
        self.mock_supabase.connect.return_value = True
        
        # Test execute_query
        with self.assertRaises(ConnectionError):
            # Set up to fail connect
            self.mock_supabase.connect.return_value = False
            self.mock_local.connect.return_value = False
            self.mock_mcp.connect.return_value = False
            
            self.manager.execute_query("SELECT * FROM test")
    
    def test_get_connection_info(self):
        """Test get_connection_info."""
        # Set up mocks
        self.manager.active_adapter = 'supabase'
        self.mock_supabase.get_connection_info.return_value = {'type': 'supabase', 'host': 'db.example.supabase.co'}
        self.mock_local.get_connection_info.return_value = {'type': 'local', 'host': 'localhost'}
        self.mock_mcp.get_connection_info.return_value = {'type': 'mcp', 'project_id': 'test-project'}
        
        # Test get_connection_info
        result = self.manager.get_connection_info()
        
        # Verify
        self.assertEqual(result['connection_mode'], 'auto')
        self.assertEqual(result['primary_adapter'], 'supabase')
        self.assertEqual(result['active_adapter'], 'supabase')
        self.assertEqual(result['adapters']['supabase'], {'type': 'supabase', 'host': 'db.example.supabase.co'})
        self.assertEqual(result['adapters']['local'], {'type': 'local', 'host': 'localhost'})
        self.assertEqual(result['adapters']['mcp'], {'type': 'mcp', 'project_id': 'test-project'})

if __name__ == '__main__':
    unittest.main()
```

For `tests/test_database_utils.py`:
```python
#!/usr/bin/env python3
"""
Tests for database utilities.
"""

import unittest
from unittest.mock import patch, MagicMock

from database.utils import (
    get_db, with_connection, with_retry, with_transaction,
    execute_query, execute_batch, get_molecule_by_id,
    get_molecule_properties, get_molecules_by_inchikey,
    insert_molecule, update_molecule, set_molecule_property,
    get_or_create_property_type, test_database_connection
)

class TestDatabaseUtils(unittest.TestCase):
    """Test the database utility functions."""
    
    def setUp(self):
        """Set up the test environment."""
        # Mock connection manager
        self.mock_manager = MagicMock()
        self.mock_manager.get_active_adapter.return_value = True
        self.mock_manager.connect.return_value = True
        
        # Create patch for get_db
        self.get_db_patcher = patch('database.utils.get_db', return_value=self.mock_manager)
        self.mock_get_db = self.get_db_patcher.start()
    
    def tearDown(self):
        """Tear down the test environment."""
        self.get_db_patcher.stop()
    
    def test_with_connection_decorator(self):
        """Test the with_connection decorator."""
        # Define test function
        @with_connection
        def test_func():
            return "success"
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, "success")
        self.mock_manager.get_active_adapter.assert_called_once()
        self.mock_manager.connect.assert_not_called()  # Already connected
    
    def test_with_connection_decorator_not_connected(self):
        """Test the with_connection decorator when not connected."""
        # Set up mock
        self.mock_manager.get_active_adapter.return_value = None
        
        # Define test function
        @with_connection
        def test_func():
            return "success"
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, "success")
        self.mock_manager.get_active_adapter.assert_called_once()
        self.mock_manager.connect.assert_called_once()
    
    def test_with_retry_decorator_success(self):
        """Test the with_retry decorator with successful execution."""
        # Define test function
        @with_retry(max_retries=3)
        def test_func():
            return "success"
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, "success")
    
    def test_with_retry_decorator_retry_success(self):
        """Test the with_retry decorator with retry success."""
        # Define test counter
        counter = [0]
        
        # Define test function
        @with_retry(max_retries=3)
        def test_func():
            counter[0] += 1
            if counter[0] < 3:
                raise ValueError("Test error")
            return "success"
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, "success")
        self.assertEqual(counter[0], 3)
    
    def test_with_retry_decorator_all_fail(self):
        """Test the with_retry decorator with all retries failing."""
        # Define test function
        @with_retry(max_retries=3)
        def test_func():
            raise ValueError("Test error")
        
        # Test function
        with self.assertRaises(ValueError):
            test_func()
    
    def test_with_transaction_decorator(self):
        """Test the with_transaction decorator."""
        # Set up mocks
        mock_transaction = MagicMock()
        self.mock_manager.begin_transaction.return_value = mock_transaction
        
        # Define test function
        @with_transaction
        def test_func(transaction=None):
            return transaction
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, mock_transaction)
        self.mock_manager.begin_transaction.assert_called_once()
        self.mock_manager.commit_transaction.assert_called_once_with(mock_transaction)
        self.mock_manager.rollback_transaction.assert_not_called()
    
    def test_with_transaction_decorator_exception(self):
        """Test the with_transaction decorator with exception."""
        # Set up mocks
        mock_transaction = MagicMock()
        self.mock_manager.begin_transaction.return_value = mock_transaction
        
        # Define test function
        @with_transaction
        def test_func(transaction=None):
            raise ValueError("Test error")
        
        # Test function
        with self.assertRaises(ValueError):
            test_func()
        
        # Verify
        self.mock_manager.begin_transaction.assert_called_once()
        self.mock_manager.commit_transaction.assert_not_called()
        self.mock_manager.rollback_transaction.assert_called_once_with(mock_transaction)
    
    def test_execute_query(self):
        """Test execute_query."""
        # Set up mock
        self.mock_manager.execute_query.return_value = [{'id': 1, 'name': 'Test'}]
        
        # Test function
        result = execute_query("SELECT * FROM test")
        
        # Verify
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])
        self.mock_manager.execute_query.assert_called_once_with("SELECT * FROM test", None)
    
    def test_get_molecule_by_id(self):
        """Test get_molecule_by_id."""
        # Set up mock
        self.mock_manager.execute_query.return_value = [{'id': '1', 'name': 'Test Molecule'}]
        
        # Test function
        result = get_molecule_by_id("1")
        
        # Verify
        self.assertEqual(result, {'id': '1', 'name': 'Test Molecule'})
        self.mock_manager.execute_query.assert_called_once()
    
    def test_insert_molecule(self):
        """Test insert_molecule."""
        # Set up mock
        self.mock_manager.execute_query.return_value = [{'id': '1', 'name': 'Test Molecule'}]
        
        # Test function
        result = insert_molecule(
            name="Test Molecule",
            formula="C10H20O",
            molecular_weight=156.27,
            smiles="CCCCCCCCCC(=O)",
            inchi="InChI=1S/C10H20O/c1-2-3-4-5-6-7-8-9-10-11/h10H,2-9H2,1H3",
            inchi_key="AXSIIYDNMFQWRZ-UHFFFAOYSA-N",
            chembl_id="CHEMBL1234",
            pubchem_cid="123456",
            data_source="test"
        )
        
        # Verify
        self.assertEqual(result, {'id': '1', 'name': 'Test Molecule'})
        self.mock_manager.execute_query.assert_called_once()

if __name__ == '__main__':
    unittest.main()
```

For `tests/test_environment_config.py`:
```python
#!/usr/bin/env python3
"""
Tests for environment configuration.
"""

import os
import unittest
from unittest.mock import patch

from config import (
    Config, DevelopmentConfig, TestingConfig, StagingConfig, ProductionConfig,
    load_environment_variables, get_config
)

class TestConfig(unittest.TestCase):
    """Test the Config classes."""
    
    def setUp(self):
        """Set up the test environment."""
        # Save original environment
        self.original_env = dict(os.environ)
        
        # Clear environment variables for testing
        os.environ.clear()
    
    def tearDown(self):
        """Tear down the test environment."""
        # Restore original environment
        os.environ.clear()
        os.environ.update(self.original_env)
    
    def test_base_config(self):
        """Test the base Config class."""
        config = Config()
        
        # Verify default values
        self.assertEqual(config.APP_NAME, 'CryoProtect v2')
        self.assertFalse(config.DEBUG)
        self.assertFalse(config.TESTING)
        self.assertEqual(config.DB_CONNECTION_MODE, 'auto')
    
    def test_development_config(self):
        """Test the DevelopmentConfig class."""
        config = DevelopmentConfig()
        
        # Verify development-specific values
        self.assertTrue(config.DEBUG)
        self.assertFalse(config.TESTING)
    
    def test_testing_config(self):
        """Test the TestingConfig class."""
        config = TestingConfig()
        
        # Verify testing-specific values
        self.assertTrue(config.DEBUG)
        self.assertTrue(config.TESTING)
        self.assertEqual(config.DB_CONNECTION_MODE, 'local')
    
    def test_staging_config(self):
        """Test the StagingConfig class."""
        config = StagingConfig()
        
        # Verify staging-specific values
        self.assertFalse(config.DEBUG)
        self.assertFalse(config.TESTING)
        self.assertEqual(config.DB_CONNECTION_MODE, 'auto')
    
    def test_production_config(self):
        """Test the ProductionConfig class."""
        config = ProductionConfig()
        
        # Verify production-specific values
        self.assertFalse(config.DEBUG)
        self.assertFalse(config.TESTING)
        self.assertEqual(config.DB_CONNECTION_MODE, 'auto')
        self.assertEqual(config.JWT_ACCESS_TOKEN_EXPIRES, 1800)  # 30 minutes
    
    def test_load_environment_variables(self):
        """Test load_environment_variables."""
        # Set some variables
        os.environ['DB_HOST'] = 'test-host'
        os.environ['SUPABASE_DB_PORT'] = '5432'
        
        # Load environment variables
        load_environment_variables()
        
        # Verify cross-filling
        self.assertEqual(os.environ['DB_HOST'], 'test-host')
        self.assertEqual(os.environ['SUPABASE_DB_HOST'], 'test-host')
        self.assertEqual(os.environ['DB_PORT'], '5432')
        self.assertEqual(os.environ['SUPABASE_DB_PORT'], '5432')
    
    def test_get_config_development(self):
        """Test get_config for development environment."""
        # Set environment
        os.environ['FLASK_ENV'] = 'development'
        
        # Get config
        config = get_config()
        
        # Verify
        self.assertIsInstance(config, DevelopmentConfig)
        self.assertTrue(config.DEBUG)
    
    def test_get_config_testing(self):
        """Test get_config for testing environment."""
        # Set environment
        os.environ['FLASK_ENV'] = 'testing'
        
        # Get config
        config = get_config()
        
        # Verify
        self.assertIsInstance(config, TestingConfig)
        self.assertTrue(config.TESTING)
    
    def test_get_config_staging(self):
        """Test get_config for staging environment."""
        # Set environment
        os.environ['FLASK_ENV'] = 'staging'
        
        # Get config
        config = get_config()
        
        # Verify
        self.assertIsInstance(config, StagingConfig)
    
    def test_get_config_production(self):
        """Test get_config for production environment."""
        # Set environment
        os.environ['FLASK_ENV'] = 'production'
        
        # Get config
        config = get_config()
        
        # Verify
        self.assertIsInstance(config, ProductionConfig)
    
    def test_get_database_config(self):
        """Test get_database_config."""
        # Set environment variables
        os.environ['DB_CONNECTION_MODE'] = 'auto'
        os.environ['LOCAL_DB_HOST'] = 'localhost'
        os.environ['SUPABASE_DB_HOST'] = 'db.example.supabase.co'
        
        # Get database config
        config = Config.get_database_config()
        
        # Verify
        self.assertIn('local', config)
        self.assertEqual(config['local']['host'], 'localhost')
        self.assertIn('supabase', config)
        self.assertEqual(config['supabase']['host'], 'db.example.supabase.co')

if __name__ == '__main__':
    unittest.main()
```

For `tests/test_integration.py`:
```python
#!/usr/bin/env python3
"""
Integration tests for database components.
"""

import os
import unittest
import psycopg2
from unittest.mock import patch

from database.connection_manager import ConnectionManager
from database.utils import (
    execute_query, get_molecule_by_id, insert_molecule,
    set_molecule_property, get_or_create_property_type
)

class TestDatabaseIntegration(unittest.TestCase):
    """Integration tests for database components."""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment."""
        # Save original environment
        cls.original_env = dict(os.environ)
        
        # Set test environment variables
        os.environ['DB_CONNECTION_MODE'] = 'local'
        os.environ['LOCAL_DB_HOST'] = 'localhost'
        os.environ['LOCAL_DB_PORT'] = '5432'
        os.environ['LOCAL_DB_NAME'] = 'cryoprotect_test'
        os.environ['LOCAL_DB_USER'] = 'postgres'
        os.environ['LOCAL_DB_PASSWORD'] = 'postgres'
        
        # Initialize test database
        try:
            # Create test database
            conn = psycopg2.connect(
                host='localhost',
                port=5432,
                database='postgres',
                user='postgres',
                password='postgres'
            )
            conn.autocommit = True
            cursor = conn.cursor()
            
            # Drop test database if it exists
            cursor.execute("DROP DATABASE IF EXISTS cryoprotect_test")
            
            # Create test database
            cursor.execute("CREATE DATABASE cryoprotect_test")
            
            # Close connection
            cursor.close()
            conn.close()
            
            # Connect to test database
            conn = psycopg2.connect(
                host='localhost',
                port=5432,
                database='cryoprotect_test',
                user='postgres',
                password='postgres'
            )
            cursor = conn.cursor()
            
            # Create test tables
            cursor.execute("""
                CREATE TABLE molecules (
                    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                    name VARCHAR NOT NULL,
                    formula VARCHAR,
                    molecular_weight DOUBLE PRECISION,
                    smiles VARCHAR,
                    inchi VARCHAR,
                    inchi_key VARCHAR,
                    chembl_id VARCHAR,
                    pubchem_cid VARCHAR,
                    data_source VARCHAR,
                    created_at TIMESTAMPTZ DEFAULT now(),
                    updated_at TIMESTAMPTZ DEFAULT now()
                )
            """)
            
            cursor.execute("""
                CREATE TABLE property_types (
                    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                    name VARCHAR NOT NULL UNIQUE,
                    description VARCHAR,
                    data_type VARCHAR NOT NULL,
                    unit VARCHAR,
                    created_at TIMESTAMPTZ DEFAULT now(),
                    updated_at TIMESTAMPTZ DEFAULT now()
                )
            """)
            
            cursor.execute("""
                CREATE TABLE molecular_properties (
                    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                    molecule_id UUID REFERENCES molecules(id),
                    property_type_id UUID REFERENCES property_types(id),
                    value DOUBLE PRECISION,
                    source VARCHAR,
                    confidence DOUBLE PRECISION,
                    created_at TIMESTAMPTZ DEFAULT now(),
                    updated_at TIMESTAMPTZ DEFAULT now(),
                    UNIQUE(molecule_id, property_type_id)
                )
            """)
            
            # Commit changes
            conn.commit()
            
            # Close connection
            cursor.close()
            conn.close()
        except Exception as e:
            print(f"Error setting up test database: {str(e)}")
            raise
    
    @classmethod
    def tearDownClass(cls):
        """Tear down the test environment."""
        # Restore original environment
        os.environ.clear()
        os.environ.update(cls.original_env)
        
        # Drop test database
        try:
            conn = psycopg2.connect(
                host='localhost',
                port=5432,
                database='postgres',
                user='postgres',
                password='postgres'
            )
            conn.autocommit = True
            cursor = conn.cursor()
            
            # Drop test database
            cursor.execute("DROP DATABASE IF EXISTS cryoprotect_test")
            
            # Close connection
            cursor.close()
            conn.close()
        except Exception as e:
            print(f"Error tearing down test database: {str(e)}")
    
    def setUp(self):
        """Set up the test case."""
        # Ensure new ConnectionManager is created
        ConnectionManager._instance = None
    
    def test_connection_manager_local(self):
        """Test connection manager with local adapter."""
        # Get connection manager
        manager = ConnectionManager.get_instance()
        
        # Connect
        result = manager.connect()
        
        # Verify
        self.assertTrue(result)
        self.assertEqual(manager.active_adapter, 'local')
    
    def test_insert_and_retrieve_molecule(self):
        """Test inserting and retrieving a molecule."""
        # Insert molecule
        molecule = insert_molecule(
            name="Test Molecule",
            formula="C10H20O",
            molecular_weight=156.27,
            smiles="CCCCCCCCCC(=O)",
            inchi="InChI=1S/C10H20O/c1-2-3-4-5-6-7-8-9-10-11/h10H,2-9H2,1H3",
            inchi_key="AXSIIYDNMFQWRZ-UHFFFAOYSA-N",
            chembl_id="CHEMBL1234",
            pubchem_cid="123456",
            data_source="test"
        )
        
        # Verify insert
        self.assertIsNotNone(molecule)
        self.assertEqual(molecule['name'], "Test Molecule")
        self.assertEqual(molecule['formula'], "C10H20O")
        
        # Get molecule by ID
        retrieved = get_molecule_by_id(molecule['id'])
        
        # Verify retrieve
        self.assertIsNotNone(retrieved)
        self.assertEqual(retrieved['id'], molecule['id'])
        self.assertEqual(retrieved['name'], "Test Molecule")
    
    def test_properties(self):
        """Test property types and molecular properties."""
        # Insert molecule
        molecule = insert_molecule(
            name="Property Test Molecule",
            inchi_key="TESTKEY123"
        )
        
        # Create property type
        prop_type = get_or_create_property_type(
            name="logP",
            description="Octanol-water partition coefficient",
            data_type="numeric",
            unit=""
        )
        
        # Verify property type
        self.assertIsNotNone(prop_type)
        self.assertEqual(prop_type['name'], "logP")
        
        # Set property
        property_result = set_molecule_property(
            molecule_id=molecule['id'],
            property_type_id=prop_type['id'],
            value=2.5,
            source="test",
            confidence=0.9
        )
        
        # Verify property
        self.assertIsNotNone(property_result)
        self.assertEqual(property_result['molecule_id'], molecule['id'])
        self.assertEqual(property_result['property_type_id'], prop_type['id'])
        self.assertEqual(property_result['value'], 2.5)

if __name__ == '__main__':
    unittest.main()
```

**Success Criteria:**
- 90%+ test coverage for database components
- All tests pass in CI environment
- Performance metrics documented
- Edge cases handled and tested

## Timeline

- Day 1: Task 3.1 (CI/CD Pipeline Implementation)
- Day 2: Task 3.2 (Docker Configuration Optimization)
- Day 3: Task 3.3 (Environment Configuration Standardization)
- Day 4: Task 4.1 (Implementation Documentation)
- Day 5: Task 4.2 (Comprehensive Testing)

## Dependency Map

```
Task 3.1 → Task 3.2 → Task 3.3 → Task 4.1 → Task 4.2
```

## Next Steps

After completing the deployment infrastructure tasks, the focus will shift to Phase 3.3: Security & Monitoring and Phase 3.4: Documentation & Knowledge Transfer.