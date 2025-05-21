# Deployment Infrastructure Enhancement Plan

## Current Status

The CryoProtect v2 project has completed:
- Phase 3.1, Stage 1: Database Foundation (RLS implementation, data population, verification)
- Toxicity data integration with Tox21

The CI/CD pipeline in `.github/workflows/deploy.yml` has a solid foundation with:
- Testing job (unit and integration tests)
- Build job with versioning
- Basic deployment to staging/production

## Next Priority: Blue/Green Deployment Implementation

Implementing a blue/green deployment strategy will provide zero-downtime deployments and better rollback capabilities, which are critical for production environment reliability.

### Task 1: Create Blue/Green Deployment Scripts

Create a deployment script that implements blue/green deployment logic:

```bash
#!/bin/bash
# scripts/blue-green-deploy.sh
#
# Blue-Green Deployment Script for CryoProtect v2
# This script implements a zero-downtime deployment using blue/green strategy

set -e

# Configuration
APP_NAME="cryoprotect"
DEPLOY_PATH="/var/www"
BLUE_DIR="$DEPLOY_PATH/blue"
GREEN_DIR="$DEPLOY_PATH/green"
LIVE_LINK="$DEPLOY_PATH/current"
PREVIOUS_LINK="$DEPLOY_PATH/previous"
VERSION="$1"
ENV_TYPE="$2"  # "production" or "staging"

# Logging
LOG_FILE="/var/log/cryoprotect-deploy-$(date +%Y%m%d-%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== Starting Blue/Green Deployment ==="
echo "Version: $VERSION"
echo "Environment: $ENV_TYPE"
echo "Timestamp: $(date)"

# Determine current live environment (blue or green)
if [[ -L "$LIVE_LINK" ]]; then
    CURRENT_ENV=$(basename $(readlink -f "$LIVE_LINK"))
    echo "Current live environment: $CURRENT_ENV"
else
    CURRENT_ENV="none"
    echo "No current live environment detected"
fi

# Determine target environment
if [[ "$CURRENT_ENV" == "blue" ]]; then
    TARGET_ENV="green"
    TARGET_DIR="$GREEN_DIR"
else
    TARGET_ENV="blue"
    TARGET_DIR="$BLUE_DIR"
fi
echo "Target environment: $TARGET_ENV ($TARGET_DIR)"

# Create target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Extract deployment package to target environment
echo "Extracting deployment package to $TARGET_DIR"
tar -xzf "deploy-$VERSION.tar.gz" -C "$TARGET_DIR"

# Set environment-specific configuration
echo "Setting up environment-specific configuration"
if [[ "$ENV_TYPE" == "production" ]]; then
    # Apply production configuration
    cp "$TARGET_DIR/config/config_production.py" "$TARGET_DIR/config.py"
else
    # Apply staging configuration
    cp "$TARGET_DIR/config/config_staging.py" "$TARGET_DIR/config.py"
fi

# Start application in target environment
echo "Starting application in $TARGET_ENV environment"
cd "$TARGET_DIR" && pm2 start app.js --name "$APP_NAME-$TARGET_ENV" || pm2 restart "$APP_NAME-$TARGET_ENV"

# Verify application is running in target environment
echo "Verifying application health in target environment"
for i in {1..10}; do
    if curl -s http://localhost:5000/api/health | grep -q "ok"; then
        echo "Application is healthy in $TARGET_ENV environment"
        HEALTHY=true
        break
    else
        echo "Waiting for application to become healthy (attempt $i/10)..."
        sleep 5
    fi
done

if [[ -z "$HEALTHY" ]]; then
    echo "ERROR: Application failed health check in target environment"
    echo "Rolling back to previous environment"
    pm2 stop "$APP_NAME-$TARGET_ENV"
    exit 1
fi

# Switch to new environment
echo "Switching live environment to $TARGET_ENV"
if [[ -L "$LIVE_LINK" ]]; then
    # Save previous environment link
    rm -f "$PREVIOUS_LINK"
    mv "$LIVE_LINK" "$PREVIOUS_LINK"
fi

# Create new live link
ln -sf "$TARGET_DIR" "$LIVE_LINK"

# Update nginx configuration to point to new environment
echo "Updating web server configuration"
if [[ -f "/etc/nginx/sites-available/$APP_NAME" ]]; then
    sed -i "s|root .*|root $LIVE_LINK;|g" "/etc/nginx/sites-available/$APP_NAME"
    nginx -t && systemctl reload nginx
fi

# Stop previous environment after grace period
if [[ "$CURRENT_ENV" != "none" ]]; then
    echo "Waiting for connections to drain from previous environment..."
    sleep 30
    echo "Stopping previous environment ($CURRENT_ENV)"
    pm2 stop "$APP_NAME-$CURRENT_ENV" || true
fi

echo "=== Blue/Green Deployment Completed Successfully ==="
echo "New live environment: $TARGET_ENV"
echo "Version deployed: $VERSION"
echo "Deployment completed at: $(date)"
```

### Task 2: Create Rollback Script

Create a rollback script for quick recovery in case of issues:

```bash
#!/bin/bash
# scripts/rollback.sh
#
# Rollback Script for CryoProtect v2
# This script rolls back to the previous deployment

set -e

# Configuration
APP_NAME="cryoprotect"
DEPLOY_PATH="/var/www"
LIVE_LINK="$DEPLOY_PATH/current"
PREVIOUS_LINK="$DEPLOY_PATH/previous"

# Logging
LOG_FILE="/var/log/cryoprotect-rollback-$(date +%Y%m%d-%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== Starting Rollback Process ==="
echo "Timestamp: $(date)"

# Check if previous deployment exists
if [[ ! -L "$PREVIOUS_LINK" ]]; then
    echo "ERROR: No previous deployment found to roll back to"
    exit 1
fi

# Determine current and previous environments
CURRENT_ENV=$(basename $(readlink -f "$LIVE_LINK"))
PREVIOUS_ENV=$(basename $(readlink -f "$PREVIOUS_LINK"))

echo "Current environment: $CURRENT_ENV"
echo "Rolling back to: $PREVIOUS_ENV"

# Start previous environment if not running
echo "Starting previous environment"
cd "$(readlink -f "$PREVIOUS_LINK")" && pm2 start app.js --name "$APP_NAME-$PREVIOUS_ENV" || pm2 restart "$APP_NAME-$PREVIOUS_ENV"

# Verify previous environment is healthy
echo "Verifying health of previous environment"
for i in {1..5}; do
    if curl -s http://localhost:5000/api/health | grep -q "ok"; then
        echo "Previous environment is healthy"
        HEALTHY=true
        break
    else
        echo "Waiting for previous environment to become healthy (attempt $i/5)..."
        sleep 5
    fi
done

if [[ -z "$HEALTHY" ]]; then
    echo "WARNING: Previous environment failed health check, but proceeding with rollback anyway"
fi

# Swap links
echo "Switching back to previous environment"
TEMP_LINK="$DEPLOY_PATH/temp"
mv "$LIVE_LINK" "$TEMP_LINK"
mv "$PREVIOUS_LINK" "$LIVE_LINK"
mv "$TEMP_LINK" "$PREVIOUS_LINK"

# Update nginx configuration
echo "Updating web server configuration"
if [[ -f "/etc/nginx/sites-available/$APP_NAME" ]]; then
    sed -i "s|root .*|root $LIVE_LINK;|g" "/etc/nginx/sites-available/$APP_NAME"
    nginx -t && systemctl reload nginx
fi

# Stop current environment after grace period
echo "Waiting for connections to drain from failed environment..."
sleep 30
echo "Stopping failed environment ($CURRENT_ENV)"
pm2 stop "$APP_NAME-$CURRENT_ENV" || true

echo "=== Rollback Completed Successfully ==="
echo "Reverted to environment: $PREVIOUS_ENV"
echo "Rollback completed at: $(date)"
```

### Task 3: Create Health Check API Endpoint

Add a health check endpoint to monitor application health:

```python
# api/system_resources.py

from flask import jsonify
from flask_restful import Resource
import datetime
import os
import psutil
import logging

logger = logging.getLogger(__name__)

class HealthCheckResource(Resource):
    """Resource for checking the health of the application."""
    
    def get(self):
        """
        Get application health status.
        
        Returns:
            JSON response with health status
        """
        try:
            # Basic health check data
            health_data = {
                "status": "ok",
                "timestamp": datetime.datetime.utcnow().isoformat(),
                "version": os.environ.get("APP_VERSION", "unknown"),
                "environment": os.environ.get("FLASK_ENV", "development"),
                "system": {
                    "cpu_usage": psutil.cpu_percent(interval=0.1),
                    "memory_usage": psutil.virtual_memory().percent,
                    "disk_usage": psutil.disk_usage('/').percent
                }
            }
            
            # Check database connection
            try:
                from api.models import BaseModel
                db_start = datetime.datetime.now()
                supabase = BaseModel.get_supabase()
                response = supabase.from_("health_check").select("*").limit(1).execute()
                db_end = datetime.datetime.now()
                db_latency = (db_end - db_start).total_seconds() * 1000  # in ms
                
                health_data["database"] = {
                    "status": "connected",
                    "latency_ms": db_latency
                }
            except Exception as e:
                logger.error(f"Database health check failed: {str(e)}")
                health_data["database"] = {
                    "status": "error",
                    "error": str(e)
                }
                health_data["status"] = "degraded"
            
            return jsonify(health_data)
            
        except Exception as e:
            logger.error(f"Health check failed: {str(e)}")
            return jsonify({
                "status": "error",
                "error": str(e),
                "timestamp": datetime.datetime.utcnow().isoformat()
            }), 500

def register_resources(api):
    """Register system resources with the API."""
    api.add_resource(HealthCheckResource, '/api/health')
```

### Task 4: Update GitHub Actions Workflow for Blue/Green Deployment

Update the CI/CD workflow to use the blue/green deployment scripts:

```yaml
# .github/workflows/deploy.yml (modifications to the deploy job)

deploy:
  name: Deploy to server
  runs-on: ubuntu-latest
  needs: build
  environment: ${{ github.event.inputs.environment || (github.ref == 'refs/heads/main' && 'production' || 'staging') }}
  
  steps:
    # ... [existing steps] ...
    
    - name: Deploy to server with Blue/Green strategy
      env:
        SSH_HOST: ${{ github.ref == 'refs/heads/main' && secrets.PRODUCTION_SSH_HOST || secrets.STAGING_SSH_HOST }}
        SSH_USER: ${{ github.ref == 'refs/heads/main' && secrets.PRODUCTION_SSH_USER || secrets.STAGING_SSH_USER }}
        TARGET_DIR: ${{ github.ref == 'refs/heads/main' && '/var/www' || '/var/www' }}
        ENV_TYPE: ${{ github.ref == 'refs/heads/main' && 'production' || 'staging' }}
        APP_VERSION: ${{ needs.build.outputs.version }}
      run: |
        # Add server to known hosts
        ssh-keyscan -H $SSH_HOST >> ~/.ssh/known_hosts
        
        # Copy the deployment archive and scripts
        scp deploy-$APP_VERSION.tar.gz $SSH_USER@$SSH_HOST:$TARGET_DIR/ || { echo "Failed to copy deployment archive"; exit 1; }
        scp scripts/blue-green-deploy.sh $SSH_USER@$SSH_HOST:$TARGET_DIR/ || { echo "Failed to copy deployment script"; exit 1; }
        scp scripts/rollback.sh $SSH_USER@$SSH_HOST:$TARGET_DIR/ || { echo "Failed to copy rollback script"; exit 1; }
        
        # Make scripts executable
        ssh $SSH_USER@$SSH_HOST "chmod +x $TARGET_DIR/blue-green-deploy.sh $TARGET_DIR/rollback.sh" || { echo "Failed to make scripts executable"; exit 1; }
        
        # Create Docker secrets for sensitive values
        # ... [existing secrets setup code] ...
        
        # Execute blue/green deployment
        ssh $SSH_USER@$SSH_HOST "cd $TARGET_DIR && ./blue-green-deploy.sh $APP_VERSION $ENV_TYPE" || { 
          echo "Deployment failed! Initiating rollback..."
          ssh $SSH_USER@$SSH_HOST "cd $TARGET_DIR && ./rollback.sh"
          exit 1
        }
        
        echo "Blue/Green deployment of version $APP_VERSION to $ENV_TYPE completed successfully!"
```

### Task 5: Docker Optimization

Enhance the Dockerfile with multi-stage builds for smaller, more secure images:

```dockerfile
# Dockerfile
# Multi-stage build for CryoProtect v2

# ===== Build Stage =====
FROM node:18-alpine AS build

# Set working directory
WORKDIR /app

# Copy package files
COPY package*.json ./

# Install dependencies
RUN npm ci

# Copy source code
COPY . .

# Build frontend
RUN npm run build

# ===== Python Dependencies Stage =====
FROM python:3.9-slim AS python-deps

# Set working directory
WORKDIR /app

# Copy requirements file
COPY requirements_updated.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements_updated.txt

# ===== Final Stage =====
FROM python:3.9-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    nginx \
    supervisor \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN groupadd -r cryoprotect && \
    useradd -r -g cryoprotect -s /bin/bash cryoprotect && \
    mkdir -p /home/cryoprotect && \
    chown -R cryoprotect:cryoprotect /home/cryoprotect

# Copy built frontend from build stage
COPY --from=build /app/dist /app/static

# Copy Python dependencies
COPY --from=python-deps /usr/local/lib/python3.9/site-packages /usr/local/lib/python3.9/site-packages
COPY --from=python-deps /usr/local/bin /usr/local/bin

# Copy application code
COPY app.py config*.py ./
COPY api ./api
COPY chemical_data ./chemical_data
COPY migrations ./migrations
COPY templates ./templates

# Copy configuration files
COPY nginx/cryoprotect.conf /etc/nginx/sites-available/default
COPY supervisor/cryoprotect.conf /etc/supervisor/conf.d/

# Set permissions
RUN chown -R cryoprotect:cryoprotect /app

# Environment variables
ENV PYTHONUNBUFFERED=1 \
    FLASK_APP=app.py \
    PORT=5000

# Expose ports
EXPOSE 80 5000

# Switch to non-root user
USER cryoprotect

# Health check
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
    CMD curl -f http://localhost:5000/api/health || exit 1

# Start supervisor to manage processes
CMD ["/usr/bin/supervisord", "-c", "/etc/supervisor/supervisord.conf"]
```

### Task 6: Create Docker Compose Configuration

Create production-optimized Docker Compose configuration:

```yaml
# docker-compose.prod.yml
version: '3.8'

services:
  cryoprotect:
    image: cryoprotect:${APP_VERSION:-latest}
    build:
      context: .
      dockerfile: Dockerfile
    restart: unless-stopped
    environment:
      - FLASK_ENV=production
      - FLASK_APP=app.py
      - USE_EXTERNAL_SECRETS=true
      - APP_VERSION=${APP_VERSION:-unknown}
    secrets:
      - cryoprotect_supabase_url
      - cryoprotect_supabase_key
      - cryoprotect_secret_key
      - cryoprotect_redis_url
    ports:
      - "80:80"
      - "5000:5000"
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:5000/api/health"]
      interval: 30s
      timeout: 5s
      retries: 3
      start_period: 30s
    deploy:
      resources:
        limits:
          cpus: '1.0'
          memory: 1G
      restart_policy:
        condition: on-failure
        delay: 5s
        max_attempts: 3
        window: 120s
    logging:
      driver: "json-file"
      options:
        max-size: "10m"
        max-file: "3"

secrets:
  cryoprotect_supabase_url:
    external: true
  cryoprotect_supabase_key:
    external: true
  cryoprotect_secret_key:
    external: true
  cryoprotect_redis_url:
    external: true

networks:
  default:
    driver: bridge
```

## Implementation Steps

1. Create the deployment scripts directory and add the scripts:
   - `scripts/blue-green-deploy.sh`
   - `scripts/rollback.sh`

2. Add the health check endpoint to the API:
   - Update `api/system_resources.py`
   - Register the endpoint in `api/__init__.py`

3. Optimize the Docker configuration:
   - Update `Dockerfile` with multi-stage builds
   - Create `docker-compose.prod.yml`

4. Update the GitHub Actions workflow:
   - Modify `.github/workflows/deploy.yml` for blue/green deployment

5. Test the deployment process in staging environment:
   - Deploy to staging with the updated workflow
   - Verify zero-downtime deployment
   - Test rollback functionality

## Success Criteria

- Zero-downtime deployments in staging and production environments
- Successful rollbacks when needed
- Reduced Docker image size and improved security
- Health check endpoint reporting system status
- All deployments completing without service interruption

## Next Steps After This Implementation

1. Add infrastructure as code (Terraform) for server provisioning
2. Implement container orchestration with Kubernetes
3. Set up advanced monitoring for deployments
4. Add database migration CI/CD automation