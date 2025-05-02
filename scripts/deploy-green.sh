#!/bin/bash
set -e

# Blue/Green Deployment Script - Green Environment
# This script deploys a new version to the green environment

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DEPLOYMENT_COLOR="green"
NGINX_CONF_DIR="$PROJECT_DIR/nginx/conf.d"
HEALTH_CHECK_RETRIES=30
HEALTH_CHECK_INTERVAL=5

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --version)
      VERSION="$2"
      shift 2
      ;;
    --auto-switch)
      AUTO_SWITCH=true
      shift
      ;;
    --help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --version VERSION    Specify the version to deploy"
      echo "  --auto-switch        Automatically switch traffic to green after successful deployment"
      echo "  --help               Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Check if version is provided
if [ -z "$VERSION" ]; then
  echo "Error: Version is required. Use --version to specify the version to deploy."
  exit 1
fi

echo "Deploying version $VERSION to $DEPLOYMENT_COLOR environment..."

# Set environment variables for deployment
export GREEN_VERSION="$VERSION"
export FLASK_ENV="production"
export USE_EXTERNAL_SECRETS=true

# Pull the new image
echo "Pulling image ghcr.io/${GITHUB_REPOSITORY:-yourusername/cryoprotect}/cryoprotect:$VERSION..."
docker pull "ghcr.io/${GITHUB_REPOSITORY:-yourusername/cryoprotect}/cryoprotect:$VERSION" || {
  echo "Error: Failed to pull image. Deployment aborted."
  exit 1
}

# Deploy the green environment
echo "Deploying to $DEPLOYMENT_COLOR environment..."
docker-compose -f "$PROJECT_DIR/docker-compose.blue-green.yml" up -d "cryoprotect-$DEPLOYMENT_COLOR" || {
  echo "Error: Failed to deploy to $DEPLOYMENT_COLOR environment. Deployment aborted."
  exit 1
}

# Wait for the service to be healthy
echo "Waiting for $DEPLOYMENT_COLOR environment to be healthy..."
for i in $(seq 1 $HEALTH_CHECK_RETRIES); do
  if docker ps --filter "name=cryoprotect-$DEPLOYMENT_COLOR" --filter "health=healthy" --quiet | grep -q .; then
    echo "$DEPLOYMENT_COLOR environment is healthy!"
    HEALTHY=true
    break
  fi
  echo "Waiting for $DEPLOYMENT_COLOR environment to be healthy... ($i/$HEALTH_CHECK_RETRIES)"
  sleep $HEALTH_CHECK_INTERVAL
done

if [ -z "$HEALTHY" ]; then
  echo "Error: $DEPLOYMENT_COLOR environment failed health checks. Deployment aborted."
  echo "Rolling back..."
  docker-compose -f "$PROJECT_DIR/docker-compose.blue-green.yml" stop "cryoprotect-$DEPLOYMENT_COLOR"
  exit 1
fi

# Apply database migrations if needed
echo "Applying database migrations..."
docker exec "$(docker ps --filter "name=cryoprotect-$DEPLOYMENT_COLOR" --quiet)" /opt/conda/envs/cryoprotect/bin/python -m database.migrations apply || {
  echo "Warning: Database migration failed. Continuing deployment..."
}

# Switch traffic to green if auto-switch is enabled
if [ "$AUTO_SWITCH" = true ]; then
  echo "Switching traffic to $DEPLOYMENT_COLOR environment..."
  "$SCRIPT_DIR/switch-to-green.sh"
fi

echo "Deployment to $DEPLOYMENT_COLOR environment completed successfully!"