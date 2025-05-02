#!/bin/bash
set -e

# Rollback Script for Blue/Green Deployment
# This script rolls back to the previous version in case of deployment failure

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
NGINX_CONF_DIR="$PROJECT_DIR/nginx/conf.d"
ACTIVE_CONF="$NGINX_CONF_DIR/active.conf"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --version)
      VERSION="$2"
      shift 2
      ;;
    --environment)
      ENVIRONMENT="$2"
      shift 2
      ;;
    --help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --version VERSION      Specify the version to rollback to"
      echo "  --environment ENV      Specify the environment to rollback (blue or green)"
      echo "  --help                 Show this help message"
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
  echo "Error: Version is required. Use --version to specify the version to rollback to."
  exit 1
fi

# Check if environment is provided
if [ -z "$ENVIRONMENT" ]; then
  echo "Error: Environment is required. Use --environment to specify the environment to rollback (blue or green)."
  exit 1
fi

# Validate environment
if [ "$ENVIRONMENT" != "blue" ] && [ "$ENVIRONMENT" != "green" ]; then
  echo "Error: Invalid environment. Must be 'blue' or 'green'."
  exit 1
fi

echo "Rolling back to version $VERSION in $ENVIRONMENT environment..."

# Set environment variables for deployment
if [ "$ENVIRONMENT" = "blue" ]; then
  export BLUE_VERSION="$VERSION"
else
  export GREEN_VERSION="$VERSION"
fi

export FLASK_ENV="production"
export USE_EXTERNAL_SECRETS=true

# Pull the rollback image
echo "Pulling image ghcr.io/${GITHUB_REPOSITORY:-yourusername/cryoprotect}/cryoprotect:$VERSION..."
docker pull "ghcr.io/${GITHUB_REPOSITORY:-yourusername/cryoprotect}/cryoprotect:$VERSION" || {
  echo "Error: Failed to pull image. Rollback aborted."
  exit 1
}

# Stop the current environment
echo "Stopping current $ENVIRONMENT environment..."
docker-compose -f "$PROJECT_DIR/docker-compose.blue-green.yml" stop "cryoprotect-$ENVIRONMENT" || {
  echo "Warning: Failed to stop current $ENVIRONMENT environment. Continuing rollback..."
}

# Deploy the rollback version
echo "Deploying rollback version to $ENVIRONMENT environment..."
docker-compose -f "$PROJECT_DIR/docker-compose.blue-green.yml" up -d "cryoprotect-$ENVIRONMENT" || {
  echo "Error: Failed to deploy rollback version to $ENVIRONMENT environment. Rollback aborted."
  exit 1
}

# Wait for the service to be healthy
echo "Waiting for $ENVIRONMENT environment to be healthy..."
HEALTH_CHECK_RETRIES=30
HEALTH_CHECK_INTERVAL=5
for i in $(seq 1 $HEALTH_CHECK_RETRIES); do
  if docker ps --filter "name=cryoprotect-$ENVIRONMENT" --filter "health=healthy" --quiet | grep -q .; then
    echo "$ENVIRONMENT environment is healthy!"
    HEALTHY=true
    break
  fi
  echo "Waiting for $ENVIRONMENT environment to be healthy... ($i/$HEALTH_CHECK_RETRIES)"
  sleep $HEALTH_CHECK_INTERVAL
done

if [ -z "$HEALTHY" ]; then
  echo "Error: $ENVIRONMENT environment failed health checks after rollback. Manual intervention required."
  exit 1
fi

# Switch traffic to the rolled back environment
echo "Switching traffic to $ENVIRONMENT environment..."
if [ "$ENVIRONMENT" = "blue" ]; then
  "$SCRIPT_DIR/switch-to-blue.sh"
else
  "$SCRIPT_DIR/switch-to-green.sh"
fi

echo "Rollback to version $VERSION in $ENVIRONMENT environment completed successfully!"