#!/bin/bash
set -e

# Initialize Blue/Green Deployment Environment
# This script sets up the initial blue/green deployment environment

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
NGINX_CONF_DIR="$PROJECT_DIR/nginx/conf.d"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --version)
      VERSION="$2"
      shift 2
      ;;
    --help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --version VERSION    Specify the initial version to deploy"
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
  echo "Error: Version is required. Use --version to specify the initial version to deploy."
  exit 1
fi

echo "Initializing blue/green deployment environment with version $VERSION..."

# Create required directories
mkdir -p "$NGINX_CONF_DIR"

# Set environment variables for deployment
export BLUE_VERSION="$VERSION"
export GREEN_VERSION="$VERSION"
export FLASK_ENV="production"
export USE_EXTERNAL_SECRETS=true

# Pull the initial image
echo "Pulling image ghcr.io/${GITHUB_REPOSITORY:-yourusername/cryoprotect}/cryoprotect:$VERSION..."
docker pull "ghcr.io/${GITHUB_REPOSITORY:-yourusername/cryoprotect}/cryoprotect:$VERSION" || {
  echo "Error: Failed to pull image. Initialization aborted."
  exit 1
}

# Start the NGINX load balancer and blue environment
echo "Starting NGINX load balancer and blue environment..."
docker-compose -f "$PROJECT_DIR/docker-compose.blue-green.yml" up -d nginx cryoprotect-blue || {
  echo "Error: Failed to start NGINX and blue environment. Initialization aborted."
  exit 1
}

# Wait for the blue environment to be healthy
echo "Waiting for blue environment to be healthy..."
HEALTH_CHECK_RETRIES=30
HEALTH_CHECK_INTERVAL=5
for i in $(seq 1 $HEALTH_CHECK_RETRIES); do
  if docker ps --filter "name=cryoprotect-blue" --filter "health=healthy" --quiet | grep -q .; then
    echo "Blue environment is healthy!"
    HEALTHY=true
    break
  fi
  echo "Waiting for blue environment to be healthy... ($i/$HEALTH_CHECK_RETRIES)"
  sleep $HEALTH_CHECK_INTERVAL
done

if [ -z "$HEALTHY" ]; then
  echo "Error: Blue environment failed health checks. Initialization aborted."
  echo "Cleaning up..."
  docker-compose -f "$PROJECT_DIR/docker-compose.blue-green.yml" down
  exit 1
fi

# Apply database migrations
echo "Applying database migrations..."
docker exec "$(docker ps --filter "name=cryoprotect-blue" --quiet)" /opt/conda/envs/cryoprotect/bin/python -m database.migrations apply || {
  echo "Warning: Database migration failed. Continuing initialization..."
}

# Ensure traffic is routed to blue
echo "Ensuring traffic is routed to blue environment..."
"$SCRIPT_DIR/switch-to-blue.sh"

echo "Blue/green deployment environment initialized successfully!"
echo "Blue environment is active with version $VERSION"
echo ""
echo "To deploy a new version to the green environment:"
echo "  $SCRIPT_DIR/deploy-green.sh --version NEW_VERSION"
echo ""
echo "To switch traffic between environments:"
echo "  $SCRIPT_DIR/switch-to-blue.sh"
echo "  $SCRIPT_DIR/switch-to-green.sh"