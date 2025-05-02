#!/bin/bash
# Blue/Green Deployment Cutover Script
# This script handles the cutover process between blue and green environments

set -e

# Default values
TARGET_COLOR=""
FORCE=false
TIMEOUT=300
CHECK_INTERVAL=5
ROLLBACK=false
CURRENT_COLOR=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --target|-t)
      TARGET_COLOR="$2"
      shift
      shift
      ;;
    --force|-f)
      FORCE=true
      shift
      ;;
    --timeout)
      TIMEOUT="$2"
      shift
      shift
      ;;
    --interval)
      CHECK_INTERVAL="$2"
      shift
      shift
      ;;
    --rollback|-r)
      ROLLBACK=true
      shift
      ;;
    --help|-h)
      echo "Usage: cutover.sh [OPTIONS]"
      echo "Options:"
      echo "  --target, -t COLOR    Target environment color (blue or green)"
      echo "  --force, -f           Force cutover without health checks"
      echo "  --timeout SECONDS     Health check timeout in seconds (default: 300)"
      echo "  --interval SECONDS    Health check interval in seconds (default: 5)"
      echo "  --rollback, -r        Rollback to previous environment"
      echo "  --help, -h            Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $key"
      exit 1
      ;;
  esac
done

# Get current active color from active.conf
CURRENT_COLOR=$(grep -oP "ACTIVE_COLOR: \K(blue|green)" /etc/nginx/conf.d/active.conf || echo "unknown")

# If rollback is requested, determine the target color
if [ "$ROLLBACK" = true ]; then
  if [ "$CURRENT_COLOR" = "blue" ]; then
    TARGET_COLOR="green"
  elif [ "$CURRENT_COLOR" = "green" ]; then
    TARGET_COLOR="blue"
  else
    echo "Error: Cannot determine current color for rollback"
    exit 1
  fi
fi

# Validate target color
if [ "$TARGET_COLOR" != "blue" ] && [ "$TARGET_COLOR" != "green" ]; then
  echo "Error: Target color must be 'blue' or 'green'"
  exit 1
fi

# Check if target is already active
if [ "$CURRENT_COLOR" = "$TARGET_COLOR" ]; then
  echo "Target color $TARGET_COLOR is already active. No action needed."
  exit 0
fi

echo "Starting cutover from $CURRENT_COLOR to $TARGET_COLOR"

# Check target environment health if not forced
if [ "$FORCE" = false ]; then
  echo "Checking $TARGET_COLOR environment health..."
  
  ELAPSED=0
  while [ $ELAPSED -lt $TIMEOUT ]; do
    HTTP_STATUS=$(curl -s -o /dev/null -w "%{http_code}" http://localhost:8080/health/$TARGET_COLOR)
    
    if [ "$HTTP_STATUS" = "200" ]; then
      echo "$TARGET_COLOR environment is healthy!"
      break
    else
      echo "$TARGET_COLOR environment is not healthy (status: $HTTP_STATUS). Retrying in $CHECK_INTERVAL seconds..."
      sleep $CHECK_INTERVAL
      ELAPSED=$((ELAPSED + CHECK_INTERVAL))
    fi
  done
  
  if [ $ELAPSED -ge $TIMEOUT ]; then
    echo "Error: $TARGET_COLOR environment health check timed out after $TIMEOUT seconds"
    exit 1
  fi
fi

# Perform the cutover
echo "Performing cutover to $TARGET_COLOR environment..."

# Create a timestamp for the deployment
TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

# Get the current version from the target config
TARGET_VERSION=$(grep -oP "ACTIVE_VERSION: \K([^\s]+)" /etc/nginx/conf.d/${TARGET_COLOR}.conf || echo "latest")

# Copy the target configuration to active.conf
cp /etc/nginx/conf.d/${TARGET_COLOR}.conf /etc/nginx/conf.d/active.conf

# Update the timestamp in active.conf
sed -i "s/DEPLOYMENT_TIMESTAMP: .*/DEPLOYMENT_TIMESTAMP: $TIMESTAMP/" /etc/nginx/conf.d/active.conf

# Reload NGINX configuration
echo "Reloading NGINX configuration..."
nginx -s reload

# Verify the cutover was successful
echo "Verifying cutover..."
sleep 2
NEW_COLOR=$(grep -oP "ACTIVE_COLOR: \K(blue|green)" /etc/nginx/conf.d/active.conf || echo "unknown")

if [ "$NEW_COLOR" = "$TARGET_COLOR" ]; then
  echo "Cutover to $TARGET_COLOR environment successful!"
  echo "Deployment timestamp: $TIMESTAMP"
  echo "Version: $TARGET_VERSION"
else
  echo "Error: Cutover verification failed. Active color is $NEW_COLOR, expected $TARGET_COLOR"
  exit 1
fi

exit 0