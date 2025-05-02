#!/bin/bash
set -e

# Health Check Script for Blue/Green Deployment
# This script checks the health of blue and green environments

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
NGINX_CONF_DIR="$PROJECT_DIR/nginx/conf.d"
ACTIVE_CONF="$NGINX_CONF_DIR/active.conf"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --environment)
      ENVIRONMENT="$2"
      shift 2
      ;;
    --verbose)
      VERBOSE=true
      shift
      ;;
    --help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --environment ENV    Specify the environment to check (blue, green, or all)"
      echo "  --verbose            Show detailed health information"
      echo "  --help               Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Default to checking all environments if not specified
if [ -z "$ENVIRONMENT" ]; then
  ENVIRONMENT="all"
fi

# Validate environment
if [ "$ENVIRONMENT" != "blue" ] && [ "$ENVIRONMENT" != "green" ] && [ "$ENVIRONMENT" != "all" ]; then
  echo "Error: Invalid environment. Must be 'blue', 'green', or 'all'."
  exit 1
fi

# Determine active environment
ACTIVE_ENV=""
if [ -f "$ACTIVE_CONF" ]; then
  if grep -q "proxy_pass http://blue" "$ACTIVE_CONF"; then
    ACTIVE_ENV="blue"
  elif grep -q "proxy_pass http://green" "$ACTIVE_CONF"; then
    ACTIVE_ENV="green"
  fi
fi

# Function to check environment health
check_environment() {
  local env=$1
  local container_name="cryoprotect-$env"
  
  echo "Checking $env environment..."
  
  # Check if container exists
  if ! docker ps --filter "name=$container_name" --quiet | grep -q .; then
    echo "  Status: Not running"
    return 1
  fi
  
  # Check container health status
  local health_status=$(docker inspect --format='{{.State.Health.Status}}' "$(docker ps --filter "name=$container_name" --quiet)")
  echo "  Status: $health_status"
  
  # Get container version
  local version=$(docker exec "$(docker ps --filter "name=$container_name" --quiet)" /opt/conda/envs/cryoprotect/bin/python -c "import os; print(os.environ.get('APP_VERSION', 'unknown'))" 2>/dev/null || echo "unknown")
  echo "  Version: $version"
  
  # Check application health endpoint
  if [ "$VERBOSE" = true ]; then
    echo "  Application health:"
    docker exec "$(docker ps --filter "name=$container_name" --quiet)" curl -s http://localhost:5000/health || echo "    Health endpoint not accessible"
  fi
  
  # Check if this is the active environment
  if [ "$env" = "$ACTIVE_ENV" ]; then
    echo "  Traffic: Receiving (ACTIVE)"
  else
    echo "  Traffic: Not receiving"
  fi
  
  # Return success if health status is healthy
  if [ "$health_status" = "healthy" ]; then
    return 0
  else
    return 1
  fi
}

# Check NGINX status
check_nginx() {
  echo "Checking NGINX load balancer..."
  
  # Check if container exists
  if ! docker ps --filter "name=nginx" --quiet | grep -q .; then
    echo "  Status: Not running"
    return 1
  fi
  
  # Check NGINX status
  local status=$(docker exec "$(docker ps --filter "name=nginx" --quiet)" nginx -t 2>&1)
  if echo "$status" | grep -q "successful"; then
    echo "  Status: Running (configuration OK)"
  else
    echo "  Status: Running (configuration ERROR)"
    if [ "$VERBOSE" = true ]; then
      echo "  Configuration test output:"
      echo "$status"
    fi
    return 1
  fi
  
  # Check NGINX health endpoint
  local health_status=$(docker exec "$(docker ps --filter "name=nginx" --quiet)" curl -s -o /dev/null -w "%{http_code}" http://localhost:8080/health 2>/dev/null || echo "failed")
  if [ "$health_status" = "200" ]; then
    echo "  Health: OK (200)"
  else
    echo "  Health: ERROR ($health_status)"
    return 1
  fi
  
  return 0
}

# Main health check logic
OVERALL_STATUS=0

# Check NGINX
check_nginx
NGINX_STATUS=$?
OVERALL_STATUS=$((OVERALL_STATUS + NGINX_STATUS))

echo ""

# Check blue environment if requested
if [ "$ENVIRONMENT" = "blue" ] || [ "$ENVIRONMENT" = "all" ]; then
  check_environment "blue"
  BLUE_STATUS=$?
  OVERALL_STATUS=$((OVERALL_STATUS + BLUE_STATUS))
  echo ""
fi

# Check green environment if requested
if [ "$ENVIRONMENT" = "green" ] || [ "$ENVIRONMENT" = "all" ]; then
  check_environment "green"
  GREEN_STATUS=$?
  OVERALL_STATUS=$((OVERALL_STATUS + GREEN_STATUS))
  echo ""
fi

# Print summary
echo "Health Check Summary:"
echo "--------------------"
if [ "$NGINX_STATUS" -eq 0 ]; then
  echo "NGINX: HEALTHY"
else
  echo "NGINX: UNHEALTHY"
fi

if [ "$ENVIRONMENT" = "blue" ] || [ "$ENVIRONMENT" = "all" ]; then
  if [ "$BLUE_STATUS" -eq 0 ]; then
    echo "Blue: HEALTHY"
  else
    echo "Blue: UNHEALTHY"
  fi
fi

if [ "$ENVIRONMENT" = "green" ] || [ "$ENVIRONMENT" = "all" ]; then
  if [ "$GREEN_STATUS" -eq 0 ]; then
    echo "Green: HEALTHY"
  else
    echo "Green: UNHEALTHY"
  fi
fi

echo ""
echo "Active environment: ${ACTIVE_ENV:-None}"

# Exit with appropriate status code
if [ "$OVERALL_STATUS" -eq 0 ]; then
  echo "Overall status: HEALTHY"
  exit 0
else
  echo "Overall status: UNHEALTHY"
  exit 1
fi