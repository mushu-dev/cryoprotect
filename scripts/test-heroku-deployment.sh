#!/bin/bash

# Test script for Heroku deployment
# This script tests the Heroku deployment configuration for CryoProtect

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if Heroku CLI is installed
check_heroku_cli() {
  if ! command -v heroku &> /dev/null; then
    echo -e "${RED}Heroku CLI is not installed. Please install it first.${NC}"
    echo "Run: npm install -g heroku or visit https://devcenter.heroku.com/articles/heroku-cli"
    exit 1
  fi

  # Check if logged in
  if ! heroku auth:whoami &> /dev/null; then
    echo -e "${YELLOW}You are not logged in to Heroku CLI. Please login first.${NC}"
    echo "Run: heroku login"
    exit 1
  fi
}

# Test Heroku environment configuration
test_heroku_env() {
  echo -e "${BLUE}Testing Heroku environment configuration...${NC}"
  
  # Ask for the app name to test
  read -p "Enter the Heroku app name to test (e.g., cryoprotect-staging): " APP_NAME
  
  if ! heroku apps:info $APP_NAME &> /dev/null; then
    echo -e "${RED}App '$APP_NAME' not found. Please check the app name.${NC}"
    exit 1
  fi
  
  # Check if essential environment variables are set
  echo -e "Checking essential environment variables..."
  
  ESSENTIAL_VARS=("FLASK_ENV" "FLASK_APP" "SUPABASE_URL" "SUPABASE_KEY" "DEPLOYMENT_COLOR")
  MISSING_VARS=()
  
  for VAR in "${ESSENTIAL_VARS[@]}"; do
    if ! heroku config:get $VAR --app $APP_NAME &> /dev/null; then
      MISSING_VARS+=("$VAR")
    else
      echo -e "${GREEN}✓ $VAR is set${NC}"
    fi
  done
  
  if [ ${#MISSING_VARS[@]} -ne 0 ]; then
    echo -e "${RED}Missing essential environment variables: ${MISSING_VARS[*]}${NC}"
    echo -e "Please set these variables using: heroku config:set VAR=value --app $APP_NAME"
    echo -e "Or run the setup-heroku-deployment.sh script"
  else
    echo -e "${GREEN}All essential environment variables are set.${NC}"
  fi
  
  # Check for add-ons
  echo -e "\nChecking for required add-ons..."
  
  if heroku addons --app $APP_NAME | grep -q "heroku-postgresql"; then
    echo -e "${GREEN}✓ PostgreSQL add-on is installed${NC}"
  else
    echo -e "${RED}× PostgreSQL add-on is not installed${NC}"
    echo -e "Install it with: heroku addons:create heroku-postgresql:hobby-dev --app $APP_NAME"
  fi
  
  if heroku addons --app $APP_NAME | grep -q "heroku-redis"; then
    echo -e "${GREEN}✓ Redis add-on is installed${NC}"
  else
    echo -e "${YELLOW}× Redis add-on is not installed${NC}"
    echo -e "Install it with: heroku addons:create heroku-redis:hobby-dev --app $APP_NAME"
  fi
}

# Test web endpoint
test_web_endpoint() {
  echo -e "\n${BLUE}Testing web endpoint...${NC}"
  
  # Get the app URL
  APP_URL=$(heroku info -s --app $APP_NAME | grep web_url | cut -d= -f2)
  
  echo -e "Testing app URL: $APP_URL"
  
  # Test the base URL
  HTTP_STATUS=$(curl -o /dev/null -s -w "%{http_code}\n" $APP_URL)
  
  if [ "$HTTP_STATUS" == "200" ]; then
    echo -e "${GREEN}✓ Base URL is responding with 200 OK${NC}"
  else
    echo -e "${RED}× Base URL is responding with $HTTP_STATUS${NC}"
  fi
  
  # Test the health check endpoint
  HEALTH_STATUS=$(curl -o /dev/null -s -w "%{http_code}\n" "${APP_URL}health")
  
  if [ "$HEALTH_STATUS" == "200" ] || [ "$HEALTH_STATUS" == "207" ]; then
    echo -e "${GREEN}✓ Health check endpoint is responding with $HEALTH_STATUS${NC}"
    
    # Get the health check response
    HEALTH_RESPONSE=$(curl -s "${APP_URL}health")
    echo -e "Health check response: $HEALTH_RESPONSE"
  else
    echo -e "${RED}× Health check endpoint is responding with $HEALTH_STATUS${NC}"
  fi
  
  # Test the liveness endpoint
  LIVENESS_STATUS=$(curl -o /dev/null -s -w "%{http_code}\n" "${APP_URL}health/liveness")
  
  if [ "$LIVENESS_STATUS" == "200" ]; then
    echo -e "${GREEN}✓ Liveness endpoint is responding with 200 OK${NC}"
  else
    echo -e "${RED}× Liveness endpoint is responding with $LIVENESS_STATUS${NC}"
  fi
}

# Check logs for errors
check_logs() {
  echo -e "\n${BLUE}Checking logs for errors...${NC}"
  
  echo -e "Last 50 log lines from Heroku:"
  heroku logs -n 50 --app $APP_NAME
  
  echo -e "\nChecking for errors in logs..."
  ERROR_COUNT=$(heroku logs -n 100 --app $APP_NAME | grep -i "error\|exception\|failed" | wc -l)
  
  if [ "$ERROR_COUNT" -gt 0 ]; then
    echo -e "${YELLOW}Found $ERROR_COUNT potential error messages in the logs${NC}"
    heroku logs -n 100 --app $APP_NAME | grep -i "error\|exception\|failed"
  else
    echo -e "${GREEN}No obvious errors found in recent logs${NC}"
  fi
}

# Check database status
check_database() {
  echo -e "\n${BLUE}Checking database status...${NC}"
  
  # Get database info
  heroku pg:info --app $APP_NAME
  
  # Check if database has tables
  echo -e "\nChecking if database has tables..."
  TABLE_COUNT=$(heroku pg:psql -c "SELECT count(*) FROM information_schema.tables WHERE table_schema = 'public';" --app $APP_NAME 2>/dev/null | grep -o '[0-9]\+' | head -1)
  
  if [ -z "$TABLE_COUNT" ]; then
    echo -e "${RED}Could not get table count from database${NC}"
  elif [ "$TABLE_COUNT" -gt 0 ]; then
    echo -e "${GREEN}Database has $TABLE_COUNT tables${NC}"
    
    # List tables
    echo -e "\nTable list:"
    heroku pg:psql -c "SELECT table_name FROM information_schema.tables WHERE table_schema = 'public' ORDER BY table_name;" --app $APP_NAME
  else
    echo -e "${RED}Database has no tables. Database setup may have failed.${NC}"
  fi
}

# Main function to run the tests
main() {
  echo -e "${BLUE}=====================================${NC}"
  echo -e "${BLUE}CryoProtect Heroku Deployment Tester${NC}"
  echo -e "${BLUE}=====================================${NC}"
  
  check_heroku_cli
  test_heroku_env
  test_web_endpoint
  check_logs
  check_database
  
  echo -e "\n${GREEN}=====================================${NC}"
  echo -e "${GREEN}Deployment test completed!${NC}"
  echo -e "${GREEN}=====================================${NC}"
  echo -e ""
  echo -e "If there are any issues, refer to the HEROKU_DEPLOYMENT_GUIDE.md file."
}

# Execute main function
main