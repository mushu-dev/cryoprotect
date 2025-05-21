#!/bin/bash

# Heroku Deployment Setup Script for CryoProtect
# This script automates the setup of Heroku applications and add-ons for CryoProtect

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

# Create Heroku applications for different environments
create_heroku_apps() {
  echo -e "${BLUE}Creating Heroku applications...${NC}"
  
  # Production app
  if heroku apps:info cryoprotect-prod &> /dev/null; then
    echo -e "${YELLOW}Production app 'cryoprotect-prod' already exists.${NC}"
  else
    echo -e "Creating production app..."
    heroku create cryoprotect-prod --region us --buildpack heroku/python
    echo -e "${GREEN}Created production app 'cryoprotect-prod'.${NC}"
  fi
  
  # Staging app
  if heroku apps:info cryoprotect-staging &> /dev/null; then
    echo -e "${YELLOW}Staging app 'cryoprotect-staging' already exists.${NC}"
  else
    echo -e "Creating staging app..."
    heroku create cryoprotect-staging --region us --buildpack heroku/python
    echo -e "${GREEN}Created staging app 'cryoprotect-staging'.${NC}"
  fi
  
  # Review app
  if heroku apps:info cryoprotect-review &> /dev/null; then
    echo -e "${YELLOW}Review app 'cryoprotect-review' already exists.${NC}"
  else
    echo -e "Creating review app..."
    heroku create cryoprotect-review --region us --buildpack heroku/python
    echo -e "${GREEN}Created review app 'cryoprotect-review'.${NC}"
  fi
}

# Set up environment variables
setup_environment_variables() {
  echo -e "${BLUE}Setting up environment variables...${NC}"
  
  # Get environment variables from user
  echo -e "${YELLOW}Please provide the following environment variables:${NC}"
  
  # Supabase Production
  read -p "Supabase Production URL: " SUPABASE_URL_PROD
  read -p "Supabase Production Key: " SUPABASE_KEY_PROD
  read -p "Supabase Production Service Role Key: " SUPABASE_SERVICE_ROLE_KEY_PROD
  
  # Supabase Staging (optional)
  read -p "Supabase Staging URL (press Enter to use same as production): " SUPABASE_URL_STAGING
  SUPABASE_URL_STAGING=${SUPABASE_URL_STAGING:-$SUPABASE_URL_PROD}
  
  read -p "Supabase Staging Key (press Enter to use same as production): " SUPABASE_KEY_STAGING
  SUPABASE_KEY_STAGING=${SUPABASE_KEY_STAGING:-$SUPABASE_KEY_PROD}
  
  read -p "Supabase Staging Service Role Key (press Enter to use same as production): " SUPABASE_SERVICE_ROLE_KEY_STAGING
  SUPABASE_SERVICE_ROLE_KEY_STAGING=${SUPABASE_SERVICE_ROLE_KEY_STAGING:-$SUPABASE_SERVICE_ROLE_KEY_PROD}
  
  # Set production environment variables
  echo -e "Setting production environment variables..."
  heroku config:set \
    FLASK_ENV=production \
    FLASK_APP=app.py \
    DEPLOYMENT_COLOR=production \
    ENABLE_HEALTH_CHECKS=true \
    SUPABASE_URL="$SUPABASE_URL_PROD" \
    SUPABASE_KEY="$SUPABASE_KEY_PROD" \
    SUPABASE_SERVICE_ROLE_KEY="$SUPABASE_SERVICE_ROLE_KEY_PROD" \
    LOG_LEVEL=INFO \
    PYTHONUNBUFFERED=1 \
    --app cryoprotect-prod
  
  # Set staging environment variables
  echo -e "Setting staging environment variables..."
  heroku config:set \
    FLASK_ENV=staging \
    FLASK_APP=app.py \
    DEPLOYMENT_COLOR=staging \
    ENABLE_HEALTH_CHECKS=true \
    SUPABASE_URL="$SUPABASE_URL_STAGING" \
    SUPABASE_KEY="$SUPABASE_KEY_STAGING" \
    SUPABASE_SERVICE_ROLE_KEY="$SUPABASE_SERVICE_ROLE_KEY_STAGING" \
    LOG_LEVEL=DEBUG \
    PYTHONUNBUFFERED=1 \
    --app cryoprotect-staging
  
  # Set review environment variables 
  echo -e "Setting review environment variables..."
  heroku config:set \
    FLASK_ENV=development \
    FLASK_APP=app.py \
    DEPLOYMENT_COLOR=review \
    ENABLE_HEALTH_CHECKS=true \
    SUPABASE_URL="$SUPABASE_URL_STAGING" \
    SUPABASE_KEY="$SUPABASE_KEY_STAGING" \
    SUPABASE_SERVICE_ROLE_KEY="$SUPABASE_SERVICE_ROLE_KEY_STAGING" \
    LOG_LEVEL=DEBUG \
    PYTHONUNBUFFERED=1 \
    --app cryoprotect-review
  
  echo -e "${GREEN}Environment variables set successfully.${NC}"
}

# Setup PostgreSQL add-ons
setup_postgresql() {
  echo -e "${BLUE}Setting up PostgreSQL add-ons...${NC}"
  
  # Production database
  echo -e "Setting up production PostgreSQL..."
  heroku addons:create heroku-postgresql:hobby-dev --app cryoprotect-prod
  
  # Staging database
  echo -e "Setting up staging PostgreSQL..."
  heroku addons:create heroku-postgresql:hobby-dev --app cryoprotect-staging
  
  # Review database
  echo -e "Setting up review PostgreSQL..."
  heroku addons:create heroku-postgresql:hobby-dev --app cryoprotect-review
  
  echo -e "${GREEN}PostgreSQL add-ons created successfully.${NC}"
}

# Setup Redis add-ons
setup_redis() {
  echo -e "${BLUE}Setting up Redis add-ons...${NC}"
  
  # Production Redis
  echo -e "Setting up production Redis..."
  heroku addons:create heroku-redis:hobby-dev --app cryoprotect-prod
  
  # Staging Redis
  echo -e "Setting up staging Redis..."
  heroku addons:create heroku-redis:hobby-dev --app cryoprotect-staging
  
  # Review Redis
  echo -e "Setting up review Redis..."
  heroku addons:create heroku-redis:hobby-dev --app cryoprotect-review
  
  echo -e "${GREEN}Redis add-ons created successfully.${NC}"
}

# Configure health checks
setup_health_checks() {
  echo -e "${BLUE}Configuring health checks...${NC}"
  
  # Production health checks
  echo -e "Setting up production health checks..."
  heroku features:enable http-health-checks --app cryoprotect-prod
  heroku health:check:update path=/health/liveness period=30 tolerance=3 --app cryoprotect-prod
  
  # Staging health checks
  echo -e "Setting up staging health checks..."
  heroku features:enable http-health-checks --app cryoprotect-staging
  heroku health:check:update path=/health/liveness period=30 tolerance=3 --app cryoprotect-staging
  
  # Review health checks
  echo -e "Setting up review health checks..."
  heroku features:enable http-health-checks --app cryoprotect-review
  heroku health:check:update path=/health/liveness period=30 tolerance=3 --app cryoprotect-review
  
  echo -e "${GREEN}Health checks configured successfully.${NC}"
}

# Setup database backups
setup_database_backups() {
  echo -e "${BLUE}Configuring database backups...${NC}"
  
  # Production backups
  echo -e "Setting up daily production database backups..."
  heroku pg:backups:schedule --at '02:00 UTC' --app cryoprotect-prod
  
  # Staging backups
  echo -e "Setting up daily staging database backups..."
  heroku pg:backups:schedule --at '03:00 UTC' --app cryoprotect-staging
  
  echo -e "${GREEN}Database backups scheduled successfully.${NC}"
}

# Main function to run the setup
main() {
  echo -e "${BLUE}=====================================${NC}"
  echo -e "${BLUE}CryoProtect Heroku Deployment Setup${NC}"
  echo -e "${BLUE}=====================================${NC}"
  
  check_heroku_cli
  
  # Ask for confirmation before proceeding
  read -p "This script will set up Heroku applications and add-ons for CryoProtect. Continue? (y/n): " confirm
  if [[ $confirm != "y" && $confirm != "Y" ]]; then
    echo -e "${YELLOW}Setup cancelled.${NC}"
    exit 0
  fi
  
  create_heroku_apps
  setup_environment_variables
  setup_postgresql
  setup_redis
  setup_health_checks
  setup_database_backups
  
  echo -e "${GREEN}=====================================${NC}"
  echo -e "${GREEN}Heroku deployment setup completed!${NC}"
  echo -e "${GREEN}=====================================${NC}"
  echo -e ""
  echo -e "Next steps:"
  echo -e "1. Set up GitHub Actions secrets for deployment"
  echo -e "2. Test the deployment workflow"
  echo -e "3. Monitor the application using Heroku dashboard"
  echo -e ""
  echo -e "For more details, refer to the HEROKU_DEPLOYMENT_GUIDE.md file."
}

# Execute main function
main