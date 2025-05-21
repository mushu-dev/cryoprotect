#!/bin/bash

# Deployment script for CryoProtect backend-frontend integration
# This script helps deploy the necessary changes to all services

set -e  # Exit on error

# Color definitions for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;36m'
NC='\033[0m' # No Color

echo -e "${BLUE}===========================================================${NC}"
echo -e "${GREEN}CryoProtect Integration Deployment Script${NC}"
echo -e "${BLUE}===========================================================${NC}"
echo ""

# Configuration variables with defaults
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect-8030e4025428}
NETLIFY_SITE_NAME=${NETLIFY_SITE_NAME:-cryoprotect}
DEPLOY_NETLIFY=${DEPLOY_NETLIFY:-true}
DEPLOY_HEROKU=${DEPLOY_HEROKU:-true}

echo -e "${YELLOW}Service Configuration:${NC}"
echo -e "  ${BLUE}Heroku App:${NC} $HEROKU_APP_NAME"
echo -e "  ${BLUE}Netlify Site:${NC} $NETLIFY_SITE_NAME"
echo -e "  ${BLUE}Deploy to Netlify:${NC} $DEPLOY_NETLIFY"
echo -e "  ${BLUE}Deploy to Heroku:${NC} $DEPLOY_HEROKU"
echo ""

# Confirm configuration
read -p "Do you want to continue with this configuration? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${RED}Deployment aborted by user.${NC}"
    exit 0
fi

# Check for required CLIs
if [[ "$DEPLOY_HEROKU" == "true" ]] && ! command -v heroku &> /dev/null; then
    echo -e "${RED}Heroku CLI not found. Please install it:${NC}"
    echo "npm install -g heroku"
    exit 1
fi

if [[ "$DEPLOY_NETLIFY" == "true" ]] && ! command -v netlify &> /dev/null; then
    echo -e "${RED}Netlify CLI not found. Please install it:${NC}"
    echo "npm install -g netlify-cli"
    exit 1
fi

# Check authentication
if [[ "$DEPLOY_HEROKU" == "true" ]]; then
    echo "Checking Heroku authentication..."
    if ! heroku auth:whoami &> /dev/null; then
        echo -e "${RED}Not authenticated with Heroku. Please run:${NC}"
        echo "heroku login"
        exit 1
    fi
    echo -e "${GREEN}Authenticated with Heroku.${NC}"
fi

if [[ "$DEPLOY_NETLIFY" == "true" ]]; then
    echo "Checking Netlify authentication..."
    if ! netlify status &> /dev/null; then
        echo -e "${RED}Not authenticated with Netlify. Please run:${NC}"
        echo "netlify login"
        exit 1
    fi
    echo -e "${GREEN}Authenticated with Netlify.${NC}"
fi

# Create deployment temp directory
echo -e "\n${YELLOW}Creating deployment files...${NC}"
DEPLOY_TEMP=$(mktemp -d)
echo "Using temporary directory: $DEPLOY_TEMP"

# Prepare Heroku deployment files
if [[ "$DEPLOY_HEROKU" == "true" ]]; then
    echo -e "\n${YELLOW}Preparing Heroku deployment...${NC}"
    
    # Create a temporary directory for Heroku files
    HEROKU_TEMP="$DEPLOY_TEMP/heroku"
    mkdir -p "$HEROKU_TEMP"
    
    # Copy the CORS config file
    cp ./api/cors_config.py "$HEROKU_TEMP/"
    
    # Copy the test_cors.py file
    cp ./api/test_cors.py "$HEROKU_TEMP/"
    
    # Create or copy the Procfile if it doesn't exist
    if [ -f ./Procfile ]; then
        cp ./Procfile "$HEROKU_TEMP/"
    else
        echo "web: python test_cors.py" > "$HEROKU_TEMP/Procfile"
    fi
    
    # Create a minimal requirements.txt if it doesn't exist
    if [ -f ./requirements.txt ]; then
        cp ./requirements.txt "$HEROKU_TEMP/"
    else
        echo "flask==2.0.1" > "$HEROKU_TEMP/requirements.txt"
        echo "flask-cors==3.0.10" >> "$HEROKU_TEMP/requirements.txt"
        echo "gunicorn==20.1.0" >> "$HEROKU_TEMP/requirements.txt"
    fi
    
    # Create a gitignore file
    echo "venv/" > "$HEROKU_TEMP/.gitignore"
    echo "__pycache__/" >> "$HEROKU_TEMP/.gitignore"
    echo "*.py[cod]" >> "$HEROKU_TEMP/.gitignore"
    
    echo -e "${GREEN}Heroku deployment files prepared in $HEROKU_TEMP${NC}"
fi

# Deploy to Heroku
if [[ "$DEPLOY_HEROKU" == "true" ]]; then
    echo -e "\n${YELLOW}Deploying to Heroku...${NC}"
    
    # Initialize Git repository in the temp directory
    cd "$HEROKU_TEMP"
    git init
    git add .
    git commit -m "Deploy CORS configuration and test endpoints"
    
    # Check if the Heroku app exists, create it if not
    if ! heroku apps:info -a "$HEROKU_APP_NAME" &> /dev/null; then
        echo "Creating Heroku app $HEROKU_APP_NAME..."
        heroku apps:create "$HEROKU_APP_NAME"
    fi
    
    # Add Heroku remote
    heroku git:remote -a "$HEROKU_APP_NAME"
    
    # Deploy to Heroku
    git push heroku master --force
    
    echo -e "${GREEN}Deployed to Heroku app $HEROKU_APP_NAME${NC}"
    
    # Set CORS environment variables
    echo "Setting CORS configuration on Heroku..."
    heroku config:set FRONTEND_URL="https://$NETLIFY_SITE_NAME.netlify.app" -a "$HEROKU_APP_NAME"
    heroku config:set RDKIT_SERVICE_URL="https://rdkit.cryoprotect.app" -a "$HEROKU_APP_NAME"
    heroku config:set CONVEX_URL="https://dynamic-mink-63.convex.cloud" -a "$HEROKU_APP_NAME"
    
    echo -e "${GREEN}Heroku environment variables configured${NC}"
    
    # Return to original directory
    cd - > /dev/null
fi

# Deploy to Netlify
if [[ "$DEPLOY_NETLIFY" == "true" ]]; then
    echo -e "\n${YELLOW}Deploying to Netlify...${NC}"
    
    # Get into the frontend directory
    if [ -d "./frontend" ]; then
        cd ./frontend
        
        # Set environment variables in Netlify
        echo "Setting environment variables in Netlify..."
        netlify env:set NEXT_PUBLIC_API_URL "https://$HEROKU_APP_NAME.herokuapp.com" -s "$NETLIFY_SITE_NAME"
        netlify env:set NEXT_PUBLIC_RDKIT_API_URL "https://rdkit.cryoprotect.app" -s "$NETLIFY_SITE_NAME"
        netlify env:set NEXT_PUBLIC_CONVEX_URL "https://dynamic-mink-63.convex.cloud" -s "$NETLIFY_SITE_NAME"
        netlify env:set NEXT_PUBLIC_USE_CONVEX "true" -s "$NETLIFY_SITE_NAME"
        netlify env:set NEXT_PUBLIC_ENVIRONMENT "production" -s "$NETLIFY_SITE_NAME"
        
        # Build with Convex
        echo "Building with Convex..."
        npm run build:with-convex
        
        # Deploy to Netlify
        echo "Deploying to Netlify..."
        netlify deploy --prod -s "$NETLIFY_SITE_NAME"
        
        echo -e "${GREEN}Deployed to Netlify site $NETLIFY_SITE_NAME${NC}"
        
        # Return to original directory
        cd - > /dev/null
    else
        echo -e "${RED}Frontend directory not found. Skipping Netlify deployment.${NC}"
    fi
fi

# Clean up
echo -e "\n${YELLOW}Cleaning up...${NC}"
rm -rf "$DEPLOY_TEMP"

# Test the deployment
echo -e "\n${YELLOW}Testing the deployment...${NC}"
echo "Running connection tests..."
./test-connection.sh

echo -e "\n${BLUE}===========================================================${NC}"
echo -e "${GREEN}Deployment Complete!${NC}"
echo -e "${BLUE}===========================================================${NC}"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "  1. Deploy the RDKit service with CORS configuration"
echo "  2. Verify Convex integration by visiting:"
echo "     https://$NETLIFY_SITE_NAME.netlify.app/convex-test"
echo "  3. If any tests are still failing, use the troubleshooting guide:"
echo "     ./test-backend-integration.sh"
echo ""
echo -e "${GREEN}Your CryoProtect application has been deployed with updated integration.${NC}"