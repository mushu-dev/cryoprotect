#!/bin/bash

# Minimal deployment script for CryoProtect backend-frontend integration
# This focuses on Heroku and Netlify which are already installed

set -e  # Exit on error

# Color definitions for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;36m'
NC='\033[0m' # No Color

echo -e "${BLUE}===========================================================${NC}"
echo -e "${GREEN}CryoProtect Minimal Integration Deployment Script${NC}"
echo -e "${BLUE}===========================================================${NC}"
echo ""

# Configuration variables with defaults
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect}
NETLIFY_SITE_NAME=${NETLIFY_SITE_NAME:-cryoprotect}
DEPLOY_NETLIFY=${DEPLOY_NETLIFY:-false}  # Default to false for now
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

# Prepare Heroku deployment without file creation
if [[ "$DEPLOY_HEROKU" == "true" ]]; then
    echo -e "\n${YELLOW}Preparing Heroku deployment...${NC}"
    
    # Set CORS environment variables
    echo "Setting CORS configuration on Heroku..."
    heroku config:set FRONTEND_URL="https://$NETLIFY_SITE_NAME.netlify.app" -a "$HEROKU_APP_NAME"
    heroku config:set ALLOWED_ORIGINS="https://$NETLIFY_SITE_NAME.netlify.app,https://rdkit.cryoprotect.app,https://dynamic-mink-63.convex.cloud" -a "$HEROKU_APP_NAME"
    
    echo -e "${GREEN}Heroku environment variables configured${NC}"
    
    # Modify existing app.py to include CORS
    echo -e "\n${YELLOW}Note: You still need to update your Heroku app code to include:${NC}"
    echo "- The cors_config.py file"
    echo "- The test_cors.py endpoints"
    echo "- CORS configuration in your main app"
    echo -e "\nExample code:"
    echo -e "${BLUE}from flask_cors import CORS"
    echo -e "app = Flask(__name__)"
    echo -e "CORS(app, resources={r\"/*\": {\"origins\": \"*\"}})"
    echo -e "${NC}"
fi

# Deploy to Netlify
if [[ "$DEPLOY_NETLIFY" == "true" ]]; then
    echo -e "\n${YELLOW}Deploying to Netlify...${NC}"
    
    # Get into the frontend directory
    if [ -d "./frontend" ]; then
        cd ./frontend
        
        # Set environment variables in Netlify
        echo "Setting environment variables in Netlify..."
        netlify env:set NEXT_PUBLIC_API_URL "https://$HEROKU_APP_NAME.herokuapp.com"
        netlify env:set NEXT_PUBLIC_RDKIT_API_URL "https://rdkit.cryoprotect.app"
        netlify env:set NEXT_PUBLIC_CONVEX_URL "https://dynamic-mink-63.convex.cloud"
        netlify env:set NEXT_PUBLIC_USE_CONVEX "true"
        netlify env:set NEXT_PUBLIC_ENVIRONMENT "production"
        
        # Just deploy the netlify.toml file
        echo "Deploying netlify.toml changes to Netlify..."
        netlify deploy --prod
        
        echo -e "${GREEN}Deployed to Netlify site $NETLIFY_SITE_NAME${NC}"
        
        # Return to original directory
        cd - > /dev/null
    else
        echo -e "${RED}Frontend directory not found. Skipping Netlify deployment.${NC}"
    fi
fi

echo -e "\n${BLUE}===========================================================${NC}"
echo -e "${GREEN}Minimal Deployment Complete!${NC}"
echo -e "${BLUE}===========================================================${NC}"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "  1. Test the API configuration:"
echo "     curl -H 'Origin: https://$NETLIFY_SITE_NAME.netlify.app' https://$HEROKU_APP_NAME.herokuapp.com/api/v1/health/connectivity"
echo "  2. For the RDKit service, you'll need to set up fly.io:"
echo "     curl -L https://fly.io/install.sh | sh"
echo "     fly auth login"
echo "     ./deploy-rdkit-service.sh"
echo "  3. Verify Convex integration by visiting:"
echo "     https://$NETLIFY_SITE_NAME.netlify.app/convex-test"
echo ""
echo -e "${GREEN}Your CryoProtect basic integration has been configured.${NC}"