#!/bin/bash

# Deployment script for RDKit service with CORS configuration
# This script prepares and deploys the RDKit service to fly.io

set -e  # Exit on error

# Color definitions for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;36m'
NC='\033[0m' # No Color

echo -e "${BLUE}===========================================================${NC}"
echo -e "${GREEN}RDKit Service Deployment Script${NC}"
echo -e "${BLUE}===========================================================${NC}"
echo ""

# Check for fly CLI
if ! command -v flyctl &> /dev/null && ! command -v fly &> /dev/null; then
    echo -e "${RED}fly.io CLI not found. Please install it:${NC}"
    echo "curl -L https://fly.io/install.sh | sh"
    exit 1
fi

# Use flyctl or fly, whichever is available
FLY_CMD="fly"
if ! command -v fly &> /dev/null; then
    FLY_CMD="flyctl"
fi

# Check authentication
echo "Checking fly.io authentication..."
if ! $FLY_CMD auth whoami &> /dev/null; then
    echo -e "${RED}Not authenticated with fly.io. Please run:${NC}"
    echo "$FLY_CMD auth login"
    exit 1
fi
echo -e "${GREEN}Authenticated with fly.io.${NC}"

# App name
APP_NAME=${APP_NAME:-cryoprotect-rdkit}
echo -e "${YELLOW}Using fly.io app:${NC} $APP_NAME"
echo ""

# If environment variable AUTO_DEPLOY is set to 1, skip confirmation
if [[ "$AUTO_DEPLOY" != "1" ]]; then
    # Confirm configuration
    read -p "Do you want to continue with this configuration? (y/n) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${RED}Deployment aborted by user.${NC}"
        exit 0
    fi
else
    echo -e "${GREEN}Auto-deployment enabled. Continuing without confirmation...${NC}"
fi

# Create deployment temp directory
echo -e "\n${YELLOW}Creating deployment files...${NC}"
DEPLOY_TEMP=$(mktemp -d)
echo "Using temporary directory: $DEPLOY_TEMP"

# Prepare fly.io deployment files
echo -e "\n${YELLOW}Preparing fly.io deployment...${NC}"

# Copy the CORS config file
cp ./rdkit-service/cors_config.py "$DEPLOY_TEMP/"

# Copy the test_cors.py file
cp ./rdkit-service/test_cors.py "$DEPLOY_TEMP/"

# Create a Dockerfile
cat > "$DEPLOY_TEMP/Dockerfile" << EOF
FROM python:3.9-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8080

CMD ["python", "test_cors.py"]
EOF

# Create a fly.toml file
cat > "$DEPLOY_TEMP/fly.toml" << EOF
app = "$APP_NAME"

[build]
  dockerfile = "Dockerfile"

[env]
  FRONTEND_URL = "https://cryoprotect.netlify.app"
  API_URL = "https://cryoprotect-8030e4025428.herokuapp.com"
  CONVEX_URL = "https://upbeat-parrot-866.convex.cloud"

[http_service]
  internal_port = 8080
  force_https = true
  auto_stop_machines = true
  auto_start_machines = true
  min_machines_running = 0

  [http_service.concurrency]
    type = "connections"
    hard_limit = 1000
    soft_limit = 500
EOF

# Create a requirements.txt file
cat > "$DEPLOY_TEMP/requirements.txt" << EOF
flask==2.0.1
flask-cors==3.0.10
gunicorn==20.1.0
werkzeug==2.0.1
EOF

echo -e "${GREEN}fly.io deployment files prepared in $DEPLOY_TEMP${NC}"

# Deploy to fly.io
echo -e "\n${YELLOW}Deploying to fly.io...${NC}"

# Change to deployment directory
cd "$DEPLOY_TEMP"

# Create app if it doesn't exist
if ! $FLY_CMD apps list | grep -q "$APP_NAME"; then
    echo "Creating fly.io app $APP_NAME..."
    $FLY_CMD apps create "$APP_NAME"
fi

# Deploy the app
echo "Deploying to fly.io..."
$FLY_CMD deploy --remote-only

echo -e "${GREEN}Deployed to fly.io app $APP_NAME${NC}"

# Set secrets
echo "Setting environment variables in fly.io..."
$FLY_CMD secrets set FRONTEND_URL="https://cryoprotect.netlify.app" -a "$APP_NAME"
$FLY_CMD secrets set API_URL="https://cryoprotect-8030e4025428.herokuapp.com" -a "$APP_NAME"
$FLY_CMD secrets set CONVEX_URL="https://upbeat-parrot-866.convex.cloud" -a "$APP_NAME"

echo -e "${GREEN}fly.io environment variables configured${NC}"

# Return to original directory
cd - > /dev/null

# Clean up
echo -e "\n${YELLOW}Cleaning up...${NC}"
rm -rf "$DEPLOY_TEMP"

echo -e "\n${BLUE}===========================================================${NC}"
echo -e "${GREEN}RDKit Service Deployment Complete!${NC}"
echo -e "${BLUE}===========================================================${NC}"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "  1. Test the RDKit service:"
echo "     curl https://$APP_NAME.fly.dev/health"
echo "  2. Test CORS configuration:"
echo "     curl -H 'Origin: https://cryoprotect.netlify.app' https://$APP_NAME.fly.dev/test-cors"
echo "  3. Update DNS if needed:"
echo "     $FLY_CMD certs create rdkit.cryoprotect.app"
echo ""
echo -e "${GREEN}Your RDKit service has been deployed with CORS configuration.${NC}"