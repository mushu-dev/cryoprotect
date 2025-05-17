#!/bin/bash
# Deploy to Netlify with analytics validation
set -e  # Exit on any error

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Starting Netlify deployment process...${NC}"

# Step 1: Validate analytics implementation
echo -e "\n${YELLOW}Step 1: Validating analytics implementation...${NC}"
./validate-analytics.js

# Check if validation was successful
if [ $? -ne 0 ]; then
  echo -e "${RED}Analytics validation failed. Please fix issues before deploying.${NC}"
  exit 1
fi

# Step 2: Navigate to frontend directory
cd frontend

# Step 3: Deploy to Netlify
echo -e "\n${YELLOW}Step 3: Deploying to Netlify (preview)...${NC}"
echo -e "Running: ${GREEN}netlify deploy${NC}"

npm run install-deps

npm run update-api-endpoints

npm run deploy:netlify

# Return to original directory
cd ..

echo -e "\n${GREEN}Deployment process completed!${NC}"
echo -e "Don't forget to enable Netlify Analytics in your Netlify dashboard:"
echo -e "${GREEN}https://app.netlify.com/sites/cryoprotect/analytics${NC}"
echo -e "It will take up to 24 hours for analytics data to start appearing."
echo -e "\nCheck ${GREEN}NETLIFY_DEPLOYMENT_WITH_ANALYTICS.md${NC} for more information and troubleshooting tips."