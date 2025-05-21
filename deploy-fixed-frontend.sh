#!/bin/bash

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Starting frontend fix deployment process...${NC}"

# 1. Navigate to the frontend directory
cd frontend

# 2. Install dependencies if needed
echo -e "${YELLOW}Installing dependencies...${NC}"
npm install

# 3. Build and export the frontend
echo -e "${YELLOW}Building and exporting the application...${NC}"
npm run build

echo -e "${YELLOW}Exporting static HTML...${NC}"
npm run export

# Check if export was successful
if [ -d "out" ]; then
  echo -e "${GREEN}Export successful! Generated files:${NC}"
  ls -la out/
  
  # Fix the HTML file extensions for Netlify
  echo -e "${YELLOW}Fixing file extensions for Netlify...${NC}"
  for file in out/*/index.html; do
    dir=$(dirname "$file")
    base=$(basename "$dir")
    if [ "$base" != "out" ]; then
      cp "$file" "out/$base.html"
      echo "Created out/$base.html from $file"
    fi
  done
else
  echo -e "${RED}Export failed! out/ directory not found.${NC}"
  exit 1
fi

# 4. Deploy to Netlify (if netlify-cli is installed)
if command -v netlify &> /dev/null; then
  echo -e "${YELLOW}Deploying to Netlify...${NC}"
  netlify deploy --dir=out --prod
else
  echo -e "${RED}Netlify CLI not found. Manual deployment required.${NC}"
  echo -e "${YELLOW}To manually deploy:${NC}"
  echo "1. Install Netlify CLI: npm install -g netlify-cli"
  echo "2. Login to Netlify: netlify login"
  echo "3. Deploy the site: netlify deploy --dir=frontend/out --prod"
fi

# 5. Verify the deployment
echo -e "${YELLOW}Waiting for deployment to complete...${NC}"
sleep 10  # Give Netlify some time to complete the deployment

echo -e "${YELLOW}Verifying routes...${NC}"
node scripts/verify-routes.js

echo -e "${GREEN}Deployment process completed!${NC}"