#!/bin/bash

# Script to deploy the frontend with Convex integration
# This script handles both Convex and frontend deployment

set -e # Exit on error

# Colors for better output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}Starting deployment with Convex integration...${NC}"

# Check for convex CLI
if ! command -v convex &> /dev/null; then
    echo -e "${RED}Error: convex CLI is not installed.${NC}"
    echo -e "Please install it with: npm install -g convex"
    exit 1
fi

# Step 1: Generate Convex types
echo -e "${YELLOW}Step 1/5: Generating Convex types...${NC}"
npm run convex:codegen
if [ $? -ne 0 ]; then
    echo -e "${RED}Error generating Convex types.${NC}"
    exit 1
fi
echo -e "${GREEN}Convex types generated successfully.${NC}"

# Step 2: Deploy Convex functions
echo -e "${YELLOW}Step 2/5: Deploying Convex functions...${NC}"
npm run convex:deploy
if [ $? -ne 0 ]; then
    echo -e "${RED}Error deploying Convex functions.${NC}"
    exit 1
fi
echo -e "${GREEN}Convex functions deployed successfully.${NC}"

# Step 3: Build Next.js project with Convex enabled
echo -e "${YELLOW}Step 3/5: Building Next.js project with Convex...${NC}"
NEXT_PUBLIC_USE_CONVEX=true npm run build
if [ $? -ne 0 ]; then
    echo -e "${RED}Error building Next.js project.${NC}"
    exit 1
fi
echo -e "${GREEN}Next.js build successful.${NC}"

# Step 4: Generate static export (if needed)
echo -e "${YELLOW}Step 4/5: Generating static export...${NC}"
npm run export
if [ $? -ne 0 ]; then
    echo -e "${RED}Error generating static export.${NC}"
    echo -e "${YELLOW}Attempting to fix issues for Netlify deployment...${NC}"
    
    # Create an empty .nojekyll file to prevent GitHub Pages from ignoring files starting with underscore
    touch out/.nojekyll
    
    # Create Netlify configuration file in the output directory if it doesn't exist
    if [ ! -f out/_redirects ]; then
        echo "/* /index.html 200" > out/_redirects
        echo "Created Netlify _redirects file for SPA routing"
    fi
    
    echo -e "${YELLOW}Manual fixes applied for static export.${NC}"
else
    echo -e "${GREEN}Static export generated successfully.${NC}"
fi

# Step 5: Deploy to Netlify using Netlify CLI (if installed)
echo -e "${YELLOW}Step 5/5: Deploying to Netlify...${NC}"
if command -v netlify &> /dev/null; then
    netlify deploy --dir=out --prod
    if [ $? -ne 0 ]; then
        echo -e "${RED}Error deploying to Netlify.${NC}"
        echo -e "${YELLOW}You may need to deploy manually from the Netlify dashboard.${NC}"
    else
        echo -e "${GREEN}Deployed to Netlify successfully.${NC}"
    fi
else
    echo -e "${YELLOW}Netlify CLI not found. Please deploy manually:${NC}"
    echo -e "${BLUE}1. Visit the Netlify dashboard${NC}"
    echo -e "${BLUE}2. Upload the 'out' directory${NC}"
    echo -e "${BLUE}3. Ensure environment variables are set:${NC}"
    echo -e "   - NEXT_PUBLIC_CONVEX_URL=https://upbeat-parrot-866.convex.cloud"
    echo -e "   - NEXT_PUBLIC_USE_CONVEX=true"
fi

echo -e "\n${GREEN}Deployment preparation complete!${NC}"
echo -e "${BLUE}Your application with Convex integration is ready to be deployed.${NC}"
echo -e "${YELLOW}To test locally with Convex integration enabled:${NC}"
echo -e "npm run dev:with-convex\n"

echo -e "${BLUE}To view your application after deployment, visit:${NC}"
echo -e "https://cryoprotect.netlify.app"
echo -e "${BLUE}To check Convex status, visit:${NC}"
echo -e "https://cryoprotect.netlify.app/convex-status"