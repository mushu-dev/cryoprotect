#!/bin/bash

# Minimal Deployment Script for Netlify
# This script deploys a minimal version to diagnose deployment issues

set -e

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}==============================================${NC}"
echo -e "${BLUE}    CryoProtect Minimal Deployment Test    ${NC}"
echo -e "${BLUE}==============================================${NC}"
echo

# Check if Netlify CLI is installed
if ! command -v netlify &> /dev/null; then
    echo -e "${YELLOW}Netlify CLI not found. Installing...${NC}"
    npm install -g netlify-cli
fi

# Clean up any existing node_modules and .next directories
echo -e "${BLUE}Cleaning previous build artifacts...${NC}"
rm -rf node_modules .next package-lock.json
echo -e "${GREEN}Build artifacts cleaned${NC}"

# Step 1: Install dependencies
echo -e "${BLUE}Installing dependencies...${NC}"
npm install
echo -e "${GREEN}Dependencies installed${NC}"

# Step 2: Install the Netlify Next.js plugin
echo -e "${BLUE}Installing Netlify Next.js plugin...${NC}"
npm install -D @netlify/plugin-nextjs
echo -e "${GREEN}Netlify Next.js plugin installed${NC}"

# Step 3: Clear Netlify cache (if any)
echo -e "${BLUE}Clearing Netlify cache...${NC}"
if command -v netlify &> /dev/null; then
    netlify build:clear-cache || echo -e "${YELLOW}No cache to clear or command not available${NC}"
else
    echo -e "${YELLOW}Netlify CLI not found, skipping cache clear${NC}"
fi
echo -e "${GREEN}Cache cleared${NC}"

# Step 4: Deploy to Netlify
echo -e "${BLUE}Deploying to Netlify...${NC}"

# Check if user is logged in to Netlify
if ! netlify status 2>&1 | grep -q "Logged in"; then
    echo -e "${YELLOW}Not logged in to Netlify. Please log in:${NC}"
    netlify login
fi

# Check if site is linked
if ! netlify status 2>&1 | grep -q "cryoprotect"; then
    echo -e "${YELLOW}Site not linked. Linking to cryoprotect site...${NC}"
    netlify unlink 2>/dev/null
    netlify link --name cryoprotect
fi

# Deploy to Netlify with full build logs
echo -e "${BLUE}Starting Netlify deployment with verbose logging...${NC}"
NETLIFY_BUILD_DEBUG=true netlify deploy --build --prod --debug

echo -e "${BLUE}==============================================${NC}"
echo -e "${GREEN}Minimal test deployment completed!${NC}"
echo -e "${BLUE}==============================================${NC}"

echo -e "${YELLOW}Note: It may take a few minutes for the changes to propagate.${NC}"
echo -e "${YELLOW}Visit your Netlify site to verify that the deployment was successful.${NC}"
echo -e "${YELLOW}Check the following URLs:${NC}"
echo -e "${YELLOW}1. Main site URL (shown above)${NC}"
echo -e "${YELLOW}2. /health endpoint to verify the Netlify Function is working${NC}"
echo -e "${YELLOW}3. /api/hello endpoint to verify API routing${NC}"

echo -e "${BLUE}If this minimal version works, the issue is in the main application configuration.${NC}"
echo -e "${BLUE}If this still doesn't work, we'll need to check Netlify's function logs for errors.${NC}"

exit 0