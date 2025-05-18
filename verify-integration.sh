#!/bin/bash
# Verification script to check the status of all CryoProtect components

set -e

# ANSI color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}===== CryoProtect Integration Verification =====${NC}"
echo

# Check Netlify
echo -e "${BLUE}Checking Netlify frontend...${NC}"
if curl -s https://cryoprotect.netlify.app/ | grep -q "CryoProtect"; then
    echo -e "${GREEN}✓ Netlify frontend is accessible${NC}"
else
    echo -e "${RED}✗ Netlify frontend check failed${NC}"
fi

# Check Heroku
echo -e "${BLUE}Checking Heroku backend...${NC}"
if curl -s https://cryoprotect.herokuapp.com/v1/health | grep -q "ok"; then
    echo -e "${GREEN}✓ Heroku backend is accessible${NC}"
else
    echo -e "${RED}✗ Heroku backend check failed${NC}"
fi

# Check RDKit Service
echo -e "${BLUE}Checking RDKit service...${NC}"
if curl -s https://rdkit.cryoprotect.app/health | grep -q "ok"; then
    echo -e "${GREEN}✓ RDKit service is accessible${NC}"
else
    echo -e "${YELLOW}⚠ RDKit service check failed - may need to be deployed${NC}"
fi

# Check Convex connection through backend
echo -e "${BLUE}Checking Convex through backend...${NC}"
if curl -s https://cryoprotect.herokuapp.com/v1/molecules?limit=1 | grep -q "data"; then
    echo -e "${GREEN}✓ Backend can successfully query Convex${NC}"
else
    echo -e "${RED}✗ Backend Convex connection check failed${NC}"
fi

# Check CORS configuration for Heroku
echo -e "${BLUE}Checking CORS configuration for Heroku...${NC}"
if curl -s -I -H "Origin: https://cryoprotect.netlify.app" \
        -H "Access-Control-Request-Method: GET" \
        https://cryoprotect.herokuapp.com/v1/health | grep -q "Access-Control-Allow-Origin"; then
    echo -e "${GREEN}✓ Heroku CORS is properly configured${NC}"
else
    echo -e "${RED}✗ Heroku CORS check failed${NC}"
fi

# Check frontend environment variables
echo -e "${BLUE}Checking frontend environment variables...${NC}"
if curl -s https://cryoprotect.netlify.app/ | grep -q "dynamic-mink-63.convex.cloud"; then
    echo -e "${GREEN}✓ Frontend includes Convex configuration${NC}"
else
    echo -e "${YELLOW}⚠ Could not verify Convex configuration in frontend${NC}"
fi

echo
echo -e "${BLUE}==========================================${NC}"
echo -e "${GREEN}Integration verification complete!${NC}"
echo -e "${BLUE}==========================================${NC}"
echo
echo "For more detailed testing, run:"
echo "  node test-full-integration.js"
echo
echo "For a comprehensive test of the adapter, run:"
echo "  python test-convex-adapter.py"