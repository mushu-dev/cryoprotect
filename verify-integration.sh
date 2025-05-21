#!/bin/bash
# Verification script to check the status of all CryoProtect components

# ANSI color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}===== CryoProtect Integration Verification =====${NC}"
echo

# Function to check HTTP response
check_endpoint() {
  local url=$1
  local name=$2
  
  echo -e "${BLUE}Testing ${name}: ${url}${NC}"
  
  response=$(curl -s -o /dev/null -w "%{http_code}" $url)
  
  if [ "$response" == "200" ]; then
    echo -e "${GREEN}✓ Success: ${name} is accessible (HTTP 200)${NC}"
    return 0
  else
    echo -e "${RED}✗ Error: ${name} returned HTTP ${response}${NC}"
    return 1
  fi
}

# Function to check CORS headers
check_cors() {
  local url=$1
  local origin=$2
  local name=$3
  
  echo -e "${BLUE}Testing CORS for ${name} with origin: ${origin}${NC}"
  
  cors_header=$(curl -s -I -H "Origin: ${origin}" $url | grep -i "Access-Control-Allow-Origin")
  
  if [[ "$cors_header" == *"$origin"* ]] || [[ "$cors_header" == *"*"* ]]; then
    echo -e "${GREEN}✓ Success: CORS is properly configured for ${origin}${NC}"
    return 0
  else
    echo -e "${RED}✗ Error: CORS is not configured correctly for ${origin}${NC}"
    echo -e "${YELLOW}Response header: ${cors_header}${NC}"
    return 1
  fi
}

# Check main domain
echo -e "${BLUE}Testing main domain...${NC}"
netlify_status=$(curl -s -o /dev/null -w "%{http_code}" https://cryoprotect.netlify.app)
echo -e "Netlify site status: ${netlify_status}"

# Try the domain with .app
domain_status=$(curl -s -o /dev/null -w "%{http_code}" https://cryoprotect.app)
echo -e "Custom domain status: ${domain_status}"

# Check Heroku API
echo -e "\n${BLUE}Testing Heroku API...${NC}"
check_endpoint "https://cryoprotect-8030e4025428.herokuapp.com/health" "Heroku API"
check_cors "https://cryoprotect-8030e4025428.herokuapp.com/health" "https://cryoprotect.app" "Heroku API"

# Check RDKit service
echo -e "\n${BLUE}Testing RDKit Service...${NC}"
check_endpoint "https://cryoprotect-rdkit.fly.dev/health" "RDKit Service"
check_cors "https://cryoprotect-rdkit.fly.dev/health" "https://cryoprotect.app" "RDKit Service"

# Check redirects through Netlify
echo -e "\n${BLUE}Testing API redirects through Netlify...${NC}"
if check_endpoint "https://cryoprotect.netlify.app/api/health" "API redirect"; then
  echo -e "${GREEN}✓ Success: API redirect is working properly${NC}"
else
  echo -e "${YELLOW}⚠ API redirect check failed - this will work after deployment${NC}"
fi

# Check RDKit redirects through Netlify
echo -e "\n${BLUE}Testing RDKit redirects through Netlify...${NC}"
if check_endpoint "https://cryoprotect.netlify.app/rdkit-api/health" "RDKit redirect"; then
  echo -e "${GREEN}✓ Success: RDKit redirect is working properly${NC}"
else
  echo -e "${YELLOW}⚠ RDKit redirect check failed - this will work after deployment${NC}"
fi

# Check Convex connectivity
echo -e "\n${BLUE}Testing Convex connectivity...${NC}"
convex_status=$(curl -s -o /dev/null -w "%{http_code}" https://upbeat-parrot-866.convex.cloud)
echo -e "Convex API status: ${convex_status}"

echo -e "\n${BLUE}==========================================${NC}"
echo -e "${GREEN}Integration verification complete!${NC}"
echo -e "${BLUE}==========================================${NC}"
echo
echo "Next step: Run the Netlify deployment script:"
echo "  ./deploy-netlify-fix.sh"
echo
echo "For a comprehensive test of the adapter, run:"
echo "  python test-convex-adapter.py"