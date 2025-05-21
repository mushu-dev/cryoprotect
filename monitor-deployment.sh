#!/bin/bash

# CryoProtect Deployment Monitoring Script
set -e

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}==============================================${NC}"
echo -e "${BLUE}    CryoProtect Deployment Monitoring    ${NC}"
echo -e "${BLUE}==============================================${NC}"
echo

check_url() {
  local url=$1
  local description=$2
  
  echo -e "${BLUE}Checking $description at $url...${NC}"
  
  # Use curl to check the URL
  status_code=$(curl -s -o /dev/null -w "%{http_code}" "$url")
  
  if [ "$status_code" -ge 200 ] && [ "$status_code" -lt 300 ]; then
    echo -e "${GREEN}✓ $description is working! (Status: $status_code)${NC}"
    return 0
  else
    echo -e "${RED}✗ $description is NOT working! (Status: $status_code)${NC}"
    return 1
  fi
}

# Check main site and key endpoints
check_url "https://cryoprotect.app" "Main Website"
check_url "https://cryoprotect.app/api/hello" "API Hello Endpoint"
check_url "https://cryoprotect.app/health" "Health Endpoint"

# Check backend service integration
check_url "https://cryoprotect-8030e4025428.herokuapp.com/health" "Heroku Backend"
check_url "https://cryoprotect-rdkit.fly.dev/health" "RDKit Service"

echo
echo -e "${BLUE}==============================================${NC}"
echo -e "${BLUE}    Netlify Function Logs    ${NC}"
echo -e "${BLUE}==============================================${NC}"
echo

echo -e "${YELLOW}Access Netlify function logs at:${NC}"
echo -e "https://app.netlify.com/projects/cryoprotect/logs/functions"
echo

echo -e "${BLUE}==============================================${NC}"
echo -e "${BLUE}    Next Steps    ${NC}"
echo -e "${BLUE}==============================================${NC}"
echo

echo -e "1. If the main website is working:${NC}"
echo -e "   - Implement Convex database integration${NC}"
echo -e "   - Restore full application functionality${NC}"
echo
echo -e "2. If the main website is still showing fallback:${NC}"
echo -e "   - Check Netlify function logs for errors${NC}"
echo -e "   - Check .next/server directory for problems${NC}"
echo -e "   - Verify NextAuth configuration${NC}"
echo
echo -e "3. For Netlify Edge Function issues:${NC}"
echo -e "   - Try disabling edge functions: NEXT_DISABLE_EDGE_IMAGES=true${NC}"
echo -e "   - Update plugins to latest version${NC}"

exit 0