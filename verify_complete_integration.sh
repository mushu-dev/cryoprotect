#!/bin/bash

# Verify Complete Integration
# This script tests all components of the CryoProtect application to verify the Convex integration

set -e

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}==============================================${NC}"
echo -e "${BLUE}    CryoProtect Convex Integration Verifier   ${NC}"
echo -e "${BLUE}==============================================${NC}"
echo

# Configuration
CONVEX_URL="https://upbeat-parrot-866.convex.cloud"
RDKIT_SERVICE_URL="https://cryoprotect-rdkit.fly.dev"
BACKEND_URL="https://cryoprotect-8030e4025428.herokuapp.com"
FRONTEND_URL="https://cryoprotect.netlify.app"

# 1. Test Convex API
echo -e "${BLUE}Testing Convex API...${NC}"
if curl -s -o /dev/null -w "%{http_code}" $CONVEX_URL | grep -q "200\|404"; then
    echo -e "  ${GREEN}✓ Convex API is accessible${NC}"
else
    echo -e "  ${RED}✗ Convex API is not accessible${NC}"
    echo -e "  ${YELLOW}(404 is actually expected here as it's not a direct HTTP endpoint)${NC}"
fi

# 2. Test RDKit Service
echo -e "\n${BLUE}Testing RDKit Service...${NC}"
HEALTH_RESPONSE=$(curl -s $RDKIT_SERVICE_URL/health)
if echo "$HEALTH_RESPONSE" | grep -q "status.*ok"; then
    echo -e "  ${GREEN}✓ RDKit Service is healthy${NC}"
    echo -e "  Response: $HEALTH_RESPONSE"
else
    echo -e "  ${RED}✗ RDKit Service health check failed${NC}"
    echo -e "  Response: $HEALTH_RESPONSE"
fi

# 3. Test CORS Configuration
echo -e "\n${BLUE}Testing CORS Configuration...${NC}"
CORS_RESPONSE=$(curl -s -H "Origin: $FRONTEND_URL" $RDKIT_SERVICE_URL/test-cors)
if echo "$CORS_RESPONSE" | grep -q "success.*true"; then
    echo -e "  ${GREEN}✓ CORS is properly configured${NC}"
    echo -e "  Allowed origins: $(echo "$CORS_RESPONSE" | grep -o "allowed_origins.*" | cut -d']' -f1)]"
else
    echo -e "  ${RED}✗ CORS test failed${NC}"
    echo -e "  Response: $CORS_RESPONSE"
fi

# 4. Test Backend API
echo -e "\n${BLUE}Testing Backend API...${NC}"
BACKEND_RESPONSE=$(curl -s $BACKEND_URL/health)
if echo "$BACKEND_RESPONSE" | grep -q "status.*ok"; then
    echo -e "  ${GREEN}✓ Backend API is healthy${NC}"
    echo -e "  Response: $BACKEND_RESPONSE"
else
    echo -e "  ${RED}✗ Backend API health check failed${NC}"
    echo -e "  Response: $BACKEND_RESPONSE"
fi

# 5. Test Frontend Accessibility
echo -e "\n${BLUE}Testing Frontend Accessibility...${NC}"
if curl -s -I $FRONTEND_URL | grep -q "HTTP/"; then
    echo -e "  ${GREEN}✓ Frontend is accessible${NC}"
else
    echo -e "  ${RED}✗ Frontend is not accessible${NC}"
    echo -e "  ${YELLOW}NOTE: This might fail if the Netlify deployment is not complete.${NC}"
fi

# 6. Test End-to-End Connectivity
echo -e "\n${BLUE}Testing End-to-End Connectivity...${NC}"
echo -e "  ${YELLOW}This would require authenticating with the application and making real requests.${NC}"
echo -e "  ${YELLOW}Manual testing is recommended for full end-to-end verification.${NC}"

# Summary
echo -e "\n${BLUE}==============================================${NC}"
echo -e "${BLUE}                 Summary                     ${NC}"
echo -e "${BLUE}==============================================${NC}"
echo -e "  Convex API: https://upbeat-parrot-866.convex.cloud"
echo -e "  RDKit Service: https://cryoprotect-rdkit.fly.dev"
echo -e "  Backend API: https://cryoprotect-8030e4025428.herokuapp.com"
echo -e "  Frontend: https://cryoprotect.netlify.app"
echo
echo -e "${GREEN}Convex integration verification complete!${NC}"
echo -e "${YELLOW}Manual testing with a real browser is recommended for final verification.${NC}"
echo -e "${BLUE}==============================================${NC}"

# Exit with success status
exit 0