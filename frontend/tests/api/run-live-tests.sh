#!/bin/bash

# This script runs integration tests against the live backend services
# (Heroku API and Fly.io RDKit service)

# Set colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}Running integration tests for CryoProtect backend services${NC}"

# Create results directory
mkdir -p ./test-results/api-live

# Allow setting custom API URLs through environment variables
if [ -n "$CUSTOM_API_URL" ]; then
  echo -e "${YELLOW}Using custom API URL: $CUSTOM_API_URL${NC}"
  export NEXT_PUBLIC_API_URL="$CUSTOM_API_URL"
fi

if [ -n "$CUSTOM_RDKIT_URL" ]; then
  echo -e "${YELLOW}Using custom RDKit URL: $CUSTOM_RDKIT_URL${NC}"
  export NEXT_PUBLIC_RDKIT_API_URL="$CUSTOM_RDKIT_URL"
fi

# Run the live API tests
echo -e "${YELLOW}Running live API tests...${NC}"
node tests/api/live-api-test.js

# Check exit status
if [ $? -eq 0 ]; then
  echo -e "${GREEN}Integration tests completed successfully.${NC}"
  exit 0
else
  echo -e "${RED}Integration tests failed.${NC}"
  exit 1
fi