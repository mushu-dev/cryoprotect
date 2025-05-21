#!/bin/bash

# This script runs manual verification tests for the Molecule and Mixture APIs

# Set colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}Running manual verification tests for Molecule and Mixture APIs${NC}"

# Create results directory
mkdir -p ./test-results/api

# Run the manual verification script
echo -e "${YELLOW}Running API verification script...${NC}"
node tests/api/verify-molecule-mixture-api.js | tee ./test-results/api/verification-results.txt

# Check exit status of the script
if [ $? -eq 0 ]; then
  echo -e "${GREEN}Manual verification complete.${NC}"
  echo -e "${GREEN}Results saved to ./test-results/api/verification-results.txt${NC}"
else
  echo -e "${RED}Error running verification script.${NC}"
  exit 1
fi

exit 0