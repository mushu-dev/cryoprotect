#!/bin/bash

# Color output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Running Protein Visualizer validation tests...${NC}"

# Check for required packages
echo -e "${YELLOW}Installing required testing dependencies...${NC}"
npm install --save-dev jest@29.5.0

# Run Jest test for the component using our config
echo -e "${YELLOW}Executing tests...${NC}"
npx jest --config=jest.config.js

# Check if tests passed
if [ $? -eq 0 ]; then
  echo -e "${GREEN}✓ Validation tests passed successfully!${NC}"
  echo -e "${YELLOW}The component passes a basic code validation.${NC}"
  echo -e "${YELLOW}For a detailed manual test, run:${NC}"
  echo -e "npm run dev"
  echo -e "Then navigate to: http://localhost:3000/protein-visualizer-demo"
else
  echo -e "${RED}✗ Validation tests failed.${NC}"
  echo -e "${YELLOW}Please check the component for issues:${NC}"
  echo -e "src/components/protein-visualizer/MolstarViewer.tsx"
fi