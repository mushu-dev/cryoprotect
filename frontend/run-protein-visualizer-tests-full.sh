#!/bin/bash

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Starting Protein Visualizer tests...${NC}"

# Create test results directory if it doesn't exist
mkdir -p ./test-results/protein-visualizer

# Check if server is already running
if nc -z localhost 3000 >/dev/null 2>&1; then
  echo -e "${YELLOW}Server already running on port 3000, using existing server${NC}"
  SERVER_ALREADY_RUNNING=true
else
  echo -e "${YELLOW}Starting Next.js development server...${NC}"
  # Start the Next.js server in the background
  npm run dev &
  SERVER_PID=$!
  
  # Store the PID in a file for cleanup later
  echo $SERVER_PID > .nextjs-server-pid
  
  # Wait for server to start (max 30 seconds)
  echo -e "${YELLOW}Waiting for server to start...${NC}"
  for i in {1..30}; do
    if nc -z localhost 3000 >/dev/null 2>&1; then
      echo -e "${GREEN}Server started successfully!${NC}"
      break
    fi
    echo -n "."
    sleep 1
    if [ $i -eq 30 ]; then
      echo -e "\n${RED}Server failed to start within 30 seconds. Aborting tests.${NC}"
      if [ -n "$SERVER_PID" ]; then
        kill $SERVER_PID
      fi
      exit 1
    fi
  done
  
  # Wait a bit more for the server to fully initialize
  echo -e "${YELLOW}Waiting a few more seconds for server to initialize...${NC}"
  sleep 5
fi

# Run the tests
echo -e "${YELLOW}Running Protein Visualizer tests...${NC}"
npx playwright test tests/protein-visualizer/visualizer.spec.js -c tests/protein-visualizer/playwright.config.js "$@"
TEST_EXIT_CODE=$?

# Cleanup if we started the server
if [ -z "$SERVER_ALREADY_RUNNING" ] && [ -f .nextjs-server-pid ]; then
  SERVER_PID=$(cat .nextjs-server-pid)
  echo -e "${YELLOW}Stopping Next.js server (PID: $SERVER_PID)...${NC}"
  kill $SERVER_PID
  rm .nextjs-server-pid
  echo -e "${GREEN}Server stopped successfully${NC}"
fi

# Print test results summary
if [ $TEST_EXIT_CODE -eq 0 ]; then
  echo -e "${GREEN}✓ All tests passed successfully!${NC}"
else
  echo -e "${RED}✗ Some tests failed. Check the test report for details.${NC}"
fi

# Exit with the test exit code
exit $TEST_EXIT_CODE