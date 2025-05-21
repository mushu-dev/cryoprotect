#!/bin/bash

# Create a temporary directory for output
mkdir -p playwright-test-output

# Run the Playwright test using our container solution
echo "Running Playwright test in container..."

# Copy the test file into the container
docker run --rm -v "$(pwd):/app" -w /app mcr.microsoft.com/playwright:latest \
  bash -c "
    # Install dependencies
    npm init -y
    npm install @playwright/test
    
    # Install browsers
    npx playwright install --with-deps chromium
    
    # Run the test
    npx playwright test test-cryoprotect.js
  "

echo "Test complete. Check the results in the playwright-test-output directory."