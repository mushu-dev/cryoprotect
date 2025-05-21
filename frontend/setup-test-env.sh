#!/bin/bash

# This script sets up the environment for running experimental data enhancement tests

echo "🔧 Setting up test environment..."

# Install npm dependencies if needed
if [ ! -d "node_modules" ]; then
  echo "Installing npm dependencies..."
  npm install
fi

# Install Playwright if needed
if [ ! -d "node_modules/@playwright" ]; then
  echo "Installing Playwright..."
  npm install --save-dev @playwright/test
  npx playwright install --with-deps
fi

# Create necessary directories
echo "Creating test directories..."
mkdir -p ./test-results
mkdir -p ./tests/playwright

# Check if test files already exist in the test directory
if [ ! -f "./tests/playwright/experiment-data-api.spec.js" ]; then
  echo "Copying test files to the proper location..."
  cp ./playwright/experiment-data-api.spec.js ./playwright/experiment-form-validation.spec.js ./playwright/experimental-data-ui.spec.js ./playwright/data-comparison.spec.js ./playwright/data-filtering.spec.js ./playwright/data-visualization.spec.js ./playwright/test-utils.js ./tests/playwright/
fi

# Make test scripts executable
chmod +x ./run-data-tests.sh ./validate-experimental-data-enhancement.sh

echo "✅ Test environment setup complete!"
echo "You can now run tests with:"
echo "  - npm run test:experimental-features  (e2e tests)"
echo "  - ./run-data-tests.sh                (specific data feature tests)" 
echo "  - ./validate-experimental-data-enhancement.sh  (full validation)"