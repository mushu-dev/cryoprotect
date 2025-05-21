#!/bin/bash

# Create the test results directory if it doesn't exist
mkdir -p test-results

# First build the app
echo "Building the app..."
npm run build

# Run the Playwright test
echo "Running UI navigation tests..."
npx playwright test tests/ui-navigation.test.js --headed

# Check the result
if [ $? -eq 0 ]; then
  echo "✅ UI tests completed successfully!"
  echo "Screenshots have been saved to the test-results directory."
else
  echo "❌ UI tests failed. Check the output above for details."
fi