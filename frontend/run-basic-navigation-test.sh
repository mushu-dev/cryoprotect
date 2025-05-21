#!/bin/bash

# Create test results directory if it doesn't exist
mkdir -p ./test-results/basic-navigation

# Run the navigation test
echo "Running basic navigation test..."
npx playwright test playwright/basic-navigation.spec.js --headed

# Check if the test completed successfully
if [ $? -eq 0 ]; then
  echo "✅ Test completed successfully!"
  echo "Screenshots saved in ./test-results/basic-navigation/"
else
  echo "❌ Test failed! Check the output above for errors."
fi