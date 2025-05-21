#!/bin/bash

# This script runs visual regression tests and prepares the directories

# Create test results directory if it doesn't exist
mkdir -p ./test-results/visual

echo "🎨 Running visual regression tests..."

# Run the visual regression tests in headed mode for accurate rendering
npx playwright test visual-regression.spec.js --headed

echo "✅ Visual regression tests completed"
echo "Screenshots saved to ./test-results/visual/"

# Optional: Open the directory with the screenshots
if [ -x "$(command -v xdg-open)" ]; then
  xdg-open ./test-results/visual/
elif [ -x "$(command -v open)" ]; then
  open ./test-results/visual/
else
  echo "Screenshots are available in ./test-results/visual/"
fi