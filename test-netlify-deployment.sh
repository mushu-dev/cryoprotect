#!/bin/bash
# Script to test the Netlify deployment with Playwright

# Get site name from argument or use default
NETLIFY_SITE_NAME=${1:-cryoprotect}
NETLIFY_URL="https://${NETLIFY_SITE_NAME}.netlify.app"

echo "Testing Netlify deployment at: $NETLIFY_URL"

# First method: Use our containerized Playwright solution
echo "Method 1: Using containerized Playwright solution..."
./mcp-playwright-final.sh browser_navigate "$NETLIFY_URL"
./mcp-playwright-final.sh browser_take_screenshot "$NETLIFY_URL" "netlify-deployment.png"

echo "Screenshot saved as netlify-deployment.png"

# Second method: Use Playwright test script
echo "Method 2: Running Playwright test script..."
cd frontend
NETLIFY_URL="$NETLIFY_URL" npx playwright test playwright/netlify-connectivity.spec.js --headed

# Check results
if [ $? -eq 0 ]; then
  echo "✅ Netlify deployment tests passed!"
else
  echo "❌ Netlify deployment tests failed!"
fi

echo "Open the screenshots to see the results:"
echo "- netlify-deployment.png (from containerized solution)"
echo "- frontend/netlify-homepage.png (from Playwright test)"
echo "- frontend/netlify-api-test.png (from Playwright test)"