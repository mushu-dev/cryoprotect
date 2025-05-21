#!/bin/bash

# This script installs Playwright and its dependencies

echo "🎭 Setting up Playwright for E2E testing"

# Install Playwright package if not already installed
if ! npm list --depth=0 | grep -q "@playwright/test"; then
  echo "Installing @playwright/test package..."
  npm install --save-dev @playwright/test
else
  echo "@playwright/test already installed."
fi

# Install browser binaries and dependencies
echo "Installing browser binaries and dependencies..."
npx playwright install --with-deps

echo "✅ Playwright setup complete!"
echo "You can now run tests using:"
echo "  ./run-e2e-tests.sh          # Run tests with default settings"
echo "  ./run-e2e-tests.sh firefox headed  # Run tests on Firefox in headed mode"
echo "  npm run test:e2e            # Run all tests with npm script"