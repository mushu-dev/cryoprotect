#!/bin/bash

# This script runs Playwright E2E tests with various options
# Usage: ./run-e2e-tests.sh [chromium|firefox|webkit|all] [headless|headed]

# Default values
BROWSER=${1:-"chromium"}
MODE=${2:-"headless"}

echo "🧪 Running E2E tests with Playwright"
echo "Browser: $BROWSER"
echo "Mode: $MODE"

# Make sure the required dependencies are installed
echo "Installing Playwright dependencies if needed..."
npx playwright install --with-deps

# Determine the headed mode flag
HEADED_FLAG=""
if [ "$MODE" == "headed" ]; then
  HEADED_FLAG="--headed"
fi

# Run tests based on browser selection
if [ "$BROWSER" == "all" ]; then
  echo "Running tests on all browsers..."
  npx playwright test $HEADED_FLAG
elif [ "$BROWSER" == "chromium" ] || [ "$BROWSER" == "firefox" ] || [ "$BROWSER" == "webkit" ]; then
  echo "Running tests on $BROWSER..."
  npx playwright test --project=$BROWSER $HEADED_FLAG
else
  echo "Invalid browser specified. Choose from: chromium, firefox, webkit, all"
  exit 1
fi

# Open the HTML report after tests complete
echo "Opening test report..."
npx playwright show-report