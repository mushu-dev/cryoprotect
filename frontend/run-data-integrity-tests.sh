#!/bin/bash
# Run all data integrity and user flow tests

# Start the development server in the background if it's not already running
SERVER_RUNNING=$(ps aux | grep "next dev" | grep -v grep | wc -l)
if [ $SERVER_RUNNING -eq 0 ]; then
  echo "Starting Next.js development server..."
  npm run dev &
  SERVER_PID=$!
  # Give the server some time to start
  echo "Waiting for server to start..."
  sleep 10
fi

echo "Running data integrity and user flow tests..."
npx playwright test tests/e2e/user-flow.spec.js tests/e2e/data-integrity.spec.js --project=chromium --headed "$@"

# Capture test result
TEST_RESULT=$?

# Check if tests ran successfully
if [ $TEST_RESULT -eq 0 ]; then
  echo "✅ All data integrity and user flow tests passed!"
else
  echo "❌ Some tests failed. Please check the report for details."
  npx playwright show-report
fi

# Shut down the dev server if we started it
if [ $SERVER_RUNNING -eq 0 ] && [ ! -z "$SERVER_PID" ]; then
  echo "Shutting down development server..."
  kill $SERVER_PID
fi

exit $TEST_RESULT