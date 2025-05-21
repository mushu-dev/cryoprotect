#!/bin/bash

# This script runs security and performance tests for the experimental data enhancement UI

echo "🔒 Running security and performance tests for the experimental data enhancement UI"

# Create directories for test results
mkdir -p ./test-results/security
mkdir -p ./test-results/performance
mkdir -p ./test-results/cross-browser

# Install dependencies if needed
if [ ! -d "node_modules/@playwright" ]; then
  echo "Installing Playwright dependencies..."
  npm install
  npx playwright install --with-deps
fi

# Run security tests
echo "Running security tests..."
npx playwright test security-testing.spec.js --reporter=list

# Run performance tests
echo "Running performance tests..."
npx playwright test performance-metrics.spec.js --reporter=list

# Run load tests with reduced settings for local usage
echo "Running load tests (with reduced settings for local testing)..."
REDUCED_LOAD=true npx playwright test load-testing.spec.js --reporter=list

# Run cross-browser tests on available browsers
echo "Running cross-browser compatibility tests..."
npx playwright test cross-browser-compatibility.spec.js --reporter=list

# Generate a timestamp for the report directory
TIMESTAMP=$(date +"%Y%m%d-%H%M%S")
REPORT_DIR="./test-results/security-performance-report-${TIMESTAMP}"

# Create summary report
mkdir -p "${REPORT_DIR}"
SUMMARY_FILE="${REPORT_DIR}/summary.md"

echo "# Security and Performance Test Summary" > "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "Tests run on: $(date)" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Security Tests" >> "${SUMMARY_FILE}"
echo "- Checked security headers" >> "${SUMMARY_FILE}"
echo "- Tested input validation" >> "${SUMMARY_FILE}"
echo "- Checked for exposed sensitive information" >> "${SUMMARY_FILE}"
echo "- Tested authentication and authorization" >> "${SUMMARY_FILE}"
echo "- Verified CSRF protection" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Performance Tests" >> "${SUMMARY_FILE}"
echo "- Measured page load times" >> "${SUMMARY_FILE}"
echo "- Tested rendering efficiency" >> "${SUMMARY_FILE}"
echo "- Verified handling of large datasets" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Load Tests" >> "${SUMMARY_FILE}"
echo "- Simulated concurrent users" >> "${SUMMARY_FILE}"
echo "- Tested form submissions under load" >> "${SUMMARY_FILE}"
echo "- Verified API response times" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Cross-Browser Tests" >> "${SUMMARY_FILE}"
echo "- Tested UI rendering across browsers" >> "${SUMMARY_FILE}"
echo "- Verified data visualization compatibility" >> "${SUMMARY_FILE}"
echo "- Checked responsive layouts" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "See the Playwright report for detailed results." >> "${SUMMARY_FILE}"

echo "✅ All tests complete!"
echo "Summary report available at: ${REPORT_DIR}/summary.md"

# Try to open the summary report
if [ -x "$(command -v xdg-open)" ]; then
  xdg-open "${SUMMARY_FILE}"
elif [ -x "$(command -v open)" ]; then
  open "${SUMMARY_FILE}"
else
  echo "To view the summary report, open ${SUMMARY_FILE} in your browser or text editor"
fi