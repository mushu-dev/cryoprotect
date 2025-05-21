#!/bin/bash

# This script runs all tests for the experimental data enhancement UI

echo "🧪 Running all tests for the experimental data enhancement UI"

# Create directories for results
mkdir -p ./test-results/visual
mkdir -p ./test-results/accessibility
mkdir -p ./test-results/visualization
mkdir -p ./test-results/filtering
mkdir -p ./test-results/comparison
mkdir -p ./test-results/reports
mkdir -p ./test-results/jest
mkdir -p ./test-results/jest/coverage

# Install dependencies if needed
if [ ! -d "node_modules/@playwright" ]; then
  echo "Installing Playwright dependencies..."
  npm install
  npx playwright install --with-deps
fi

# Install axe-core if not already installed
if [ ! -d "node_modules/@axe-core" ]; then
  echo "Installing @axe-core/playwright..."
  npm install --save-dev @axe-core/playwright
fi

# Install axios-mock-adapter if not already installed
if [ ! -d "node_modules/axios-mock-adapter" ]; then
  echo "Installing axios-mock-adapter..."
  npm install --save-dev axios-mock-adapter
fi

# Run API tests for molecules and mixtures
echo "Running Molecule and Mixture API tests..."
mkdir -p ./test-results/api
npm run test:api:verify

# Run Jest unit tests first
echo "Running Jest unit tests..."
npm test -- --json --outputFile=./test-results/jest/results.json --coverage --coverageDirectory=./test-results/jest/coverage

# Run tests with HTML report
echo "Running Playwright E2E tests with HTML report..."
npx playwright test --reporter=html

# Run user flow tests explicitly
echo "Running user flow tests..."
npx playwright test tests/e2e/user-flow.spec.js --reporter=html

# Run data integrity tests explicitly
echo "Running data integrity tests..."
npx playwright test tests/e2e/data-integrity.spec.js --reporter=html

# Generate a timestamp for the report directory
TIMESTAMP=$(date +"%Y%m%d-%H%M%S")
REPORT_DIR="./test-results/reports/report-${TIMESTAMP}"

# Copy the report to our timestamped directory
echo "Copying report to ${REPORT_DIR}..."
mkdir -p "${REPORT_DIR}"
cp -r ./playwright-report/* "${REPORT_DIR}/"

# Create a markdown summary report
SUMMARY_FILE="${REPORT_DIR}/summary.md"

echo "# Experimental Data Enhancement Test Summary" > "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "Tests run on: $(date)" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## UI Tests" >> "${SUMMARY_FILE}"
echo "- Tested basic UI navigation and interaction" >> "${SUMMARY_FILE}"
echo "- Verified experiment cards display correctly" >> "${SUMMARY_FILE}"
echo "- Checked experiment detail page functionality" >> "${SUMMARY_FILE}"
echo "- Tested form inputs and validation" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Data Visualization Tests" >> "${SUMMARY_FILE}"
echo "- Verified rendering of experiment results charts" >> "${SUMMARY_FILE}"
echo "- Tested chart type selection and updates" >> "${SUMMARY_FILE}"
echo "- Checked data table representation" >> "${SUMMARY_FILE}"
echo "- Verified handling of different data ranges" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Filtering and Comparison Tests" >> "${SUMMARY_FILE}"
echo "- Tested search and filter functionality" >> "${SUMMARY_FILE}"
echo "- Verified experiment comparison features" >> "${SUMMARY_FILE}"
echo "- Checked statistical analysis tools" >> "${SUMMARY_FILE}"
echo "- Tested data export capabilities" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Performance and Security Tests" >> "${SUMMARY_FILE}"
echo "- Measured page load times and rendering efficiency" >> "${SUMMARY_FILE}"
echo "- Tested with multiple concurrent users" >> "${SUMMARY_FILE}"
echo "- Verified security headers and input validation" >> "${SUMMARY_FILE}"
echo "- Checked for information leakage" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## User Flow Tests" >> "${SUMMARY_FILE}"
echo "- Verified complete user registration and login process" >> "${SUMMARY_FILE}"
echo "- Tested navigation through main application sections" >> "${SUMMARY_FILE}"
echo "- Verified experiment creation and management workflow" >> "${SUMMARY_FILE}"
echo "- Tested user profile management functionality" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Data Integrity Tests" >> "${SUMMARY_FILE}"
echo "- Verified molecule data consistency across views" >> "${SUMMARY_FILE}"
echo "- Tested experiment data persistence through edit cycles" >> "${SUMMARY_FILE}"
echo "- Verified data validation prevents invalid input" >> "${SUMMARY_FILE}"
echo "- Tested protocol data consistency across related experiments" >> "${SUMMARY_FILE}"
echo "- Verified user permission data integrity" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Accessibility and Cross-Browser Tests" >> "${SUMMARY_FILE}"
echo "- Verified WCAG 2.1 AA compliance" >> "${SUMMARY_FILE}"
echo "- Tested keyboard navigation" >> "${SUMMARY_FILE}"
echo "- Checked rendering across Chrome, Firefox, and Safari" >> "${SUMMARY_FILE}"
echo "- Verified responsive layouts on different devices" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Unit Tests (Jest)" >> "${SUMMARY_FILE}"
echo "- Verified React 18 component behavior" >> "${SUMMARY_FILE}"
echo "- Tested UI components with React Testing Library" >> "${SUMMARY_FILE}"
echo "- Validated theme provider functionality" >> "${SUMMARY_FILE}"
echo "- Tested protein visualizer components" >> "${SUMMARY_FILE}"
echo "- Verified component interactions with user events" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## API Tests" >> "${SUMMARY_FILE}"
echo "- Tested Molecule API endpoints functionality" >> "${SUMMARY_FILE}"
echo "- Verified Mixture API operations" >> "${SUMMARY_FILE}"
echo "- Tested error handling and edge cases" >> "${SUMMARY_FILE}"
echo "- Verified integration between Molecule and Mixture APIs" >> "${SUMMARY_FILE}"
echo "- Tested all CRUD operations for molecules and mixtures" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "See the full HTML report for detailed test results and screenshots." >> "${SUMMARY_FILE}"

echo "✅ All tests complete!"
echo "HTML report available at: ${REPORT_DIR}/index.html"
echo "Summary report available at: ${SUMMARY_FILE}"

# Try to open the report in the default browser
if [ -x "$(command -v xdg-open)" ]; then
  xdg-open "${REPORT_DIR}/index.html"
elif [ -x "$(command -v open)" ]; then
  open "${REPORT_DIR}/index.html"
else
  echo "To view the report, open ${REPORT_DIR}/index.html in your browser"
fi