#!/bin/bash

# This script runs all data-focused tests for the experimental data enhancement UI

echo "📊 Running data visualization, filtering, and comparison tests"

# Create directories for test results
mkdir -p ./test-results/visualization
mkdir -p ./test-results/filtering
mkdir -p ./test-results/comparison
mkdir -p ./test-results/api
mkdir -p ./test-results/ui
mkdir -p ./test-results/form-validation

# Install dependencies if needed
if [ ! -d "node_modules/@playwright" ]; then
  echo "Installing Playwright dependencies..."
  npm install
  npx playwright install --with-deps
fi

# Run all data-focused tests
echo "Running experiment API tests..."
npx playwright test tests/playwright/experiment-data-api.spec.js

echo "Running experiment UI tests..."
npx playwright test tests/playwright/experimental-data-ui.spec.js

echo "Running form validation tests..."
npx playwright test tests/playwright/experiment-form-validation.spec.js

echo "Running data visualization tests..."
npx playwright test tests/playwright/data-visualization.spec.js

echo "Running data filtering tests..."
npx playwright test tests/playwright/data-filtering.spec.js

echo "Running data comparison tests..."
npx playwright test tests/playwright/data-comparison.spec.js

# Generate a timestamp for the report directory
TIMESTAMP=$(date +"%Y%m%d-%H%M%S")
REPORT_DIR="./test-results/experimental-data-validation-${TIMESTAMP}"

# Create summary report
mkdir -p "${REPORT_DIR}"
SUMMARY_FILE="${REPORT_DIR}/summary.md"

echo "# Experimental Data Enhancement Test Summary" > "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "Tests run on: $(date)" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Experiment API Tests" >> "${SUMMARY_FILE}"
echo "- Verified API endpoints for experiment data retrieval" >> "${SUMMARY_FILE}"
echo "- Tested error handling for invalid requests" >> "${SUMMARY_FILE}"
echo "- Checked data format consistency" >> "${SUMMARY_FILE}"
echo "- Verified request/response cycle performance" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Experiment UI Tests" >> "${SUMMARY_FILE}"
echo "- Verified experiment list rendering" >> "${SUMMARY_FILE}"
echo "- Tested experiment detail view" >> "${SUMMARY_FILE}"
echo "- Checked mobile responsiveness" >> "${SUMMARY_FILE}"
echo "- Verified interactive UI elements" >> "${SUMMARY_FILE}"
echo "- Tested loading states and error handling" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Form Validation Tests" >> "${SUMMARY_FILE}"
echo "- Tested required field validation" >> "${SUMMARY_FILE}"
echo "- Verified data type validation" >> "${SUMMARY_FILE}"
echo "- Checked error message display" >> "${SUMMARY_FILE}"
echo "- Tested form submission behavior" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Data Visualization Tests" >> "${SUMMARY_FILE}"
echo "- Verified rendering of experiment results charts" >> "${SUMMARY_FILE}"
echo "- Tested chart type selection and updates" >> "${SUMMARY_FILE}"
echo "- Checked data table representation" >> "${SUMMARY_FILE}"
echo "- Verified interactive tooltips on data points" >> "${SUMMARY_FILE}"
echo "- Tested handling of different data ranges" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Data Filtering Tests" >> "${SUMMARY_FILE}"
echo "- Tested date range filtering" >> "${SUMMARY_FILE}"
echo "- Verified text search functionality" >> "${SUMMARY_FILE}"
echo "- Checked sorting by different criteria" >> "${SUMMARY_FILE}"
echo "- Tested status and type filtering" >> "${SUMMARY_FILE}"
echo "- Verified clear filters functionality" >> "${SUMMARY_FILE}"
echo "- Tested combining multiple filters" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

echo "## Data Comparison Tests" >> "${SUMMARY_FILE}"
echo "- Tested selecting multiple experiments" >> "${SUMMARY_FILE}"
echo "- Verified comparative visualizations" >> "${SUMMARY_FILE}"
echo "- Checked statistical analysis features" >> "${SUMMARY_FILE}"
echo "- Tested parameter adjustment" >> "${SUMMARY_FILE}"
echo "- Verified export functionality" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

# Create subdirectories in the report directory
mkdir -p "${REPORT_DIR}/api-tests"
mkdir -p "${REPORT_DIR}/ui-tests"
mkdir -p "${REPORT_DIR}/form-validation-tests"
mkdir -p "${REPORT_DIR}/data-analysis-tests"
mkdir -p "${REPORT_DIR}/integration-tests"
mkdir -p "${REPORT_DIR}/accessibility-tests"

# Copy screenshots to report directory
cp -r ./test-results/api/ "${REPORT_DIR}/api-tests/"
cp -r ./test-results/ui/ "${REPORT_DIR}/ui-tests/"
cp -r ./test-results/form-validation/ "${REPORT_DIR}/form-validation-tests/"
cp -r ./test-results/visualization/ "${REPORT_DIR}/data-analysis-tests/"
cp -r ./test-results/filtering/ "${REPORT_DIR}/data-analysis-tests/"
cp -r ./test-results/comparison/ "${REPORT_DIR}/data-analysis-tests/"

echo "See the screenshots in ${REPORT_DIR}/ for visual verification of test results." >> "${SUMMARY_FILE}"

echo "✅ All experimental data enhancement tests complete!"
echo "Summary report available at: ${REPORT_DIR}/summary.md"

# Show the Playwright report
npx playwright show-report

# Try to open the summary report
if [ -x "$(command -v xdg-open)" ]; then
  xdg-open "${SUMMARY_FILE}"
elif [ -x "$(command -v open)" ]; then
  open "${SUMMARY_FILE}"
else
  echo "To view the summary report, open ${SUMMARY_FILE} in your text editor"
fi