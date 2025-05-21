# Experimental Data Enhancement Testing

This document provides a quick reference for testing the Experimental Data Enhancement feature.

## Quick Start

1. Setup test environment:
```bash
./setup-test-env.sh
```

2. Run all tests:
```bash
./validate-experimental-data-enhancement.sh
```

## Available Scripts

- `./setup-test-env.sh` - Sets up the test environment
- `./run-data-tests.sh` - Runs specialized data feature tests
- `./validate-experimental-data-enhancement.sh` - Runs complete validation suite
- `./run-experimental-tests.sh` - Runs existing experimental data tests

## Test Files

### End-to-End Tests
- `tests/e2e/experimental-data.spec.js` - Complete workflow tests
- `tests/e2e/data-analysis.spec.js` - Analysis workflow tests

### Specialized Feature Tests
- `tests/playwright/experiment-data-api.spec.js` - API tests
- `tests/playwright/experimental-data-ui.spec.js` - UI component tests
- `tests/playwright/experiment-form-validation.spec.js` - Form validation tests
- `tests/playwright/data-visualization.spec.js` - Visualization tests
- `tests/playwright/data-filtering.spec.js` - Filtering tests
- `tests/playwright/data-comparison.spec.js` - Comparison tests

## Test Coverage

### UI Tests
- Experiment list page rendering
- Individual experiment details view
- Mobile responsiveness
- Interactive UI elements
- Loading states and error handling

### API Tests
- API endpoints for experiment data retrieval
- Error handling for invalid requests
- Data format consistency
- Request/response cycle performance

### Form Validation Tests
- Required field validation
- Data type validation
- Error message display
- Form submission behavior

### Data Feature Tests
- Chart rendering and visualization types
- Filtering and sorting capabilities
- Data comparison tools
- Statistical analysis

## Documentation

For detailed information on the test strategy and implementation, see:
- [Experimental Data Validation Guide](/docs/EXPERIMENTAL_DATA_VALIDATION.md)

## Troubleshooting

If tests fail, try:

1. Reinstalling browser dependencies:
```bash
npx playwright install --with-deps
```

2. Running a specific test with debugging:
```bash
npx playwright test tests/playwright/experimental-data-ui.spec.js --debug
```

3. Checking the test report:
```bash
npx playwright show-report
```

4. Reviewing test results in the test-results directory:
```bash
ls -la ./test-results/experimental-data-validation-*
```