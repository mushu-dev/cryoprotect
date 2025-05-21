#!/bin/bash

# Create test results directory if it doesn't exist
mkdir -p ./test-results/protein-visualizer

# Run protein visualizer tests with custom configuration
npx playwright test tests/protein-visualizer/visualizer.spec.js -c tests/protein-visualizer/playwright.config.js "$@"