# Protein Visualizer Component Testing Guide

## Testing Overview

The Protein Visualizer component includes both unit tests and end-to-end (E2E) tests to verify its functionality. Due to the complex nature of 3D visualization and browser compatibility, we've provided multiple testing approaches.

## Testing Approaches

### 1. Unit Tests (Recommended)

Unit tests use Jest and React Testing Library to verify the component's structure and behavior without requiring a browser or server. This approach is fast, reliable, and doesn't depend on external factors.

### 2. End-to-End Tests (Advanced)

E2E tests use Playwright to test the component in a real browser. These tests are more comprehensive but require additional setup and may be affected by environment-specific issues.

#### WebKit Compatibility Issues

Some E2E tests may fail in WebKit due to missing dependencies (`libicudata.so.66`). This is a known issue with WebKit on some Linux distributions. The tests have been configured to run only on Chromium to avoid these issues.

#### Test Server Configuration

The E2E tests require a running Next.js development server. We've provided two approaches:

1. **Automated Server Management**: The `run-protein-visualizer-tests-full.sh` script will start a development server, run the tests, and clean up afterward.

2. **Manual Server Management**: If you prefer to manage the server yourself, you can start it manually with `npm run dev` and then use the `test:protein-visualizer:manual` script.

## Running the Tests

### Unit Tests (Recommended)

```bash
npm run test:protein-visualizer
```

This will:
1. Run the Jest unit tests for the MolstarViewer component
2. Verify the component's structure and behavior in a controlled environment
3. Provide fast and reliable test results

### End-to-End Tests (Advanced)

#### Fully Automated E2E Test

```bash
npm run test:protein-visualizer:e2e
```

This will:
1. Check if a server is already running on port 3000
2. If not, start a Next.js development server
3. Wait for the server to fully initialize
4. Run the tests using Chromium browser in headed mode
5. Generate screenshots in the `test-results/protein-visualizer` directory
6. Stop the server (only if it was started by the script)

#### Manual Server Management

If you prefer to have more control over the server, you can run:

```bash
# First, start the server in a separate terminal
npm run dev

# Then, in another terminal, run the tests
npm run test:protein-visualizer:manual
```

### Additional Options for E2E Tests

You can pass additional options to the E2E test runner:

```bash
# Run with debugger
npm run test:protein-visualizer:e2e -- --debug

# Run headless
npm run test:protein-visualizer:e2e -- --headless

# Run with UI mode
npm run test:protein-visualizer:e2e -- --ui
```

## Test Details

### Unit Tests

The unit test suite includes the following checks:

1. **Component Rendering**: Verifies that the component renders without crashing.

2. **Loading State**: Verifies that the loading indicator is displayed initially.

3. **Error Handling**: Verifies that error states are properly handled and displayed.

4. **Props Handling**: Verifies that the component correctly applies props like dimensions, styles, and behaviors.

5. **UI Elements**: Verifies that control elements are conditionally rendered based on props.

### E2E Tests

The E2E test suite includes the following checks (if you choose to run them):

1. **Server Connectivity**: Verifies that the Next.js server is running and responding to requests.

2. **Demo Page**: Verifies that the protein visualizer demo page loads correctly.

3. **Component Features**: Tests for specific component features like style changes, rotations, etc.

## Manual Testing

The Protein Visualizer component can also be manually tested by:

1. Starting the development server: `npm run dev`
2. Navigating to `/protein-visualizer-demo` in your browser
3. Interacting with the various visualization options and controls

## Troubleshooting

### Unit Tests

If you encounter issues with unit tests:

1. **Missing Dependencies**: Ensure all Node.js dependencies are installed with `npm install`.

2. **Jest Configuration**: Verify that the Jest configuration in package.json is correct.

3. **Module Resolution**: If you encounter module resolution errors, check the moduleNameMapper in the Jest config.

### E2E Tests

If you encounter issues with E2E tests:

1. **Server Not Starting**: Check if port 3000 is already in use by another process.

2. **Missing Dependencies**: Some browser engines might require additional system dependencies.

3. **Browser Issues**: If tests fail due to browser issues, try running with just Chromium:
   ```bash
   npm run test:protein-visualizer:e2e -- --project=chromium
   ```

4. **Screenshot Location**: All test screenshots are saved in `test-results/protein-visualizer` and can help diagnose rendering issues.

5. **Debugging**: For detailed debugging, run the tests with the `--debug` flag:
   ```bash
   npm run test:protein-visualizer:e2e -- --debug
   ```