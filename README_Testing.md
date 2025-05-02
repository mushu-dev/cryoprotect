We need to create the openapi.yaml file that's mentioned in docs/api/README.md
  but doesn't yet exist.# CryoProtect v2 Testing Framework

This document provides instructions for using the CryoProtect v2 testing framework to validate the system.

## Overview

The testing framework consists of several test scripts that validate different components of the system:

1. **Database Schema and Data Integrity Tests** (`test_database_schema.py`): Verifies that the database schema matches the expected structure and that data integrity constraints are enforced.

2. **RLS Policy Tests** (`test_rls_policies.py`): Verifies that Row Level Security (RLS) policies are correctly implemented and enforced.

3. **API Integration Tests** (`test_api_integration.py`): Verifies that all API endpoints are functioning correctly.

4. **RDKit Integration Tests** (`test_rdkit_integration.py`): Verifies that RDKit is correctly integrated and functioning properly.

5. **End-to-End Workflow Tests** (`test_end_to_end_workflows.py`): Verifies that complete user workflows function correctly.

The test runner script (`run_all_tests.py`) runs all these tests and generates a comprehensive test report.

## Prerequisites

Before running the tests, ensure that you have the following:

1. Python 3.9 or higher installed
2. The CryoProtect v2 application installed and configured
3. A Supabase database set up with the correct schema
4. RDKit installed and configured
5. The API server running

## Running the Tests

### On Windows

1. Open a command prompt
2. Navigate to the CryoProtect v2 directory
3. Run the `run_tests.bat` script:

```
run_tests.bat
```

### On Linux/macOS

1. Open a terminal
2. Navigate to the CryoProtect v2 directory
3. Make the `run_tests.sh` script executable:

```
chmod +x run_tests.sh
```

4. Run the script:

```
./run_tests.sh
```

## Test Results

After running the tests, a comprehensive test report will be generated in the file `TEST_RESULTS_REPORT_<timestamp>.md`. This report includes:

- A summary of all test results
- Detailed results for each component
- Issues found during testing
- Recommendations for addressing any issues

Additionally, each test script generates a JSON file with detailed test results:

- `database_schema_test_results.json`
- `rls_policy_test_results.json`
- `api_integration_test_results.json`
- `rdkit_integration_test_results.json`
- `end_to_end_test_results.json`

## Test Coverage Reporting

### Generating Coverage Reports Locally

To measure and visualize code coverage for the entire CryoProtect v2 codebase:

1. Ensure your environment is activated and all dependencies are installed (pytest and pytest-cov are required).
2. Run the following command from the project root:

   ```bash
   python run_coverage.py
   ```

3. This will:
   - Execute all tests with coverage measurement
   - Generate a unified HTML coverage report in `reports/htmlcov/index.html`
   - Print a coverage summary by module to the terminal

4. Open `reports/htmlcov/index.html` in your browser to view detailed coverage metrics by module and file. The report highlights which lines are covered, missing, or excluded.

If any tests fail, the script will exit with a nonzero status code.

### Coverage Target and Interpretation

The project has a **minimum required coverage target of 70%**. This target applies to:
- Overall project coverage
- Individual modules where possible

When reviewing coverage reports:
1. Focus first on modules with less than 70% coverage
2. Prioritize critical modules and core functionality
3. Pay special attention to complex code paths and error handling
4. Look for untested public methods and functions

Low coverage often indicates:
- Insufficient test cases
- Code that's difficult to test (may need refactoring)
- Dead code that could be removed

### CI/CD Integration

The project's CI/CD workflow automatically runs coverage analysis on every push:

1. The workflow executes `run_coverage.py` to generate coverage reports
2. The HTML report is uploaded as an artifact that can be downloaded from the GitHub Actions page
   - Navigate to the specific workflow run
   - Click on "Artifacts"
   - Download the "coverage-report" artifact
   - Extract and open `index.html` in your browser

3. Coverage data is also uploaded to Codecov for tracking and visualization:
   - Visit the project's Codecov page to view historical coverage trends
   - Coverage badges are available for inclusion in documentation
   - Codecov provides PR comments showing coverage impact of changes

This integration ensures that coverage is continuously monitored and that coverage regressions can be quickly identified and addressed.

## Running Individual Tests

You can also run individual test scripts if you want to focus on a specific component:

```
python test_database_schema.py
python test_rls_policies.py
python test_api_integration.py
python test_rdkit_integration.py
python test_end_to_end_workflows.py
```

Each script will generate its own log file and JSON results file.

## Customizing the Tests

### Database Schema Tests

The database schema tests are configured in the `DatabaseSchemaTest` class in `test_database_schema.py`. You can modify the `expected_tables` and `expected_foreign_keys` lists to match your specific database schema.

### RLS Policy Tests

The RLS policy tests are configured in the `RLSPolicyTest` class in `test_rls_policies.py`. You can modify the test methods to match your specific RLS policies.

### API Integration Tests

The API integration tests are configured in the `APIIntegrationTest` class in `test_api_integration.py`. You can modify the test methods to match your specific API endpoints.

### RDKit Integration Tests

The RDKit integration tests are configured in the `RDKitIntegrationTest` class in `test_rdkit_integration.py`. You can modify the test methods to match your specific RDKit integration.

### End-to-End Workflow Tests

The end-to-end workflow tests are configured in the `EndToEndWorkflowTest` class in `test_end_to_end_workflows.py`. You can modify the test methods to match your specific workflows.

## Troubleshooting

### Database Connection Issues

If you encounter database connection issues, check the following:

1. Ensure that the Supabase URL and API keys are correctly configured in the `config.py` file
2. Verify that the Supabase database is running and accessible
3. Check the network connectivity between the test machine and the Supabase database

### API Connection Issues

If you encounter API connection issues, check the following:

1. Ensure that the API server is running and accessible
2. Verify that the API base URL is correctly configured in the test scripts
3. Check the network connectivity between the test machine and the API server

### RDKit Issues

If you encounter RDKit issues, check the following:

1. Ensure that RDKit is correctly installed and configured
2. Verify that the Python environment has access to the RDKit libraries
3. Check the RDKit version compatibility with the CryoProtect v2 application

## Extending the Tests

To add new tests, you can:

1. Add new test methods to the existing test classes
2. Create new test scripts for additional components
3. Update the `run_all_tests.py` script to include your new test scripts

When adding new tests, follow the existing pattern of:

1. Setting up any necessary test data
2. Performing the test operations
3. Verifying the results
4. Cleaning up any test data
5. Adding the test results to the test report

## Conclusion

The CryoProtect v2 testing framework provides a comprehensive set of tests to validate the system. By running these tests regularly, you can ensure that the system continues to function correctly as changes are made.