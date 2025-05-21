# Experiment API Integration Testing

This document describes how to run the integration tests for the Experiment API. The tests cover both the backend API and frontend service integration.

## Overview

The test suite includes:

1. Python script (`test_experiment_api_integration.py`) for testing the backend API directly
2. JavaScript script (`frontend/test-experiment-integration.js`) for testing the frontend service integration

These tests validate the full functionality of the experiment API including:
- Creating experiments
- Adding experiment results
- Retrieving experiment details
- Updating experiment data
- Running experiment analysis
- Searching and listing experiments

## Prerequisites

### For backend testing:

- Python 3.7+
- `requests` package

Install dependencies:
```bash
pip install requests
```

### For frontend testing:

- Node.js 14+
- npm or yarn

## Running the Backend Tests

The backend tests communicate directly with the API endpoints. Make sure your API server is running before executing the tests.

```bash
# Make sure the API server is running first
python test_experiment_api_integration.py
```

By default, the script assumes the API is running at `http://localhost:5000/api`. You can modify the `BASE_URL` variable in the script if your API is at a different location.

## Running the Frontend Tests

The frontend tests validate the frontend service layer that communicates with the API. These tests can run in two modes:

1. **Mock mode** (default): Uses mock implementations to simulate API responses
2. **Real API mode**: Makes actual API calls to the backend

```bash
# From the project root directory
cd frontend
node test-experiment-integration.js
```

To switch between mock mode and real API mode, modify the `useMockMode` variable in the test script.

## Test Outcomes

The tests will print detailed information about each step, including:

- ✅ Success messages for passed tests
- ❌ Error messages for failed tests
- Details about created or retrieved experiments

## Troubleshooting

If the tests fail, check the following:

1. **API server running**: Ensure the API server is running and accessible
2. **Authentication**: If your API requires authentication, make sure the test scripts include the proper authentication headers
3. **Database connectivity**: Verify the API server can connect to the database
4. **Network issues**: Check for firewall or network configuration problems
5. **Data validation**: Ensure the test data meets all validation requirements

## Extending the Tests

To add more test cases:

1. Add new test functions to the test scripts
2. Call these functions from the main test execution function
3. Update the test data generators if needed

## Continuous Integration

These tests can be incorporated into a CI/CD pipeline by:

1. Running the backend tests after deploying the API
2. Running the frontend tests after building the frontend application
3. Including both tests in pre-merge checks for relevant code changes

---

For any questions or issues, please contact the development team.