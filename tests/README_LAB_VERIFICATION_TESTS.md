# Lab Verification Workflow Tests

This document provides an overview of the test suite for the lab verification workflow in the CryoProtect v2 application.

## Overview

The lab verification workflow allows laboratory technicians to verify experimental results by recording verification data, including:
- Verification status (Verified, Pending, Rejected)
- Verifier information
- Equipment used
- Additional comments

The test suite covers all aspects of this workflow:
1. Backend model methods
2. API endpoints
3. Frontend UI components
4. End-to-end integration tests

## Test Files

### Backend Tests

- **Model Tests**: `tests/unit/api/test_lab_verification_model.py`
  - Tests the `LabVerification` model methods
  - Covers: record_verification, get_verification, update_verification_status
  - Tests validation, error handling, and successful operations

- **API Endpoint Tests**: `tests/unit/api/test_lab_verification_resources.py`
  - Tests the REST API endpoints for lab verification
  - Covers: GET, POST, and PUT operations
  - Tests validation, error handling, and successful responses

### Frontend Tests

- **UI Tests**: `tests/unit/frontend/test_lab_verification_ui.js`
  - Tests the JavaScript functionality in `static/js/lab-verification.js`
  - Covers: loading verification data, rendering UI components, form validation, API interactions
  - Uses JSDOM to simulate the browser environment

### Integration Tests

- **Workflow Tests**: `tests/integration/test_lab_verification_workflow.py`
  - Tests the complete lab verification workflow
  - Covers: creating, retrieving, and updating verifications
  - Tests the interaction between model methods and API endpoints

## Running the Tests

### Backend Tests

Run the Python tests using pytest:

```bash
# Run all lab verification tests
pytest tests/unit/api/test_lab_verification_model.py tests/unit/api/test_lab_verification_resources.py tests/integration/test_lab_verification_workflow.py

# Run specific test file
pytest tests/unit/api/test_lab_verification_model.py
```

### Frontend Tests

Run the JavaScript tests using Mocha:

```bash
# Install dependencies if not already installed
npm install --save-dev mocha chai jsdom sinon

# Run the tests
npx mocha tests/unit/frontend/test_lab_verification_ui.js
```

## Test Coverage

The test suite covers the following scenarios:

### Success Cases
- Creating a new verification record
- Retrieving an existing verification
- Updating a verification's status
- Rendering verification data in the UI
- Form validation and submission

### Error Cases
- Invalid verification status
- Missing required fields
- Non-existent experiment or verification
- API errors
- Form validation errors

### Edge Cases
- Different verification statuses (Verified, Pending, Rejected)
- Empty comments field
- Status updates with and without comments

## Mocking Strategy

The tests use mocking to isolate components and avoid external dependencies:

- **Supabase**: The database interactions are mocked using the `MockSupabaseBaseTestCase`
- **API Requests**: Frontend API calls are mocked using Sinon.js
- **DOM**: The browser DOM is simulated using JSDOM
- **Bootstrap**: Modal functionality is mocked for UI tests

## Extending the Tests

When adding new features to the lab verification workflow, extend the tests by:

1. Adding new test methods to the existing test classes
2. Creating new test files for entirely new components
3. Updating the integration tests to cover new workflow scenarios

Follow the existing patterns for consistency and maintain high test coverage for all new functionality.