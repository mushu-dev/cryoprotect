# CryoProtect v2 API Testing with Supabase Integration

This document provides instructions for testing the CryoProtect v2 API endpoints that interact with Supabase.

## Overview

The API testing framework focuses on verifying that the API endpoints correctly interact with the Supabase database, with a particular emphasis on:

1. **Authentication scenarios**:
   - Unauthenticated access
   - Authenticated user accessing own data
   - Authenticated user attempting to access others' data
   - Admin user access
   - Service role access

2. **RLS (Row Level Security) policy enforcement**:
   - Verifying that RLS policies correctly restrict access based on user roles
   - Ensuring that admin users can access any data
   - Confirming that service roles can bypass RLS policies

3. **CRUD operations for key data types**:
   - Molecules
   - Mixtures
   - Predictions
   - Experiments

4. **Property calculation endpoints**:
   - Testing RDKit integration for molecular property calculations

5. **Search endpoints**:
   - Testing search functionality for molecules and mixtures

## Prerequisites

Before running the tests, ensure that you have the following:

1. Python 3.9 or higher installed
2. The CryoProtect v2 application installed and configured
3. A Supabase project set up with the correct schema (project ID: tsdlmynydfuypiugmkev)
4. The required Python packages installed:
   ```
   pip install -r requirements.txt
   ```

## Running the Tests

### On Windows

1. Open a command prompt
2. Navigate to the CryoProtect v2 directory
3. Run the test script:
   ```
   run_supabase_api_tests.bat
   ```

### On Linux/macOS

1. Open a terminal
2. Navigate to the CryoProtect v2 directory
3. Make the test script executable:
   ```
   chmod +x run_supabase_api_tests.sh
   ```
4. Run the script:
   ```
   ./run_supabase_api_tests.sh
   ```

## Test Results

After running the tests, a comprehensive test report will be generated in the `reports` directory:

- `SUPABASE_API_TEST_REPORT_<timestamp>.md`: A Markdown report with detailed test results
- `supabase_api_test_results_<timestamp>.json`: A JSON file with raw test results

The report includes:
- A summary of test results (passed, failed, errors)
- Details of any failures or errors
- Authentication test results
- RLS policy enforcement evaluation
- Recommendations for fixing any issues

## Test Structure

The tests are organized into the following categories:

### Unauthenticated Access Tests

These tests verify that:
- Public endpoints are accessible without authentication
- Protected endpoints require authentication

### Authenticated User Tests - Own Data

These tests verify that authenticated users can:
- Create their own data
- Read their own data
- Update their own data
- Delete their own data

### Authenticated User Tests - Other Users' Data

These tests verify that authenticated users cannot:
- Read other users' private data
- Update other users' data
- Delete other users' data

### Admin User Tests

These tests verify that admin users can:
- Access any data in the system
- Update any data in the system
- Delete any data in the system

### Service Role Tests

These tests verify that service roles can:
- Bypass RLS policies
- Access any data in the system

### Property Calculation Tests

These tests verify that:
- Molecular property calculation endpoints work correctly
- RDKit integration functions properly

### Search Tests

These tests verify that:
- Search endpoints return correct results
- Search parameters are properly validated

### Error Handling Tests

These tests verify that:
- 404 errors are properly handled
- Validation errors are properly handled
- Server errors are properly handled

## Customizing the Tests

### Adding New Tests

To add new tests, you can:

1. Add new test methods to the `TestSupabaseAPIIntegration` class in `tests/test_supabase_api_integration.py`
2. Follow the existing pattern of:
   - Setting up test data
   - Mocking Supabase responses
   - Making API requests
   - Asserting expected results

### Modifying Existing Tests

To modify existing tests:

1. Open `tests/test_supabase_api_integration.py`
2. Locate the test method you want to modify
3. Update the test as needed

## Troubleshooting

### Common Issues

1. **Authentication Failures**:
   - Check that the Supabase URL and API keys are correctly configured in the `.env` file
   - Verify that the authentication headers are correctly formatted

2. **RLS Policy Issues**:
   - Check that the RLS policies are correctly defined in Supabase
   - Verify that the user roles are correctly assigned

3. **Mock Supabase Issues**:
   - Check that the mock Supabase client is correctly configured
   - Verify that the mock responses match the expected format

### Debugging

To enable more detailed debugging:

1. Add the following to the top of `tests/run_supabase_api_tests.py`:
   ```python
   import logging
   logging.basicConfig(level=logging.DEBUG)
   ```

2. Run the tests again to see more detailed logs

## Conclusion

The CryoProtect v2 API testing framework provides a comprehensive set of tests to validate the API endpoints that interact with Supabase. By running these tests regularly, you can ensure that the API continues to function correctly as changes are made to the system.