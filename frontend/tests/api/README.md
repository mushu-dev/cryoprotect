# CryoProtect API Tests

This directory contains tests for the CryoProtect API, focusing on the Molecule and Mixture endpoints.

## Overview

We have three different approaches to testing the API:

1. **Mock API Tests** - `molecule-mixture-api.test.js`
   - Uses Jest and axios-mock-adapter
   - Tests API client code against mock responses
   - No actual API calls are made
   - Run with `npm run test:api:molecule-mixture`

2. **API Verification** - `verify-molecule-mixture-api.js`
   - Standalone script that verifies API structure
   - Uses mock data to simulate responses
   - Performs validation on expected data structure
   - Run with `npm run test:api:verify`

3. **Live API Tests** - `live-api-test.js`
   - Integration tests that connect to real backend services
   - Makes actual API calls to Heroku and Fly.io services
   - Validates real responses against expected formats
   - Run with `npm run test:api:live`

## Testing Strategy

Our API testing strategy follows a layered approach:

- **Unit Tests** - Test individual API client functions in isolation
- **Integration Tests** - Test integration between frontend and backend
- **Live Tests** - Test against real backend services

This approach allows us to identify issues at different levels of the application.

## Running the Tests

### Mock API Tests

```bash
npm run test:api:molecule-mixture
```

These tests use Jest and axios-mock-adapter to test the API client code without making actual API calls.

### API Verification

```bash
npm run test:api:verify
```

This script verifies the expected structure of API responses using mock data.

### Live API Tests

```bash
npm run test:api:live
```

These tests make actual API calls to the backend services and validate the responses.

#### Using Custom API URLs

You can specify custom API URLs for the live tests:

```bash
CUSTOM_API_URL=http://localhost:5000/api/v1 CUSTOM_RDKIT_URL=http://localhost:5001 npm run test:api:live
```

## Test Files

- `molecule-mixture-api.test.js` - Jest tests for API client
- `basic-api-test.js` - Simplified Jest tests for troubleshooting
- `verify-molecule-mixture-api.js` - API structure verification
- `live-api-test.js` - Live API integration tests

## Backend Services

The tests interact with the following backend services:

- **Main API** - Heroku: `https://cryoprotect-8030e4025428.herokuapp.com/v1`
- **RDKit Service** - Fly.io: `https://cryoprotect-rdkit.fly.dev`

## Test Results

Test results are saved in the following directories:

- Mock API test results: `test-results/jest/`
- API verification results: `test-results/api/`
- Live API test results: `test-results/api-live/`