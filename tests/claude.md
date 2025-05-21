# Testing Framework Reference

## Overview

The testing framework for CryoProtect provides comprehensive tools for unit, integration, and end-to-end testing, ensuring code quality and reliability across the application.

## Test Categories

### Unit Tests
- **API Tests**: Testing individual API endpoints and resources
- **Database Tests**: Testing database adapters and operations
- **Utility Tests**: Testing helper functions and utilities
- **Model Tests**: Testing data models and validation

### Integration Tests
- **API Integration**: End-to-end API workflow testing
- **Database Integration**: Testing database migrations and complex queries
- **External Service Integration**: Testing ChEMBL and PubChem client integration
- **Authentication Flow**: Testing JWT and service role authentication

### Performance Tests
- **Query Performance**: Database query execution time testing
- **API Response Time**: Endpoint response time benchmarking
- **Concurrency Tests**: Testing under parallel load
- **Connection Pool Tests**: Testing connection pooling under load

## Test Infrastructure

### Framework Components
- **Test Runner**: Custom test runner in `run_tests.py`
- **Mock Supabase**: Mock implementation in `tests/mock_supabase/`
- **Test Data**: Fixture data in `tests/test_data/`
- **Test Utilities**: Helper functions in `tests/utils/`

### Mock Implementations
- **Mock Database**: In-memory database simulation
- **Mock ChEMBL API**: Cached response simulation
- **Mock PubChem API**: Fixture-based response simulation
- **Mock Authentication**: Test JWT generation and validation

## Key Test Files

- **`test_api_endpoints.py`**: Comprehensive API endpoint testing
- **`test_database_performance.py`**: Database performance benchmarks
- **`test_connection_pool.py`**: Connection pool behavior verification
- **`test_rdkit_integration.py`**: RDKit integration validation
- **`test_end_to_end_workflows.py`**: Complete workflow validation

## Testing Patterns

### Fixture Management
- **Test Data Creation**: Dynamic test data generation
- **Database Seeding**: Reproducible database state initialization
- **Cleanup Operations**: Proper test cleanup to prevent state leakage
- **Context Managers**: Environment setup and teardown

### Mocking Strategies
- **Service Mocking**: External service simulation
- **Database Mocking**: Database adapter substitution
- **Request Mocking**: HTTP request interception
- **Environment Mocking**: Environment variable overrides

## Running Tests

### Basic Test Execution
```bash
# Run all tests
python tests/run_tests.py

# Run specific test file
python tests/run_tests.py -t test_file_name.py

# Run specific test pattern
python -m unittest tests/test_specific_file.py::TestClass::test_method
```

### Test Options
- **Verbosity**: `-v` for verbose output
- **Coverage**: `--coverage` for code coverage report
- **Skip Categories**: `--skip db` to skip database tests
- **Performance Only**: `--performance` for performance tests only
- **Parallel Execution**: `--parallel` to run tests in parallel

## Test Writing Guidelines

1. **Test Isolation**: Tests should not depend on each other
2. **Clean Setup/Teardown**: Initialize and clean test environment properly
3. **Clear Assertions**: Use descriptive assertion messages
4. **Comprehensive Coverage**: Test happy path and error cases
5. **Mock External Dependencies**: Don't rely on external services in unit tests
6. **Performance Awareness**: Include performance assertions where relevant

## Common Pitfalls

1. **Database State Leakage**: Tests affecting each other through database state
2. **External Dependencies**: Relying on external services in unit tests
3. **Hardcoded Expectations**: Tests with hardcoded values that may change
4. **Missing Cleanup**: Not properly cleaning up created resources
5. **Incomplete Mocking**: Not mocking all external dependencies
6. **Brittle Tests**: Tests that fail due to implementation details rather than behavior

## Best Practices

1. **Test-Driven Development**: Write tests before implementation
2. **Continuous Testing**: Run tests automatically on code changes
3. **Separating Concerns**: Unit tests should test one thing only
4. **Integration Test Coverage**: Ensure key workflows have integration tests
5. **Performance Baseline**: Establish performance baselines for critical operations
6. **Test Data Management**: Use generators instead of hardcoded test data