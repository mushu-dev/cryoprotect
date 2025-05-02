# CryoProtect v2 - Database Connection Mocking System

This document provides an overview of the database connection mocking system for offline testing of the CryoProtect v2 application.

## Overview

The database connection mocking system allows you to test the CryoProtect v2 application without requiring an actual Supabase connection. It provides realistic mock responses for common database queries, supports CRUD operations on all relevant tables, and allows for easy configuration of test scenarios.

## Features

- **Complete Supabase Client Mocking**: Simulates the behavior of the real Supabase client
- **Realistic Test Data**: Includes predefined test data based on the application's schema
- **CRUD Operation Support**: Supports all common database operations (select, insert, update, delete)
- **RPC Function Mocking**: Simulates remote procedure calls used in the application
- **Authentication Mocking**: Simulates Supabase authentication
- **Configurable Test Scenarios**: Supports success, failure, and timeout scenarios
- **Flexible Integration**: Multiple approaches to integrate with your tests

## Directory Structure

```
tests/mock_supabase/
├── __init__.py         # Package initialization and exports
├── client.py           # Mock Supabase client implementation
├── data.py             # Mock data storage and management
├── query.py            # Mock query builders and response handling
├── rpc.py              # Mock RPC implementation
├── auth.py             # Mock authentication implementation
├── helpers.py          # Helper functions for using the mocks in tests
├── patch.py            # Functions to patch the Supabase client
└── README.md           # Detailed documentation
```

## Usage Examples

### Using the Decorator Approach

```python
from tests.mock_supabase.helpers import patch_supabase

class TestYourFeature(unittest.TestCase):
    
    @patch_supabase(load_data=True)
    def test_your_function(self, mock_client):
        # Your test code here
        # The mock_client is automatically passed to your test
        
        # Example: Test getting all molecules
        molecules = Molecule.get_all()
        self.assertGreater(len(molecules), 0)
```

### Using the Base Test Case Approach

```python
from tests.mock_supabase.helpers import MockSupabaseTestCase

class TestYourFeature(MockSupabaseTestCase):
    
    def test_your_function(self):
        # Your test code here
        # self.mock_client is available in your test
        
        # Example: Test getting all molecules
        molecules = Molecule.get_all()
        self.assertGreater(len(molecules), 0)
```

### Configuring Test Scenarios

```python
from tests.mock_supabase.helpers import configure_test_scenario

# Configure success scenario (default)
configure_test_scenario('success')

# Configure failure scenario
configure_test_scenario('failure')

# Configure timeout scenario
configure_test_scenario('timeout')
```

## Test Data

The mocking system includes predefined test data for:

- **Molecules**: DMSO, Glycerol, Ethylene glycol, Propylene glycol, Trehalose, Water
- **Molecular Properties**: Molecular Weight, LogP, H-Bond Donors, H-Bond Acceptors
- **Mixtures**: DMSO/EG Mixture
- **Mixture Components**: Components for the DMSO/EG Mixture

This test data is loaded automatically when using the mocking system, providing a realistic dataset for your tests.

## Integration with Existing Tests

The mocking system is designed to integrate seamlessly with the existing test framework. You can use it in your existing tests by:

1. Using the `patch_supabase` decorator
2. Extending the `MockSupabaseTestCase` class
3. Using the `SupabasePatchMixin` mixin
4. Directly patching the Supabase client using the functions in `patch.py`

## Example Test File

A complete example test file is provided in `tests/test_mock_supabase_example.py`. This file demonstrates how to use the mocking system in different scenarios and can be used as a reference for your own tests.

## Detailed Documentation

For more detailed information about the mocking system, including advanced usage and customization options, see the `README.md` file in the `tests/mock_supabase` directory.