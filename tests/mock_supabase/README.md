# CryoProtect v2 - Supabase Mocking System

This directory contains a comprehensive mocking system for the Supabase database connection used in the CryoProtect v2 application. It allows for offline testing without requiring an actual Supabase connection.

## Overview

The mocking system provides:

1. A mock Supabase client that mimics the behavior of the real Supabase client
2. Mock implementations of common database operations (select, insert, update, delete)
3. Mock implementations of Supabase-specific features (RPC, auth)
4. Predefined test data based on the application's schema
5. Helper functions for using the mocks in tests
6. Support for configuring test scenarios (success, failure, timeout)

## Directory Structure

- `__init__.py` - Package initialization and exports
- `client.py` - Mock Supabase client implementation
- `data.py` - Mock data storage and management
- `query.py` - Mock query builders and response handling
- `rpc.py` - Mock RPC (Remote Procedure Call) implementation
- `auth.py` - Mock authentication implementation
- `helpers.py` - Helper functions for using the mocks in tests
- `README.md` - This documentation file

## Usage

There are two main approaches to using the mocking system in your tests:

### 1. Using the Decorator Approach

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

### 2. Using the Base Test Case Approach

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

## Test Data

The mocking system comes with predefined test data that matches the schema used in the CryoProtect v2 application. The test data includes:

- Molecules (DMSO, Glycerol, Ethylene glycol, etc.)
- Molecular properties (Molecular Weight, LogP, H-Bond Donors, etc.)
- Mixtures (DMSO/EG Mixture)
- Mixture components

You can load this test data using the `load_test_data()` function or by setting `load_data=True` in the `patch_supabase` decorator.

## Customizing Mock Behavior

### Custom RPC Functions

You can register custom RPC functions to mock specific behavior:

```python
from tests.mock_supabase.helpers import mock_rpc_function

# Register a custom RPC function
mock_rpc_function('your_function_name', result={'success': True, 'data': 'test'})

# Or with an error
mock_rpc_function('your_function_name', error='Something went wrong')
```

### Test Scenarios

You can configure different test scenarios to test how your code handles different situations:

```python
from tests.mock_supabase.helpers import configure_test_scenario

# Configure success scenario (default)
configure_test_scenario('success')

# Configure failure scenario
configure_test_scenario('failure')

# Configure timeout scenario
configure_test_scenario('timeout')
```

## Example

Here's a complete example of how to use the mocking system:

```python
import unittest
from tests.mock_supabase.helpers import MockSupabaseTestCase

class TestMixtures(MockSupabaseTestCase):
    
    def test_create_mixture(self):
        # Create a new mixture
        mixture_data = {
            'name': 'Test Mixture',
            'description': 'A test mixture',
            'created_by': '00000000-0000-0000-0000-000000000001'
        }
        
        mixture = Mixture.create(mixture_data)
        
        # Verify the mixture was created
        self.assertIsNotNone(mixture)
        self.assertEqual(mixture['name'], 'Test Mixture')
        
        # Get the mixture
        retrieved = Mixture.get(mixture['id'])
        self.assertEqual(retrieved['name'], 'Test Mixture')
```

For more examples, see the `test_mock_supabase_example.py` file.

## Extending the Mocking System

If you need to extend the mocking system to support additional features:

1. Add new mock data structures in `data.py`
2. Add new mock implementations in the appropriate module
3. Update the `__init__.py` file to export any new functions or classes

## Limitations

The mocking system has some limitations compared to a real Supabase connection:

1. It does not support all Supabase features (only the ones used in the CryoProtect v2 application)
2. It does not validate SQL syntax in RPC calls
3. It does not enforce foreign key constraints
4. It does not support full-text search

These limitations should not affect most tests, but if you need to test these specific features, you may need to use a real Supabase connection.