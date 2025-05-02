# Test Data Fixtures for CryoProtect v2

This directory contains test data fixtures for the CryoProtect v2 application. These fixtures provide standardized test data and data generation utilities for all entity types in the system.

## Overview

The Test Data Fixtures component is designed to:

1. Provide standardized test data sets for all entity types
2. Offer flexible data generation utilities for creating custom test data
3. Support both simple and complex test scenarios
4. Ensure consistent data relationships between entities
5. Facilitate integration with other test components (Database Fixtures, Mock Objects, API Fixtures)

## Directory Structure

```
tests/fixtures/data/
├── __init__.py              # Package exports
├── generators.py            # Data generation utilities
├── loaders.py               # Data loading utilities
├── README.md                # This documentation file
├── molecules.json           # Standard molecule test data
├── property_types.json      # Standard property type test data
├── molecular_properties.json # Standard molecular property test data
├── mixtures.json            # Standard mixture test data
├── mixture_components.json  # Standard mixture component test data
├── calculation_methods.json # Standard calculation method test data
├── predictions.json         # Standard prediction test data
├── experiments.json         # Standard experiment test data
├── experiment_properties.json # Standard experiment property test data
├── teams.json               # Standard team test data
├── projects.json            # Standard project test data
└── user_profiles.json       # Standard user profile test data
```

## Available Fixtures

### Standard Test Data Sets

The following standard test data sets are available:

- **Molecules**: Basic cryoprotectant molecules with properties
- **Property Types**: Common property types for molecules and mixtures
- **Molecular Properties**: Properties of standard test molecules
- **Mixtures**: Predefined mixtures of cryoprotectant molecules
- **Mixture Components**: Components of standard test mixtures
- **Calculation Methods**: Methods for property predictions
- **Predictions**: Predicted properties for molecules and mixtures
- **Experiments**: Experimental data for mixtures
- **Experiment Properties**: Properties measured in experiments
- **Teams**: Research teams
- **Projects**: Research projects
- **User Profiles**: User accounts and metadata

### Data Generation Utilities

The following data generation functions are available:

- `generate_molecule()`: Generate a test molecule
- `generate_property_type()`: Generate a test property type
- `generate_molecular_property()`: Generate a test molecular property
- `generate_mixture()`: Generate a test mixture
- `generate_mixture_component()`: Generate a test mixture component
- `generate_calculation_method()`: Generate a test calculation method
- `generate_prediction()`: Generate a test prediction
- `generate_experiment()`: Generate a test experiment
- `generate_experiment_property()`: Generate a test experiment property
- `generate_project()`: Generate a test project
- `generate_team()`: Generate a test team
- `generate_user_profile()`: Generate a test user profile

### Data Loading Utilities

The following data loading functions are available:

- `load_molecules()`: Load standard test molecules
- `load_property_types()`: Load standard test property types
- `load_mixtures()`: Load standard test mixtures
- `load_experiments()`: Load standard test experiments
- `load_calculation_methods()`: Load standard test calculation methods
- `load_teams()`: Load standard test teams
- `load_projects()`: Load standard test projects
- `load_related_data()`: Load related data for a specific entity

## Usage Examples

### Loading Standard Test Data

```python
from tests.fixtures.data import load_molecules, load_mixtures

# Load standard test molecules
molecules = load_molecules()

# Load standard test mixtures
mixtures = load_mixtures()

# Use the data in tests
assert len(molecules) > 0
assert len(mixtures) > 0
```

### Generating Custom Test Data

```python
from tests.fixtures.data import (
    generate_molecule,
    generate_mixture,
    generate_experiment
)

# Generate a custom molecule
molecule = generate_molecule(
    name="Custom Molecule",
    formula="C3H8O3",
    smiles="C(C(CO)O)O",
    molecular_weight=92.09
)

# Generate a custom mixture with components
mixture = generate_mixture(
    name="Custom Mixture",
    description="A custom test mixture",
    with_components=True,
    component_count=3
)

# Generate a custom experiment with properties
experiment = generate_experiment(
    name="Custom Experiment",
    mixture_id=mixture['id'],
    with_properties=True,
    property_count=2
)

# Use the generated data in tests
assert molecule['name'] == "Custom Molecule"
assert len(mixture['components']) == 3
assert len(experiment['properties']) == 2
```

### Testing with Related Data

```python
from tests.fixtures.data import load_mixtures, load_related_data

# Load standard test mixtures
mixtures = load_mixtures()
mixture_id = mixtures[0]['id']

# Load related data for a mixture
related_data = load_related_data('mixture', mixture_id)

# Access components, experiments, and predictions
components = related_data['components']
experiments = related_data['experiments']
predictions = related_data['predictions']

# Use the related data in tests
assert len(components) > 0
assert len(experiments) > 0
assert len(predictions) > 0
```

### Integration with Database Fixtures

```python
from tests.fixtures.data import load_molecules, generate_molecular_property

def test_with_database_fixtures(mock_db):
    # Load standard test molecules
    molecules = load_molecules()
    
    # Add molecules to the mock database
    for molecule in molecules:
        mock_db.table('molecules').insert(molecule).execute()
    
    # Generate and add a custom property
    molecule_id = molecules[0]['id']
    property_type_id = "prop-type-1"  # Molecular Weight
    
    property = generate_molecular_property(
        molecule_id=molecule_id,
        property_type_id=property_type_id,
        value=180.16
    )
    
    mock_db.table('molecular_properties').insert(property).execute()
    
    # Query the database
    result = mock_db.table('molecular_properties')\
                    .select()\
                    .eq('molecule_id', molecule_id)\
                    .execute()
    
    # Verify the result
    assert len(result.data) == 1
    assert result.data[0]['value'] == 180.16
```

### Integration with API Fixtures

```python
from tests.fixtures.data import load_molecules

def test_with_api_fixtures(authenticated_client, mock_db):
    # Load standard test molecules
    molecules = load_molecules()
    
    # Add molecules to the mock database
    for molecule in molecules:
        mock_db.table('molecules').insert(molecule).execute()
    
    # Make API request
    response = authenticated_client.get('/api/molecules')
    
    # Verify the response
    assert response.status_code == 200
    data = authenticated_client.parse_json(response)
    assert len(data) == len(molecules)
```

### Integration with Mock Objects

```python
from tests.fixtures.data import load_molecules
from tests.fixtures.mocks.rdkit import patch_rdkit

def test_with_mock_objects():
    # Load standard test molecules
    molecules = load_molecules()
    molecule = molecules[0]
    
    # Use mock RDKit
    with patch_rdkit() as rdkit:
        # Create a molecule from SMILES
        mol = rdkit.Chem.MolFromSmiles(molecule['smiles'])
        
        # Calculate properties
        weight = rdkit.Descriptors.MolWt(mol)
        
        # Verify the result
        assert weight == 180.16
```

## Test Data Schema

### Molecules

```json
{
  "id": "UUID",
  "name": "String",
  "formula": "String",
  "smiles": "String",
  "molecular_weight": "Number",
  "cid": "Number",
  "inchikey": "String",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Property Types

```json
{
  "id": "UUID",
  "name": "String",
  "data_type": "String (numeric, text, boolean)",
  "description": "String",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Molecular Properties

```json
{
  "id": "UUID",
  "molecule_id": "UUID",
  "property_type_id": "UUID",
  "value": "Any (Number, String, Boolean)",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Mixtures

```json
{
  "id": "UUID",
  "name": "String",
  "description": "String",
  "created_by": "UUID (optional)",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Mixture Components

```json
{
  "id": "UUID",
  "mixture_id": "UUID",
  "molecule_id": "UUID",
  "concentration": "Number",
  "concentration_unit": "String",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Calculation Methods

```json
{
  "id": "UUID",
  "name": "String",
  "description": "String",
  "version": "String",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Predictions

```json
{
  "id": "UUID",
  "molecule_id": "UUID (optional)",
  "mixture_id": "UUID (optional)",
  "property_type_id": "UUID",
  "calculation_method_id": "UUID",
  "predicted_value": "Any (Number, String, Boolean)",
  "confidence": "Number",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Experiments

```json
{
  "id": "UUID",
  "mixture_id": "UUID",
  "name": "String",
  "date_performed": "ISO DateTime",
  "conditions": "Object",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Experiment Properties

```json
{
  "id": "UUID",
  "experiment_id": "UUID",
  "property_type_id": "UUID",
  "value": "Any (Number, String, Boolean)",
  "unit": "String",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Teams

```json
{
  "id": "UUID",
  "name": "String",
  "description": "String",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### Projects

```json
{
  "id": "UUID",
  "name": "String",
  "description": "String",
  "team_id": "UUID",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

### User Profiles

```json
{
  "id": "UUID",
  "user_id": "UUID",
  "name": "String",
  "email": "String",
  "team_id": "UUID (optional)",
  "role": "String",
  "created_at": "ISO DateTime",
  "updated_at": "ISO DateTime"
}
```

## Best Practices

1. **Use standard test data** when possible to ensure consistency across tests.

2. **Generate custom data** when you need specific test scenarios or edge cases.

3. **Maintain relationships** between entities to ensure realistic test data.

4. **Combine with other fixtures** (Database, API, Mock Objects) for comprehensive testing.

5. **Verify data integrity** by checking that generated data meets expected constraints.

6. **Use descriptive names** for custom test data to make tests more readable.

7. **Document test data** used in tests to make them easier to understand and maintain.

## Extending the Fixtures

### Adding New Entity Types

To add a new entity type:

1. Create a new JSON file with standard test data
2. Add a generator function in `generators.py`
3. Add a loader function in `loaders.py`
4. Update the `__init__.py` exports
5. Update this README with documentation and examples

### Modifying Existing Entity Types

To modify an existing entity type:

1. Update the JSON file with new standard test data
2. Update the generator function in `generators.py` if needed
3. Update the loader function in `loaders.py` if needed
4. Update this README with documentation and examples

## Troubleshooting

### Missing Data

If test data is missing:

1. Check that the JSON file exists and is properly formatted
2. Verify that the loader function is correctly implemented
3. Ensure that the `__init__.py` exports the loader function

### Inconsistent Data

If test data is inconsistent:

1. Check the relationships between entities
2. Verify that IDs match between related entities
3. Ensure that data types are consistent

### Integration Issues

If integration with other fixtures fails:

1. Check that the data format matches the expected format
2. Verify that the other fixtures are correctly implemented
3. Ensure that the integration code is correct