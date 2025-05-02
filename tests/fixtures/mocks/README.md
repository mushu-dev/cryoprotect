# Mock Objects

This directory contains mock implementations of external dependencies for testing.

## Available Mocks

### Supabase Mock
- `MockSupabase`: Mock implementation of the Supabase client
- `patch_supabase_client`: Context manager to patch Supabase client imports

### RDKit Mock
- `MockRDKit`: Mock implementation of RDKit
- `MockMolecule`: Mock molecule implementation
- `patch_rdkit`: Context manager to patch RDKit imports
- `mock_rdkit`: Fixture providing a mock RDKit instance
- `molecule_factory`: Factory function for creating mock molecules

## Usage Examples

### Supabase Mock

```python
from tests.fixtures.mocks.supabase import patch_supabase_client

def test_with_mock_supabase():
    with patch_supabase_client() as mock_db:
        # Add test data
        mock_db.add_test_data('molecules', [
            {'id': 'mol-1', 'name': 'Test Molecule'}
        ])
        
        # Use the mock in your test
        result = mock_db.table('molecules').select().execute()
        assert len(result.data) == 1
```

### RDKit Mock

```python
from tests.fixtures.mocks.rdkit import patch_rdkit

def test_with_mock_rdkit():
    with patch_rdkit() as rdkit:
        # Create a molecule
        mol = rdkit.Chem.MolFromSmiles('CCO')
        
        # Calculate properties
        weight = rdkit.Descriptors.MolWt(mol)
        logp = rdkit.Descriptors.MolLogP(mol)
        
        # Use the properties in your test
        assert weight == 180.16
        assert logp == -0.5
```

### Using Fixtures

```python
def test_with_mock_fixtures(mock_rdkit, molecule_factory):
    # Use the mock_rdkit fixture
    mol = mock_rdkit.Chem.MolFromSmiles('CCO')
    
    # Create a custom molecule with the factory
    custom_mol = molecule_factory(
        smiles='CCC',
        mol_weight=44.1,
        formula='C2H6'
    )
    
    # Use the custom molecule
    assert custom_mol.formula == 'C2H6'
```

### Combined Usage

```python
from tests.fixtures.mocks.supabase import patch_supabase_client
from tests.fixtures.mocks.rdkit import patch_rdkit

def test_with_both_mocks():
    with patch_supabase_client() as mock_db, patch_rdkit() as rdkit:
        # Add test data
        mock_db.add_test_data('molecules', [
            {'id': 'mol-1', 'name': 'Ethanol', 'smiles': 'CCO'}
        ])
        
        # Query database
        result = mock_db.table('molecules').select().eq('name', 'Ethanol').execute()
        molecule_data = result.data[0]
        
        # Process with RDKit
        mol = rdkit.Chem.MolFromSmiles(molecule_data['smiles'])
        weight = rdkit.Descriptors.MolWt(mol)
        
        # Use the result
        assert weight == 180.16
```

## Extending the Mocks

### Adding Custom RPC Handlers

```python
def test_with_custom_rpc():
    with patch_supabase_client() as mock_db:
        # Register a custom RPC handler
        mock_db.register_rpc_handler(
            'calculate_similarity',
            lambda params: {'similarity': 0.85}
        )
        
        # Use the custom RPC
        result = mock_db.rpc('calculate_similarity', {
            'smiles1': 'CCO',
            'smiles2': 'CCCO'
        }).execute()
        
        assert result.data[0]['similarity'] == 0.85
```

### Creating Custom Molecules

```python
from tests.fixtures.mocks.rdkit import MockMolecule

def test_with_custom_molecule():
    # Create a custom molecule directly
    mol = MockMolecule(
        smiles='C1CCCCC1',
        mol_weight=84.16,
        logp=1.7,
        h_donors=0,
        h_acceptors=0,
        formula='C6H12'
    )
    
    # Use the custom molecule
    assert mol.GetFormula() == 'C6H12'
    assert mol.GetMolWt() == 84.16
```

## Best Practices

1. **Use context managers to ensure proper cleanup**:
   ```python
   with patch_supabase_client() as mock_db:
       # Test code here
   ```

2. **Add custom test data for each test**:
   ```python
   mock_db.add_test_data('molecules', [
       {'id': 'mol-1', 'name': 'Test Molecule', 'smiles': 'CCO'}
   ])
   ```

3. **Use the molecule factory for creating molecules with custom properties**:
   ```python
   custom_mol = molecule_factory(
       smiles='CCC',
       mol_weight=44.1,
       formula='C2H6'
   )
   ```

4. **Register custom RPC handlers for testing specific functionality**:
   ```python
   mock_db.register_rpc_handler(
       'my_function',
       lambda params: {'result': 'success'}
   )
   ```

5. **Verify mock state after operations**:
   ```python
   # After an insert
   result = mock_db.table('molecules').select().execute()
   assert len(result.data) == expected_count