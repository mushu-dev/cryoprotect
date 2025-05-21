# CryoProtect RDKit Integration

This README provides a comprehensive guide to the RDKit integration in the CryoProtect project, explaining how to use the wrapper module, run the containerized environment, and extend the functionality.

## Quick Start

### Using the RDKit Wrapper

```python
# Import the wrapper
import rdkit_wrapper

# Create a molecule
molecule = rdkit_wrapper.create_molecule_from_smiles("CCO")

# Calculate properties
properties = rdkit_wrapper.calculate_properties(molecule)
print(f"Molecular Weight: {properties['molecular_weight']}")
print(f"LogP: {properties['logp']}")

# Generate SVG visualization
svg = rdkit_wrapper.generate_molecule_svg(molecule)

# Check if RDKit is available
if rdkit_wrapper.RDKIT_AVAILABLE:
    # Perform operations that require full RDKit
    fingerprint = rdkit_wrapper.generate_fingerprint(molecule)
    mol_3d = rdkit_wrapper.generate_molecule_3d_coordinates(molecule)
```

### Running with RDKit Container

```bash
# Run a script with the RDKit container
./run_with_rdkit_container.sh your_script.py

# Build or rebuild the container
./build_rdkit_container.sh
```

## Components

1. **rdkit_wrapper.py**: Unified interface to RDKit functionality with fallback to mock implementation
2. **mock_rdkit_formula.py**: Mock implementation for basic property calculations
3. **build_rdkit_container.sh**: Script to build the RDKit container
4. **run_with_rdkit_container.sh**: Script to run Python scripts with RDKit container
5. **rdkit_integration_test.py**: Comprehensive integration test

## Setup

### Local Development

1. Install RDKit via pip or conda:
   ```bash
   pip install rdkit-pypi
   # OR with conda
   conda install -c conda-forge rdkit
   ```

2. Make sure the wrapper modules are in your Python path:
   ```bash
   cp rdkit_wrapper.py mock_rdkit_formula.py /path/to/your/project/
   ```

### Container Deployment

1. Clone the repository:
   ```bash
   git clone https://github.com/your-org/cryoprotect.git
   cd cryoprotect
   ```

2. Build the container:
   ```bash
   ./build_rdkit_container.sh
   ```

3. Run scripts with the container:
   ```bash
   ./run_with_rdkit_container.sh path/to/your/script.py
   ```

## Wrapper Functions

| Function | Description |
|----------|-------------|
| `create_molecule_from_smiles(smiles)` | Creates a molecule from SMILES string |
| `calculate_properties(mol_or_smiles)` | Calculates molecular properties |
| `generate_fingerprint(mol_or_smiles)` | Generates Morgan fingerprint |
| `calculate_similarity(mol1, mol2)` | Calculates similarity between molecules |
| `perform_substructure_search(query, target)` | Performs substructure search |
| `generate_molecule_svg(mol_or_smiles)` | Generates SVG visualization |
| `generate_molecule_3d_coordinates(mol_or_smiles)` | Generates 3D coordinates |
| `get_rdkit_status()` | Returns RDKit status information |

## Property Calculations

The `calculate_properties()` function returns a dictionary with the following properties:

- molecular_formula
- molecular_weight
- logp
- tpsa
- h_donors
- h_acceptors
- rotatable_bonds
- ring_count
- aromatic_ring_count
- heavy_atom_count
- rdkit_available
- rdkit_version
- calculation_method

When RDKit is available, additional properties are calculated:
- fraction_csp3
- num_stereocenters
- qed (Quantitative Estimate of Drug-likeness)

## Testing

Run the integration test to verify RDKit functionality:

```bash
# Test in host environment
python rdkit_integration_test.py

# Test in container
./run_with_rdkit_container.sh rdkit_integration_test.py
```

## Extending the Wrapper

To add new functionality to the wrapper:

1. Add a new function to `rdkit_wrapper.py`:
   ```python
   def my_new_function(mol_or_smiles: Union[str, Any]) -> Any:
       """
       My new function documentation
       """
       if not RDKIT_AVAILABLE:
           logger.warning("Function not available in mock mode")
           return None
           
       # Convert SMILES to molecule if needed
       if isinstance(mol_or_smiles, str):
           mol = create_molecule_from_smiles(mol_or_smiles)
       else:
           mol = mol_or_smiles
           
       if not mol:
           return None
           
       # Implement RDKit functionality
       result = ...
       
       return result
   ```

2. If applicable, add a mock version to `mock_rdkit_formula.py`

3. Add tests for the new function to `rdkit_integration_test.py`

## Troubleshooting

### Common Issues

1. **ImportError: No module named 'rdkit'**
   - Verify RDKit is installed: `pip list | grep rdkit`
   - Make sure the wrapper files are in your Python path

2. **Container fails to build**
   - Check Docker/Podman installation: `podman version`
   - Verify you have sufficient permissions
   - Check internet connection for downloading dependencies

3. **Property calculations differ between environments**
   - This is expected - mock implementation provides approximations
   - For critical calculations, always use the RDKit container

## References

- [RDKit Documentation](https://www.rdkit.org/docs/index.html)
- [RDKit on GitHub](https://github.com/rdkit/rdkit)
- [Detailed Integration Guide](RDKIT_INTEGRATION_GUIDE.md)
- [Deployment Checklist](RDKIT_DEPLOYMENT_CHECKLIST.md)