# Molecule Transformation Module

The molecule transform module provides utilities for transforming molecular data between different formats, calculating properties, and handling cross-references between databases.

## MoleculeTransformer Class

The main class `MoleculeTransformer` handles all molecule-related transformations with these key features:

### Standardization

Converts between chemical identifiers and ensures all common formats are available:
- SMILES (canonical and isomeric)
- InChI and InChIKey
- Molecular formula
- Molecular weight

### Property Calculation

Calculates and standardizes molecular properties including:
- LogP (lipophilicity)
- TPSA (topological polar surface area)
- Hydrogen bond donors and acceptors
- Rotatable bond count
- Ring count
- Aromatic ring count
- Complexity metrics
- Druglikeness scores

### Cross-Reference Resolution

Resolves identifiers between different chemical databases:
- ChEMBL IDs to PubChem CIDs
- PubChem CIDs to ChEMBL IDs
- Uses multiple resolution strategies (direct API lookup, InChIKey matching, structure similarity)
- Supports batch resolution for efficiency

### Mixture Handling

Detects and handles molecular mixtures:
- Identifies disconnected structures in SMILES strings
- Extracts individual components from mixtures
- Calculates properties for each component
- Supports component naming and identification

### Molecular Scaffold Analysis

Identifies core structural features:
- Bemis-Murcko scaffolds
- Generic scaffolds for structure classification
- Enables grouping of molecules by common structural elements

### Similarity Comparison

Provides utilities for comparing molecular structures:
- Morgan fingerprints (ECFP)
- MACCS keys
- Topological fingerprints
- Tanimoto similarity calculations

## Usage Examples

### Basic Standardization

```python
transformer = MoleculeTransformer()
molecule_data = {
    'name': 'Aspirin',
    'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'
}
standardized = await transformer.standardize_molecule(molecule_data)
```

### Batch Processing

```python
molecules = [molecule1, molecule2, molecule3]
results = await transformer.standardize_molecules_batch(molecules)
```

### Cross-Reference Resolution

```python
# Enable cross-reference resolution
config = {'resolve_cross_references': True}
transformer = MoleculeTransformer(config=config)

# Add PubChem CID to a ChEMBL molecule
chembl_mol = {'name': 'Aspirin', 'chembl_id': 'CHEMBL25'}
result = await transformer.standardize_molecule(chembl_mol)
print(f"PubChem CID: {result.get('pubchem_cid')}")
```

### Mixture Analysis

```python
mixture = {'name': 'Sodium Chloride Solution', 'smiles': '[Na+].[Cl-].[H]O[H]'}
result = await transformer.standardize_molecule(mixture)

if result['is_mixture']:
    print(f"Components: {len(result['mixture_components'])}")
    for component in result['mixture_components']:
        print(f"- {component['name']}: {component['formula']}")
```

## Utility Functions

The module also provides standalone functions for common operations:

- `normalize_smiles(smiles)`: Standardize a SMILES string for consistent representation
- `get_molecule_fingerprint(smiles)`: Calculate molecular fingerprint for similarity comparisons
- `calculate_similarity(smiles1, smiles2, fingerprint_type)`: Calculate similarity between molecules

## Configuration Options

The `MoleculeTransformer` accepts several configuration options:

- `resolve_cross_references`: Whether to resolve IDs between databases (default: True)
- `handle_mixtures`: Whether to analyze mixtures (default: True)
- `separate_mixtures`: Whether to split mixtures into separate entries (default: False)
- `api_delay`: Delay between API calls to avoid rate limiting (default: 0.5 seconds)
- `batch_size`: Number of molecules to process in a batch (default: 50)
- `add_standardized_structure`: Whether to add standardized structures (default: False)

## Dependencies

This module requires:
- RDKit for chemical structure manipulation and property calculation
- aiohttp for asynchronous API calls
- asyncio for concurrent processing