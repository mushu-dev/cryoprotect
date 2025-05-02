# CryoProtect Analyzer - RDKit Integration

This document provides an overview of the RDKit integration in the CryoProtect Analyzer application, including the enhanced functionality for cryoprotectant analysis.

## Overview

The CryoProtect Analyzer uses RDKit, a powerful cheminformatics toolkit, to perform molecular operations, property calculations, and chemical searches. Since direct RDKit integration in Supabase PostgreSQL is not available, we've implemented RDKit functionality at the application level.

The integration consists of two main components:

1. **Core RDKit Utilities** (`api/rdkit_utils.py`): Basic molecular property calculations and operations
2. **Enhanced RDKit Module** (`api/rdkit_enhanced.py`): Advanced functionality specific to cryoprotectant analysis

## Features

### Core RDKit Functionality

- Molecular property calculation (LogP, TPSA, molecular weight, etc.)
- Hydrogen bonding analysis
- Functional group identification
- Permeability estimation
- Molecular visualization
- Substructure search
- Similarity calculation

### Enhanced RDKit Functionality

- **Property Caching**: Avoids redundant calculations for the same molecules
- **Cryoprotectant-Specific Properties**:
  - Glass transition temperature estimation
  - Vitrification tendency scoring
  - Cell membrane permeability for cryoprotection
  - Toxicity risk assessment
  - Overall cryoprotectant score
- **Batch Processing**: Efficient handling of multiple molecules
- **Advanced Search Capabilities**: Find similar molecules or substructures across multiple targets
- **Scaffold Analysis**: Identify and analyze molecular scaffolds
- **Molecule Grid Visualization**: Generate grids of multiple molecules

## API Endpoints

### Core RDKit Endpoints

- `/api/v1/rdkit/properties`: Calculate molecular properties
- `/api/v1/rdkit/visualization`: Generate molecule visualization
- `/api/v1/rdkit/substructure`: Perform substructure search
- `/api/v1/rdkit/similarity`: Calculate molecular similarity
- `/api/v1/molecules/<molecule_id>/calculate-properties`: Calculate and store properties for a specific molecule

### Enhanced RDKit Endpoints

- `/api/v1/rdkit/cryoprotectant-properties`: Calculate cryoprotectant-specific properties
- `/api/v1/rdkit/batch-properties`: Calculate properties for multiple molecules
- `/api/v1/rdkit/batch-similarity`: Find similar molecules from a list of targets
- `/api/v1/rdkit/batch-substructure`: Perform substructure search on multiple targets
- `/api/v1/rdkit/molecule-grid`: Generate a grid of molecule visualizations
- `/api/v1/rdkit/scaffold-analysis`: Analyze molecular scaffolds
- `/api/v1/rdkit/property-cache`: Manage the property cache
- `/api/v1/molecules/batch-calculate-properties`: Calculate and store properties for multiple molecules

## Usage Examples

### Calculating Cryoprotectant Properties

```python
import requests
import json

# API endpoint
url = "http://localhost:5000/api/v1/rdkit/cryoprotectant-properties"

# Glycerol (a common cryoprotectant)
data = {
    "molecule_data": "C(C(CO)O)O",
    "input_format": "smiles"
}

# Send request
response = requests.post(url, json=data)
properties = response.json()

# Access cryoprotectant-specific properties
cryo_props = properties["cryoprotectant_properties"]
print(f"Vitrification tendency: {cryo_props['vitrification_tendency']}")
print(f"Overall cryoprotectant score: {cryo_props['overall_cryoprotectant_score']}")
```

### Batch Property Calculation

```python
import requests
import json

# API endpoint
url = "http://localhost:5000/api/v1/rdkit/batch-properties"

# List of molecules to analyze
data = {
    "molecules": [
        {"data": "CCO", "format": "smiles"},  # Ethanol
        {"data": "C(C(CO)O)O", "format": "smiles"},  # Glycerol
        {"data": "CS(=O)C", "format": "smiles"}  # DMSO
    ]
}

# Send request
response = requests.post(url, json=data)
results = response.json()

# Process results
for result in results:
    molecule_index = result["molecule_index"]
    smiles = result["smiles"]
    logp = result["logp"]
    print(f"Molecule {molecule_index} ({smiles}): LogP = {logp}")
```

### Finding Similar Molecules

```python
import requests
import json

# API endpoint
url = "http://localhost:5000/api/v1/rdkit/batch-similarity"

# Find molecules similar to glycerol
data = {
    "query_mol_data": "C(C(CO)O)O",  # Glycerol
    "target_molecules": [
        "CCO",  # Ethanol
        "C(C(CO)O)O",  # Glycerol
        "CS(=O)C",  # DMSO
        "OCCO",  # Ethylene glycol
        "C(CO)(CO)CO"  # Glycerol
    ],
    "fingerprint_type": "morgan",
    "similarity_threshold": 0.5
}

# Send request
response = requests.post(url, json=data)
results = response.json()

# Process results
for result in results:
    smiles = result["smiles"]
    tanimoto = result["tanimoto"]
    print(f"Molecule {smiles}: Tanimoto similarity = {tanimoto}")
```

## Implementation Details

### Property Caching

The property cache uses a combination of in-memory LRU cache and disk-based storage to avoid redundant calculations:

1. When a property calculation is requested, the system first checks if the result is already in the cache
2. If found, the cached result is returned immediately
3. If not found, the calculation is performed and the result is stored in the cache for future use

The cache key is based on a hash of the canonical SMILES string, ensuring consistent caching regardless of SMILES representation variations.

### Cryoprotectant Property Models

The cryoprotectant property models are based on molecular features known to be important for cryoprotection:

1. **Glass Transition Temperature**: Estimated based on molecular weight, hydrogen bonding capacity, and LogP
2. **Vitrification Tendency**: Scored based on factors that promote vitrification (high molecular weight, many hydrogen bonds, low LogP, high TPSA)
3. **Cell Membrane Permeability**: Optimized for cryoprotection, balancing permeability with stability
4. **Toxicity Risk**: Assessed based on structural features and physicochemical properties
5. **Overall Cryoprotectant Score**: Weighted combination of the above factors

These models provide a computational approach to screening potential cryoprotectants before experimental validation.

## Testing

A comprehensive test suite (`test_rdkit_enhanced.py`) is provided to verify the functionality of the RDKit integration. To run the tests:

```bash
python test_rdkit_enhanced.py
```

The tests cover all major features, including:
- Property calculation and caching
- Cryoprotectant-specific properties
- Batch processing
- Similarity and substructure searching
- Visualization and scaffold analysis

## Troubleshooting

If you encounter issues with the RDKit integration, please refer to the `README_RDKit_Troubleshooting.md` document for detailed troubleshooting steps.

Common issues include:
- RDKit installation problems
- Missing dependencies
- Memory errors with large molecules or batch operations

## Future Enhancements

Planned enhancements for the RDKit integration include:

1. **Machine Learning Models**: Train models on experimental cryoprotectant data to improve property predictions
2. **3D Conformer Generation**: Generate and analyze 3D conformers for more accurate property calculations
3. **Molecular Dynamics Integration**: Interface with molecular dynamics simulations for deeper analysis
4. **Interactive Visualizations**: Enhance the web interface with interactive molecular visualizations
5. **Database Integration**: Explore options for direct RDKit integration with Supabase PostgreSQL

## References

1. RDKit Documentation: [https://www.rdkit.org/docs/](https://www.rdkit.org/docs/)
2. Lipinski, C. A., et al. (1997). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. Advanced Drug Delivery Reviews, 23(1-3), 3-25.
3. Veber, D. F., et al. (2002). Molecular properties that influence the oral bioavailability of drug candidates. Journal of Medicinal Chemistry, 45(12), 2615-2623.
4. Daina, A., & Zoete, V. (2016). A BOILED-Egg to predict gastrointestinal absorption and brain penetration of small molecules. ChemMedChem, 11(11), 1117-1121.