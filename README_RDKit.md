# CryoProtect Analyzer - RDKit Integration

This document provides information about the RDKit integration in the CryoProtect Analyzer application, including setup instructions, API endpoints, and usage examples.

## Overview

The RDKit integration adds powerful molecular property calculation, visualization, and searching capabilities to the CryoProtect Analyzer. It enables:

- Calculation of hydrogen bonding capacity (donors/acceptors)
- Computation of XLogP (partition coefficient)
- Determination of topological polar surface area
- Calculation of molecular weight and volume
- Identification of functional groups
- Estimation of permeability coefficients
- Molecular structure visualization
- Substructure and similarity searching

## Installation

RDKit is a cheminformatics and machine learning toolkit that must be installed separately from the main application dependencies. The recommended installation method is via conda:

```bash
# Create a new conda environment
conda create -n cryoprotect python=3.9
conda activate cryoprotect

# Install RDKit
conda install -c conda-forge rdkit

# Install other dependencies
pip install -r requirements.txt
```

Note: If you're not using conda, you can install RDKit via pip, but it may be more complex:

```bash
pip install rdkit
```

## API Endpoints

The RDKit integration adds the following API endpoints:

### Calculate Molecular Properties

**Endpoint:** `/api/v1/rdkit/properties`  
**Method:** POST  
**Description:** Calculate molecular properties using RDKit  
**Request Body:**
```json
{
  "molecule_data": "CCO",
  "input_format": "smiles"
}
```
**Response:**
```json
{
  "hydrogen_bonding": {
    "donors": 1,
    "acceptors": 1,
    "total": 2
  },
  "logp": -0.14,
  "tpsa": 20.23,
  "molecular_properties": {
    "molecular_weight": 46.07,
    "exact_mass": 46.042,
    "heavy_atom_count": 3,
    "atom_count": 9,
    "rotatable_bond_count": 1,
    "ring_count": 0,
    "aromatic_ring_count": 0,
    "fraction_csp3": 1.0,
    "molecular_volume": 53.12
  },
  "functional_groups": {
    "alcohol": 1,
    "hydroxyl": 1
  },
  "permeability": {
    "rule_of_5_violations": 0,
    "veber_violations": 0,
    "bbb_permeant": true,
    "intestinal_absorption": true,
    "estimated_log_papp": -4.5
  },
  "smiles": "CCO",
  "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
  "inchi_key": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
}
```

### Generate Molecular Visualization

**Endpoint:** `/api/v1/rdkit/visualization`  
**Method:** POST  
**Description:** Generate a visualization of a molecule as SVG  
**Request Body:**
```json
{
  "molecule_data": "CCO",
  "input_format": "smiles",
  "width": 400,
  "height": 300,
  "highlight_atoms": [0, 1]
}
```
**Response:**
```json
{
  "svg": "<svg version='1.1' baseProfile='full' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' xml:space='preserve' width='400' height='300' viewBox='0 0 400 300'>...</svg>",
  "width": 400,
  "height": 300
}
```

### Perform Substructure Search

**Endpoint:** `/api/v1/rdkit/substructure`  
**Method:** POST  
**Description:** Search for a substructure within a molecule  
**Request Body:**
```json
{
  "query_mol_data": "[OH]",
  "target_mol_data": "CCO",
  "query_format": "smarts",
  "target_format": "smiles"
}
```
**Response:**
```json
{
  "match": true,
  "match_count": 1,
  "matches": [[2]],
  "visualization": "<svg version='1.1' baseProfile='full' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' xml:space='preserve' width='400' height='300' viewBox='0 0 400 300'>...</svg>"
}
```

### Calculate Molecular Similarity

**Endpoint:** `/api/v1/rdkit/similarity`  
**Method:** POST  
**Description:** Calculate similarity between two molecules  
**Request Body:**
```json
{
  "mol1_data": "CCO",
  "mol2_data": "CC(=O)O",
  "mol1_format": "smiles",
  "mol2_format": "smiles",
  "fingerprint_type": "morgan"
}
```
**Response:**
```json
{
  "tanimoto": 0.42,
  "dice": 0.59,
  "fingerprint_type": "morgan"
}
```

### Calculate Properties for Database Molecule

**Endpoint:** `/api/v1/molecules/{molecule_id}/calculate-properties`  
**Method:** POST  
**Description:** Calculate and store properties for a molecule in the database  
**Response:**
```json
{
  "message": "Calculated and stored 15 properties for molecule 12345"
}
```

## Web Interface

The RDKit integration includes a dedicated web interface for molecular visualization and property calculation, accessible at `/molecules/rdkit`. This interface provides:

1. **Structure Viewer**: Visualize molecular structures and view basic properties
2. **Property Calculator**: Calculate detailed molecular properties
3. **Structure Search**: Perform substructure and similarity searches
4. **Upload Structure**: Upload new molecular structures to the database

## JavaScript API

The RDKit integration includes a JavaScript API for client-side interaction with the RDKit endpoints. The API is available as `window.RDKitAPI` and includes the following functions:

- `calculateProperties(moleculeData, inputFormat)`: Calculate molecular properties
- `generateVisualization(moleculeData, inputFormat, width, height, highlightAtoms)`: Generate a molecular visualization
- `performSubstructureSearch(queryMolData, targetMolData, queryFormat, targetFormat)`: Perform a substructure search
- `calculateSimilarity(mol1Data, mol2Data, mol1Format, mol2Format, fingerprintType)`: Calculate molecular similarity
- `calculateMoleculeProperties(moleculeId)`: Calculate properties for a molecule in the database
- `displayMolecularProperties(properties, container)`: Display molecular properties in a container element
- `displayMolecularVisualization(svgData, container)`: Display a molecular visualization in a container element
- `initMoleculeSearchForm(formId, resultsContainerId)`: Initialize a molecule search form

## Usage Examples

### Python Example: Calculate Properties

```python
from api.rdkit_utils import calculate_all_properties

# Calculate properties for ethanol
properties = calculate_all_properties("CCO")
print(f"LogP: {properties['logp']}")
print(f"TPSA: {properties['tpsa']}")
print(f"H-Bond Donors: {properties['hydrogen_bonding']['donors']}")
print(f"H-Bond Acceptors: {properties['hydrogen_bonding']['acceptors']}")
```

### JavaScript Example: Visualize Molecule

```javascript
// Get a reference to the container element
const container = document.getElementById('moleculeVisualization');

// Generate and display a visualization
RDKitAPI.generateVisualization("CCO", "smiles", 400, 300)
  .then(result => {
    RDKitAPI.displayMolecularVisualization(result.svg, container);
  })
  .catch(error => {
    console.error("Error visualizing molecule:", error);
  });
```

### JavaScript Example: Search for a Substructure

```javascript
// Search for a hydroxyl group in ethanol
RDKitAPI.performSubstructureSearch("[OH]", "CCO", "smarts", "smiles")
  .then(result => {
    if (result.match) {
      console.log(`Found ${result.match_count} matches!`);
      // Display the visualization with highlighted matches
      document.getElementById('visualization').innerHTML = result.visualization;
    } else {
      console.log("No matches found");
    }
  })
  .catch(error => {
    console.error("Error performing substructure search:", error);
  });
```

## Troubleshooting

### Common Issues

1. **RDKit Import Error**: If you see an error like "No module named 'rdkit'", make sure you have installed RDKit correctly and are using the correct Python environment.

2. **Molecule Parsing Error**: If you see an error when parsing molecules, check that your SMILES, MOL, or SDF data is valid. You can validate SMILES strings using online tools like the [SMILES Translator](https://cactus.nci.nih.gov/translate/).

3. **SVG Rendering Issues**: If molecular visualizations don't appear correctly, check that your browser supports SVG rendering and that the SVG data is valid.

### Getting Help

If you encounter issues with the RDKit integration, check the following resources:

- [RDKit Documentation](https://www.rdkit.org/docs/index.html)
- [RDKit GitHub Repository](https://github.com/rdkit/rdkit)
- [RDKit Google Group](https://groups.google.com/g/rdkit-discuss)

## References

- [RDKit: Open-source cheminformatics](https://www.rdkit.org/)
- [Lipinski's Rule of Five](https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five)
- [SMILES Notation](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
- [InChI](https://en.wikipedia.org/wiki/International_Chemical_Identifier)