# Using RDKit with CryoProtect

This guide explains how to use RDKit with the CryoProtect project, which is essential for molecular property calculations and visualizations.

## Overview

RDKit is a collection of cheminformatics and machine-learning tools that provide functionality for:
- Working with molecular structures
- Calculating molecular properties and descriptors
- Generating molecular fingerprints
- Visualizing molecules
- Performing substructure searching

In the CryoProtect application, RDKit is used for:
1. Calculating physical properties of cryoprotectants
2. Generating visualizations of molecules
3. Performing similarity calculations
4. Creating and analyzing molecular scaffolds

## Installation Options

There are three ways to use RDKit with CryoProtect:

### 1. Using the Podman Container (Recommended)

The simplest way to use RDKit is through the pre-configured Podman container:

```bash
# Set up the container with RDKit
./quick_conda_container.sh

# Run a script with RDKit support
./run_with_rdkit.sh your_script.py
```

### 2. Local Conda Installation

If you prefer a local installation:

```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create environment with RDKit
conda create -n cryoprotect python=3.10
conda activate cryoprotect
conda install -c conda-forge rdkit
```

### 3. Using Mock RDKit (For Testing Only)

For simple testing without full RDKit functionality:

```bash
# Set up mock RDKit
python mock_rdkit.py
export PYTHONPATH="/tmp/mock_modules:$PYTHONPATH"
```

## Wrapper Scripts

The repository includes wrapper scripts to simplify RDKit usage:

1. `run_with_rdkit.sh` - Run any Python script with RDKit support
2. `run_unified_import_in_container.sh` - Run the unified ChEMBL importer with RDKit support
3. `verify_rdkit_standalone.py` - Verify RDKit installation and functionality

## Usage Examples

### Calculating Properties with RDKit

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

# Create a molecule from SMILES
mol = Chem.MolFromSmiles("CCO")  # Ethanol

# Calculate properties
molecular_weight = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)
tpsa = Descriptors.TPSA(mol)

print(f"Molecular Weight: {molecular_weight:.2f}")
print(f"LogP: {logp:.2f}")
print(f"TPSA: {tpsa:.2f}")
```

### Generating Molecule Visualizations

```python
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

# Create a molecule from SMILES
mol = Chem.MolFromSmiles("CCO")  # Ethanol

# Generate SVG visualization
drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
svg = drawer.GetDrawingText()

# Save to file
with open("molecule.svg", "w") as f:
    f.write(svg)
```

## Troubleshooting

### Common Issues

1. **ImportError: No module named 'rdkit'**
   - Solution: Use the container with `./run_with_rdkit.sh` or check your conda environment

2. **Error in molecule parsing**
   - Solution: Verify your SMILES strings are valid

3. **3D coordinate generation fails**
   - Solution: This requires the full RDKit installation, not the mock version

### Verifying RDKit Installation

To verify your RDKit installation:

```bash
# Run the standalone verification script
./run_with_rdkit.sh verify_rdkit_standalone.py --all
```

## Best Practices

1. Always use the container for consistent RDKit availability
2. Add error handling around RDKit operations
3. Include fallback strategies when RDKit is unavailable
4. For production, ensure your environment has the correct RDKit version

## Additional Resources

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [RDKit GitHub Repository](https://github.com/rdkit/rdkit)
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)