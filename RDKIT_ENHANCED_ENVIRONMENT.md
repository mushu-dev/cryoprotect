# Enhanced RDKit Environment Setup

This document outlines the setup and usage of the enhanced RDKit environment for CryoProtect. The enhanced environment provides advanced molecular modeling capabilities beyond the basic RDKit functionality, including 3D conformer generation, sophisticated descriptor calculations, and integration with machine learning models.

## Components

1. **Docker Container**
   - Custom container based on Fedora with full RDKit installation
   - Additional Python packages for machine learning and molecular modeling
   - GPU acceleration support for conformer generation and calculations

2. **Advanced Molecular Descriptors**
   - 3D pharmacophore features
   - Electronic property calculations
   - Topological descriptors
   - Molecular interaction field descriptors
   - Fragment-based descriptors

3. **3D Conformer Generation**
   - MMFF94 force field optimization
   - Conformer ensemble generation
   - Energy minimization
   - Conformer clustering

4. **Integration with Scientific Models**
   - Python interface for all models
   - Database storage of calculated descriptors
   - Caching system for expensive calculations
   - Batch processing capabilities

## Setup Instructions

### 1. Building the Container

```bash
# Build the enhanced RDKit container
cd /home/mushu/Projects/CryoProtect
docker build -f Dockerfile.rdkit-enhanced -t cryoprotect-rdkit-enhanced .
```

### 2. Running the Container

```bash
# Run the container with GPU support (if available)
docker run --gpus all -d -p 5001:5000 --name cryoprotect-rdkit-enhanced cryoprotect-rdkit-enhanced
```

### 3. Using the Python API

```python
from scientific_models.rdkit_enhanced import EnhancedRDKitCalculator

# Initialize the calculator
calculator = EnhancedRDKitCalculator()

# Calculate descriptors for a molecule
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
descriptors = calculator.calculate_descriptors(smiles, include_3d=True)

# Generate conformers
conformers = calculator.generate_conformers(smiles, n_conformers=10)

# Predict property using 3D information
predicted_value = calculator.predict_property("solubility", smiles)
```

## API Endpoints

The enhanced RDKit environment provides the following API endpoints:

- `POST /api/v1/rdkit/descriptors` - Calculate molecular descriptors
- `POST /api/v1/rdkit/conformers` - Generate 3D conformers
- `POST /api/v1/rdkit/minimize` - Energy minimize a molecule
- `POST /api/v1/rdkit/similarity/3d` - Calculate 3D similarity between molecules
- `POST /api/v1/rdkit/pharmacophore` - Generate pharmacophore model
- `POST /api/v1/rdkit/property/predict` - Predict properties using ML models

## Integration with Scientific Models

The enhanced RDKit environment is integrated with our scientific models:

1. Temperature models use electronic properties to predict temperature dependence
2. Concentration models use topological descriptors for concentration effects
3. Mixture optimization uses 3D similarity and interaction potential
4. Glass transition temperature prediction uses fragment contributions

## Troubleshooting

Common issues and their solutions:

1. **GPU Acceleration Issues**
   - Ensure NVIDIA drivers are installed
   - Check Docker GPU support with `docker info | grep -i gpu`

2. **Memory Limitations**
   - For large molecules, increase container memory allocation
   - Use batch processing for multiple molecules

3. **Integration Issues**
   - Check API URL configuration in `config.py`
   - Verify network connectivity between services

## References

1. RDKit Documentation: https://www.rdkit.org/docs/index.html
2. Scientific Models Integration Guide: See `scientific_models/README.md`
3. Docker GPU Support: https://docs.docker.com/config/containers/resource_constraints/#gpu