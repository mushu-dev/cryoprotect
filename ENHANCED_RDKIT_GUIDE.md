# Enhanced RDKit Environment Guide

This guide explains how to use the enhanced RDKit environment for CryoProtect, which provides advanced molecular modeling capabilities beyond the basic RDKit functionality.

## Overview

The enhanced RDKit environment consists of:

1. A dedicated Docker container with a full RDKit installation and additional tools
2. A standalone service that provides advanced molecular modeling capabilities
3. Python interfaces for integrating with the CryoProtect scientific models
4. API resources for exposing these capabilities through the CryoProtect API

## Getting Started

### Starting the Enhanced RDKit Service

```bash
# Build and start the enhanced RDKit service
cd /home/mushu/Projects/CryoProtect
./run_rdkit_enhanced.sh
```

This script builds the Docker container and starts the service on port 5001.

### Verifying the Service

```bash
# Check if the service is running
curl http://localhost:5001/health

# Get information about the service
curl http://localhost:5001/api/v1/rdkit/info
```

## Using the Python Interface

The `scientific_models/rdkit_enhanced.py` module provides a Python interface for interacting with the enhanced RDKit service:

```python
from scientific_models.rdkit_enhanced import EnhancedRDKitCalculator

# Initialize the calculator
calculator = EnhancedRDKitCalculator()

# Calculate descriptors for a molecule
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
descriptors = calculator.calculate_descriptors(smiles, include_3d=True)

# Generate conformers
conformers = calculator.generate_conformers(smiles, n_conformers=10)

# Calculate 3D similarity between molecules
similarity = calculator.calculate_3d_similarity(smiles1, smiles2)

# Get pharmacophore features
features = calculator.get_pharmacophore(smiles)

# Predict properties
solubility = calculator.predict_property("solubility", smiles)
glass_transition = calculator.predict_property("glass_transition", smiles)
```

## Using the Scientific Models

The enhanced RDKit environment provides several scientific models that use the RDKit service:

```python
from scientific_models.rdkit_enhanced import (
    DescriptorBasedModel,
    MolecularSimilarityModel,
    ConformerAnalysisModel,
    PharmacophoreModel
)

# Analyze conformers
model = ConformerAnalysisModel({"n_conformers": 10})
result = model.calculate({"smiles": smiles})
print(f"Flexibility score: {result['flexibility_score']}")
print(f"Stability score: {result['stability_score']}")

# Analyze pharmacophore
model = PharmacophoreModel()
result = model.calculate({"smiles": smiles})
print(f"Number of features: {result['total_features']}")
print(f"Feature types: {result['feature_counts']}")

# Calculate similarity
model = MolecularSimilarityModel({"similarity_mode": "3d"})
result = model.calculate({
    "query_smiles": smiles1,
    "target_smiles_list": [smiles2, smiles3, smiles4]
})
for mol in result["results"]:
    print(f"Similarity to {mol['target_smiles']}: {mol['similarity_score']}")
```

## Using the API

The enhanced RDKit environment provides the following API endpoints:

### Calculate Descriptors

```
POST /api/v1/rdkit/enhanced/descriptors

Request:
{
    "molecule": "SMILES string",
    "include_3d": true/false
}

Response:
{
    "molecule": "SMILES string",
    "include_3d": true/false,
    "descriptors": {
        "MolWt": 180.159,
        "MolLogP": 1.43,
        ...
    }
}
```

### Generate Conformers

```
POST /api/v1/rdkit/enhanced/conformers

Request:
{
    "molecule": "SMILES string",
    "n_conformers": 10
}

Response:
{
    "molecule": "SMILES string",
    "n_conformers": 10,
    "conformers": [
        {
            "id": 0,
            "energy": 25.6,
            "molblock": "...",
            "relative_energy": 0.0
        },
        ...
    ]
}
```

### Calculate Similarity

```
POST /api/v1/rdkit/enhanced/similarity

Request:
{
    "query_molecule": "SMILES string",
    "target_molecules": ["SMILES1", "SMILES2", ...],
    "similarity_mode": "3d"
}

Response:
{
    "query_molecule": "SMILES string",
    "similarity_mode": "3d",
    "results": [
        {
            "target_smiles": "SMILES1",
            "similarity_score": 0.85
        },
        ...
    ]
}
```

### Analyze Pharmacophore

```
POST /api/v1/rdkit/enhanced/pharmacophore

Request:
{
    "molecule": "SMILES string"
}

Response:
{
    "molecule": "SMILES string",
    "total_features": 5,
    "feature_counts": {
        "Donor": 1,
        "Acceptor": 2,
        "Hydrophobe": 1,
        "Aromatic": 1
    },
    "feature_percentages": {
        "Donor": 20.0,
        "Acceptor": 40.0,
        "Hydrophobe": 20.0,
        "Aromatic": 20.0
    },
    "interaction_potential": 5.5,
    "features": [
        {
            "type": "Donor",
            "atoms": [11],
            "position": [1.2, 2.3, 3.4]
        },
        ...
    ]
}
```

### Analyze Conformers

```
POST /api/v1/rdkit/enhanced/conformer-analysis

Request:
{
    "molecule": "SMILES string",
    "n_conformers": 10
}

Response:
{
    "molecule": "SMILES string",
    "n_conformers": 10,
    "energy_min": 10.0,
    "energy_max": 20.0,
    "energy_range": 10.0,
    "energy_avg": 15.0,
    "flexibility_score": 3.3,
    "stability_score": 90.9,
    "conformer_details": [...]
}
```

### Predict Properties

```
POST /api/v1/rdkit/enhanced/predict-property

Request:
{
    "molecule": "SMILES string",
    "properties": ["solubility", "glass_transition", ...],
    "use_3d": true/false
}

Response:
{
    "molecule": "SMILES string",
    "use_3d": true/false,
    "predictions": {
        "solubility": 1.5,
        "glass_transition": 125.0,
        ...
    }
}
```

### Get Service Info

```
GET /api/v1/rdkit/enhanced/info

Response:
{
    "version": "1.0.0",
    "rdkit_version": "2023.09.1",
    "descriptors_count": {
        "2d": 150,
        "3d": 20,
        "total": 170
    },
    "api": {
        "endpoints": [...],
        "description": "Enhanced RDKit service API for CryoProtect"
    }
}
```

## Advanced Usage

### Customizing the Environment

You can customize the enhanced RDKit environment by modifying the following files:

- `Dockerfile.rdkit-enhanced`: Docker container definition
- `requirements_rdkit_enhanced.txt`: Python dependencies
- `rdkit_enhanced_service.py`: Service implementation
- `scientific_models/rdkit_enhanced.py`: Python interface

### Adding New Models

To add new scientific models that use the enhanced RDKit service:

1. Create a new class that inherits from `DescriptorBasedModel` in `scientific_models/rdkit_enhanced.py`
2. Implement the `calculate` method to interact with the RDKit service
3. Add a new API resource in `api/rdkit_enhanced_resources.py`
4. Register the new resource in the `register_resources` function

### Troubleshooting

#### Service Not Starting

If the enhanced RDKit service fails to start:

```bash
# Check the container logs
docker logs cryoprotect-rdkit-enhanced

# Check if the port is already in use
sudo lsof -i :5001
```

#### Service Connection Issues

If the Python interface can't connect to the service:

```python
# Set the service URL explicitly
calculator = EnhancedRDKitCalculator(service_url="http://localhost:5001")

# Test the connection
calculator._check_connection()
```

#### Container Performance Issues

If the container is using too many resources:

```bash
# Limit CPU and memory usage
docker run --cpus=2 --memory=4g -d -p 5001:5000 --name cryoprotect-rdkit-enhanced cryoprotect-rdkit-enhanced
```

## References

- [RDKit Documentation](https://www.rdkit.org/docs/index.html)
- [Docker Documentation](https://docs.docker.com/)
- [Flask-RESTful Documentation](https://flask-restful.readthedocs.io/en/latest/)