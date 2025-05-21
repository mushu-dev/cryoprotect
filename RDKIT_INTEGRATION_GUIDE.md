# RDKit Integration Guide for CryoProtect

This guide provides comprehensive recommendations for properly integrating RDKit into the CryoProtect project for long-term stability, maintainability, and reliability.

## Current State Analysis

Based on our investigation, we found:

1. **Host Environment:**
   - RDKit 2025.3.1 is installed via pip on the host system
   - Python 3.13.3 is used on the host system
   - RDKit is functional on the host with basic functionality

2. **Container Environment:**
   - RDKit 2022.9.5 is now installed in the Docker/Podman container
   - Container uses Python 3.9.12
   - Container allows for isolated testing with consistent environment

3. **Existing Integration Approaches:**
   - `mock_rdkit_formula.py` - For fallback when RDKit isn't available
   - `verify_rdkit.py` - Tests RDKit functionality
   - `rdkit_test_environment.sh` - Manages RDKit test environment
   - `run_with_rdkit.sh` - Wrapper to run scripts with RDKit in container

## Recommendations for Optimal RDKit Integration

### 1. Standardize on Container-Based Approach

The most reliable approach is to standardize on a container-based deployment for RDKit:

```
Host System → Container with RDKit → Production Environment
```

**Benefits:**
- Consistent environment across development, testing, and production
- Avoids "works on my machine" problems
- Easier reproducibility and deployment

### 2. Create a Layered Architecture for RDKit Integration

```
┌─────────────────────────────┐
│   High-Level API Functions  │
│                             │
│ calculate_properties(smiles)│
│ generate_svg(molecule)      │
│ etc.                        │
└───────────────┬─────────────┘
                │
                ▼
┌─────────────────────────────┐
│   Integration Layer         │
│                             │
│ try:                        │
│   import rdkit              │
│ except ImportError:         │
│   import mock_rdkit_formula │
│                             │
└───────────────┬─────────────┘
                │
                ▼
┌─────────────────┬─────────────┐
│  Real RDKit     │  Mock RDKit │
└─────────────────┴─────────────┘
```

### 3. Implementation Approach

1. **Create a unified RDKit wrapper module (`rdkit_wrapper.py`):**

```python
"""
RDKit Wrapper - Provides a unified interface to RDKit functionality
with fallback to mock implementation when needed.
"""

import logging
import os
import sys
from typing import Dict, Any, Optional, List, Tuple, Union

logger = logging.getLogger(__name__)

# Try to import RDKit
RDKIT_AVAILABLE = False
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, MolSurf, AllChem, rdMolDescriptors
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    
    RDKIT_AVAILABLE = True
    RDKIT_VERSION = getattr(Chem, '__version__', 'Unknown')
    logger.info(f"Using RDKit version {RDKIT_VERSION}")
except ImportError:
    logger.warning("RDKit not available, falling back to mock implementation")
    
    # Import mock implementation
    try:
        from mock_rdkit_formula import (
            calculate_molecular_formula,
            calculate_molecular_weight,
            calculate_logp,
            calculate_tpsa,
            calculate_h_donors,
            calculate_h_acceptors,
            calculate_rotatable_bonds,
            calculate_ring_count,
            calculate_aromatic_ring_count,
            calculate_heavy_atom_count
        )
    except ImportError:
        logger.error("Failed to import mock_rdkit_formula, RDKit functionality will be limited")

# Unified interface functions that work with both real and mock RDKit

def create_molecule_from_smiles(smiles: str) -> Optional[Any]:
    """Create a molecule object from SMILES string."""
    if not smiles:
        return None
        
    if RDKIT_AVAILABLE:
        return Chem.MolFromSmiles(smiles)
    else:
        # For mock implementation, just return the SMILES as a placeholder
        # This allows functions to work with both real and mock implementations
        return smiles

def calculate_properties(mol_or_smiles: Union[str, Any]) -> Dict[str, Any]:
    """Calculate molecular properties for a molecule or SMILES string."""
    properties = {}
    
    # Convert SMILES to molecule if needed
    if isinstance(mol_or_smiles, str):
        mol = create_molecule_from_smiles(mol_or_smiles)
        smiles = mol_or_smiles
    else:
        mol = mol_or_smiles
        smiles = Chem.MolToSmiles(mol) if RDKIT_AVAILABLE else str(mol)
    
    if not mol:
        return properties
    
    # Calculate properties using either real or mock RDKit
    if RDKIT_AVAILABLE:
        properties["molecular_formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
        properties["molecular_weight"] = round(Descriptors.MolWt(mol), 2)
        properties["logp"] = round(Descriptors.MolLogP(mol), 2)
        properties["tpsa"] = round(MolSurf.TPSA(mol), 2)
        properties["h_donors"] = Lipinski.NumHDonors(mol)
        properties["h_acceptors"] = Lipinski.NumHAcceptors(mol)
        properties["rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
        properties["ring_count"] = Chem.rdMolDescriptors.CalcNumRings(mol)
        properties["aromatic_ring_count"] = Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
        properties["heavy_atom_count"] = mol.GetNumHeavyAtoms()
    else:
        # Use mock implementation
        properties["molecular_formula"] = calculate_molecular_formula(smiles)
        properties["molecular_weight"] = calculate_molecular_weight(smiles)
        properties["logp"] = calculate_logp(smiles)
        properties["tpsa"] = calculate_tpsa(smiles)
        properties["h_donors"] = calculate_h_donors(smiles)
        properties["h_acceptors"] = calculate_h_acceptors(smiles)
        properties["rotatable_bonds"] = calculate_rotatable_bonds(smiles)
        properties["ring_count"] = calculate_ring_count(smiles)
        properties["aromatic_ring_count"] = calculate_aromatic_ring_count(smiles)
        properties["heavy_atom_count"] = calculate_heavy_atom_count(smiles)
    
    return properties

def generate_molecule_svg(mol_or_smiles: Union[str, Any], width: int = 300, height: int = 200) -> Optional[str]:
    """Generate SVG visualization of a molecule."""
    if not RDKIT_AVAILABLE:
        logger.warning("RDKit not available, SVG generation is not supported in mock mode")
        return None
    
    # Convert SMILES to molecule if needed
    if isinstance(mol_or_smiles, str):
        mol = create_molecule_from_smiles(mol_or_smiles)
    else:
        mol = mol_or_smiles
    
    if not mol:
        return None
    
    # Generate SVG
    try:
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        logger.error(f"Error generating SVG: {e}")
        return None

# Add more unified interface functions as needed
```

2. **Update all code to use this wrapper instead of direct RDKit imports:**

```python
# Before:
from rdkit import Chem
from rdkit.Chem import Descriptors

mol = Chem.MolFromSmiles(smiles)
mw = Descriptors.MolWt(mol)

# After:
from rdkit_wrapper import create_molecule_from_smiles, calculate_properties

mol = create_molecule_from_smiles(smiles)
properties = calculate_properties(mol)
mw = properties["molecular_weight"]
```

### 4. Container Configuration Improvements

1. **Create a dedicated `Dockerfile.rdkit` for RDKit-specific setup:**

```dockerfile
FROM continuumio/miniconda3:4.12.0

WORKDIR /app

# Install mamba for faster conda installs and essential tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    ca-certificates && \
    conda install -y -c conda-forge mamba && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy environment file
COPY environment.yml .

# Create the conda environment with RDKit
RUN mamba env create -f environment.yml && \
    conda clean -afy && \
    rm -rf /root/.cache

# Set up environment activation in the entrypoint
COPY docker-entrypoint-rdkit.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Set the default command
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["python", "-c", "import rdkit; print(f'RDKit {rdkit.__version__} is available')"]
```

2. **Create a dedicated entry point script (`docker-entrypoint-rdkit.sh`):**

```bash
#!/bin/bash
set -e

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate cryoprotect

# Run the specified command
exec "$@"
```

### 5. Deployment Strategy

1. **Development Mode:**
   - Local developers can use `run_with_rdkit.sh` to test specific scripts
   - CI pipeline should use the dedicated RDKit container for all tests

2. **Production Mode:**
   - Deploy using the fixed RDKit container
   - Consider a microservices approach where RDKit operations run in a separate container

3. **Continuous Testing:**
   - Regularly run `verify_rdkit.py` in CI/CD pipeline
   - Set up automated tests for the wrapper functions

### 6. Upgrade Strategy

1. **Container Version Management:**
   - Pin RDKit version in `environment.yml` (e.g., `rdkit=2024.03.4`)
   - Create a process for regular updates (quarterly or with major releases)

2. **Test Coverage for Backward Compatibility:**
   - Create tests with known expected results
   - Run these tests on each update to ensure API compatibility

3. **Version-specific Adaptations:**
   - Add version check logic for any breaking changes
   ```python
   if RDKIT_AVAILABLE and RDKIT_VERSION >= "2024.03":
       # Use newer API
   else:
       # Use compatible API
   ```

## Implementation Roadmap

1. **Phase 1: Wrapper Creation** (1-2 days)
   - Create the unified `rdkit_wrapper.py` module
   - Write tests for all wrapper functions

2. **Phase 2: Integration** (2-3 days)
   - Update existing code to use the wrapper
   - Ensure all tests pass with both real and mock implementations

3. **Phase 3: Containerization** (1-2 days)
   - Create or update the dedicated RDKit container
   - Update deployment scripts

4. **Phase 4: Documentation and Training** (1 day)
   - Update project documentation with new standards
   - Create usage examples for developers

## Conclusion

By implementing these recommendations, the CryoProtect project will have a robust, maintainable, and consistent approach to RDKit integration. This will ensure that both development and production environments can reliably perform molecular calculations, visualization, and analysis, while maintaining backward compatibility and allowing for future upgrades.

The layered architecture provides the flexibility to handle different environments (with or without RDKit) while the containerization strategy ensures consistency across all environments.