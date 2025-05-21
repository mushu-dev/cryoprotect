# Property-Based PubChem Import for CryoProtect

This document provides an overview of the property-based PubChem import functionality for the CryoProtect application.

## Overview

The property-based PubChem import is an enhanced approach to identifying and importing potential cryoprotectant molecules from the PubChem database. Unlike the earlier name-based import, this approach leverages molecular properties that are characteristic of effective cryoprotectants to find a more diverse and comprehensive set of candidate molecules.

## Key Features

- **Multi-stage search approach**: Combines name-based, property-based, and similarity-based searches
- **Property filters**: Uses optimized property ranges to target different classes of cryoprotectants
- **Resilient API handling**: Implements proper error handling, rate limiting, and retry logic
- **Batch processing**: Efficiently processes compounds in configurable batches
- **Checkpointing**: Supports resumable imports through checkpoint files
- **Detailed reporting**: Generates comprehensive reports on import results

## Fixed Database Schema Issues

The fixed implementation (`property_based_pubchem_import_fixed.py`) resolves several important issues:

1. **Property Type Handling**: Properly handles the database schema requirement for non-null `data_type` values in the `property_types` table
2. **Individual Property Searches**: Uses individual property searches due to limitations in PubChem's combined property search API
3. **API Error Handling**: Improved handling of PubChem API errors and rate limits
4. **Database Connection Resilience**: Better error handling for database connections and transactions

## Running the Import

There are several ways to run the property-based PubChem import:

### Standalone Import

To run just the property-based PubChem import:

```bash
./run_property_pubchem_import_fixed.sh
```

This script sets up the necessary container environment and runs the import with reasonable default parameters.

### Full Database Population

To run a complete database population including both ChEMBL and PubChem imports:

```bash
./populate_full_database_improved.sh
```

This script:
1. Runs the ChEMBL import first
2. Then runs the property-based PubChem import
3. Generates a comprehensive report with results from both imports

## Configuration Options

The import script accepts several configuration options:

- `--target`: Number of compounds to import (default: 100)
- `--api-delay`: Base delay between API calls in seconds (default: 0.5)
- `--batch-size`: Number of compounds to process in each batch (default: 10)
- `--similarity-threshold`: Similarity threshold for finding similar compounds (default: 0.8)
- `--checkpoint`: Path to a checkpoint file to resume a previous import

## Property Filtering Approach

The property-based filtering targets three main classes of potential cryoprotectants:

1. **Small polar molecules with H-bonding capability**:
   - Molecular Weight: 32-200 Da
   - XLogP: -4 to 3
   - At least 1 H-bond donor
   - At least 2 H-bond acceptors

2. **Sugars and polyols**:
   - Molecular Weight: 90-500 Da
   - XLogP: -4 to 0 (very hydrophilic)
   - Multiple H-bond donors (3+)
   - Multiple H-bond acceptors (4+)

3. **Amino acids and derivatives**:
   - Molecular Weight: 70-200 Da
   - XLogP: -3 to 1
   - At least 2 H-bond donors
   - At least 3 H-bond acceptors

## Database Schema Integration

The implementation properly integrates with the CryoProtect database schema:

- **Molecule Table**: Stores basic information about each compound
- **Property Types Table**: Defines the types of molecular properties
- **Molecular Properties Table**: Stores the property values for each molecule

## Property Types

The script handles the following property types:

- Molecular Weight
- XLogP3
- Hydrogen Bond Donor Count
- Hydrogen Bond Acceptor Count
- Rotatable Bond Count
- Topological Polar Surface Area
- Heavy Atom Count
- Complexity

Each property type is created with the appropriate data_type, description, and units.

## Future Enhancements

Potential future enhancements include:

1. **Refined Property Filters**: Fine-tuning property ranges based on effectiveness of discovered compounds
2. **Machine Learning Integration**: Adding ML models to evaluate and score candidate compounds
3. **Structural Classification**: Adding substructure analysis to identify chemical classes
4. **Experimental Data Integration**: Linking to experimental cryoprotectant effectiveness data

## Troubleshooting

If you encounter issues:

1. **API Rate Limiting**: The script implements exponential backoff, but you might need to increase the --api-delay parameter
2. **Database Connection**: Ensure the database credentials and connection parameters are correct
3. **Container Issues**: If the container fails, try running `podman logs cryoprotect-rdkit-minimal` to see the error messages