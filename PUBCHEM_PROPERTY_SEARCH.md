# PubChem Property-Based Search

This document describes the property-based approach to searching for cryoprotectant molecules in PubChem, similar to what we've implemented for ChEMBL data import.

## Overview

Cryoprotectants typically exhibit certain molecular properties that make them effective for their purpose, such as appropriate molecular weight, polarity, and hydrogen bonding capabilities. Our property-based search approach leverages these characteristics to identify potential cryoprotectants in the PubChem database.

## Property Filters

We use several sets of property filters to capture different classes of potential cryoprotectants:

### Filter Set 1: Small Polar Molecules with H-bonding Capability

These properties target small molecules with moderate hydrophilicity and good hydrogen bonding capability:

```python
{
    "XLogP": "-3:3",          # Range from -3 to 3 (relatively hydrophilic)
    "TPSA": "40:200",         # Topological Polar Surface Area: moderate to high polarity
    "MolecularWeight": "32:500", # Small to medium molecules
    "HBondDonorCount": "1:",  # At least 1 H-bond donor
    "HBondAcceptorCount": "2:"  # At least 2 H-bond acceptors
}
```

### Filter Set 2: Focus on Sugars and Polyols

These properties target sugar-like molecules and polyols, which are common cryoprotectants:

```python
{
    "XLogP": "-4:0",          # Very hydrophilic
    "TPSA": "60:200",         # High polarity
    "MolecularWeight": "90:500", # Typical sugar/polyol size
    "HBondDonorCount": "3:",  # Multiple H-bond donors
    "HBondAcceptorCount": "4:"  # Multiple H-bond acceptors
}
```

### Filter Set 3: Amino Acids and Derivatives

These properties target amino acids and similar compounds, which can function as cryoprotectants:

```python
{
    "XLogP": "-3:1",          # Hydrophilic to slightly hydrophobic
    "MolecularWeight": "70:200", # Amino acid range
    "HBondDonorCount": "2:",  # At least 2 H-bond donors (NH2, COOH)
    "HBondAcceptorCount": "3:"  # At least 3 H-bond acceptors
}
```

## Multi-Stage Search Approach

Our property-based PubChem search follows a multi-stage approach:

1. **Seed Collection**: Collect known cryoprotectants by name as seed compounds
2. **Property Filtering**: Apply multiple property filter sets to find compounds with similar properties
3. **Similarity Expansion**: Find compounds structurally similar to the seed compounds
4. **Data Retrieval**: Fetch comprehensive data for all identified compounds
5. **Database Storage**: Store the compounds and their properties in the database

## Implementation

The implementation is in the `property_based_pubchem_import.py` script, which:

1. Uses the PubChem API for property-based searches
2. Handles property type creation and mapping in the database
3. Implements batch processing for efficient data import
4. Provides checkpointing for resumable imports
5. Generates detailed reports on the import process

## Property Types in Database

A key improvement in this implementation is proper handling of property types in the database:

1. The script checks for existing property types in the database
2. If property types don't exist, it creates them
3. It uses the property type IDs when inserting properties, resolving the issue we encountered with the earlier PubChem import

## Running the Import

To run the property-based PubChem import:

```bash
./run_property_pubchem_import.sh
```

Or to run a full database population including both ChEMBL and property-based PubChem imports:

```bash
./populate_full_database_updated.sh
```

## Results

The property-based approach significantly improves the quality and diversity of cryoprotectant candidates identified from PubChem compared to the name-based approach alone. By focusing on molecular properties characteristic of cryoprotectants, we capture compounds that might not be explicitly labeled as cryoprotectants but have the potential to function as such.

## Future Enhancements

Possible enhancements to the property-based search include:

1. Fine-tuning property ranges based on effectiveness of discovered compounds
2. Adding machine learning models trained on known cryoprotectants to evaluate candidates
3. Clustering compounds to identify new chemical classes with cryoprotectant potential
4. Integrating experimental data to validate and refine the property filters