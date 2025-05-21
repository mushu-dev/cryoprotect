# Unified Molecular Data Import Plan

This document outlines a plan for unifying the ChEMBL and PubChem data import processes into a single streamlined system.

## Current State

Currently, there are multiple import scripts with overlapping functionality:

- `ChEMBL_Integrated_Import.py`: Imports from ChEMBL with structured logging, checkpointing, and direct Supabase connections
- `import_pubchem_data_direct.py`: Imports from PubChem with parallel processing, batch operations, and direct database connections
- Additional variants in the codebase with different optimization strategies

## Proposed Architecture

A modular, component-based architecture with a unified core:

```
UnifiedMolecularImporter/
├── core/
│   ├── config.py             # Unified configuration system
│   ├── logging.py            # Enhanced logging with structured JSON logs
│   ├── checkpoint.py         # Robust checkpoint system
│   ├── database.py           # Database connection and operations
│   ├── progress.py           # Progress tracking and reporting
│   └── validation.py         # Molecule validation and filtering
├── sources/
│   ├── chembl_source.py      # ChEMBL specific import logic
│   ├── pubchem_source.py     # PubChem specific import logic
│   └── source_base.py        # Base class for data sources
├── transforms/
│   ├── molecule_transform.py # Transform to common molecule format
│   └── property_transform.py # Transform to common property format
└── main.py                   # Main entry point
```

## Key Features

### 1. Modular Data Sources

- Implement a base `MolecularDataSource` class
- Source-specific implementations for ChEMBL and PubChem
- Common interface for fetching compounds, validation, and transformations

### 2. Enhanced Error Handling

- Progressive backoff strategy for API requests
- Transaction-based database operations with retry logic
- Error categorization and appropriate recovery strategies

### 3. Robust Checkpointing

- Regular checkpoints at configurable intervals
- Automatically resume from last successful batch
- Store detailed state for accurate resumption
- Backup checkpoints to prevent data loss

### 4. Parallel Processing

- Configurable worker pool with thread and async hybrid approach
- Semaphore-based concurrency control
- Asynchronous API requests with batched database operations

### 5. Progress Tracking

- Real-time progress metrics and ETA calculations
- Separate tracking for different error types (API vs. DB)
- Comprehensive completion reports
- Detailed logging of skipped compounds with reasons

## Implementation Plan

### Phase 1: Core Framework (1-2 days)

1. Create the directory structure
2. Implement base classes and interfaces
3. Set up the configuration and logging systems
4. Develop the checkpoint mechanism

### Phase 2: First Data Source (2-3 days)

1. Implement the ChEMBL data source
2. Create molecule and property transformations
3. Implement database operations with transaction support
4. Test with a small subset of data

### Phase 3: Second Data Source (1-2 days)

1. Implement the PubChem data source
2. Adapt transformations for PubChem data
3. Test with a small subset of data

### Phase 4: Refinement (2-3 days)

1. Optimize performance parameters
2. Add comprehensive error handling
3. Improve reporting and monitoring
4. Create detailed documentation

## Command-Line Interface

The unified importer will support:

```
python unified_import.py [options]

Options:
  --source SOURCE         Data source (chembl, pubchem, all)
  --limit LIMIT           Maximum compounds to import
  --batch-size BATCH_SIZE Batch size for processing
  --workers WORKERS       Number of worker threads
  --api-delay API_DELAY   Delay between API calls
  --resume                Resume from last checkpoint
  --dry-run               Simulate without inserting data
```

## Technical Details

### Data Transformation

Both ChEMBL and PubChem data will be transformed into a common molecule format:

```python
{
    "name": "Compound Name",
    "chembl_id": "CHEMBLXXXX",  # For ChEMBL compounds
    "pubchem_cid": "XXXXX",     # For PubChem compounds
    "formula": "C6H12O6",
    "molecular_weight": 180.156,
    "smiles": "C1C(C(C(C(C1O)O)O)O)O",
    "inchi": "InChI=1S/C6H12O6/...",
    "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
    "data_source": "ChEMBL" or "PubChem"
}
```

Properties will be transformed to a standardized format with typed values:

```python
{
    "molecule_id": "molecule-uuid",
    "property_name": "LogP",
    "property_type": "physicochemical",
    "numeric_value": 1.45,     # Type-specific value fields
    "text_value": null,
    "boolean_value": null,
    "unit": "log units",
    "source": "ChEMBL" or "PubChem"
}
```

## Benefits of Unified Approach

1. **Code Reusability**: Common operations consolidated into base classes
2. **Source Extensibility**: New sources can be added by implementing the interface
3. **Error Resilience**: Improved handling at both API and database levels
4. **Progress Tracking**: Detailed statistics reported in real-time
5. **Optimized Database Operations**: Batch operations with transaction support

## Future Enhancements

1. Add support for additional data sources
2. Implement more advanced filtering and validation
3. Add data enrichment from secondary sources
4. Create a web dashboard for monitoring imports
5. Support incremental updates of existing compounds