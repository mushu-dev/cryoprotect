# Batch Processor for Molecules with None Names

This guide explains the functionality and usage of the batch processor for molecules with None names.

## Overview

During data import and consolidation, some molecules may end up with NULL values in their `name` field. The batch processor identifies these molecules and assigns appropriate names using various strategies:

1. Structural information (SMILES)
2. External identifiers (PubChem CID, ChEMBL ID)
3. Chemical properties (formula, molecular weight)
4. Fallback to generated unique names

The processor works in batches to efficiently handle large datasets and includes checkpointing to allow resuming interrupted operations.

## Usage

### Basic Usage

```bash
./batch_process_none_names.py
```

This will process all molecules with None names using the default batch size (100).

### Command Line Options

The script supports several command line options:

```bash
./batch_process_none_names.py [--batch-size=100] [--dry-run] [--checkpoint-file=filename.json] [--load-checkpoint=previous.json]
```

- `--batch-size`: Number of molecules to process in each batch (default: 100, max: 1000)
- `--dry-run`: Run in dry-run mode without making any changes to the database
- `--checkpoint-file`: Specify a custom filename for checkpoint data
- `--load-checkpoint`: Resume processing from a previous checkpoint file

### Examples

Process with a larger batch size:
```bash
./batch_process_none_names.py --batch-size=500
```

Run in dry-run mode to see what would be changed:
```bash
./batch_process_none_names.py --dry-run
```

Resume from a previous checkpoint:
```bash
./batch_process_none_names.py --load-checkpoint=none_name_process_checkpoint_20250513_120000.json
```

## Name Generation Strategies

The processor uses multiple strategies to generate meaningful names:

### 1. SMILES-based Name Generation

If the molecule has a valid SMILES string, the processor attempts to use RDKit (if available) to generate a meaningful name based on the molecular structure. If RDKit is not available, it falls back to using a simplified version of the SMILES string.

### 2. External ID-based Name Generation

If the molecule has identifiers from external databases, the processor uses these to generate a name:
- PubChem molecules: "PubChem-{CID}"
- ChEMBL molecules: "ChEMBL-{ID}"

### 3. Formula-based Name Generation

If the molecule has a chemical formula, the processor generates a name based on this:
- With molecular weight: "{Formula} (MW: {weight})"
- Without molecular weight: "Formula-{Formula}"

### 4. Fallback to Generated Names

If none of the above strategies yield a name, the processor generates a unique name using a UUID.

## Checkpointing System

The batch processor includes a robust checkpointing system that records:
- Last molecule ID processed
- Total number of molecules processed
- Number of successful updates
- Timestamp of the checkpoint

This allows operations to be resumed if they are interrupted, which is particularly important for large datasets.

Checkpoint files are JSON-formatted and named with a timestamp by default:
```
none_name_process_checkpoint_YYYYMMDD_HHMMSS.json
```

## Logging

The processor logs detailed information about its operation to both the console and a log file:
```
none_name_batch_process_YYYYMMDD_HHMMSS.log
```

The log includes:
- Initial count of molecules with None names
- Progress updates for each batch
- Error messages for any failures
- Final summary of results

## Testing

The processor includes a comprehensive test suite in `tests/test_batch_process_none_names.py`. Run the tests with:

```bash
python -m unittest tests/test_batch_process_none_names.py
```

The tests verify:
- Molecule retrieval functions
- Name generation strategies
- Database update operations
- Checkpoint functionality
- Overall batch processing flow

## Verification

A verification script is provided to check that all components of the batch processor have been correctly implemented:

```bash
./verify_batch_processor.py
```

This script checks for the existence of required files, functions, and test methods to ensure complete implementation.

## Implementation Details

The main components of the batch processor are:

1. **Molecule Retrieval**: `get_molecules_with_none_names` retrieves molecules with NULL names in batches
2. **Name Generation**: Various functions to generate names based on available data
3. **Database Update**: `update_molecule_name` updates the molecule name in the database
4. **Batch Processing**: `process_molecules` handles the overall flow of retrieving, naming, and updating molecules
5. **Checkpointing**: `save_checkpoint` and `load_checkpoint` manage the checkpoint data

## Performance Considerations

- The processor uses batching to avoid overwhelming the database
- A small delay between batches prevents excessive load
- Database operations are performed individually to minimize transaction size
- Progress logging provides visibility into long-running operations

## Conclusion

The batch processor provides an efficient, resumable solution for resolving None names in the molecule database, enhancing data quality and usability in the CryoProtect system.