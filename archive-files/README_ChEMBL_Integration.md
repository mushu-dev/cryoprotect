# ChEMBL Integration for CryoProtect v2

This document describes the integration of ChEMBL as an additional chemical data source for the CryoProtect v2 system, as part of TASK_2.2 - Implement Integration for Additional Chemical Data Source.

## Overview

The ChEMBL integration allows CryoProtect v2 to fetch and store chemical data from the ChEMBL database, a large bioactivity database containing information on drug-like bioactive compounds. This integration complements the existing PubChem data source, providing a more comprehensive dataset for cryoprotectant analysis.

## Components

The integration consists of the following components:

1. **ChEMBL Client** - A resilient client for interacting with the ChEMBL API
   - Located in the `chembl/` directory
   - Features caching, rate limiting, retry logic, and circuit breaking
   - Handles API errors gracefully with fallbacks to cached data

2. **ChEMBL Data Importer** - A script to import cryoprotectant data from ChEMBL into Supabase
   - Located at `ChEMBL_CryoProtectants_Supabase.py`
   - Filters molecules based on cryoprotectant criteria
   - Scores molecules based on cryoprotective properties
   - Handles deduplication based on InChIKey
   - Provides source attribution for provenance tracking

3. **Test Script** - A script to test the ChEMBL client functionality
   - Located at `test_chembl_client.py`
   - Tests basic client operations like fetching molecules and searching

## Architecture

The ChEMBL integration follows the same modular architecture as the existing PubChem integration:

```
chembl/
├── __init__.py        # Package initialization
├── cache.py           # Caching system for ChEMBL API responses
├── client.py          # Main client for ChEMBL API
├── rate_limiter.py    # Rate limiting for ChEMBL API requests
└── utils.py           # Utility functions

ChEMBL_CryoProtectants_Supabase.py  # Main data importer script
test_chembl_client.py               # Test script
```

## Data Flow

1. The ChEMBL client fetches molecule data from the ChEMBL API
2. The data is normalized to match the CryoProtect v2 schema
3. Molecules are filtered based on cryoprotectant criteria
4. Filtered molecules are scored based on cryoprotective properties
5. Molecules are checked for duplicates based on InChIKey
6. New molecules are inserted into the Supabase database with source attribution
7. Molecular properties are inserted into the database

## Features

- **Resilient API Client**
  - Adaptive rate limiting to respect ChEMBL API limits
  - Two-level caching (memory and disk) for improved performance
  - Exponential backoff retry logic for transient errors
  - Circuit breaker to prevent repeated failures
  - Fallback to cached data when API is unavailable

- **Data Processing**
  - Batch processing with checkpointing for resumable operation
  - Deduplication based on InChIKey to avoid duplicate entries
  - Source attribution for provenance tracking
  - Comprehensive logging for monitoring and debugging

- **Integration with Existing System**
  - Uses the same database schema as the PubChem integration
  - Compatible with existing filtering and scoring logic
  - Provides the same data fields for consistent analysis

## Usage

### Testing the ChEMBL Client

```bash
python test_chembl_client.py
```

### Importing Data from ChEMBL

```bash
python ChEMBL_CryoProtectants_Supabase.py --batch-size 10
```

Additional options:
- `--batch-size <size>` - Set the batch size for processing (default: 10)
- `--checkpoint <path>` - Specify the checkpoint file path (default: chembl_checkpoint.json)
- `--resume` - Resume from the last checkpoint
- `--reset` - Reset checkpoint and start from scratch
- `--skipped-log <path>` - Specify the path for the skipped molecules log

## Configuration

The ChEMBL integration uses the following environment variables:

- `SUPABASE_URL` - Supabase project URL
- `SUPABASE_KEY` - Supabase API key
- `SUPABASE_USER` - Supabase user email (optional)
- `SUPABASE_PASSWORD` - Supabase user password (optional)
- `CHEMBL_API_DELAY` - Delay between ChEMBL API requests in seconds (default: 0.3)
- `CHEMBL_ID_FILE` - Path to a file containing ChEMBL IDs to process (optional)

These can be set in a `.env` file in the project root directory.

## Deduplication Strategy

To avoid duplicate entries when importing from multiple sources, the integration:

1. Uses InChIKey as a unique identifier for molecules
2. Checks if a molecule with the same InChIKey already exists in the database
3. Skips insertion if a duplicate is found
4. Logs the skipped molecule with the reason "Duplicate"

## Source Attribution

For provenance tracking, each molecule is tagged with its data source:

- The `data_source` field is set to "ChEMBL" for molecules from ChEMBL
- The `modification_history` field includes the source script name
- The ChEMBL ID is stored as a property for reference

## Future Improvements

- Add support for more ChEMBL-specific properties
- Implement more sophisticated scoring based on ChEMBL bioactivity data
- Add support for ChEMBL's compound hierarchy and relationships
- Improve search capabilities to find more potential cryoprotectants