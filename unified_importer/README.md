# Unified Molecular Data Importer

A comprehensive framework for importing molecular data from multiple sources (PubChem, ChEMBL) into the CryoProtect database.

## Features

- **Multi-Source Support**: Import from PubChem, ChEMBL, and expandable to other sources
- **Advanced Molecule Transformation**: Standardization, property calculation, cross-reference resolution
- **Efficient Batch Processing**: Parallel async operations with rate limiting and connection pooling
- **Enhanced Connection Pool**: Health checking, dynamic sizing, circuit breaker pattern, and detailed metrics
- **Adaptive Batch Processing**: Dynamic batch sizing, parallel processing, and memory-aware execution
- **Comprehensive Caching System**: In-memory and disk caching with smart policies and automatic expiration
- **PubChem Property Filters**: Configurable filtering of compounds based on properties and search terms
- **Checkpointing**: Resume interrupted imports from where they left off
- **Progress Tracking**: Detailed statistics and progress reporting
- **Mixture Handling**: Automatic detection and processing of molecular mixtures
- **Comprehensive Test Suite**: Unit, integration, and functional tests

## Directory Structure

```
unified_importer/
├── core/                  # Core functionality
│   ├── config.py          # Configuration management
│   ├── database.py        # Database operations with connection pooling
│   ├── checkpoint.py      # Checkpoint system for resumable imports
│   ├── progress.py        # Progress tracking and reporting
│   └── validation.py      # Molecule validation
├── sources/               # Data source implementations
│   ├── source_base.py     # Base class for data sources
│   ├── chembl_source.py   # ChEMBL-specific implementation
│   └── pubchem_source.py  # PubChem-specific implementation
├── transforms/            # Data transformation utilities
│   ├── molecule_transform.py # Molecule standardization and conversion
│   └── property_transform.py # Property standardization
├── migrations/            # Database migrations
├── config/                # Configuration files
├── examples/              # Usage examples
├── tests/                 # Unit and integration tests
└── docs/                  # Documentation
```

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/cryoprotect.git
cd cryoprotect

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

### 1. Apply Database Migration

```bash
# Apply the required database migration
python -m unified_importer.migrations.apply_migration
```

### 2. Basic Import Usage

```python
import asyncio
from unified_importer.sources.chembl_source import ChEMBLDataSource
from unified_importer.sources.pubchem_source import PubChemDataSource
from unified_importer.core.database import DatabaseOperations

async def import_compounds():
    # Initialize database operations
    db = DatabaseOperations('your_supabase_url', 'your_api_key')
    
    # Create ChEMBL data source
    chembl_source = ChEMBLDataSource(db)
    
    # Import specific ChEMBL compounds
    results = await chembl_source.import_compounds(['CHEMBL25', 'CHEMBL1201'])
    print(f"Imported {results[0]} ChEMBL compounds")
    
    # Create PubChem data source
    pubchem_source = PubChemDataSource(db)
    
    # Search and import PubChem compounds
    results = await pubchem_source.search_and_import('glycerol', max_results=10)
    print(f"Imported {results[0]} PubChem compounds")

# Run the async function
asyncio.run(import_compounds())
```

### 3. Using Configuration Files

```python
from unified_importer.main import MolecularImporter

# Create importer with configuration file
importer = MolecularImporter(config_file='unified_importer/config/config_combined.json')

# Run import
importer.run_import(
    sources=['chembl', 'pubchem'],
    query='cryoprotectant',
    limit=100
)
```

## Advanced Usage

### Transformation Example

```python
from unified_importer.transforms.molecule_transform import MoleculeTransformer

# Create transformer
transformer = MoleculeTransformer()

# Process a molecule
molecule_data = {'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O', 'name': 'Aspirin'}
result = await transformer.standardize_molecule(molecule_data)

# Print calculated properties
for prop_name, value in result['properties'].items():
    print(f"{prop_name}: {value}")
```

### Progress Tracking

```python
from unified_importer.core.progress import ProgressTracker, ConsoleProgressReporter

# Create tracker and reporter
tracker = ProgressTracker(total_items=1000)
reporter = ConsoleProgressReporter(tracker, update_interval=5.0)

# Use in import process
for i in range(100):
    # Process items...
    tracker.update(processed=10, successful=8, failed=2)
    # Sleep to simulate work...
    
# Generate report
report = tracker.generate_report()
tracker.save_report('import_report.json')
```

### Checkpointing

```python
from unified_importer.core.checkpoint import CheckpointManager

# Create checkpoint manager
checkpoint = CheckpointManager('checkpoint.json')

# Use in import process
for item_id in item_ids:
    if checkpoint.is_processed(item_id):
        continue  # Skip already processed items
        
    # Process item...
    checkpoint.mark_processed(item_id)
    
    # Save checkpoint every 10 items
    if i % 10 == 0:
        checkpoint.save()
```

## Configuration

See the [Configuration Guide](config/README.md) for detailed configuration options and examples.

## Documentation

- [Molecule Transformation](docs/molecule_transform.md)
- [Property Standardization](docs/property_transform.md)
- [ChEMBL Source](docs/chembl_source.md)
- [PubChem Source](docs/pubchem_source.md)
- [Property Filters](docs/property_filters.md)
- [Enhanced Connection Pool](docs/enhanced_connection_pool.md)
- [Enhanced Batch Processing](docs/enhanced_batch_processing.md)
- [Caching Mechanism](docs/caching_mechanism.md)
- [Checkpoint System](docs/checkpoint.md)
- [Progress Tracking](docs/progress.md)
- [Phase 4 Progress Report](docs/PHASE4_PROGRESS_REPORT.md)
- [Batch Processing Report](docs/PHASE4_BATCH_PROCESSING_REPORT.md)
- [Caching Report](docs/PHASE4_CACHING_REPORT.md)

## Testing

```bash
# Run unit tests
pytest unified_importer/tests/

# Run specific test module
pytest unified_importer/tests/test_molecule_transform.py

# Run integration tests
pytest unified_importer/tests/test_checkpoint_integration.py
```

## License

MIT