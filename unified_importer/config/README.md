# Configuration Files for Unified Molecular Importer

This directory contains configuration files for the unified molecular data importer. These configurations control various aspects of the import process, including database connections, API settings, and transformation options.

## Available Configuration Files

- `config_example.json`: Example configuration with documentation for all options
- `config_pubchem.json`: Configuration optimized for PubChem imports
- `config_chembl.json`: Configuration optimized for ChEMBL imports 
- `config_combined.json`: Configuration for importing from both sources

## Configuration Format

The configuration files use a structured JSON format with the following main sections:

### Database Configuration

Controls database connections and behavior:

```json
"database": {
    "url": "YOUR_SUPABASE_URL",
    "key": "YOUR_SUPABASE_KEY",
    "use_connection_pool": true,
    "pool_size": 10,
    "max_retries": 3,
    "retry_delay": 2.0
}
```

### Checkpoint Configuration

Controls the resumable import checkpointing system:

```json
"checkpoints": {
    "directory": "checkpoints",
    "backup_interval": 5,
    "enabled": true
}
```

### Logging Configuration

Controls logging behavior:

```json
"logging": {
    "level": "INFO",
    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    "file": null,
    "structured_output": true
}
```

### Data Sources Configuration

Controls settings for each data source:

```json
"sources": {
    "pubchem": {
        "api_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
        "api_delay": 0.5,
        "batch_size": 20,
        // Additional source-specific settings...
    },
    "chembl": {
        "api_url": "https://www.ebi.ac.uk/chembl/api/data",
        "api_delay": 0.5,
        // Additional source-specific settings...
    }
}
```

### Transform Configuration

Controls data transformation behavior:

```json
"transforms": {
    "molecule_transform": {
        "resolve_cross_references": true,
        "handle_mixtures": true
        // Additional transform settings...
    },
    "property_transform": {
        "standardize_units": true,
        "add_missing_properties": true
    }
}
```

### Import Configuration

Controls the overall import process:

```json
"import": {
    "batch_size": 50,
    "worker_count": 4,
    "max_retries": 3,
    "retry_delay": 2.0,
    "checkpoint_frequency": 10
}
```

## Usage

To use these configuration files, load them when initializing the importer:

```python
from unified_importer.main import MolecularImporter

# Load a configuration file
importer = MolecularImporter(config_file='config/config_pubchem.json')

# Or provide configuration programmatically
importer = MolecularImporter(config={
    'database': {
        'url': 'https://your-project.supabase.co',
        'key': 'your-api-key'
    },
    # Other configuration sections...
})
```

## Creating Custom Configurations

You can create custom configuration files by copying the example file and modifying it for your specific needs. Save your custom configurations in this directory with descriptive names.