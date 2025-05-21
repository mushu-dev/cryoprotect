# Unified Molecular Importer Configuration

This document explains how to configure the Unified Molecular Importer system for importing chemical data from various sources.

## Configuration File Format

The configuration is stored in JSON format and can be provided in several ways:

1. As a configuration file path when initializing the `MolecularImporter` class
2. As a dictionary when initializing the `MolecularImporter` class
3. Through environment variables with the prefix `CRYOPROTECT_IMPORT_`

## Core Configuration Sections

### Database Configuration

Controls the database connection and pooling settings.

```json
"database": {
    "url": "https://your-project.supabase.co",
    "key": "your-supabase-key",
    "use_connection_pool": true,
    "pool_size": 10,
    "pool_min_size": 2,
    "pool_max_size": 20,
    "max_retries": 3,
    "retry_delay": 2.0,
    "batch_size": 100
}
```

### Data Sources Configuration

Controls the data sources to import from and their specific settings.

```json
"sources": {
    "pubchem": {
        "api_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
        "api_delay": 0.5,
        "max_compounds_per_request": 50,
        "timeout": 30,
        "batch_size": 20,
        "property_filters": [
            {
                "name": "cryoprotectants", 
                "description": "Known cryoprotectants and related compounds",
                "terms": ["glycerol", "dmso", "ethylene glycol", "trehalose"]
            },
            {
                "name": "small_molecules",
                "description": "Small molecules suitable for cryoprotection",
                "molecular_weight_max": 200,
                "logp_max": 1.0,
                "rotatable_bonds_max": 5
            }
        ]
    },
    "chembl": {
        "api_url": "https://www.ebi.ac.uk/chembl/api/data",
        "api_delay": 0.5,
        "timeout": 30,
        "batch_size": 20,
        "use_client": true,
        "property_filters": [
            {
                "full_mwt__lte": 200.0,
                "hba__gte": 2,
                "hbd__gte": 1
            }
        ],
        "reference_compounds": [
            "CHEMBL1234",
            "CHEMBL230130"
        ]
    }
}
```

### Transforms Configuration

Controls how molecules are processed and standardized during import.

```json
"transforms": {
    "molecule_transform": {
        "resolve_cross_references": true,
        "handle_mixtures": true,
        "separate_mixtures": false,
        "add_standardized_structure": true
    },
    "property_transform": {
        "standardize_units": true,
        "add_missing_properties": true
    }
}
```

### Checkpoints Configuration

Controls the checkpointing system for resumable imports.

```json
"checkpoints": {
    "directory": "checkpoints",
    "backup_interval": 5,
    "enabled": true
}
```

### Logging Configuration

Controls the logging behavior.

```json
"logging": {
    "level": "INFO",
    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    "file": "import.log",
    "structured_output": true
}
```

### Import Configuration

Controls general import behavior.

```json
"import": {
    "batch_size": 50,
    "worker_count": 4, 
    "max_retries": 3,
    "retry_delay": 2.0,
    "checkpoint_frequency": 10
}
```

## Sample Configurations

We provide sample configuration files in the `config` directory:

- `config_example.json`: Example configuration with all options documented
- `config_pubchem.json`: Configuration optimized for PubChem imports
- `config_chembl.json`: Configuration optimized for ChEMBL imports
- `config_combined.json`: Configuration for importing from both sources

## Command-Line Arguments

When running the importer from the command line, the following arguments can be used to override configuration settings:

```
usage: python -m unified_importer.main [options]

options:
  --config CONFIG         Path to configuration file
  --sources SOURCES       Comma-separated list of sources (pubchem,chembl)
  --query QUERY           Search query
  --identifiers IDS       Comma-separated list of identifiers to import
  --limit LIMIT           Maximum number of results per source
  --checkpoint-dir DIR    Directory for checkpoint files
  --db-url URL            Database URL
  --db-key KEY            Database API key
  --log-level LEVEL       Logging level (DEBUG, INFO, WARNING, ERROR)
  --log-file FILE         Log file path
```

## Environment Variables

Configuration can also be provided through environment variables with the prefix `CRYOPROTECT_IMPORT_`.
For example, to set the batch size, you could use:

```
export CRYOPROTECT_IMPORT_BATCH_SIZE=100
```

## Complete Example

Here is a complete configuration example:

```json
{
    "database": {
        "url": "https://your-project.supabase.co",
        "key": "your-supabase-key",
        "use_connection_pool": true,
        "pool_size": 10,
        "max_retries": 3,
        "retry_delay": 2.0
    },
    "checkpoints": {
        "directory": "checkpoints",
        "backup_interval": 5,
        "enabled": true
    },
    "logging": {
        "level": "INFO",
        "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        "file": "import.log",
        "structured_output": true
    },
    "sources": {
        "pubchem": {
            "api_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
            "api_delay": 0.5,
            "max_compounds_per_request": 50,
            "timeout": 30,
            "batch_size": 20,
            "property_filters": [
                {
                    "name": "cryoprotectants", 
                    "description": "Known cryoprotectants and related compounds",
                    "terms": ["glycerol", "dmso", "ethylene glycol", "trehalose"]
                }
            ]
        },
        "chembl": {
            "api_url": "https://www.ebi.ac.uk/chembl/api/data",
            "api_delay": 0.5,
            "timeout": 30,
            "batch_size": 20,
            "use_client": true,
            "property_filters": [
                {
                    "full_mwt__lte": 200.0,
                    "hba__gte": 2,
                    "hbd__gte": 1
                }
            ]
        }
    },
    "transforms": {
        "molecule_transform": {
            "resolve_cross_references": true,
            "handle_mixtures": true,
            "separate_mixtures": false,
            "add_standardized_structure": true
        },
        "property_transform": {
            "standardize_units": true,
            "add_missing_properties": true
        }
    },
    "import": {
        "batch_size": 50,
        "worker_count": 4, 
        "max_retries": 3,
        "retry_delay": 2.0,
        "checkpoint_frequency": 10
    }
}
```