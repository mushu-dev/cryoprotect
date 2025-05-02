# Task: Create Database Operations Module

## Objective
Consolidate multiple database population scripts into a single, modular database operations module.

## Context
The codebase currently contains 12+ separate database population scripts with redundant functionality. These scripts should be consolidated into a modular system to reduce complexity and improve maintainability.

## Acceptance Criteria
- A structured `database` package is created with appropriate submodules
- All functionality from existing population scripts is preserved
- A unified CLI interface for database operations is implemented
- Tests verify that the consolidated module performs the same functions
- Documentation is updated to reflect the new structure

## Implementation Steps

1. Create the modular database directory structure:
   ```bash
   mkdir -p database/population
   mkdir -p database/migrations
   mkdir -p database/models
   mkdir -p database/utils
   mkdir -p database/verification
   ```

2. Create the package initialization files:
   ```bash
   touch database/__init__.py
   touch database/population/__init__.py
   touch database/migrations/__init__.py
   touch database/models/__init__.py
   touch database/utils/__init__.py
   touch database/verification/__init__.py
   ```

3. Create the core population modules:

   First, create a utilities module for shared functionality:
   ```python
   # database/utils/connection.py
   
   import os
   from supabase import create_client
   
   def get_supabase_client():
       """Get a Supabase client for database operations."""
       url = os.environ.get("SUPABASE_URL")
       key = os.environ.get("SUPABASE_KEY")
       
       if not url or not key:
           raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in environment variables")
       
       return create_client(url, key)
   
   def execute_sql(sql, params=None):
       """Execute SQL directly on the database."""
       client = get_supabase_client()
       return client.rpc("execute_sql", {"sql": sql, "params": params}).execute()
   ```

4. Create the population modules for different data types:

   ```python
   # database/population/molecules.py
   
   """
   Module for populating molecule data in the database.
   
   This module consolidates functionality from:
   - populate_molecules.py
   - populate_molecules_production.py
   - populate_molecules_remediation.py
   """
   
   import os
   import json
   import logging
   from database.utils.connection import get_supabase_client, execute_sql
   
   logger = logging.getLogger(__name__)
   
   def populate(config=None):
       """
       Populate the molecules table from data sources.
       
       Args:
           config: Configuration dictionary with options
       
       Returns:
           dict: Results of the population operation
       """
       if config is None:
           config = {}
       
       # Get data source - allow configuration to override
       data_source = config.get('data_source', os.environ.get('MOLECULE_DATA_SOURCE', 'default'))
       
       # Load molecule data based on source
       if data_source == 'pubchem':
           molecule_data = _load_from_pubchem(config)
       elif data_source == 'json':
           molecule_data = _load_from_json(config)
       elif data_source == 'production':
           molecule_data = _load_production_data(config)
       else:
           molecule_data = _load_default_data()
       
       # Insert molecules into database
       return _insert_molecules(molecule_data, config)
   
   def _load_from_pubchem(config):
       """Load molecule data from PubChem."""
       # Implementation from populate_molecules.py
       logger.info("Loading molecule data from PubChem")
       # ...
   
   def _load_from_json(config):
       """Load molecule data from JSON file."""
       # Implementation from populate_molecules.py
       json_path = config.get('json_path', 'data/molecules.json')
       logger.info(f"Loading molecule data from {json_path}")
       # ...
   
   def _load_production_data(config):
       """Load production molecule data."""
       # Implementation from populate_molecules_production.py
       logger.info("Loading production molecule data")
       # ...
   
   def _load_default_data():
       """Load default molecule data for development."""
       # Basic development data
       logger.info("Loading default molecule data")
       # ...
   
   def _insert_molecules(molecule_data, config):
       """Insert molecules into the database."""
       # Implementation from various populate scripts
       supabase = get_supabase_client()
       logger.info(f"Inserting {len(molecule_data)} molecules into database")
       
       results = {
           'inserted': 0,
           'updated': 0,
           'skipped': 0,
           'errors': 0
       }
       
       for molecule in molecule_data:
           try:
               # Check if molecule exists
               response = supabase.from_("molecules").select("id").eq("inchikey", molecule['inchikey']).execute()
               
               if response.data:
                   # Update existing molecule if configured to do so
                   if config.get('update_existing', False):
                       supabase.from_("molecules").update(molecule).eq("id", response.data[0]['id']).execute()
                       results['updated'] += 1
                   else:
                       results['skipped'] += 1
               else:
                   # Insert new molecule
                   supabase.from_("molecules").insert(molecule).execute()
                   results['inserted'] += 1
           except Exception as e:
               logger.error(f"Error inserting molecule {molecule.get('name', 'unknown')}: {str(e)}")
               results['errors'] += 1
       
       logger.info(f"Molecule population complete: {results}")
       return results
   ```

5. Create similar modules for other data types:
   - `database/population/mixtures.py`
   - `database/population/experiments.py`
   - `database/population/predictions.py`

6. Create a unified runner module:

   ```python
   # database/population/runner.py
   
   """
   Main entry point for database population operations.
   
   This module provides a unified interface for populating 
   different types of data in the database.
   """
   
   import argparse
   import logging
   import json
   from pathlib import Path
   
   from database.population import (
       molecules, mixtures, experiments, predictions
   )
   
   logger = logging.getLogger(__name__)
   
   def populate_all(config=None):
       """
       Populate all database tables.
       
       Args:
           config: Configuration dictionary
       
       Returns:
           dict: Results of the population operation
       """
       if config is None:
           config = {}
       
       results = {}
       
       logger.info("Starting full database population")
       
       # Populate in dependency order
       results['molecules'] = molecules.populate(config)
       results['mixtures'] = mixtures.populate(config)
       results['experiments'] = experiments.populate(config)
       results['predictions'] = predictions.populate(config)
       
       logger.info("Database population complete")
       return results
   
   def run_from_cli():
       """Run database population from command line."""
       parser = argparse.ArgumentParser(description="Database Population Tool")
       parser.add_argument('--all', action='store_true', help='Populate all tables')
       parser.add_argument('--molecules', action='store_true', help='Populate molecules table')
       parser.add_argument('--mixtures', action='store_true', help='Populate mixtures table')
       parser.add_argument('--experiments', action='store_true', help='Populate experiments table')
       parser.add_argument('--predictions', action='store_true', help='Populate predictions table')
       parser.add_argument('--config', type=str, help='Configuration file path')
       parser.add_argument('--env', type=str, default='development', help='Environment (development, production)')
       parser.add_argument('--report', action='store_true', help='Generate report')
       
       args = parser.parse_args()
       
       # Load configuration
       config = {}
       if args.config:
           with open(args.config, 'r') as f:
               config = json.load(f)
       
       # Add environment to config
       config['environment'] = args.env
       
       results = {}
       
       # Determine what to populate
       if args.all:
           results = populate_all(config)
       else:
           if args.molecules:
               results['molecules'] = molecules.populate(config)
           if args.mixtures:
               results['mixtures'] = mixtures.populate(config)
           if args.experiments:
               results['experiments'] = experiments.populate(config)
           if args.predictions:
               results['predictions'] = predictions.populate(config)
       
       # Generate report if requested
       if args.report:
           report_path = Path(f"reports/database/population_{args.env}.json")
           report_path.parent.mkdir(parents=True, exist_ok=True)
           with open(report_path, 'w') as f:
               json.dump(results, f, indent=2)
           logger.info(f"Report generated at {report_path}")
       
       return 0
   
   if __name__ == "__main__":
       # Configure logging
       logging.basicConfig(
           level=logging.INFO,
           format="%(asctime)s [%(levelname)s] %(message)s"
       )
       run_from_cli()
   ```

7. Create a main module script at the project root level:

   ```python
   # database_ops.py
   
   #!/usr/bin/env python3
   """
   Database Operations CLI
   
   This script provides a command-line interface for database operations.
   """
   
   import sys
   from database.population.runner import run_from_cli
   
   if __name__ == "__main__":
       sys.exit(run_from_cli())
   ```

8. Make the script executable:
   ```bash
   chmod +x database_ops.py
   ```

9. Create a simple test for the new module:
   ```python
   # tests/database/test_population.py
   
   """Tests for the database population module."""
   
   import unittest
   from unittest.mock import patch, MagicMock
   
   from database.population import molecules
   from database.population.runner import populate_all
   
   class TestDatabasePopulation(unittest.TestCase):
       """Test the database population functionality."""
       
       @patch('database.utils.connection.get_supabase_client')
       def test_molecule_population(self, mock_get_client):
           """Test populating molecules."""
           # Mock Supabase client
           mock_client = MagicMock()
           mock_client.from_.return_value.select.return_value.eq.return_value.execute.return_value.data = []
           mock_client.from_.return_value.insert.return_value.execute.return_value = None
           mock_get_client.return_value = mock_client
           
           # Call the function
           results = molecules.populate({
               'data_source': 'json',
               'json_path': 'tests/data/test_molecules.json'
           })
           
           # Assertions
           self.assertIn('inserted', results)
           self.assertIn('errors', results)
           
       @patch('database.population.molecules.populate')
       @patch('database.population.mixtures.populate')
       @patch('database.population.experiments.populate')
       @patch('database.population.predictions.populate')
       def test_populate_all(self, mock_predictions, mock_experiments, mock_mixtures, mock_molecules):
           """Test populating all tables."""
           # Set up mocks
           mock_molecules.return_value = {'inserted': 10}
           mock_mixtures.return_value = {'inserted': 5}
           mock_experiments.return_value = {'inserted': 2}
           mock_predictions.return_value = {'inserted': 3}
           
           # Call the function
           results = populate_all()
           
           # Assertions
           self.assertIn('molecules', results)
           self.assertIn('mixtures', results)
           self.assertIn('experiments', results)
           self.assertIn('predictions', results)
           
           # Verify the mocks were called
           mock_molecules.assert_called_once()
           mock_mixtures.assert_called_once()
           mock_experiments.assert_called_once()
           mock_predictions.assert_called_once()
   ```

10. Update the documentation with the new module:
    ```markdown
    # Database Operations
    
    This module provides utilities for database operations, including:
    
    - Population of initial data
    - Verification of database structure
    - Database migrations
    
    ## Usage
    
    To populate the database with initial data:
    
    ```bash
    # Populate all tables
    python database_ops.py --all
    
    # Populate specific tables
    python database_ops.py --molecules --mixtures
    
    # Use production data source
    python database_ops.py --all --env production
    
    # Generate a report
    python database_ops.py --all --report
    ```
    
    See the code for more options and details.
    ```

## Files to Modify
- Create new directory structure and files as described above
- No modification of existing files is needed yet, as this is an additive change

## Verification
1. Run the tests to verify the new module:
   ```bash
   python -m unittest tests/database/test_population.py
   ```

2. Run the CLI tool to populate the database:
   ```bash
   python database_ops.py --all --report
   ```

3. Verify that the data was inserted correctly:
   ```bash
   # Query the database to check the data
   python -c "from database.utils.connection import get_supabase_client; print(get_supabase_client().from_('molecules').select('*').execute().data)"
   ```

## Notes for Roo Code Agent
- This implementation preserves all functionality from the existing scripts
- The modular design allows for easy extension with new data types
- Configuration options make it flexible for different environments
- After verification, the old scripts can be deprecated and eventually removed