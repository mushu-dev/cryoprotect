# Task 1.2: Implement Database Population Module

## Objective
Consolidate multiple population scripts into a coherent, modular database population system that follows consistent patterns and reduces code duplication.

## Context
The project currently has 10+ separate database population scripts (populate_molecules.py, populate_mixtures.py, etc.) with significant code duplication and inconsistent patterns. This makes maintenance difficult and increases the risk of bugs. By creating a unified population module, we'll improve code organization, reduce duplication, and enable better dependency management.

## Acceptance Criteria
- Create a functional database population package with clear modules for different entity types
- Implement a clean, consistent API for all population operations
- Provide a unified CLI entry point for database population
- Support different environments (development, staging, production)
- Include proper error handling and logging
- Maintain compatibility with existing data formats

## Implementation Steps

1. Create the core population module structure:
   ```
   database/population/
   ├── __init__.py          # Package exports
   ├── runner.py            # Main entry point
   ├── molecules.py         # Molecule population
   ├── mixtures.py          # Mixture population 
   └── utils.py             # Shared utilities
   ```

2. Implement the population module `__init__.py`:
   ```python
   """
   Database population module for CryoProtect v2.
   
   This module provides tools for populating the database with
   molecules, mixtures, experiments, and other scientific data.
   """
   
   from database.population.runner import (
       populate_all,
       populate_specific,
       populate_from_file
   )
   
   from database.population.molecules import populate_molecules
   from database.population.mixtures import populate_mixtures
   
   __all__ = [
       'populate_all',
       'populate_specific',
       'populate_from_file',
       'populate_molecules',
       'populate_mixtures'
   ]
   ```

3. Implement the main runner module:
   ```python
   """
   Main entry point for database population operations.
   
   This module provides high-level functions for populating the database
   with various types of data.
   """
   
   import logging
   import os
   from typing import Dict, List, Optional, Union
   from database.utils.connection import create_connection
   from database.population.molecules import populate_molecules
   from database.population.mixtures import populate_mixtures
   
   logger = logging.getLogger(__name__)
   
   def populate_all(
       environment: str = 'development',
       config: Optional[Dict] = None
   ) -> Dict[str, int]:
       """
       Populate all database tables based on environment.
       
       Args:
           environment: Target environment ('development', 'staging', 'production')
           config: Optional configuration dictionary
           
       Returns:
           Dictionary with table names as keys and count of inserted records as values
       """
       logger.info(f"Populating all tables in {environment} environment")
       
       conn = create_connection(config)
       results = {}
       
       # Populate in dependency order
       results['molecules'] = populate_molecules(conn, environment)
       results['mixtures'] = populate_mixtures(conn, environment)
       
       logger.info(f"Population complete. Results: {results}")
       return results
       
   def populate_specific(
       tables: List[str],
       environment: str = 'development',
       config: Optional[Dict] = None
   ) -> Dict[str, int]:
       """
       Populate specific tables based on environment.
       
       Args:
           tables: List of tables to populate
           environment: Target environment ('development', 'staging', 'production')
           config: Optional configuration dictionary
           
       Returns:
           Dictionary with table names as keys and count of inserted records as values
       """
       logger.info(f"Populating specific tables {tables} in {environment} environment")
       
       conn = create_connection(config)
       results = {}
       
       if 'molecules' in tables:
           results['molecules'] = populate_molecules(conn, environment)
           
       if 'mixtures' in tables:
           results['mixtures'] = populate_mixtures(conn, environment)
       
       logger.info(f"Population complete. Results: {results}")
       return results
       
   def populate_from_file(
       file_path: str, 
       table: str,
       environment: str = 'development',
       config: Optional[Dict] = None
   ) -> int:
       """
       Populate a table from a data file.
       
       Args:
           file_path: Path to the data file (JSON or CSV)
           table: Name of the table to populate
           environment: Target environment ('development', 'staging', 'production')
           config: Optional configuration dictionary
           
       Returns:
           Number of records inserted
       """
       logger.info(f"Populating {table} from file {file_path}")
       
       if not os.path.exists(file_path):
           logger.error(f"File not found: {file_path}")
           return 0
           
       conn = create_connection(config)
       
       # Dispatch to appropriate handler based on table name
       if table == 'molecules':
           return populate_molecules(conn, environment, file_path)
       elif table == 'mixtures':
           return populate_mixtures(conn, environment, file_path)
       else:
           logger.error(f"Unsupported table: {table}")
           return 0
       
   def main():
       """
       CLI entry point for database population.
       """
       import argparse
       
       parser = argparse.ArgumentParser(description='Populate CryoProtect database tables')
       parser.add_argument(
           '--env', 
           choices=['development', 'staging', 'production'],
           default='development',
           help='Target environment'
       )
       parser.add_argument(
           '--tables', 
           nargs='+',
           help='Specific tables to populate'
       )
       parser.add_argument(
           '--file',
           help='Path to data file (for single table population)'
       )
       parser.add_argument(
           '--table',
           help='Table to populate from file (used with --file)'
       )
       
       args = parser.parse_args()
       
       if args.file and args.table:
           count = populate_from_file(args.file, args.table, args.env)
           print(f"Populated {args.table} with {count} records from {args.file}")
       elif args.tables:
           results = populate_specific(args.tables, args.env)
           for table, count in results.items():
               print(f"Populated {table} with {count} records")
       else:
           results = populate_all(args.env)
           for table, count in results.items():
               print(f"Populated {table} with {count} records")
       
   if __name__ == '__main__':
       main()
   ```

4. Implement the molecules population module:
   ```python
   """
   Module for populating the molecules table.
   
   This module provides functions for populating the molecules table
   from various data sources.
   """
   
   import json
   import logging
   import os
   from typing import Any, Dict, List, Optional, Union
   
   logger = logging.getLogger(__name__)
   
   def _load_molecule_data(
       environment: str,
       file_path: Optional[str] = None
   ) -> List[Dict]:
       """
       Load molecule data from appropriate source.
       
       Args:
           environment: Target environment
           file_path: Optional path to data file
           
       Returns:
           List of molecule records
       """
       # If a specific file is provided, use it
       if file_path and os.path.exists(file_path):
           with open(file_path, 'r') as f:
               data = json.load(f)
               return data
               
       # Otherwise, use environment-specific data
       base_dir = os.path.abspath(os.path.dirname(os.path.dirname(
           os.path.dirname(os.path.dirname(__file__)))))
           
       if environment == 'production':
           data_path = os.path.join(base_dir, 'data', 'production', 'molecules.json')
       elif environment == 'staging':
           data_path = os.path.join(base_dir, 'data', 'staging', 'molecules.json')
       else:  # development
           data_path = os.path.join(base_dir, 'data', 'development', 'molecules.json')
           
       if not os.path.exists(data_path):
           logger.warning(f"Data file not found: {data_path}")
           return []
           
       with open(data_path, 'r') as f:
           data = json.load(f)
           return data
   
   def populate_molecules(
       conn: Any,
       environment: str = 'development',
       file_path: Optional[str] = None
   ) -> int:
       """
       Populate the molecules table.
       
       Args:
           conn: Database connection
           environment: Target environment
           file_path: Optional path to data file
           
       Returns:
           Number of records inserted
       """
       logger.info(f"Populating molecules table in {environment} environment")
       
       # Load data
       molecule_data = _load_molecule_data(environment, file_path)
       if not molecule_data:
           logger.warning("No molecule data found")
           return 0
           
       # Insert data
       count = 0
       for molecule in molecule_data:
           try:
               # Example using Supabase client
               result = conn.table('molecules').insert(molecule).execute()
               count += 1
           except Exception as e:
               logger.error(f"Error inserting molecule {molecule.get('name', 'unknown')}: {str(e)}")
               
       logger.info(f"Inserted {count} molecules")
       return count
   ```

5. Implement the mixtures population module:
   ```python
   """
   Module for populating the mixtures table.
   
   This module provides functions for populating the mixtures table
   from various data sources.
   """
   
   import json
   import logging
   import os
   from typing import Any, Dict, List, Optional, Union
   
   logger = logging.getLogger(__name__)
   
   def _load_mixture_data(
       environment: str,
       file_path: Optional[str] = None
   ) -> List[Dict]:
       """
       Load mixture data from appropriate source.
       
       Args:
           environment: Target environment
           file_path: Optional path to data file
           
       Returns:
           List of mixture records
       """
       # If a specific file is provided, use it
       if file_path and os.path.exists(file_path):
           with open(file_path, 'r') as f:
               data = json.load(f)
               return data
               
       # Otherwise, use environment-specific data
       base_dir = os.path.abspath(os.path.dirname(os.path.dirname(
           os.path.dirname(os.path.dirname(__file__)))))
           
       if environment == 'production':
           data_path = os.path.join(base_dir, 'data', 'production', 'mixtures.json')
       elif environment == 'staging':
           data_path = os.path.join(base_dir, 'data', 'staging', 'mixtures.json')
       else:  # development
           data_path = os.path.join(base_dir, 'data', 'development', 'mixtures.json')
           
       if not os.path.exists(data_path):
           logger.warning(f"Data file not found: {data_path}")
           return []
           
       with open(data_path, 'r') as f:
           data = json.load(f)
           return data
   
   def populate_mixtures(
       conn: Any,
       environment: str = 'development',
       file_path: Optional[str] = None
   ) -> int:
       """
       Populate the mixtures table.
       
       Args:
           conn: Database connection
           environment: Target environment
           file_path: Optional path to data file
           
       Returns:
           Number of records inserted
       """
       logger.info(f"Populating mixtures table in {environment} environment")
       
       # Load data
       mixture_data = _load_mixture_data(environment, file_path)
       if not mixture_data:
           logger.warning("No mixture data found")
           return 0
           
       # Insert data
       count = 0
       for mixture in mixture_data:
           try:
               # Example using Supabase client
               result = conn.table('mixtures').insert(mixture).execute()
               
               # Handle mixture components if present
               if 'components' in mixture and mixture['components']:
                   for component in mixture['components']:
                       component['mixture_id'] = result.data[0]['id']
                       conn.table('mixture_components').insert(component).execute()
                       
               count += 1
           except Exception as e:
               logger.error(f"Error inserting mixture {mixture.get('name', 'unknown')}: {str(e)}")
               
       logger.info(f"Inserted {count} mixtures")
       return count
   ```

6. Update the main database `__init__.py` to include the population module:
   ```python
   """
   Database package for CryoProtect v2.
   
   This package centralizes all database-related functionality including
   models, operations, migrations, and verification tools.
   """
   
   from database.main import (
       initialize_database,
       populate_database,
       verify_database,
       backup_database,
       restore_database
   )
   
   # Include the population module exports
   from database.population import (
       populate_all,
       populate_specific,
       populate_from_file
   )
   
   __all__ = [
       'initialize_database',
       'populate_database',
       'verify_database',
       'backup_database',
       'restore_database',
       'populate_all',
       'populate_specific',
       'populate_from_file'
   ]
   ```

7. Update the main entry point to use the new population module:
   ```python
   """
   Main entry points for the database package.
   
   This module provides high-level functions for common database operations.
   """
   
   import logging
   from typing import Dict, List, Optional, Union
   
   # Import from the population module
   from database.population.runner import populate_all, populate_specific
   
   logger = logging.getLogger(__name__)
   
   def initialize_database(config: Dict = None) -> bool:
       """
       Initialize the database with schema and required tables.
       
       Args:
           config: Configuration dictionary with connection details
           
       Returns:
           True if initialization was successful, False otherwise
       """
       logger.info("Initializing database")
       # Implementation to be added
       return True
       
   def populate_database(
       tables: Optional[List[str]] = None,
       environment: str = 'development',
       config: Optional[Dict] = None
   ) -> Dict[str, int]:
       """
       Populate the database with data.
       
       Args:
           tables: List of tables to populate, or None for all tables
           environment: Environment to use ('development', 'staging', 'production')
           config: Optional configuration dictionary
           
       Returns:
           Dictionary with table names as keys and count of inserted records as values
       """
       logger.info(f"Populating database in {environment} environment")
       
       if tables:
           return populate_specific(tables, environment, config)
       else:
           return populate_all(environment, config)
       
   # ... rest of the file unchanged
   ```

## Files to Modify
- Create new directory: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/population/`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/population/__init__.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/population/runner.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/population/molecules.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/population/mixtures.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/population/README.md`
- Update existing file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/__init__.py`
- Update existing file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/main.py`

## Verification
1. Verify all modules import properly without errors
2. Check that population functions can be called from the main package
3. Verify CLI entry point works with all parameters
4. Ensure error handling works for invalid input
5. Verify environment-specific data loading works as expected
6. Ensure code follows project style guidelines