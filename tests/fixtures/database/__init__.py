"""
Database fixtures for testing.

This module provides fixtures for database testing, including
connections, data generation, and migration testing.
"""

from tests.fixtures.database.connection import (
    db_config,
    mock_db,
    real_db_connection,
    populated_mock_db
)

from tests.fixtures.database.migrations import (
    migration_db,
    temp_migration_script
)

__all__ = [
    'db_config',
    'mock_db',
    'real_db_connection',
    'populated_mock_db',
    'migration_db',
    'temp_migration_script'
]