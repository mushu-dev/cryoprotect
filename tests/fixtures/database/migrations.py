"""
Migration testing fixtures.

This module provides fixtures for testing database migrations.
"""

import os
import pytest
import tempfile
from typing import Dict, Any, Generator, List

from database.migrations import (
    initialize_migration_tracking,
    apply_migrations,
    rollback_migrations
)

@pytest.fixture
def migration_db(mock_db) -> Generator[Any, None, None]:
    """
    Provide a mock database prepared for migration testing.
    
    Yields:
        MockSupabase instance set up for migrations
    """
    # Initialize migration tracking
    initialize_migration_tracking(mock_db)
    
    # Register additional RPC handlers for migration testing
    def has_column_handler(params):
        table = params.get('table_name', '')
        column = params.get('column_name', '')
        
        if table in mock_db.tables:
            # Simplified mock - assume columns exist in test tables
            return {'data': [True]}
        return {'data': [False]}
        
    mock_db.register_rpc_handler('has_column', has_column_handler)
    
    yield mock_db

@pytest.fixture
def temp_migration_script() -> Generator[str, None, None]:
    """
    Create a temporary migration script for testing.
    
    Yields:
        Path to the temporary migration script
    """
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a migration script
        script_path = os.path.join(temp_dir, '999_test_migration.py')
        with open(script_path, 'w') as f:
            f.write('''
def apply(conn, environment):
    """Apply the migration."""
    conn.sql("""
        CREATE TABLE IF NOT EXISTS test_migration_table (
            id SERIAL PRIMARY KEY,
            name TEXT NOT NULL,
            created_at TIMESTAMP DEFAULT NOW()
        )
    """).execute()
    
def rollback(conn, environment):
    """Roll back the migration."""
    conn.sql("DROP TABLE IF EXISTS test_migration_table").execute()
''')
        yield script_path