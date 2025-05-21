"""
Tests for Convex migrations.

This module provides tests for the Convex migration system.
"""

import os
import unittest
from unittest.mock import patch, MagicMock, Mock

from database.migrations.convex_bridge import (
    _is_convex_enabled,
    _run_convex_cli,
    get_convex_migration_status,
    apply_convex_migrations,
    rollback_convex_migrations,
    create_convex_migration,
    sync_supabase_to_convex,
    run_universal_migrations
)

class TestConvexMigrations(unittest.TestCase):
    """Tests for Convex migrations."""
    
    def setUp(self):
        """Set up test environment."""
        # Set environment variables
        os.environ['USE_CONVEX'] = 'true'
        os.environ['CONVEX_URL'] = 'https://test-123.convex.cloud'
        
        # Mock subprocess.Popen
        self.popen_mock = Mock()
        self.popen_mock.returncode = 0
        self.popen_mock.communicate.return_value = ('{"success": true}', '')
        
    @patch('subprocess.Popen')
    def test_run_convex_cli(self, popen_mock):
        """Test running Convex CLI."""
        popen_mock.return_value = self.popen_mock
        
        # Test basic command
        result = _run_convex_cli('status')
        
        popen_mock.assert_called_once()
        self.assertEqual(result, {'success': True})
        
        # Reset mock
        popen_mock.reset_mock()
        
        # Test with arguments
        result = _run_convex_cli('apply', {'dry-run': True, 'target': '002'})
        
        popen_mock.assert_called_once()
        self.assertEqual(result, {'success': True})
    
    def test_is_convex_enabled(self):
        """Test checking if Convex is enabled."""
        # Test with USE_CONVEX=true
        self.assertTrue(_is_convex_enabled())
        
        # Test with USE_CONVEX=false
        os.environ['USE_CONVEX'] = 'false'
        self.assertFalse(_is_convex_enabled())
        
        # Reset for other tests
        os.environ['USE_CONVEX'] = 'true'
    
    @patch('database.migrations.convex_bridge._run_convex_cli')
    def test_get_convex_migration_status(self, run_cli_mock):
        """Test getting Convex migration status."""
        run_cli_mock.return_value = {
            'success': True,
            'migrations': {
                '001': {
                    'version': '001',
                    'name': 'initial_schema',
                    'applied': True,
                    'appliedAt': 1620000000000
                }
            }
        }
        
        result = get_convex_migration_status()
        
        run_cli_mock.assert_called_once_with('status', {'format': 'json'})
        self.assertEqual(result, run_cli_mock.return_value)
    
    @patch('database.migrations.convex_bridge._run_convex_cli')
    def test_apply_convex_migrations(self, run_cli_mock):
        """Test applying Convex migrations."""
        run_cli_mock.return_value = {
            'success': True,
            'applied': [
                {
                    'version': '001',
                    'name': 'initial_schema',
                    'status': 'applied'
                }
            ],
            'message': 'Applied 1 migrations'
        }
        
        # Test without target version
        result = apply_convex_migrations()
        
        run_cli_mock.assert_called_once_with('apply', {'dry-run': False})
        self.assertEqual(result, run_cli_mock.return_value)
        
        # Reset mock
        run_cli_mock.reset_mock()
        
        # Test with target version and dry run
        result = apply_convex_migrations(target_version='002', dry_run=True)
        
        run_cli_mock.assert_called_once_with('apply', {'dry-run': True, 'target': '002'})
        self.assertEqual(result, run_cli_mock.return_value)
    
    @patch('database.migrations.convex_bridge._run_convex_cli')
    def test_rollback_convex_migrations(self, run_cli_mock):
        """Test rolling back Convex migrations."""
        run_cli_mock.return_value = {
            'success': True,
            'rolledBack': [
                {
                    'version': '002',
                    'name': 'add_teams',
                    'status': 'rolled_back'
                }
            ],
            'message': 'Rolled back 1 migrations'
        }
        
        # Test without target version
        result = rollback_convex_migrations()
        
        run_cli_mock.assert_called_once_with('rollback', {'dry-run': False})
        self.assertEqual(result, run_cli_mock.return_value)
        
        # Reset mock
        run_cli_mock.reset_mock()
        
        # Test with target version and dry run
        result = rollback_convex_migrations(target_version='001', dry_run=True)
        
        run_cli_mock.assert_called_once_with('rollback', {'dry-run': True, 'target': '001'})
        self.assertEqual(result, run_cli_mock.return_value)
    
    @patch('database.migrations.convex_bridge._run_convex_cli')
    def test_create_convex_migration(self, run_cli_mock):
        """Test creating a Convex migration."""
        run_cli_mock.return_value = {
            'success': True,
            'message': 'Created migration file: /path/to/003_add_features.json'
        }
        
        result = create_convex_migration('add_features')
        
        run_cli_mock.assert_called_once_with('create', {'name': 'add_features'})
        self.assertEqual(result, run_cli_mock.return_value)
    
    @patch('database.migrations.convex_bridge.get_db_client')
    def test_sync_supabase_to_convex(self, get_db_client_mock):
        """Test syncing data from Supabase to Convex."""
        # Mock Supabase client
        supabase_mock = MagicMock()
        supabase_mock.table.return_value.select.return_value.limit.return_value.execute.return_value = MagicMock(
            data=[{'id': 1, 'name': 'Test'}],
            error=None
        )
        
        # Mock Convex client
        convex_mock = MagicMock()
        convex_mock.table.return_value.insert.return_value.execute.return_value = MagicMock(
            error=None
        )
        
        # Configure get_db_client to return our mocks
        def get_db_client_side_effect(force_convex):
            return convex_mock if force_convex else supabase_mock
        
        get_db_client_mock.side_effect = get_db_client_side_effect
        
        # Test syncing
        result = sync_supabase_to_convex('users')
        
        # Verify calls
        supabase_mock.table.assert_called_once_with('users')
        supabase_mock.table.return_value.select.assert_called_once_with('*')
        convex_mock.table.assert_called_once_with('users')
        
        # Verify result
        self.assertTrue(result['success'])
        self.assertEqual(result['count'], 1)
    
    @patch('database.migrations.convex_bridge.apply_convex_migrations')
    @patch('database.migrations.convex_bridge.get_db_client')
    @patch('database.migrations.runner.apply_migrations')
    def test_run_universal_migrations(
        self, apply_migrations_mock, get_db_client_mock, apply_convex_migrations_mock
    ):
        """Test running universal migrations."""
        # Mock responses
        apply_migrations_mock.return_value = [
            {'version': '001', 'name': 'initial_schema', 'status': 'applied'}
        ]
        
        apply_convex_migrations_mock.return_value = {
            'success': True,
            'applied': [
                {'version': '001', 'name': 'initial_schema', 'status': 'applied'}
            ]
        }
        
        # Test running migrations
        result = run_universal_migrations(target_version='001', dry_run=True)
        
        # Verify calls
        get_db_client_mock.assert_called_once_with(force_convex=False)
        apply_migrations_mock.assert_called_once()
        apply_convex_migrations_mock.assert_called_once_with(
            target_version='001',
            dry_run=True
        )
        
        # Verify result
        self.assertIn('supabase', result)
        self.assertIn('convex', result)
        self.assertTrue(result['supabase']['success'])
        self.assertEqual(result['convex'], apply_convex_migrations_mock.return_value)

if __name__ == '__main__':
    unittest.main()