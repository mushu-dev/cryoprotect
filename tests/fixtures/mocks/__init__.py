"""
Mock implementations for testing.

This module provides mock implementations of external dependencies
for testing without requiring the actual dependencies.
"""

from tests.fixtures.mocks.supabase import (
    MockSupabase,
    patch_supabase_client
)

from tests.fixtures.mocks.rdkit import (
    MockRDKit,
    MockMolecule,
    patch_rdkit,
    mock_rdkit,
    molecule_factory
)

__all__ = [
    'MockSupabase',
    'patch_supabase_client',
    'MockRDKit',
    'MockMolecule',
    'patch_rdkit',
    'mock_rdkit',
    'molecule_factory'
]