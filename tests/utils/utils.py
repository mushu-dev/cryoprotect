"""
Testing utilities for CryoProtect test suite.
Provides helper functions and fixtures for test modules.
"""

import os
import tempfile
import shutil

def create_temp_dir():
    """
    Create a temporary directory for use in tests.
    Returns the path to the created directory.
    """
    return tempfile.mkdtemp()

def remove_temp_dir(path):
    """
    Remove a temporary directory and all its contents.
    Args:
        path (str): Path to the directory to remove.
    """
    shutil.rmtree(path, ignore_errors=True)

def get_test_file_path(filename):
    """
    Get the absolute path to a test file in the tests/data directory.
    Args:
        filename (str): Name of the test file.
    Returns:
        str: Absolute path to the test file.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base_dir, '..', 'data', filename)
