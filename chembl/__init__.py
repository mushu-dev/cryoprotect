"""
ChEMBL API client for CryoProtect v2.

This package provides a robust client for interacting with the ChEMBL API,
featuring adaptive rate limiting, caching, retry logic, and circuit breaking.
"""

from .client import ResilientChEMBLClient

__all__ = ["ResilientChEMBLClient"]