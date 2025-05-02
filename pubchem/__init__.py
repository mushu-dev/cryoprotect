"""
PubChem API client package for CryoProtect v2.

This package provides a resilient client for interacting with the PubChem API,
featuring adaptive rate limiting, caching, retry logic, and circuit breaking.
It also includes simplified versions of the client and rate limiter for use
with the enhanced PubChem importer.
"""

from .client import ResilientPubChemClient
from .simple_client import PubChemClient
from .rate_limiter import AdaptiveRateLimiter
from .simple_rate_limiter import RateLimiter
from .cache import PubChemCache

__all__ = [
    "ResilientPubChemClient",
    "PubChemClient",
    "AdaptiveRateLimiter",
    "RateLimiter",
    "PubChemCache"
]