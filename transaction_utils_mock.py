#!/usr/bin/env python3
"""
Simplified mock version of transaction_utils.py for testing purposes.
"""

import logging
import functools
from typing import Any, Callable, TypeVar
from contextlib import contextmanager

from db_connection_utils_mock import get_db_connection, safe_transaction

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('transaction_utils_mock')

# Type variable for generic function return type
T = TypeVar('T')

def with_transaction_retry(max_retries: int = 3, retry_delay: float = 0.5):
    """
    Mock decorator to execute a function within a transaction with retry logic.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for attempt in range(max_retries):
                try:
                    with safe_transaction():
                        # Pass None as the connection since we're mocking
                        return func(None, *args, **kwargs)
                except Exception as e:
                    logger.warning(f"Mock transaction attempt {attempt+1}/{max_retries} failed: {str(e)}")
                    if attempt == max_retries - 1:
                        raise
            return None
        return wrapper
    return decorator

def execute_in_transaction(func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
    """
    Mock function to execute a function within a transaction.
    """
    with safe_transaction():
        return func(None, *args, **kwargs)

def is_transaction_active(connection: Any) -> bool:
    """
    Mock function to check if a transaction is active.
    """
    return True