#!/usr/bin/env python3
"""
Transaction Utilities for CryoProtect v2

This module provides robust transaction management utilities:
1. Safe transaction context manager with proper cleanup
2. Transaction retry logic with exponential backoff
3. Detection and handling of aborted transactions
4. Verification of transaction cleanup

Based on specifications in DATABASE_POPULATION_ISSUES.md (Section 2.3.1)
"""

import logging
import time
import functools
from typing import Any, Optional, Callable, TypeVar, cast
from contextlib import contextmanager

from db_connection_utils import get_db_connection

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('transaction_utils')

# Type variable for generic function return type
T = TypeVar('T')

@contextmanager
def safe_transaction():
    """
    Safe transaction context manager with cleanup.
    
    This context manager ensures proper transaction handling with:
    1. Detection of existing aborted transactions
    2. Automatic rollback in case of errors
    3. Verification of transaction cleanup
    4. Proper connection release
    
    Yields:
        Database connection with active transaction
    """
    connection = None
    transaction = None
    try:
        # Get a database connection
        connection = get_db_connection()
        if not connection:
            raise ConnectionError("Failed to establish database connection")
        
        # Check if already in transaction
        status_query = "SELECT current_setting('transaction_isolation')"
        try:
            if hasattr(connection, 'execute_query'):
                connection.execute_query(status_query)
            else:
                with connection.cursor() as cursor:
                    cursor.execute(status_query)
        except Exception as e:
            if "current transaction is aborted" in str(e):
                logger.warning("Detected aborted transaction, rolling back")
                if hasattr(connection, 'rollback'):
                    connection.rollback()
        
        # Start new transaction
        if hasattr(connection, 'begin_transaction'):
            transaction = connection.begin_transaction()
        else:
            connection.autocommit = False
            transaction = connection
        
        yield connection
        
        # Commit transaction
        if hasattr(connection, 'commit_transaction'):
            connection.commit_transaction(transaction)
        else:
            connection.commit()
            
    except Exception as e:
        # Rollback transaction on error
        logger.error(f"Transaction error: {str(e)}")
        if connection:
            if hasattr(connection, 'rollback_transaction') and transaction:
                try:
                    connection.rollback_transaction(transaction)
                except Exception as rollback_error:
                    logger.error(f"Rollback failed: {rollback_error}")
            elif hasattr(connection, 'rollback'):
                try:
                    connection.rollback()
                except Exception as rollback_error:
                    logger.error(f"Rollback failed: {rollback_error}")
        raise
    finally:
        # Verify transaction cleanup
        if connection:
            try:
                # Execute a simple query to verify the connection is in a good state
                if hasattr(connection, 'execute_query'):
                    connection.execute_query("SELECT 1")
                else:
                    with connection.cursor() as cursor:
                        cursor.execute("SELECT 1")
            except Exception:
                # If query fails, attempt emergency rollback
                logger.warning("Connection may be in an invalid state, attempting emergency rollback")
                try:
                    if hasattr(connection, 'rollback'):
                        connection.rollback()
                except Exception:
                    logger.error("Emergency rollback failed")
            
            # Release connection
            if hasattr(connection, 'close'):
                connection.close()

def with_transaction_retry(max_retries: int = 3, retry_delay: float = 0.5, 
                          backoff_factor: float = 2.0) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator to execute a function within a transaction with retry logic.
    
    Args:
        max_retries: Maximum number of retry attempts
        retry_delay: Initial delay between retries in seconds
        backoff_factor: Multiplier for the delay after each retry
        
    Returns:
        Decorator function
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> T:
            last_exception = None
            current_delay = retry_delay
            
            for attempt in range(max_retries):
                try:
                    with safe_transaction() as conn:
                        # Pass the connection as the first argument
                        return func(conn, *args, **kwargs)
                except Exception as e:
                    last_exception = e
                    logger.warning(f"Transaction attempt {attempt+1}/{max_retries} failed: {str(e)}")
                    
                    # Check if we should retry
                    if attempt < max_retries - 1:
                        logger.info(f"Retrying in {current_delay:.2f} seconds...")
                        time.sleep(current_delay)
                        current_delay *= backoff_factor
                    else:
                        logger.error(f"All {max_retries} transaction attempts failed")
                        break
            
            # If we get here, all retries failed
            if last_exception:
                raise last_exception
            
            # This should never happen, but just in case
            raise RuntimeError("Transaction failed for unknown reason")
        
        return cast(Callable[..., T], wrapper)
    
    return decorator

def execute_in_transaction(func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
    """
    Execute a function within a transaction.
    
    Args:
        func: Function to execute
        *args: Arguments to pass to the function
        **kwargs: Keyword arguments to pass to the function
        
    Returns:
        Result of the function
    """
    with safe_transaction() as conn:
        return func(conn, *args, **kwargs)

def is_transaction_active(connection: Any) -> bool:
    """
    Check if a transaction is currently active on the connection.
    
    Args:
        connection: Database connection
        
    Returns:
        True if a transaction is active, False otherwise
    """
    try:
        # Try to execute a query that would fail if in an aborted transaction
        if hasattr(connection, 'execute_query'):
            connection.execute_query("SELECT 1")
        else:
            with connection.cursor() as cursor:
                cursor.execute("SELECT 1")
        
        # If we get here, we're in a good state, but we need to check if autocommit is off
        if hasattr(connection, 'autocommit'):
            return not connection.autocommit
        
        # For connection wrappers, check if they have a transaction_active attribute
        if hasattr(connection, 'transaction_active'):
            return connection.transaction_active
        
        # If we can't determine, assume a transaction is active (safer)
        return True
    except Exception as e:
        # If we get an error about aborted transaction, then a transaction is active (but aborted)
        if "current transaction is aborted" in str(e):
            return True
        
        # For other errors, log and assume no transaction is active
        logger.warning(f"Error checking transaction status: {str(e)}")
        return False