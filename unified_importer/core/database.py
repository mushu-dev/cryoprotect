"""
Database operations for the unified molecular importer.

This module provides a unified interface for database operations,
supporting both direct database connections and Supabase.
"""

import os
import time
import asyncio
import logging
from typing import Dict, List, Any, Optional, Tuple, Union, Callable, AsyncIterator
import contextlib
import json
import uuid
import queue
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from unittest.mock import MagicMock

# Optional imports based on available backends
try:
    import psycopg2
    import psycopg2.extras
    PSYCOPG2_AVAILABLE = True
except ImportError:
    PSYCOPG2_AVAILABLE = False

try:
    from supabase import create_client, Client
    SUPABASE_AVAILABLE = True
except ImportError:
    SUPABASE_AVAILABLE = False


class DatabaseError(Exception):
    """Base class for database-related exceptions."""
    pass


class ConnectionError(DatabaseError):
    """Exception raised for connection-related errors."""
    pass


class TransactionError(DatabaseError):
    """Exception raised for transaction-related errors."""
    pass


class QueryError(DatabaseError):
    """Exception raised for query execution errors."""
    pass


@dataclass
class ConnectionInfo:
    """Information about a database connection."""
    host: str = "localhost"
    port: int = 5432
    database: str = "postgres"
    user: str = "postgres"
    password: str = ""
    connection_string: str = ""
    connected_at: float = 0.0
    last_used: float = 0.0
    is_healthy: bool = False
    connection_type: str = "direct"
    extra: Dict[str, Any] = field(default_factory=dict)


class ConnectionPool:
    """
    Connection pool for database connections.

    This class manages a pool of database connections that can be
    reused across multiple operations. It supports both synchronous
    and asynchronous access patterns.
    """

    def __init__(
        self,
        min_size: int = 2,
        max_size: int = 10,
        timeout: float = 30.0,
        max_lifetime: float = 3600.0,
        connection_type: str = "direct",
        connection_params: Optional[Dict[str, Any]] = None,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the connection pool.

        Args:
            min_size: Minimum number of connections to maintain
            max_size: Maximum number of connections allowed
            timeout: Connection timeout in seconds
            max_lifetime: Maximum lifetime of a connection in seconds
            connection_type: Type of connection ("direct", "supabase")
            connection_params: Parameters for creating connections
            logger: Logger instance
        """
        self.min_size = min_size
        self.max_size = max_size
        self.timeout = timeout
        self.max_lifetime = max_lifetime
        self.connection_type = connection_type
        self.connection_params = connection_params or {}
        self.logger = logger or logging.getLogger(__name__)

        # Pool state
        self.pool = queue.Queue(maxsize=max_size)
        self.active_connections = 0
        self.lock = threading.RLock()
        self.stats = {
            "created": 0,
            "acquired": 0,
            "released": 0,
            "discarded": 0,
            "errors": 0,
            "timeouts": 0,
            "last_error": None,
            "peak_active": 0
        }

        # Async pool state
        self.async_pool = asyncio.Queue(maxsize=max_size)
        self.async_lock = asyncio.Lock()
        self.executor = ThreadPoolExecutor(max_workers=max_size)

        # Connection tracking
        self.connection_info = {}

        # Initialize the pool
        self._initialize_pool()

    def _initialize_pool(self) -> None:
        """Initialize the pool with minimum connections."""
        self.logger.info(
            f"Initializing connection pool with {self.min_size} connections"
        )

        for _ in range(self.min_size):
            try:
                conn = self._create_connection()
                conn_id = str(uuid.uuid4())

                # Store connection info
                self.connection_info[conn_id] = ConnectionInfo(
                    host=self.connection_params.get("host", "localhost"),
                    port=self.connection_params.get("port", 5432),
                    database=self.connection_params.get("database", "postgres"),
                    user=self.connection_params.get("user", "postgres"),
                    connected_at=time.time(),
                    last_used=time.time(),
                    is_healthy=True,
                    connection_type=self.connection_type
                )

                # Put connection in the pool with its ID
                self.pool.put((conn, conn_id))

                # Also initialize async pool
                self.async_pool.put_nowait((conn, conn_id))
            except Exception as e:
                self.logger.error(f"Error initializing connection: {str(e)}")
                self.stats["errors"] += 1

    def _create_connection(self):
        """Create a new database connection based on connection type."""
        with self.lock:
            self.active_connections += 1
            self.stats["created"] += 1

            if self.active_connections > self.stats["peak_active"]:
                self.stats["peak_active"] = self.active_connections

        if self.connection_type == "direct":
            # Direct PostgreSQL connection
            import psycopg2
            import psycopg2.extras

            conn_params = {
                'dbname': self.connection_params.get('database', 'postgres'),
                'user': self.connection_params.get('user', 'postgres'),
                'password': self.connection_params.get('password', ''),
                'host': self.connection_params.get('host', 'localhost'),
                'port': self.connection_params.get('port', 5432)
            }

            conn = psycopg2.connect(**conn_params)
            return conn
        elif self.connection_type == "supabase":
            # Supabase connection
            from supabase import create_client

            supabase_url = self.connection_params.get('url') or os.environ.get('SUPABASE_URL')
            supabase_key = self.connection_params.get('key') or os.environ.get('SUPABASE_KEY')

            if not supabase_url or not supabase_key:
                raise ValueError("Supabase URL and key are required")

            conn = create_client(supabase_url, supabase_key)
            return conn
        else:
            raise ValueError(f"Unsupported connection type: {self.connection_type}")

    def get_connection(self):
        """Get a connection from the pool, creating a new one if needed."""
        conn = None
        conn_id = None
        is_new = False

        try:
            # Try to get a connection from the pool
            try:
                conn, conn_id = self.pool.get(block=False)

                # Check connection age
                conn_info = self.connection_info.get(conn_id)
                if conn_info and time.time() - conn_info.connected_at > self.max_lifetime:
                    # Connection is too old, discard it
                    self._close_connection(conn)
                    with self.lock:
                        self.active_connections -= 1
                        self.stats["discarded"] += 1
                        del self.connection_info[conn_id]

                    conn = None
                    conn_id = None
                else:
                    # Update last used time
                    if conn_info:
                        conn_info.last_used = time.time()

                    with self.lock:
                        self.stats["acquired"] += 1
            except queue.Empty:
                # Pool is empty, we'll create a new connection
                pass

            # Create a new connection if needed
            if conn is None:
                with self.lock:
                    if self.active_connections < self.max_size:
                        conn = self._create_connection()
                        conn_id = str(uuid.uuid4())
                        is_new = True

                        # Store connection info
                        self.connection_info[conn_id] = ConnectionInfo(
                            host=self.connection_params.get("host", "localhost"),
                            port=self.connection_params.get("port", 5432),
                            database=self.connection_params.get("database", "postgres"),
                            user=self.connection_params.get("user", "postgres"),
                            connected_at=time.time(),
                            last_used=time.time(),
                            is_healthy=True,
                            connection_type=self.connection_type
                        )

                        with self.lock:
                            self.stats["acquired"] += 1
                    else:
                        # We've reached max_size, wait for a connection
                        try:
                            conn, conn_id = self.pool.get(block=True, timeout=self.timeout)

                            # Update last used time
                            conn_info = self.connection_info.get(conn_id)
                            if conn_info:
                                conn_info.last_used = time.time()

                            with self.lock:
                                self.stats["acquired"] += 1
                        except queue.Empty:
                            with self.lock:
                                self.stats["timeouts"] += 1
                                self.stats["errors"] += 1
                                self.stats["last_error"] = "Timed out waiting for connection"

                            raise TimeoutError("Timed out waiting for database connection")
        except Exception as e:
            with self.lock:
                self.stats["errors"] += 1
                self.stats["last_error"] = str(e)

            self.logger.error(f"Error getting connection: {str(e)}")
            raise

        return conn, conn_id, is_new

    def release_connection(self, conn, conn_id):
        """Release a connection back to the pool."""
        try:
            # Update last used time
            conn_info = self.connection_info.get(conn_id)
            if conn_info:
                conn_info.last_used = time.time()

            # Check connection age
            if conn_info and time.time() - conn_info.connected_at > self.max_lifetime:
                # Connection is too old, discard it
                self._close_connection(conn)
                with self.lock:
                    self.active_connections -= 1
                    self.stats["discarded"] += 1
                    if conn_id in self.connection_info:
                        del self.connection_info[conn_id]
            else:
                # Put the connection back in the pool
                self.pool.put((conn, conn_id))
                with self.lock:
                    self.stats["released"] += 1
        except Exception as e:
            with self.lock:
                self.stats["errors"] += 1
                self.stats["last_error"] = str(e)

            # Ensure we don't leak connections
            try:
                self._close_connection(conn)
            except Exception:
                pass

            with self.lock:
                self.active_connections -= 1

            self.logger.error(f"Error releasing connection: {str(e)}")

    def _close_connection(self, conn):
        """Close a database connection."""
        try:
            if self.connection_type == "direct":
                conn.close()
            # Supabase connections don't need explicit closing
        except Exception as e:
            self.logger.error(f"Error closing connection: {str(e)}")

    def close_all(self):
        """Close all connections in the pool."""
        self.logger.info("Closing all connections in the pool")

        # Clear the pool and close each connection
        closed_count = 0
        error_count = 0

        while not self.pool.empty():
            try:
                conn, conn_id = self.pool.get(block=False)
                try:
                    self._close_connection(conn)
                    closed_count += 1
                except Exception as e:
                    self.logger.error(f"Error closing connection: {str(e)}")
                    error_count += 1
                finally:
                    with self.lock:
                        self.active_connections -= 1
                        if conn_id in self.connection_info:
                            del self.connection_info[conn_id]
            except queue.Empty:
                break

        # Clear async pool as well
        while not self.async_pool.empty():
            try:
                self.async_pool.get_nowait()
            except asyncio.QueueEmpty:
                break

        # Shutdown executor
        self.executor.shutdown(wait=False)

        with self.lock:
            self.active_connections = 0

        self.logger.info(f"Closed {closed_count} connections with {error_count} errors")

    def get_stats(self):
        """Get statistics about the connection pool."""
        with self.lock:
            stats = self.stats.copy()
            stats["pool_size"] = self.pool.qsize()
            stats["active_connections"] = self.active_connections
            stats["min_size"] = self.min_size
            stats["max_size"] = self.max_size
            stats["connection_age_stats"] = self._get_connection_age_stats()
            return stats

    def _get_connection_age_stats(self):
        """Get statistics about connection ages."""
        now = time.time()
        connection_ages = [now - info.connected_at for info in self.connection_info.values()]

        if not connection_ages:
            return {
                "count": 0,
                "min": 0,
                "max": 0,
                "avg": 0
            }

        return {
            "count": len(connection_ages),
            "min": min(connection_ages),
            "max": max(connection_ages),
            "avg": sum(connection_ages) / len(connection_ages)
        }

    # Async methods

    async def get_connection_async(self):
        """Get a connection asynchronously from the pool."""
        conn = None
        conn_id = None
        is_new = False

        try:
            # Try to get a connection from the pool
            try:
                conn, conn_id = await self.async_pool.get_nowait()

                # Check connection age
                conn_info = self.connection_info.get(conn_id)
                if conn_info and time.time() - conn_info.connected_at > self.max_lifetime:
                    # Connection is too old, discard it
                    await asyncio.to_thread(self._close_connection, conn)
                    async with self.async_lock:
                        self.active_connections -= 1
                        self.stats["discarded"] += 1
                        if conn_id in self.connection_info:
                            del self.connection_info[conn_id]

                    conn = None
                    conn_id = None
                else:
                    # Update last used time
                    if conn_info:
                        conn_info.last_used = time.time()

                    async with self.async_lock:
                        self.stats["acquired"] += 1
            except asyncio.QueueEmpty:
                # Pool is empty, we'll create a new connection
                pass

            # Create a new connection if needed
            if conn is None:
                async with self.async_lock:
                    if self.active_connections < self.max_size:
                        # Use ThreadPoolExecutor to run blocking connection creation
                        conn = await asyncio.to_thread(self._create_connection)
                        conn_id = str(uuid.uuid4())
                        is_new = True

                        # Store connection info
                        self.connection_info[conn_id] = ConnectionInfo(
                            host=self.connection_params.get("host", "localhost"),
                            port=self.connection_params.get("port", 5432),
                            database=self.connection_params.get("database", "postgres"),
                            user=self.connection_params.get("user", "postgres"),
                            connected_at=time.time(),
                            last_used=time.time(),
                            is_healthy=True,
                            connection_type=self.connection_type
                        )

                        self.stats["acquired"] += 1
                    else:
                        # We've reached max_size, wait for a connection
                        try:
                            # Use asyncio.wait_for to implement timeout
                            conn, conn_id = await asyncio.wait_for(
                                self.async_pool.get(),
                                timeout=self.timeout
                            )

                            # Update last used time
                            conn_info = self.connection_info.get(conn_id)
                            if conn_info:
                                conn_info.last_used = time.time()

                            async with self.async_lock:
                                self.stats["acquired"] += 1
                        except asyncio.TimeoutError:
                            async with self.async_lock:
                                self.stats["timeouts"] += 1
                                self.stats["errors"] += 1
                                self.stats["last_error"] = "Timed out waiting for connection"

                            raise TimeoutError("Timed out waiting for database connection")
        except Exception as e:
            async with self.async_lock:
                self.stats["errors"] += 1
                self.stats["last_error"] = str(e)

            self.logger.error(f"Error getting connection asynchronously: {str(e)}")
            raise

        return conn, conn_id, is_new

    async def release_connection_async(self, conn, conn_id):
        """Release a connection back to the async pool."""
        try:
            # Update last used time
            conn_info = self.connection_info.get(conn_id)
            if conn_info:
                conn_info.last_used = time.time()

            # Check connection age
            if conn_info and time.time() - conn_info.connected_at > self.max_lifetime:
                # Connection is too old, discard it
                await asyncio.to_thread(self._close_connection, conn)
                async with self.async_lock:
                    self.active_connections -= 1
                    self.stats["discarded"] += 1
                    if conn_id in self.connection_info:
                        del self.connection_info[conn_id]
            else:
                # Put the connection back in the pool
                await self.async_pool.put((conn, conn_id))
                async with self.async_lock:
                    self.stats["released"] += 1
        except Exception as e:
            async with self.async_lock:
                self.stats["errors"] += 1
                self.stats["last_error"] = str(e)

            # Ensure we don't leak connections
            try:
                await asyncio.to_thread(self._close_connection, conn)
            except Exception:
                pass

            async with self.async_lock:
                self.active_connections -= 1

            self.logger.error(f"Error releasing connection asynchronously: {str(e)}")

    async def close_all_async(self):
        """Close all connections in the pool asynchronously."""
        self.logger.info("Closing all connections in the pool asynchronously")

        # Clear the async pool and close each connection
        closed_count = 0
        error_count = 0

        while not self.async_pool.empty():
            try:
                conn, conn_id = await self.async_pool.get_nowait()
                try:
                    await asyncio.to_thread(self._close_connection, conn)
                    closed_count += 1
                except Exception as e:
                    self.logger.error(f"Error closing connection asynchronously: {str(e)}")
                    error_count += 1
                finally:
                    async with self.async_lock:
                        self.active_connections -= 1
                        if conn_id in self.connection_info:
                            del self.connection_info[conn_id]
            except asyncio.QueueEmpty:
                break

        # Clear sync pool as well
        while not self.pool.empty():
            try:
                conn, conn_id = self.pool.get(block=False)
                try:
                    self._close_connection(conn)
                    closed_count += 1
                except Exception as e:
                    self.logger.error(f"Error closing connection: {str(e)}")
                    error_count += 1
                finally:
                    with self.lock:
                        self.active_connections -= 1
                        if conn_id in self.connection_info:
                            del self.connection_info[conn_id]
            except queue.Empty:
                break

        # Shutdown executor
        self.executor.shutdown(wait=False)

        async with self.async_lock:
            self.active_connections = 0

        self.logger.info(f"Closed {closed_count} connections with {error_count} errors")

    async def get_stats_async(self):
        """Get statistics about the connection pool asynchronously."""
        async with self.async_lock:
            stats = self.stats.copy()
            stats["pool_size"] = self.async_pool.qsize()
            stats["active_connections"] = self.active_connections
            stats["min_size"] = self.min_size
            stats["max_size"] = self.max_size
            stats["connection_age_stats"] = self._get_connection_age_stats()
            return stats


class DatabaseOperations:
    """
    Unified interface for database operations across different backends.

    This class provides methods for common database operations like
    inserting molecules, properties, and synonyms. It supports both
    direct PostgreSQL connections and Supabase with advanced connection
    pooling, transaction management, and retry capabilities.
    """
    
    def __init__(
        self,
        connection_type: str = "direct",
        connection_params: Optional[Dict[str, Any]] = None,
        batch_size: int = 100,
        max_retries: int = 3,
        retry_delay: float = 2.0,
        pool_min_size: int = 2,
        pool_max_size: int = 10,
        pool_timeout: float = 30.0,
        max_pool_connections: int = 20,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize database operations with the specified connection.

        Args:
            connection_type: Type of connection ("direct", "supabase")
            connection_params: Connection parameters
            batch_size: Number of items to include in batch operations
            max_retries: Maximum number of retry attempts
            retry_delay: Delay between retries (in seconds)
            pool_min_size: Minimum pool size
            pool_max_size: Maximum pool size
            pool_timeout: Pool connection timeout
            max_pool_connections: Maximum number of connections in all pools
            logger: Logger instance
        """
        self.connection_type = connection_type.lower()
        self.connection_params = connection_params or {}
        self.batch_size = batch_size
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.pool_min_size = pool_min_size
        self.pool_max_size = pool_max_size
        self.pool_timeout = pool_timeout
        self.logger = logger or logging.getLogger(__name__)

        # Connection state
        self._connection = None
        self._transaction_count = 0

        # Track active transactions
        self._active_transactions = {}
        self._active_transactions_lock = threading.RLock()
        self._transaction_id_counter = 0

        # Connection pool
        self._pool = ConnectionPool(
            min_size=pool_min_size,
            max_size=pool_max_size,
            timeout=pool_timeout,
            connection_type=connection_type,
            connection_params=connection_params,
            logger=logger
        )

        # Initialize an async lock for transaction management
        self._async_lock = asyncio.Lock()

        # Initialize the connection if needed
        self._init_connection()
    
    def _init_connection(self) -> None:
        """Initialize database connection - primarily for backward compatibility."""
        # The actual connections are managed by the connection pool now
        if self.connection_type == "direct":
            if not PSYCOPG2_AVAILABLE:
                raise ImportError("psycopg2 is required for direct PostgreSQL connections")
        elif self.connection_type == "supabase":
            if not SUPABASE_AVAILABLE:
                raise ImportError("supabase-py is required for Supabase connections")
        elif self.connection_type == "mock":
            # Mock connection for testing without a database
            self.logger.info("Using mock database connection for testing")
            self._connection = MagicMock()
        else:
            raise ValueError(f"Unsupported connection type: {self.connection_type}")

        # Log initialization of the connection pool
        self.logger.info(f"Connection pool initialized with {self.pool_min_size} connections")
    
    @contextlib.contextmanager
    def transaction(self):
        """
        Context manager for database transactions.

        Use this to ensure operations are atomic and rolled back on error.

        Example:
            with db.transaction():
                db.insert_molecule(molecule_data)
                db.insert_properties(molecule_id, properties)
        """
        # Get connection from the pool
        conn, conn_id, is_new = self._pool.get_connection()
        transaction_id = self._get_next_transaction_id()

        # Track the transaction
        with self._active_transactions_lock:
            self._active_transactions[transaction_id] = {
                "connection": conn,
                "connection_id": conn_id,
                "started_at": time.time(),
                "operations": 0
            }
            self._transaction_count += 1

        if self.connection_type == "direct":
            # For direct PostgreSQL connections
            try:
                self.logger.debug(f"Starting transaction {transaction_id}")
                yield transaction_id
                conn.commit()
                self.logger.debug(f"Transaction {transaction_id} committed")
            except Exception as e:
                conn.rollback()
                self.logger.error(f"Transaction {transaction_id} rolled back: {str(e)}")
                raise TransactionError(f"Transaction failed: {str(e)}")
            finally:
                # Clean up transaction tracking
                with self._active_transactions_lock:
                    if transaction_id in self._active_transactions:
                        del self._active_transactions[transaction_id]
                    self._transaction_count -= 1

                # Release connection back to the pool
                self._pool.release_connection(conn, conn_id)
        elif self.connection_type == "supabase":
            # Supabase doesn't have explicit transactions, so we just yield
            try:
                self.logger.debug(f"Starting operations block {transaction_id}")
                yield transaction_id
                self.logger.debug(f"Operations block {transaction_id} completed")
            except Exception as e:
                self.logger.error(f"Operations block {transaction_id} failed: {str(e)}")
                raise TransactionError(f"Operations failed: {str(e)}")
            finally:
                # Clean up transaction tracking
                with self._active_transactions_lock:
                    if transaction_id in self._active_transactions:
                        del self._active_transactions[transaction_id]
                    self._transaction_count -= 1

                # Release connection back to the pool
                self._pool.release_connection(conn, conn_id)
    
    def _get_next_transaction_id(self) -> str:
        """Get the next transaction ID."""
        with self._active_transactions_lock:
            self._transaction_id_counter += 1
            return f"tx-{self._transaction_id_counter}"

    async def begin_transaction_async(self) -> str:
        """Begin a database transaction asynchronously.

        Returns:
            Transaction ID as string
        """
        # Get connection from the pool
        conn, conn_id, is_new = await self._pool.get_connection_async()
        transaction_id = self._get_next_transaction_id()

        # Track the transaction
        async with self._async_lock:
            self._active_transactions[transaction_id] = {
                "connection": conn,
                "connection_id": conn_id,
                "started_at": time.time(),
                "operations": 0
            }
            self._transaction_count += 1

        if self.connection_type == "direct":
            # Begin transaction
            self.logger.debug(f"Starting async transaction {transaction_id}")
        else:
            # Supabase doesn't have explicit transactions
            self.logger.debug(f"Starting async operations block {transaction_id}")

        return transaction_id

    async def commit_transaction_async(self, transaction_id: str) -> bool:
        """Commit a database transaction asynchronously.

        Args:
            transaction_id: Transaction ID returned by begin_transaction_async

        Returns:
            True if commit was successful, False otherwise
        """
        # Get transaction info
        transaction_info = self._active_transactions.get(transaction_id)
        if not transaction_info:
            self.logger.error(f"Transaction {transaction_id} not found")
            return False

        conn = transaction_info["connection"]
        conn_id = transaction_info["connection_id"]

        try:
            if self.connection_type == "direct":
                # Commit transaction
                await asyncio.to_thread(conn.commit)
                self.logger.debug(f"Transaction {transaction_id} committed")
            else:
                # Supabase doesn't have explicit transactions
                self.logger.debug(f"Operations block {transaction_id} completed")

            return True
        except Exception as e:
            self.logger.error(f"Error committing transaction {transaction_id}: {str(e)}")
            return False
        finally:
            # Clean up transaction tracking
            async with self._async_lock:
                if transaction_id in self._active_transactions:
                    del self._active_transactions[transaction_id]
                self._transaction_count -= 1

            # Release connection back to the pool
            await self._pool.release_connection_async(conn, conn_id)

    async def rollback_transaction_async(self, transaction_id: str) -> bool:
        """Rollback a database transaction asynchronously.

        Args:
            transaction_id: Transaction ID returned by begin_transaction_async

        Returns:
            True if rollback was successful, False otherwise
        """
        # Get transaction info
        transaction_info = self._active_transactions.get(transaction_id)
        if not transaction_info:
            self.logger.error(f"Transaction {transaction_id} not found")
            return False

        conn = transaction_info["connection"]
        conn_id = transaction_info["connection_id"]

        try:
            if self.connection_type == "direct":
                # Rollback transaction
                await asyncio.to_thread(conn.rollback)
                self.logger.debug(f"Transaction {transaction_id} rolled back")
            else:
                # Supabase doesn't have explicit transactions
                self.logger.debug(f"Operations block {transaction_id} aborted")

            return True
        except Exception as e:
            self.logger.error(f"Error rolling back transaction {transaction_id}: {str(e)}")
            return False
        finally:
            # Clean up transaction tracking
            async with self._async_lock:
                if transaction_id in self._active_transactions:
                    del self._active_transactions[transaction_id]
                self._transaction_count -= 1

            # Release connection back to the pool
            await self._pool.release_connection_async(conn, conn_id)

    @contextlib.asynccontextmanager
    async def transaction_async(self):
        """Asynchronous context manager for database transactions.

        Example:
            async with db.transaction_async() as tx_id:
                await db.insert_molecule(tx_id, molecule_data)
                await db.insert_properties(tx_id, molecule_id, properties)
        """
        transaction_id = await self.begin_transaction_async()

        try:
            yield transaction_id
            await self.commit_transaction_async(transaction_id)
        except Exception as e:
            await self.rollback_transaction_async(transaction_id)
            self.logger.error(f"Transaction {transaction_id} rolled back: {str(e)}")
            raise TransactionError(f"Transaction failed: {str(e)}")

    async def insert_molecule(self, molecule_data: Dict[str, Any], transaction_id: Optional[str] = None) -> str:
        """
        Insert a molecule into the database.

        Args:
            molecule_data: Dictionary containing molecule data
            transaction_id: Optional transaction ID for transaction context

        Returns:
            ID of the inserted molecule
        """
        molecule_id = molecule_data.get('id') or str(uuid.uuid4())

        # Add standard fields if not present
        data = {
            'id': molecule_id,
            'created_at': time.time(),
            'updated_at': time.time(),
            **molecule_data
        }

        try:
            # Get connection - either from transaction or pool
            if transaction_id and transaction_id in self._active_transactions:
                # Use the transaction's connection
                conn = self._active_transactions[transaction_id]["connection"]
                conn_id = self._active_transactions[transaction_id]["connection_id"]
                need_release = False

                # Update transaction stats
                with self._active_transactions_lock:
                    self._active_transactions[transaction_id]["operations"] += 1
            else:
                # Get a new connection from the pool
                conn, conn_id, _ = await self._pool.get_connection_async()
                need_release = True

            # Execute query based on connection type
            if self.connection_type == "direct":
                # Use psycopg2 connection
                cursor = await asyncio.to_thread(conn.cursor)

                # Build query
                columns = list(data.keys())
                placeholders = [f"%({col})s" for col in columns]

                query = f"""
                INSERT INTO molecules ({', '.join(columns)})
                VALUES ({', '.join(placeholders)})
                RETURNING id
                """

                # Execute query
                await asyncio.to_thread(cursor.execute, query, data)
                result = await asyncio.to_thread(cursor.fetchone)

                # Close cursor
                await asyncio.to_thread(cursor.close)

                # Return result
                inserted_id = result[0] if result else molecule_id
            else:
                # Use Supabase connection
                result = await asyncio.to_thread(
                    lambda: conn.table('molecules').insert(data).execute()
                )
                inserted_id = result.data[0]['id'] if result.data else molecule_id

            # Log the operation
            self.logger.debug(f"Inserted molecule with ID {inserted_id}")

            # Release connection if not part of a transaction
            if need_release:
                await self._pool.release_connection_async(conn, conn_id)

            return inserted_id
        except Exception as e:
            self.logger.error(f"Failed to insert molecule: {str(e)}")

            # Release connection if not part of a transaction and an error occurred
            if 'need_release' in locals() and need_release and 'conn' in locals() and 'conn_id' in locals():
                await self._pool.release_connection_async(conn, conn_id)

            raise QueryError(f"Failed to insert molecule: {str(e)}")
    
    async def insert_properties(
        self,
        molecule_id: str,
        properties: List[Dict[str, Any]],
        transaction_id: Optional[str] = None
    ) -> List[str]:
        """
        Insert molecular properties into the database.
        
        Args:
            molecule_id: ID of the molecule
            properties: List of property dictionaries
            
        Returns:
            List of inserted property IDs
        """
        property_ids = []
        
        try:
            for prop in properties:
                prop_id = str(uuid.uuid4())
                
                # Add standard fields
                data = {
                    'id': prop_id,
                    'molecule_id': molecule_id,
                    'created_at': time.time(),
                    **prop
                }
                
                if self.connection_type == "direct":
                    with self._connection.cursor() as cursor:
                        columns = list(data.keys())
                        placeholders = [f"%({col})s" for col in columns]
                        
                        query = f"""
                        INSERT INTO molecular_properties ({', '.join(columns)})
                        VALUES ({', '.join(placeholders)})
                        RETURNING id
                        """
                        
                        cursor.execute(query, data)
                        result = cursor.fetchone()
                        property_ids.append(result[0] if result else prop_id)
                else:
                    # Supabase
                    result = self._connection.table('molecular_properties').insert(data).execute()
                    property_ids.append(result.data[0]['id'] if result.data else prop_id)
            
            return property_ids
        except Exception as e:
            self.logger.error(f"Failed to insert properties: {str(e)}")
            raise QueryError(f"Failed to insert properties: {str(e)}")
    
    async def insert_synonyms(
        self,
        molecule_id: str,
        synonyms: List[str],
        transaction_id: Optional[str] = None
    ) -> List[str]:
        """
        Insert molecular synonyms into the database.
        
        Args:
            molecule_id: ID of the molecule
            synonyms: List of synonym strings
            
        Returns:
            List of inserted synonym IDs
        """
        synonym_ids = []
        
        try:
            for name in synonyms:
                synonym_id = str(uuid.uuid4())
                
                # Add standard fields
                data = {
                    'id': synonym_id,
                    'molecule_id': molecule_id,
                    'name': name,
                    'created_at': time.time()
                }
                
                if self.connection_type == "direct":
                    with self._connection.cursor() as cursor:
                        columns = list(data.keys())
                        placeholders = [f"%({col})s" for col in columns]
                        
                        query = f"""
                        INSERT INTO synonyms ({', '.join(columns)})
                        VALUES ({', '.join(placeholders)})
                        RETURNING id
                        """
                        
                        cursor.execute(query, data)
                        result = cursor.fetchone()
                        synonym_ids.append(result[0] if result else synonym_id)
                else:
                    # Supabase
                    result = self._connection.table('synonyms').insert(data).execute()
                    synonym_ids.append(result.data[0]['id'] if result.data else synonym_id)
            
            return synonym_ids
        except Exception as e:
            self.logger.error(f"Failed to insert synonyms: {str(e)}")
            raise QueryError(f"Failed to insert synonyms: {str(e)}")
    
    async def molecule_exists(
        self,
        identifier: str,
        identifier_type: str = 'pubchem_cid'
    ) -> Tuple[bool, Optional[str]]:
        """
        Check if a molecule exists in the database.
        
        Args:
            identifier: Molecule identifier (e.g., PubChem CID, ChEMBL ID)
            identifier_type: Type of identifier
            
        Returns:
            Tuple of (exists, molecule_id)
        """
        try:
            if self.connection_type == "direct":
                with self._connection.cursor() as cursor:
                    query = f"""
                    SELECT id FROM molecules
                    WHERE {identifier_type} = %s
                    """
                    
                    cursor.execute(query, (identifier,))
                    result = cursor.fetchone()
                    
                    return (True, result[0]) if result else (False, None)
            else:
                # Supabase
                result = self._connection.table('molecules').select('id').eq(identifier_type, identifier).execute()
                
                return (True, result.data[0]['id']) if result.data else (False, None)
        except Exception as e:
            self.logger.error(f"Failed to check if molecule exists: {str(e)}")
            raise QueryError(f"Failed to check if molecule exists: {str(e)}")
    
    async def execute_query(
        self,
        query: str,
        params: Optional[Dict[str, Any]] = None
    ) -> List[Dict[str, Any]]:
        """
        Execute a custom query on the database.
        
        Args:
            query: SQL query string
            params: Query parameters
            
        Returns:
            List of query result rows
        """
        params = params or {}
        
        try:
            if self.connection_type == "direct":
                with self._connection.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
                    cursor.execute(query, params)
                    
                    if cursor.description:
                        return cursor.fetchall()
                    return []
            else:
                # For Supabase, we can only execute raw queries with RPC
                # This is a simplified approach and may need adjustment
                function_name = "execute_query"
                
                # Ensure the function exists
                self._ensure_rpc_function_exists()
                
                result = self._connection.rpc(
                    function_name,
                    {"query_text": query, "query_params": json.dumps(params)}
                ).execute()
                
                return result.data if result.data else []
        except Exception as e:
            self.logger.error(f"Failed to execute query: {str(e)}")
            raise QueryError(f"Failed to execute query: {str(e)}")
    
    def _ensure_rpc_function_exists(self) -> None:
        """
        Ensure the RPC function for executing queries exists.
        
        This is only needed for Supabase connections.
        """
        if self.connection_type != "supabase":
            return
            
        # Check if function exists and create it if not
        try:
            # This is a simplified check and may need adjustment
            result = self._connection.rpc(
                "function_exists",
                {"function_name": "execute_query"}
            ).execute()
            
            function_exists = result.data and result.data[0].get('exists', False)
            
            if not function_exists:
                # Create the function
                # This would typically be done through migrations
                self.logger.warning(
                    "The execute_query RPC function does not exist. "
                    "Please create it using migrations."
                )
        except Exception as e:
            self.logger.warning(
                f"Could not check if execute_query function exists: {str(e)}"
            )
    
    async def batch_insert(
        self,
        table: str,
        records: List[Dict[str, Any]],
        returning: Optional[str] = 'id',
        transaction_id: Optional[str] = None
    ) -> List[Any]:
        """
        Insert multiple records in a batch operation.
        
        Args:
            table: Table name
            records: List of record dictionaries
            returning: Column to return from inserted records
            
        Returns:
            List of values from the specified returning column
        """
        if not records:
            return []
            
        try:
            if self.connection_type == "direct":
                # For direct PostgreSQL
                results = []
                
                # Get columns from the first record
                columns = list(records[0].keys())
                
                # Split into batches
                for i in range(0, len(records), self.batch_size):
                    batch = records[i:i + self.batch_size]
                    
                    with self._connection.cursor() as cursor:
                        # Prepare values part of the query
                        values_template = []
                        values_params = []
                        
                        for record in batch:
                            record_placeholders = []
                            
                            for col in columns:
                                values_params.append(record.get(col))
                                record_placeholders.append('%s')
                                
                            values_template.append(f"({', '.join(record_placeholders)})")
                            
                        query = f"""
                        INSERT INTO {table} ({', '.join(columns)})
                        VALUES {', '.join(values_template)}
                        """
                        
                        if returning:
                            query += f" RETURNING {returning}"
                            
                        cursor.execute(query, values_params)
                        
                        if returning:
                            batch_results = cursor.fetchall()
                            results.extend([row[0] for row in batch_results])
                
                return results
            else:
                # For Supabase
                results = []
                
                # Split into batches
                for i in range(0, len(records), self.batch_size):
                    batch = records[i:i + self.batch_size]
                    
                    result = self._connection.table(table).insert(batch).execute()
                    
                    if returning and result.data:
                        results.extend([row[returning] for row in result.data])
                
                return results
        except Exception as e:
            self.logger.error(f"Failed to batch insert: {str(e)}")
            raise QueryError(f"Failed to batch insert: {str(e)}")
    
    async def with_retry(
        self,
        operation: Callable,
        *args,
        **kwargs
    ) -> Any:
        """
        Execute an operation with retry logic.
        
        Args:
            operation: Function to execute
            *args: Arguments for the operation
            **kwargs: Keyword arguments for the operation
            
        Returns:
            Result from the operation
        """
        last_error = None
        
        for attempt in range(1, self.max_retries + 1):
            try:
                return await operation(*args, **kwargs)
            except (ConnectionError, QueryError, TransactionError) as e:
                last_error = e
                
                if attempt < self.max_retries:
                    delay = self.retry_delay * (2 ** (attempt - 1))  # Exponential backoff
                    self.logger.warning(
                        f"Operation failed, retrying in {delay:.2f}s "
                        f"(attempt {attempt}/{self.max_retries}): {str(e)}"
                    )
                    await asyncio.sleep(delay)
                else:
                    self.logger.error(
                        f"Operation failed after {self.max_retries} attempts: {str(e)}"
                    )
        
        # If we get here, all retries failed
        raise last_error