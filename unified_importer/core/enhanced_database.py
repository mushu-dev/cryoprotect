"""
Enhanced database operations with improved connection pooling and error handling.

This module extends the unified database operations interface with advanced
features for reliability, performance optimization, and monitoring.
"""

import os
import time
import asyncio
import logging
import uuid
import json
import contextlib
from typing import Dict, List, Any, Optional, Tuple, Union, Callable, AsyncIterator

# Import the base database operations
from .database import DatabaseOperations, ConnectionError, TransactionError, QueryError

# Import the enhanced connection pool
from .enhanced_connection_pool import EnhancedConnectionPool


class EnhancedDatabaseOperations(DatabaseOperations):
    """
    Enhanced database operations with improved connection pooling and error handling.
    
    This class extends the base DatabaseOperations with advanced features:
    - Enhanced connection pooling with health checking
    - Dynamic pool sizing based on workload
    - Circuit breaker pattern for error handling
    - Detailed metrics and monitoring
    - Connection validation
    """
    
    def __init__(
        self,
        connection_type: str = "direct",
        connection_params: Optional[Dict[str, Any]] = None,
        batch_size: int = 100,
        max_retries: int = 3,
        retry_delay: float = 2.0,
        pool_min_size: int = 5,
        pool_max_size: int = 20,
        pool_target_utilization: float = 0.7,
        pool_timeout: float = 30.0,
        pool_max_lifetime: float = 3600.0,
        pool_validation_interval: float = 300.0,
        pool_health_check_interval: float = 60.0,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize enhanced database operations.
        
        Args:
            connection_type: Type of connection ("direct", "supabase")
            connection_params: Connection parameters
            batch_size: Number of items to include in batch operations
            max_retries: Maximum number of retry attempts
            retry_delay: Delay between retries (in seconds)
            pool_min_size: Minimum pool size
            pool_max_size: Maximum pool size
            pool_target_utilization: Target pool utilization (0.0-1.0)
            pool_timeout: Pool connection timeout
            pool_max_lifetime: Maximum connection lifetime
            pool_validation_interval: Interval between connection validations
            pool_health_check_interval: Interval between pool health checks
            logger: Logger instance
        """
        # Initialize with minimal configuration
        super().__init__(
            connection_type=connection_type,
            connection_params=connection_params,
            batch_size=batch_size,
            max_retries=max_retries,
            retry_delay=retry_delay,
            pool_min_size=pool_min_size,
            pool_max_size=pool_max_size,
            pool_timeout=pool_timeout,
            logger=logger
        )
        
        # Replace the basic connection pool with the enhanced one
        self._pool = EnhancedConnectionPool(
            min_size=pool_min_size,
            max_size=pool_max_size,
            target_utilization=pool_target_utilization,
            timeout=pool_timeout,
            max_lifetime=pool_max_lifetime,
            validation_interval=pool_validation_interval,
            health_check_interval=pool_health_check_interval,
            connection_type=connection_type,
            connection_params=connection_params,
            logger=logger
        )
        
        self.logger.info(
            f"Enhanced database operations initialized with "
            f"connection type: {connection_type}, "
            f"pool size: {pool_min_size}-{pool_max_size}"
        )
    
    def get_pool_stats(self) -> Dict[str, Any]:
        """
        Get detailed statistics about the connection pool.
        
        Returns:
            Dictionary of pool statistics
        """
        return self._pool.get_stats()
    
    async def get_pool_stats_async(self) -> Dict[str, Any]:
        """
        Get detailed statistics about the connection pool asynchronously.
        
        Returns:
            Dictionary of pool statistics
        """
        return await self._pool.get_stats_async()
    
    def run_health_check(self) -> None:
        """Run a manual health check on the connection pool."""
        self._pool._run_health_check()
    
    async def run_health_check_async(self) -> Dict[str, Any]:
        """
        Run a manual health check on the connection pool asynchronously.
        
        Returns:
            Dictionary with health check results
        """
        return await self._pool.run_health_check_async()
    
    @contextlib.contextmanager
    def transaction(self):
        """
        Context manager for database transactions with enhanced error handling.
        
        Use this to ensure operations are atomic and rolled back on error.
        
        Example:
            with db.transaction() as tx_id:
                db.insert_molecule(molecule_data, transaction_id=tx_id)
                db.insert_properties(molecule_id, properties, transaction_id=tx_id)
        """
        # Get connection from the pool with enhanced error handling
        conn, conn_id, is_new = self._pool.get_connection()
        transaction_id = self._get_next_transaction_id()
        
        # Track the transaction
        with self._active_transactions_lock:
            self._active_transactions[transaction_id] = {
                "connection": conn,
                "connection_id": conn_id,
                "started_at": time.time(),
                "operations": 0,
                "is_new_connection": is_new
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
                
                # Track error in connection info
                _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
                if conn_info:
                    conn_info.error_count += 1
                    conn_info.consecutive_errors += 1
                
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
                
                # Track error in connection info
                _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
                if conn_info:
                    conn_info.error_count += 1
                    conn_info.consecutive_errors += 1
                
                raise TransactionError(f"Operations failed: {str(e)}")
            finally:
                # Clean up transaction tracking
                with self._active_transactions_lock:
                    if transaction_id in self._active_transactions:
                        del self._active_transactions[transaction_id]
                    self._transaction_count -= 1
                
                # Release connection back to the pool
                self._pool.release_connection(conn, conn_id)
    
    @contextlib.asynccontextmanager
    async def transaction_async(self):
        """
        Asynchronous context manager for database transactions with enhanced error handling.
        
        Example:
            async with db.transaction_async() as tx_id:
                await db.insert_molecule(molecule_data, transaction_id=tx_id)
                await db.insert_properties(molecule_id, properties, transaction_id=tx_id)
        """
        transaction_id = await self.begin_transaction_async()
        
        try:
            yield transaction_id
            await self.commit_transaction_async(transaction_id)
        except Exception as e:
            await self.rollback_transaction_async(transaction_id)
            self.logger.error(f"Transaction {transaction_id} rolled back: {str(e)}")
            
            # Error already tracked in rollback method
            raise TransactionError(f"Transaction failed: {str(e)}")
    
    async def begin_transaction_async(self) -> str:
        """
        Begin a database transaction asynchronously with enhanced error handling.
        
        Returns:
            Transaction ID as string
        """
        # Get connection from the pool with enhanced error handling
        conn, conn_id, is_new = await self._pool.get_connection_async()
        transaction_id = self._get_next_transaction_id()
        
        # Track the transaction
        async with self._async_lock:
            self._active_transactions[transaction_id] = {
                "connection": conn,
                "connection_id": conn_id,
                "started_at": time.time(),
                "operations": 0,
                "is_new_connection": is_new
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
        """
        Commit a database transaction asynchronously with enhanced error handling.
        
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
            
            # Reset consecutive errors on success
            _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
            if conn_info:
                conn_info.consecutive_errors = 0
            
            return True
        except Exception as e:
            self.logger.error(f"Error committing transaction {transaction_id}: {str(e)}")
            
            # Track error in connection info
            _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
            if conn_info:
                conn_info.error_count += 1
                conn_info.consecutive_errors += 1
            
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
        """
        Rollback a database transaction asynchronously with enhanced error handling.
        
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
            
            # Track error in connection info
            _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
            if conn_info:
                conn_info.error_count += 1
                conn_info.consecutive_errors += 1
            
            return True
        except Exception as e:
            self.logger.error(f"Error rolling back transaction {transaction_id}: {str(e)}")
            
            # Still track error in connection info
            _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
            if conn_info:
                conn_info.error_count += 1
                conn_info.consecutive_errors += 1
            
            return False
        finally:
            # Clean up transaction tracking
            async with self._async_lock:
                if transaction_id in self._active_transactions:
                    del self._active_transactions[transaction_id]
                self._transaction_count -= 1
            
            # Release connection back to the pool
            await self._pool.release_connection_async(conn, conn_id)
    
    async def insert_molecule(
        self,
        molecule_data: Dict[str, Any],
        transaction_id: Optional[str] = None
    ) -> str:
        """
        Insert a molecule into the database with enhanced error handling.
        
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
            
            # Track successful query in connection info
            _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
            if conn_info:
                conn_info.query_count += 1
                conn_info.consecutive_errors = 0
            
            # Release connection if not part of a transaction
            if need_release:
                await self._pool.release_connection_async(conn, conn_id)
            
            return inserted_id
        except Exception as e:
            self.logger.error(f"Failed to insert molecule: {str(e)}")
            
            # Track error in connection info
            _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
            if conn_info:
                conn_info.error_count += 1
                conn_info.consecutive_errors += 1
            
            # Release connection if not part of a transaction and an error occurred
            if 'need_release' in locals() and need_release and 'conn' in locals() and 'conn_id' in locals():
                await self._pool.release_connection_async(conn, conn_id)
            
            raise QueryError(f"Failed to insert molecule: {str(e)}")
    
    async def batch_insert(
        self,
        table: str,
        records: List[Dict[str, Any]],
        returning: Optional[str] = 'id',
        transaction_id: Optional[str] = None
    ) -> List[Any]:
        """
        Insert multiple records in a batch operation with enhanced error handling.
        
        Args:
            table: Table name
            records: List of record dictionaries
            returning: Column to return from inserted records
            transaction_id: Optional transaction ID for transaction context
            
        Returns:
            List of values from the specified returning column
        """
        if not records:
            return []
        
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
            
            results = []
            
            if self.connection_type == "direct":
                # For direct PostgreSQL
                # Get columns from the first record
                columns = list(records[0].keys())
                
                # Split into batches
                for i in range(0, len(records), self.batch_size):
                    batch = records[i:i + self.batch_size]
                    
                    # Execute in a thread to avoid blocking
                    batch_results = await asyncio.to_thread(self._execute_batch_insert_direct, 
                                                           conn, table, columns, batch, returning)
                    results.extend(batch_results)
            else:
                # For Supabase
                # Split into batches
                for i in range(0, len(records), self.batch_size):
                    batch = records[i:i + self.batch_size]
                    
                    # Execute in a thread to avoid blocking
                    result = await asyncio.to_thread(
                        lambda: conn.table(table).insert(batch).execute()
                    )
                    
                    if returning and result.data:
                        results.extend([row[returning] for row in result.data])
            
            # Track successful query in connection info
            _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
            if conn_info:
                conn_info.query_count += len(records)
                conn_info.consecutive_errors = 0
            
            # Release connection if not part of a transaction
            if need_release:
                await self._pool.release_connection_async(conn, conn_id)
            
            return results
        except Exception as e:
            self.logger.error(f"Failed to batch insert: {str(e)}")
            
            # Track error in connection info
            _, conn_info = self._pool.all_connections.get(conn_id, (None, None))
            if conn_info:
                conn_info.error_count += 1
                conn_info.consecutive_errors += 1
            
            # Release connection if not part of a transaction and an error occurred
            if 'need_release' in locals() and need_release and 'conn' in locals() and 'conn_id' in locals():
                await self._pool.release_connection_async(conn, conn_id)
            
            raise QueryError(f"Failed to batch insert: {str(e)}")
    
    def _execute_batch_insert_direct(
        self,
        conn: Any,
        table: str,
        columns: List[str],
        batch: List[Dict[str, Any]],
        returning: Optional[str]
    ) -> List[Any]:
        """
        Execute a batch insert on a direct PostgreSQL connection.
        
        Args:
            conn: Database connection
            table: Table name
            columns: Column names
            batch: Batch of records
            returning: Column to return
            
        Returns:
            List of values from the specified returning column
        """
        with conn.cursor() as cursor:
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
                return [row[0] for row in batch_results]
            
            return []
    
    def close(self) -> None:
        """Close all database connections."""
        self._pool.close_all()
    
    async def close_async(self) -> None:
        """Close all database connections asynchronously."""
        await self._pool.close_all_async()
    
    def health_check(self) -> Dict[str, Any]:
        """
        Run a health check on the database.
        
        Returns:
            Dictionary with health check results
        """
        pool_stats = self.get_pool_stats()
        
        # Run a simple query to verify database connectivity
        try:
            with self.transaction() as tx_id:
                # Use the transaction's connection
                conn = self._active_transactions[tx_id]["connection"]
                
                if self.connection_type == "direct":
                    with conn.cursor() as cursor:
                        cursor.execute("SELECT 1 as health_check")
                        result = cursor.fetchone()
                        db_healthy = result and result[0] == 1
                else:
                    # For Supabase, just check that the connection exists
                    db_healthy = conn is not None
        except Exception as e:
            self.logger.error(f"Database health check failed: {str(e)}")
            db_healthy = False
        
        # Collect health data
        health_data = {
            "database_healthy": db_healthy,
            "pool_healthy": pool_stats["idle_connections"] > 0,
            "connection_errors": pool_stats["errors"],
            "pool_size": pool_stats["pool_size"],
            "idle_connections": pool_stats["idle_connections"],
            "active_connections": pool_stats["active_connections"],
            "circuit_breaker_state": pool_stats.get("circuit_breaker_state", "unknown"),
            "validation_failures": pool_stats.get("validation_failures", 0)
        }
        
        # Determine overall health
        health_data["healthy"] = health_data["database_healthy"] and health_data["pool_healthy"]
        
        return health_data
    
    async def health_check_async(self) -> Dict[str, Any]:
        """
        Run a health check on the database asynchronously.
        
        Returns:
            Dictionary with health check results
        """
        pool_stats = await self.get_pool_stats_async()
        
        # Run a simple query to verify database connectivity
        try:
            async with self.transaction_async() as tx_id:
                # Use the transaction's connection
                conn = self._active_transactions[tx_id]["connection"]
                
                if self.connection_type == "direct":
                    cursor = await asyncio.to_thread(conn.cursor)
                    await asyncio.to_thread(cursor.execute, "SELECT 1 as health_check")
                    result = await asyncio.to_thread(cursor.fetchone)
                    db_healthy = result and result[0] == 1
                else:
                    # For Supabase, just check that the connection exists
                    db_healthy = conn is not None
        except Exception as e:
            self.logger.error(f"Database health check failed: {str(e)}")
            db_healthy = False
        
        # Collect health data
        health_data = {
            "database_healthy": db_healthy,
            "pool_healthy": pool_stats["idle_connections"] > 0,
            "connection_errors": pool_stats["errors"],
            "pool_size": pool_stats["pool_size"],
            "idle_connections": pool_stats["idle_connections"],
            "active_connections": pool_stats["active_connections"],
            "circuit_breaker_state": pool_stats.get("circuit_breaker_state", "unknown"),
            "validation_failures": pool_stats.get("validation_failures", 0)
        }
        
        # Determine overall health
        health_data["healthy"] = health_data["database_healthy"] and health_data["pool_healthy"]
        
        return health_data