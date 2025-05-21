# Direct Supabase Connection Plan for ChEMBL Integration

## Overview

This document outlines the implementation of a direct Supabase connection for SQL execution to replace the slower MCP approach. This will significantly improve database operation performance during the ChEMBL data integration process.

## Implementation Strategy

### 1. Direct Connection Components

**Connection Module**
- File: `supabase_direct.py`
- Purpose: Establish and manage direct PostgreSQL connections to Supabase
- Key Features:
  - Connection pooling
  - Retry logic
  - Transaction management
  - Service role authentication

**SQL Execution Utilities**
- File: `sql_executor.py`
- Purpose: Provide optimized SQL execution functions
- Key Features:
  - Batch execution
  - Parameterized queries (prevent SQL injection)
  - Transaction management
  - Result processing

### 2. Database Configuration

**Environment Setup**
- Store Supabase PostgreSQL connection details in environment variables:
  - `SUPABASE_DB_HOST` - Supabase PostgreSQL host
  - `SUPABASE_DB_PORT` - PostgreSQL port (typically 5432 or 6543)
  - `SUPABASE_DB_NAME` - Database name
  - `SUPABASE_DB_USER` - Service role user
  - `SUPABASE_DB_PASSWORD` - Service role password
  - `SUPABASE_DB_MAX_CONNECTIONS` - Connection pool size (default: 10)

**Service Role Authentication**
- Use service role credentials for direct database access
- Implement proper security measures for credential management

### 3. Connection Pooling Architecture

**Pool Configuration**
- Default pool size: 10 connections (configurable)
- Connection timeout: 30 seconds
- Connection recycling: Every 10 minutes
- Idle timeout: 60 seconds

**Pool Management**
- Automatic connection health checks
- Connection reset on errors
- Graceful shutdown process

## Implementation Details

### Direct Connection Module

```python
# supabase_direct.py

import os
import time
import logging
import threading
from typing import Dict, Any, Optional, List, Union, Tuple
import psycopg2
from psycopg2 import pool
from psycopg2.extras import RealDictCursor

logger = logging.getLogger(__name__)

class SupabaseDirectConnection:
    """Direct PostgreSQL connection manager for Supabase."""
    
    _instance = None
    _lock = threading.Lock()
    
    @classmethod
    def get_instance(cls):
        """Singleton pattern to ensure one connection pool per application."""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance
    
    def __init__(self):
        """Initialize connection pool."""
        # Extract connection parameters from environment
        self.db_params = {
            'host': os.getenv('SUPABASE_DB_HOST'),
            'port': os.getenv('SUPABASE_DB_PORT', '5432'),
            'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
            'user': os.getenv('SUPABASE_DB_USER'),
            'password': os.getenv('SUPABASE_DB_PASSWORD')
        }
        
        max_connections = int(os.getenv('SUPABASE_DB_MAX_CONNECTIONS', '10'))
        
        # Validate required parameters
        missing = [k for k, v in self.db_params.items() if not v]
        if missing:
            raise ValueError(f"Missing required database parameters: {', '.join(missing)}")
        
        # Create connection pool
        self.connection_pool = psycopg2.pool.ThreadedConnectionPool(
            minconn=1,
            maxconn=max_connections,
            **self.db_params
        )
        
        # Set up connection management
        self.active_connections = 0
        self._lock = threading.Lock()
        
        logger.info(f"Initialized Supabase direct connection pool with {max_connections} max connections")
    
    def get_connection(self):
        """Get a connection from the pool."""
        with self._lock:
            try:
                conn = self.connection_pool.getconn()
                self.active_connections += 1
                logger.debug(f"Acquired connection (active: {self.active_connections})")
                return conn
            except Exception as e:
                logger.error(f"Error getting database connection: {str(e)}")
                raise
    
    def release_connection(self, conn):
        """Return a connection to the pool."""
        with self._lock:
            try:
                self.connection_pool.putconn(conn)
                self.active_connections -= 1
                logger.debug(f"Released connection (active: {self.active_connections})")
            except Exception as e:
                logger.error(f"Error releasing database connection: {str(e)}")
                raise
    
    def execute_query(self, query: str, params: Optional[tuple] = None, 
                     fetch_one: bool = False) -> Union[List[Dict[str, Any]], Dict[str, Any], None]:
        """
        Execute SQL query with proper connection management.
        
        Args:
            query: SQL query to execute
            params: Query parameters for parameterized query
            fetch_one: If True, return only the first result
            
        Returns:
            Query results as dictionaries, or None for non-SELECT queries
        """
        conn = None
        try:
            conn = self.get_connection()
            with conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute(query, params)
                
                # Check if query returns results
                if query.strip().upper().startswith('SELECT') or 'RETURNING' in query.upper():
                    if fetch_one:
                        result = cursor.fetchone()
                    else:
                        result = cursor.fetchall()
                    return result
                else:
                    # For non-SELECT queries like INSERT/UPDATE/DELETE
                    conn.commit()
                    return None
        except Exception as e:
            if conn:
                conn.rollback()
            logger.error(f"Database error executing query: {str(e)}")
            raise
        finally:
            if conn:
                self.release_connection(conn)
    
    def execute_batch(self, queries: List[str], transaction: bool = True) -> List[Any]:
        """
        Execute multiple SQL queries, optionally in a transaction.
        
        Args:
            queries: List of SQL queries to execute
            transaction: If True, execute as a single transaction
            
        Returns:
            List of results for each query
        """
        conn = None
        results = []
        
        try:
            conn = self.get_connection()
            
            # Start transaction if requested
            if transaction:
                conn.autocommit = False
            else:
                conn.autocommit = True
                
            with conn.cursor(cursor_factory=RealDictCursor) as cursor:
                for query in queries:
                    cursor.execute(query)
                    
                    # Only collect results for SELECT queries
                    if query.strip().upper().startswith('SELECT') or 'RETURNING' in query.upper():
                        results.append(cursor.fetchall())
                    else:
                        results.append(None)
                
                # Commit transaction if applicable
                if transaction:
                    conn.commit()
                    
            return results
        except Exception as e:
            if conn and transaction:
                conn.rollback()
            logger.error(f"Database error executing batch: {str(e)}")
            raise
        finally:
            if conn:
                self.release_connection(conn)
    
    def close_all(self):
        """Close all connections in the pool."""
        if hasattr(self, 'connection_pool'):
            self.connection_pool.closeall()
            logger.info("Closed all database connections")
```

### SQL Execution Utility

```python
# sql_executor.py

import json
import logging
from typing import Dict, Any, Optional, List, Union, Tuple
from supabase_direct import SupabaseDirectConnection

logger = logging.getLogger(__name__)

def execute_sql(query: str, params: Optional[tuple] = None) -> Union[List[Dict[str, Any]], Dict[str, Any], None]:
    """
    Execute a SQL query using direct Supabase connection.
    
    Args:
        query: SQL query to execute
        params: Query parameters for parameterized query
        
    Returns:
        Query results as dictionaries, or None for non-SELECT queries
    """
    db = SupabaseDirectConnection.get_instance()
    return db.execute_query(query, params)

def execute_transaction(queries: List[str]) -> List[Any]:
    """
    Execute multiple SQL queries as a single transaction.
    
    Args:
        queries: List of SQL queries to execute
        
    Returns:
        List of results for each query
    """
    db = SupabaseDirectConnection.get_instance()
    return db.execute_batch(queries, transaction=True)

def batch_insert(table: str, data_list: List[Dict[str, Any]], 
                 returning: bool = True) -> Optional[List[Dict[str, Any]]]:
    """
    Insert multiple rows into a table efficiently.
    
    Args:
        table: Target table name
        data_list: List of dictionaries with column values
        returning: Whether to return inserted IDs
        
    Returns:
        List of inserted records with IDs if returning=True
    """
    if not data_list:
        return []
    
    db = SupabaseDirectConnection.get_instance()
    
    # Get column names from first data item
    columns = list(data_list[0].keys())
    
    # Build values statement for each data item
    values_parts = []
    for data in data_list:
        # Convert each value to SQL string representation
        values = []
        for col in columns:
            val = data.get(col)
            if val is None:
                values.append("NULL")
            elif isinstance(val, (int, float)):
                values.append(str(val))
            elif isinstance(val, bool):
                values.append("TRUE" if val else "FALSE")
            elif isinstance(val, (dict, list)):
                # Escape and encode JSON objects
                json_str = json.dumps(val).replace("'", "''")
                values.append(f"'{json_str}'")
            else:
                # Escape string values
                str_val = str(val).replace("'", "''")
                values.append(f"'{str_val}'")
        
        values_parts.append(f"({', '.join(values)})")
    
    # Build the complete query
    query = f"INSERT INTO {table} ({', '.join(columns)})\nVALUES\n"
    query += ",\n".join(values_parts)
    
    if returning:
        query += "\nRETURNING id"
    
    # Execute query
    try:
        result = db.execute_query(query)
        return result
    except Exception as e:
        logger.error(f"Error in batch insert: {str(e)}")
        raise
```

### ChEMBL Integration Updates

Update the ChEMBL integration script to use the direct connection instead of MCP:

```python
# ChEMBL_Integrated_Import.py

# Replace MCP imports
# from use_mcp_tool import execute_sql, get_project_id
from sql_executor import execute_sql, execute_transaction, batch_insert

# Then replace all MCP-based SQL execution calls
# For example, replace:
# result = execute_sql(sql, get_project_id_for_mcp())
# With:
# result = execute_sql(sql)

# And replace batch operations:
# transaction_queries = ["BEGIN;", ..., "COMMIT;"]
# combined_query = "\n".join(transaction_queries)
# return execute_sql_through_mcp(combined_query)
# With:
# return execute_transaction(transaction_queries)
```

## Performance Optimizations

### 1. Bulk Insert/Update

- Use `COPY` command for very large datasets
- Use multi-row inserts for batches (100-1000 rows)
- Use `ON CONFLICT` instead of separate existence checks

### 2. Transaction Management

- Group related operations into transactions
- Use appropriate transaction isolation levels
- Minimize transaction duration

### 3. Connection Management

- Implement connection pooling
- Reuse connections for sequential operations
- Clean shutdown of connections

### 4. Query Optimization

- Use indexed columns in WHERE clauses
- Minimize nested queries
- Use EXPLAIN ANALYZE to identify bottlenecks

## Security Considerations

### 1. Credential Management

- Store credentials in environment variables
- Never hardcode passwords
- Rotate credentials periodically

### 2. SQL Injection Prevention

- Use parameterized queries
- Validate and sanitize all inputs
- Escape string values properly

### 3. Service Role Permissions

- Use least privilege principle
- Create specific service role for data imports
- Audit service role operations

## Testing Plan

### 1. Connection Testing

- Test connection establishment
- Test connection pooling
- Test connection recovery after network issues

### 2. Performance Benchmarking

- Compare MCP vs direct connection performance
- Measure throughput (operations/second)
- Measure latency for different operation types

### 3. Integration Testing

- Test ChEMBL data import with direct connection
- Verify data consistency and integrity
- Test failover and recovery scenarios

## Implementation Steps

1. Set up environment variables for database connection
2. Implement `supabase_direct.py` for connection management
3. Implement `sql_executor.py` for optimized SQL operations
4. Update ChEMBL integration script to use direct connection
5. Test and benchmark performance
6. Deploy and monitor

## Troubleshooting

### Common Issues

1. **Connection Errors**
   - Verify network connectivity
   - Check credentials and connection parameters
   - Verify firewall and security group settings

2. **Performance Issues**
   - Check transaction isolation levels
   - Monitor connection pool utilization
   - Analyze query plans for bottlenecks

3. **Security Restrictions**
   - Verify service role permissions
   - Check RLS policies and their impact