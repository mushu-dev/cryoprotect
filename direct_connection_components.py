"""
Direct Database Connection Components

Key code snippets for implementing direct PostgreSQL connections to Supabase.
These components are designed to be used as building blocks for the Roo agents
to assemble complete solutions.
"""

# 1. Basic connection setup with DNS fallback
def _resolve_host(host):
    """
    Resolve hostname to IP address to work around DNS issues.
    Falls back to the original hostname if resolution fails.
    """
    import socket
    import logging
    
    try:
        # Try to resolve the hostname to an IP address
        ip_address = socket.gethostbyname(host)
        logging.info(f"Resolved {host} to IP: {ip_address}")
        return ip_address
    except socket.gaierror as e:
        # If resolution fails, log the error and use the original hostname
        logging.warning(f"Could not resolve hostname {host}: {str(e)}")
        return host

# 2. Connection pool configuration
def initialize_connection_pool(min_connections=1, max_connections=10):
    """
    Initialize the PostgreSQL connection pool with proper configuration.
    """
    import os
    import psycopg2.pool
    from dotenv import load_dotenv
    
    load_dotenv()
    
    host = os.getenv('SUPABASE_DB_HOST')
    port = os.getenv('SUPABASE_DB_PORT', '5432')
    dbname = os.getenv('SUPABASE_DB_NAME', 'postgres')
    user = os.getenv('SUPABASE_DB_USER', 'postgres')
    password = os.getenv('SUPABASE_DB_PASSWORD')
    
    # Resolve hostname to IP if needed
    resolved_host = _resolve_host(host)
    
    # Create connection parameters
    conn_params = {
        'host': resolved_host,
        'port': port,
        'dbname': dbname,
        'user': user,
        'password': password,
        'connect_timeout': 10
    }
    
    # Create the connection pool
    pool = psycopg2.pool.ThreadedConnectionPool(
        min_connections,
        max_connections,
        **conn_params
    )
    
    return pool

# 3. Context manager for connection handling
def get_connection_context(pool):
    """Context manager for getting and returning a connection from the pool."""
    from contextlib import contextmanager
    
    @contextmanager
    def _get_connection():
        """Get a connection from the pool with automatic return."""
        connection = None
        try:
            connection = pool.getconn()
            yield connection
        finally:
            if connection:
                pool.putconn(connection)
    
    return _get_connection

# 4. Optimized bulk insert function
def bulk_insert(pool, table, data, columns=None, return_ids=False):
    """
    Efficiently insert multiple rows into a table.
    
    Args:
        pool: Connection pool
        table: Table name
        data: List of dictionaries with column:value pairs
        columns: Specific columns to insert (if None, use all keys from first data item)
        return_ids: If True, return the inserted IDs
    """
    import psycopg2.extras
    
    if not data:
        return None
    
    # If columns not specified, use keys from first data item
    if columns is None:
        columns = list(data[0].keys())
    
    # Prepare values and parameters
    values_placeholder = ', '.join(['%s'] * len(columns))
    columns_str = ', '.join([f'"{col}"' for col in columns])
    
    # Create the base query
    query = f'INSERT INTO {table} ({columns_str}) VALUES '
    
    # Prepare values for each row
    rows_values = []
    all_params = []
    
    for row in data:
        row_params = [row.get(col) for col in columns]
        all_params.extend(row_params)
        rows_values.append(f'({values_placeholder})')
    
    # Complete the query with all rows
    query += ', '.join(rows_values)
    
    # Add returning clause if needed
    if return_ids:
        query += ' RETURNING id'
    
    # Execute the query
    with pool.getconn() as conn:
        try:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
                cursor.execute(query, all_params)
                conn.commit()
                
                if return_ids and cursor.description:
                    return [row['id'] for row in cursor.fetchall()]
        finally:
            pool.putconn(conn)
    
    return None

# 5. Transaction management with context manager
def transaction_context(pool):
    """Context manager for database transactions."""
    from contextlib import contextmanager
    
    @contextmanager
    def _transaction():
        conn = None
        try:
            conn = pool.getconn()
            yield conn
            conn.commit()
        except Exception:
            if conn:
                conn.rollback()
            raise
        finally:
            if conn:
                pool.putconn(conn)
    
    return _transaction

# 6. Batch processing for large datasets
def process_in_batches(items, batch_size=1000, process_func=None):
    """
    Process a large list of items in batches.
    
    Args:
        items: List of items to process
        batch_size: Number of items to process in each batch
        process_func: Function to call for each batch
    """
    import logging
    
    total_batches = (len(items) + batch_size - 1) // batch_size
    logging.info(f"Processing {len(items)} items in {total_batches} batches")
    
    results = []
    
    for i in range(0, len(items), batch_size):
        batch = items[i:i+batch_size]
        batch_num = i // batch_size + 1
        logging.info(f"Processing batch {batch_num}/{total_batches}: {len(batch)} items")
        
        if process_func:
            result = process_func(batch)
            if result is not None:
                results.append(result)
    
    return results

# 7. Retry logic for resilient database operations
def with_retry(func, max_retries=3, retry_delay=2):
    """
    Execute a function with retry logic.
    
    Args:
        func: Function to execute
        max_retries: Maximum number of retry attempts
        retry_delay: Delay between retries in seconds
    """
    import time
    import logging
    
    def wrapper(*args, **kwargs):
        attempt = 0
        while attempt < max_retries:
            try:
                return func(*args, **kwargs)
            except Exception as e:
                attempt += 1
                if attempt < max_retries:
                    logging.warning(f"Operation failed (attempt {attempt}/{max_retries}): {str(e)}")
                    time.sleep(retry_delay)
                else:
                    logging.error(f"Operation failed after {max_retries} attempts: {str(e)}")
                    raise
    
    return wrapper