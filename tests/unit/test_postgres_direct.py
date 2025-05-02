"""
Unit tests for PostgreSQL direct connection and SQL executor utilities.
"""

import os
import unittest
from unittest.mock import patch, MagicMock, call
import psycopg2
import psycopg2.pool
import psycopg2.extras

# Import the modules to test
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from postgres_direct import PostgresDirectConnection
import sql_executor

class TestPostgresDirectConnection(unittest.TestCase):
    """Test cases for the PostgresDirectConnection class."""
    
    @patch('postgres_direct.psycopg2.pool.ThreadedConnectionPool')
    @patch('postgres_direct.load_dotenv')
    @patch('postgres_direct.os.getenv')
    def test_singleton_pattern(self, mock_getenv, mock_load_dotenv, mock_pool):
        """Test that PostgresDirectConnection follows the singleton pattern."""
        # Mock environment variables
        mock_getenv.side_effect = lambda key, default=None: {
            'SUPABASE_DB_HOST': 'test-host',
            'SUPABASE_DB_PORT': '5432',
            'SUPABASE_DB_NAME': 'test-db',
            'SUPABASE_DB_USER': 'test-user',
            'SUPABASE_DB_PASSWORD': 'test-password'
        }.get(key, default)
        
        # Create two instances
        instance1 = PostgresDirectConnection()
        instance2 = PostgresDirectConnection()
        
        # Verify they are the same instance
        self.assertIs(instance1, instance2)
        
        # Verify the pool was only created once
        mock_pool.assert_called_once()

    @patch('postgres_direct.psycopg2.pool.ThreadedConnectionPool')
    @patch('postgres_direct.load_dotenv')
    @patch('postgres_direct.os.getenv')
    @patch('postgres_direct.socket.gethostbyname')
    def test_host_resolution(self, mock_gethostbyname, mock_getenv, mock_load_dotenv, mock_pool):
        """Test that hostname resolution works correctly."""
        # Mock environment variables
        mock_getenv.side_effect = lambda key, default=None: {
            'SUPABASE_DB_HOST': 'test-host',
            'SUPABASE_DB_PORT': '5432',
            'SUPABASE_DB_NAME': 'test-db',
            'SUPABASE_DB_USER': 'test-user',
            'SUPABASE_DB_PASSWORD': 'test-password'
        }.get(key, default)
        
        # Mock hostname resolution
        mock_gethostbyname.return_value = '192.168.1.1'
        
        # Create an instance
        instance = PostgresDirectConnection()
        
        # Verify hostname resolution was called
        mock_gethostbyname.assert_called_once_with('test-host')
        
        # Verify the pool was created with the resolved IP
        mock_pool.assert_called_once()
        args, kwargs = mock_pool.call_args
        self.assertEqual(kwargs['host'], '192.168.1.1')

    @patch('postgres_direct.psycopg2.pool.ThreadedConnectionPool')
    @patch('postgres_direct.load_dotenv')
    @patch('postgres_direct.os.getenv')
    def test_connection_context_manager(self, mock_getenv, mock_load_dotenv, mock_pool):
        """Test the connection context manager."""
        # Mock environment variables
        mock_getenv.side_effect = lambda key, default=None: {
            'SUPABASE_DB_HOST': 'test-host',
            'SUPABASE_DB_PORT': '5432',
            'SUPABASE_DB_NAME': 'test-db',
            'SUPABASE_DB_USER': 'test-user',
            'SUPABASE_DB_PASSWORD': 'test-password'
        }.get(key, default)
        
        # Mock the connection pool
        mock_pool_instance = MagicMock()
        mock_pool.return_value = mock_pool_instance
        
        # Mock a connection
        mock_conn = MagicMock()
        mock_pool_instance.getconn.return_value = mock_conn
        
        # Create an instance
        instance = PostgresDirectConnection()
        
        # Use the connection context manager
        with instance.get_connection() as conn:
            self.assertEqual(conn, mock_conn)
        
        # Verify the connection was returned to the pool
        mock_pool_instance.getconn.assert_called_once()
        mock_pool_instance.putconn.assert_called_once_with(mock_conn)

    @patch('postgres_direct.psycopg2.pool.ThreadedConnectionPool')
    @patch('postgres_direct.load_dotenv')
    @patch('postgres_direct.os.getenv')
    def test_transaction_context_manager(self, mock_getenv, mock_load_dotenv, mock_pool):
        """Test the transaction context manager."""
        # Mock environment variables
        mock_getenv.side_effect = lambda key, default=None: {
            'SUPABASE_DB_HOST': 'test-host',
            'SUPABASE_DB_PORT': '5432',
            'SUPABASE_DB_NAME': 'test-db',
            'SUPABASE_DB_USER': 'test-user',
            'SUPABASE_DB_PASSWORD': 'test-password'
        }.get(key, default)
        
        # Mock the connection pool
        mock_pool_instance = MagicMock()
        mock_pool.return_value = mock_pool_instance
        
        # Mock a connection
        mock_conn = MagicMock()
        mock_pool_instance.getconn.return_value = mock_conn
        
        # Create an instance
        instance = PostgresDirectConnection()
        
        # Use the transaction context manager
        with instance.transaction() as conn:
            self.assertEqual(conn, mock_conn)
        
        # Verify the transaction was committed
        mock_conn.commit.assert_called_once()
        mock_pool_instance.putconn.assert_called_once_with(mock_conn)
        
        # Test transaction rollback on exception
        mock_pool_instance.reset_mock()
        mock_conn.reset_mock()
        
        try:
            with instance.transaction() as conn:
                raise ValueError("Test exception")
        except ValueError:
            pass
        
        # Verify the transaction was rolled back
        mock_conn.rollback.assert_called_once()
        mock_pool_instance.putconn.assert_called_once_with(mock_conn)

    @patch('postgres_direct.psycopg2.pool.ThreadedConnectionPool')
    @patch('postgres_direct.load_dotenv')
    @patch('postgres_direct.os.getenv')
    def test_execute_query(self, mock_getenv, mock_load_dotenv, mock_pool):
        """Test the execute_query method."""
        # Mock environment variables
        mock_getenv.side_effect = lambda key, default=None: {
            'SUPABASE_DB_HOST': 'test-host',
            'SUPABASE_DB_PORT': '5432',
            'SUPABASE_DB_NAME': 'test-db',
            'SUPABASE_DB_USER': 'test-user',
            'SUPABASE_DB_PASSWORD': 'test-password'
        }.get(key, default)
        
        # Mock the connection pool
        mock_pool_instance = MagicMock()
        mock_pool.return_value = mock_pool_instance
        
        # Mock a connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_pool_instance.getconn.return_value = mock_conn
        mock_conn.cursor.return_value.__enter__.return_value = mock_cursor
        
        # Mock cursor results
        mock_cursor.description = ['id', 'name']
        mock_cursor.fetchall.return_value = [{'id': 1, 'name': 'Test'}]
        
        # Create an instance
        instance = PostgresDirectConnection()
        
        # Execute a query
        result = instance.execute_query("SELECT * FROM test", {'param': 'value'})
        
        # Verify the query was executed
        mock_cursor.execute.assert_called_once_with("SELECT * FROM test", {'param': 'value'})
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])

    @patch('postgres_direct.psycopg2.pool.ThreadedConnectionPool')
    @patch('postgres_direct.load_dotenv')
    @patch('postgres_direct.os.getenv')
    def test_bulk_insert(self, mock_getenv, mock_load_dotenv, mock_pool):
        """Test the bulk_insert method."""
        # Mock environment variables
        mock_getenv.side_effect = lambda key, default=None: {
            'SUPABASE_DB_HOST': 'test-host',
            'SUPABASE_DB_PORT': '5432',
            'SUPABASE_DB_NAME': 'test-db',
            'SUPABASE_DB_USER': 'test-user',
            'SUPABASE_DB_PASSWORD': 'test-password'
        }.get(key, default)
        
        # Create a mock instance with a mocked execute_query method
        with patch('postgres_direct.PostgresDirectConnection.execute_query') as mock_execute:
            # Mock the execute_query result for returning IDs
            mock_execute.return_value = [{'id': 1}, {'id': 2}]
            
            # Create an instance
            instance = PostgresDirectConnection()
            
            # Test bulk insert
            data = [
                {'name': 'Test 1', 'value': 100},
                {'name': 'Test 2', 'value': 200}
            ]
            
            result = instance.bulk_insert('test_table', data, return_ids=True)
            
            # Verify the query was executed
            mock_execute.assert_called_once()
            args, kwargs = mock_execute.call_args
            
            # Check that the query is an INSERT statement
            self.assertTrue(args[0].startswith('INSERT INTO test_table'))
            
            # Check that the parameters include all values
            self.assertEqual(len(args[1]), 4)  # 2 rows * 2 columns
            
            # Check the returned IDs
            self.assertEqual(result, [1, 2])


class TestSQLExecutor(unittest.TestCase):
    """Test cases for the SQL executor utilities."""
    
    @patch('sql_executor.PostgresDirectConnection')
    def test_execute_query(self, mock_db_class):
        """Test the execute_query function."""
        # Mock the database instance
        mock_db = MagicMock()
        mock_db_class.return_value = mock_db
        
        # Mock the execute_query method
        mock_db.execute_query.return_value = [{'id': 1, 'name': 'Test'}]
        
        # Execute a query
        result = sql_executor.execute_query("SELECT * FROM test", {'param': 'value'})
        
        # Verify the query was executed
        mock_db.execute_query.assert_called_once_with(
            "SELECT * FROM test", {'param': 'value'}, False, True
        )
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])

    @patch('sql_executor.PostgresDirectConnection')
    def test_bulk_insert(self, mock_db_class):
        """Test the bulk_insert function."""
        # Mock the database instance
        mock_db = MagicMock()
        mock_db_class.return_value = mock_db
        
        # Mock the bulk_insert method
        mock_db.bulk_insert.return_value = [1, 2]
        
        # Test bulk insert
        data = [
            {'name': 'Test 1', 'value': 100},
            {'name': 'Test 2', 'value': 200}
        ]
        
        result = sql_executor.bulk_insert('test_table', data, return_ids=True)
        
        # Verify the method was called
        mock_db.bulk_insert.assert_called_once_with(
            'test_table', data, None, True, 1000
        )
        self.assertEqual(result, [1, 2])

    @patch('sql_executor.PostgresDirectConnection')
    def test_with_transaction(self, mock_db_class):
        """Test the with_transaction decorator."""
        # Mock the database instance
        mock_db = MagicMock()
        mock_db_class.return_value = mock_db
        
        # Mock the transaction methods
        mock_transaction = MagicMock()
        mock_db.begin_transaction.return_value = mock_transaction
        
        # Define a test function
        @sql_executor.with_transaction
        def test_func(transaction, arg1, arg2=None):
            return transaction, arg1, arg2
        
        # Call the decorated function
        result = test_func('value1', arg2='value2')
        
        # Verify the transaction methods were used
        mock_db.begin_transaction.assert_called_once()
        mock_db.commit_transaction.assert_called_once_with(mock_transaction)
        self.assertEqual(result, (mock_transaction, 'value1', 'value2'))

    @patch('sql_executor.time.sleep')
    def test_with_retry(self, mock_sleep):
        """Test the with_retry decorator."""
        # Define a test function that fails twice then succeeds
        attempt_count = 0
        
        def test_func():
            nonlocal attempt_count
            attempt_count += 1
            if attempt_count < 3:
                raise ValueError("Test error")
            return "Success"
        
        # Decorate the function
        decorated_func = sql_executor.with_retry(test_func, max_retries=3, retry_delay=0.1)
        
        # Call the decorated function
        result = decorated_func()
        
        # Verify the function was retried
        self.assertEqual(attempt_count, 3)
        self.assertEqual(result, "Success")
        self.assertEqual(mock_sleep.call_count, 2)

    @patch('sql_executor.logger')
    def test_process_in_batches(self, mock_logger):
        """Test the process_in_batches function."""
        # Define a test function
        processed_items = []
        
        def process_func(batch):
            processed_items.extend(batch)
            return len(batch)
        
        # Process items in batches
        items = list(range(25))
        results = sql_executor.process_in_batches(items, batch_size=10, process_func=process_func)
        
        # Verify the items were processed
        self.assertEqual(processed_items, items)
        self.assertEqual(results, [10, 10, 5])
        self.assertEqual(mock_logger.info.call_count, 9)  # 3 batches * 3 log messages per batch

    @patch('sql_executor.execute_batch')
    def test_upsert(self, mock_execute_batch):
        """Test the upsert function."""
        # Test data
        data = [
            {'id': 1, 'name': 'Test 1', 'value': 100},
            {'id': 2, 'name': 'Test 2', 'value': 200}
        ]
        
        # Call upsert
        sql_executor.upsert('test_table', data, key_columns=['id'], update_columns=['name', 'value'])
        
        # Verify execute_batch was called
        mock_execute_batch.assert_called_once()
        args, kwargs = mock_execute_batch.call_args
        
        # Check that the query is an INSERT ... ON CONFLICT statement
        self.assertTrue('INSERT INTO test_table' in args[0])
        self.assertTrue('ON CONFLICT ("id")' in args[0])
        self.assertTrue('DO UPDATE SET' in args[0])


if __name__ == '__main__':
    unittest.main()