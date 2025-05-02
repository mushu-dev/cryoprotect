#!/usr/bin/env python3
"""
Unit tests for property_utils.py locking mechanism

Tests the lock-based property type creation to prevent race conditions,
focusing on the thread-safety of the get_or_create_property_type method.
"""

import unittest
from unittest.mock import MagicMock, patch, call
import threading
import time
import uuid
import queue
from concurrent.futures import ThreadPoolExecutor

from property_utils import PropertyManager
from transaction_utils import safe_transaction

class TestPropertyUtilsLocking(unittest.TestCase):
    """Test cases for property_utils.py locking mechanism."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a mock connection
        self.mock_connection = MagicMock()
        
        # Create a mock execute_query function that simulates database operations
        self.property_types = {}
        self.property_type_id_counter = 0
        self.query_delay = 0.05  # Simulate database query delay
        
        # Clear the class-level locks to ensure tests are isolated
        PropertyManager._property_type_locks = {}
        PropertyManager._property_type_locks_mutex = threading.Lock()
        
        # Create a patch for the execute_query function
        self.execute_query_patcher = patch('property_utils.execute_query')
        self.mock_execute_query = self.execute_query_patcher.start()
        self.mock_execute_query.side_effect = self.mock_execute_query_func
        
        # Create a patch for the safe_transaction context manager
        self.safe_transaction_patcher = patch('property_utils.safe_transaction')
        self.mock_safe_transaction = self.safe_transaction_patcher.start()
        self.mock_safe_transaction.return_value.__enter__.return_value = self.mock_connection
        self.mock_safe_transaction.return_value.__exit__.return_value = None
        
        # Create a patch for the get_db_connection function
        self.get_db_connection_patcher = patch('property_utils.get_db_connection')
        self.mock_get_db_connection = self.get_db_connection_patcher.start()
        self.mock_get_db_connection.return_value = self.mock_connection
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.execute_query_patcher.stop()
        self.safe_transaction_patcher.stop()
        self.get_db_connection_patcher.stop()
    
    def mock_execute_query_func(self, query, params=None, fetch_one=False, dict_cursor=True):
        """Mock implementation of execute_query that simulates database operations."""
        # Simulate query delay
        time.sleep(self.query_delay)
        
        # Handle SELECT query for property types
        if query.strip().upper().startswith('SELECT') and 'property_types' in query:
            if 'WHERE name =' in query and params and params[0]:
                property_name = params[0]
                if property_name in self.property_types:
                    result = {
                        'id': self.property_types[property_name]['id'],
                        'data_type': self.property_types[property_name]['data_type']
                    }
                    return result if fetch_one else [result]
                return None if fetch_one else []
            
            elif 'WHERE id =' in query and params and params[0]:
                property_id = params[0]
                for prop_name, prop_data in self.property_types.items():
                    if prop_data['id'] == property_id:
                        result = {
                            'id': prop_data['id'],
                            'data_type': prop_data['data_type']
                        }
                        return result if fetch_one else [result]
                return None if fetch_one else []
            
            else:
                # Return all property types
                results = []
                for prop_name, prop_data in self.property_types.items():
                    results.append({
                        'id': prop_data['id'],
                        'name': prop_name,
                        'data_type': prop_data['data_type']
                    })
                return results
        
        # Handle INSERT query for property types
        elif query.strip().upper().startswith('INSERT INTO property_types'):
            property_name = params[0]
            data_type = params[1]
            
            # Generate a new UUID for the property type
            self.property_type_id_counter += 1
            new_id = uuid.UUID(f'00000000-0000-0000-0000-{self.property_type_id_counter:012d}')
            
            # Store the property type
            self.property_types[property_name] = {
                'id': new_id,
                'data_type': data_type
            }
            
            return [{'id': new_id}]
        
        # Handle other queries
        return []
    
    def test_get_property_type_id_single_thread(self):
        """Test get_property_type_id method in a single thread."""
        # Create a PropertyManager instance
        property_manager = PropertyManager(connection=self.mock_connection)
        
        # Get a property type ID for a new property
        property_name = 'test_property'
        data_type = 'numeric'
        
        property_id = property_manager.get_property_type_id(property_name, data_type)
        
        # Verify the property type was created
        self.assertIn(property_name, self.property_types)
        self.assertEqual(self.property_types[property_name]['id'], property_id)
        self.assertEqual(self.property_types[property_name]['data_type'], data_type)
        
        # Get the same property type ID again
        property_id2 = property_manager.get_property_type_id(property_name)
        
        # Verify the same ID is returned
        self.assertEqual(property_id, property_id2)
    
    def test_get_property_type_id_concurrent(self):
        """Test get_property_type_id method with concurrent threads."""
        # Create a PropertyManager instance
        property_manager = PropertyManager(connection=self.mock_connection)
        
        # Number of concurrent threads
        num_threads = 10
        
        # Property name to use
        property_name = 'concurrent_property'
        data_type = 'numeric'
        
        # Queue to collect results from threads
        results_queue = queue.Queue()
        
        def worker():
            """Worker function to get property type ID."""
            try:
                property_id = property_manager.get_property_type_id(property_name, data_type)
                results_queue.put(property_id)
            except Exception as e:
                results_queue.put(e)
        
        # Create and start threads
        threads = []
        for _ in range(num_threads):
            thread = threading.Thread(target=worker)
            threads.append(thread)
            thread.start()
        
        # Wait for all threads to complete
        for thread in threads:
            thread.join()
        
        # Collect results
        results = []
        while not results_queue.empty():
            results.append(results_queue.get())
        
        # Verify all threads got the same property type ID
        self.assertEqual(len(results), num_threads)
        first_id = results[0]
        for result in results:
            self.assertEqual(result, first_id)
        
        # Verify only one property type was created
        self.assertEqual(len(self.property_types), 1)
        self.assertIn(property_name, self.property_types)
        self.assertEqual(self.property_types[property_name]['id'], first_id)
    
    def test_get_property_type_id_different_properties_concurrent(self):
        """Test get_property_type_id method with different properties concurrently."""
        # Create a PropertyManager instance
        property_manager = PropertyManager(connection=self.mock_connection)
        
        # Number of different properties
        num_properties = 10
        
        # Queue to collect results from threads
        results_queue = queue.Queue()
        
        def worker(prop_name):
            """Worker function to get property type ID."""
            try:
                property_id = property_manager.get_property_type_id(prop_name, 'numeric')
                results_queue.put((prop_name, property_id))
            except Exception as e:
                results_queue.put((prop_name, e))
        
        # Create and start threads for different properties
        with ThreadPoolExecutor(max_workers=num_properties) as executor:
            for i in range(num_properties):
                prop_name = f'property_{i}'
                executor.submit(worker, prop_name)
        
        # Collect results
        results = {}
        while not results_queue.empty():
            prop_name, prop_id = results_queue.get()
            results[prop_name] = prop_id
        
        # Verify all properties were created
        self.assertEqual(len(results), num_properties)
        self.assertEqual(len(self.property_types), num_properties)
        
        # Verify each property has a unique ID
        for prop_name, prop_id in results.items():
            self.assertIn(prop_name, self.property_types)
            self.assertEqual(self.property_types[prop_name]['id'], prop_id)
    
    def test_get_property_type_id_with_cache(self):
        """Test get_property_type_id method with caching."""
        # Create a PropertyManager instance
        property_manager = PropertyManager(connection=self.mock_connection)
        
        # Get a property type ID for a new property
        property_name = 'cached_property'
        data_type = 'numeric'
        
        # First call should create the property
        property_id = property_manager.get_property_type_id(property_name, data_type)
        
        # Reset the mock to track new calls
        self.mock_execute_query.reset_mock()
        
        # Second call should use the cache
        property_id2 = property_manager.get_property_type_id(property_name)
        
        # Verify the same ID is returned
        self.assertEqual(property_id, property_id2)
        
        # Verify no database queries were made for the second call
        self.mock_execute_query.assert_not_called()
    
    def test_get_property_type_id_with_lock_contention(self):
        """Test get_property_type_id method with lock contention."""
        # Create a PropertyManager instance
        property_manager = PropertyManager(connection=self.mock_connection)
        
        # Property name to use
        property_name = 'contention_property'
        data_type = 'numeric'
        
        # Create a lock to control thread execution
        control_lock = threading.Lock()
        control_lock.acquire()  # Lock initially
        
        # Event to signal when first thread has acquired the property lock
        first_thread_has_lock = threading.Event()
        
        # Event to signal when second thread is waiting for the lock
        second_thread_waiting = threading.Event()
        
        # Results from threads
        results = {'thread1': None, 'thread2': None}
        
        def thread1_func():
            """First thread function."""
            try:
                # Acquire the property lock
                with PropertyManager._property_type_locks_mutex:
                    if property_name not in PropertyManager._property_type_locks:
                        PropertyManager._property_type_locks[property_name] = threading.Lock()
                    property_lock = PropertyManager._property_type_locks[property_name]
                
                # Acquire the property lock
                property_lock.acquire()
                
                # Signal that we have the lock
                first_thread_has_lock.set()
                
                # Wait for the second thread to be waiting
                second_thread_waiting.wait()
                
                # Wait for the control lock to be released
                control_lock.acquire()
                control_lock.release()
                
                # Now proceed with getting the property type ID
                results['thread1'] = property_manager.get_property_type_id(property_name, data_type)
                
                # Release the property lock
                property_lock.release()
                
            except Exception as e:
                results['thread1'] = e
        
        def thread2_func():
            """Second thread function."""
            try:
                # Wait for the first thread to acquire the lock
                first_thread_has_lock.wait()
                
                # Signal that we're about to wait for the lock
                second_thread_waiting.set()
                
                # Try to get the property type ID (will wait for the lock)
                results['thread2'] = property_manager.get_property_type_id(property_name, data_type)
                
            except Exception as e:
                results['thread2'] = e
        
        # Create and start threads
        thread1 = threading.Thread(target=thread1_func)
        thread2 = threading.Thread(target=thread2_func)
        
        thread1.start()
        thread2.start()
        
        # Wait for both threads to be in position
        first_thread_has_lock.wait()
        second_thread_waiting.wait()
        
        # Release the control lock to allow thread1 to proceed
        control_lock.release()
        
        # Wait for both threads to complete
        thread1.join()
        thread2.join()
        
        # Verify both threads got the same property type ID
        self.assertIsNotNone(results['thread1'])
        self.assertIsNotNone(results['thread2'])
        self.assertEqual(results['thread1'], results['thread2'])
        
        # Verify only one property type was created
        self.assertEqual(len(self.property_types), 1)
        self.assertIn(property_name, self.property_types)
        self.assertEqual(self.property_types[property_name]['id'], results['thread1'])


if __name__ == '__main__':
    unittest.main()