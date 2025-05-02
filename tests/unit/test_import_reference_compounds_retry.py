#!/usr/bin/env python3
"""
Test for the with_retry decorator in import_reference_compounds.py
"""

import unittest
from unittest.mock import patch, MagicMock
import sys
import os

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Import the module to test
from sql_executor import with_retry
import import_reference_compounds

class TestImportReferenceCompoundsRetry(unittest.TestCase):
    """Test the with_retry decorator in import_reference_compounds.py"""

    @patch('sql_executor.time.sleep')
    def test_with_retry_parameter(self, mock_sleep):
        """Test that the with_retry decorator accepts max_retries parameter"""
        
        # Create a mock function that fails twice then succeeds
        mock_func = MagicMock(side_effect=[Exception("First failure"),
                                          Exception("Second failure"),
                                          "Success"])
        
        # Apply the decorator with max_retries parameter
        decorated_func = with_retry(max_retries=3, retry_delay=0.1)(mock_func)
        
        # Call the decorated function
        result = decorated_func()
        
        # Verify the function was called the expected number of times
        self.assertEqual(mock_func.call_count, 3)
        self.assertEqual(result, "Success")
        
        # Verify sleep was called with the expected delays
        mock_sleep.assert_any_call(0.1)
        
        # Check that the second sleep call used a value close to 0.15 (with backoff of 1.5)
        second_call_value = mock_sleep.call_args_list[1][0][0]
        self.assertAlmostEqual(second_call_value, 0.15, places=1)
    
    @patch('sql_executor.time.sleep')
    def test_with_retry_backward_compatibility(self, mock_sleep):
        """Test that the with_retry decorator accepts max_attempts parameter for backward compatibility"""
        
        # Create a mock function that fails twice then succeeds
        mock_func = MagicMock(side_effect=[Exception("First failure"),
                                          Exception("Second failure"),
                                          "Success"])
        
        # Apply the decorator with max_attempts parameter (backward compatibility)
        decorated_func = with_retry(max_attempts=4, retry_delay=0.1)(mock_func)
        
        # Call the decorated function
        result = decorated_func()
        
        # Verify the function was called the expected number of times
        self.assertEqual(mock_func.call_count, 3)
        self.assertEqual(result, "Success")
        
        # Verify sleep was called with the expected delays
        mock_sleep.assert_any_call(0.1)
        
        # Check that the second sleep call used a value close to 0.15 (with backoff of 1.5)
        second_call_value = mock_sleep.call_args_list[1][0][0]
        self.assertAlmostEqual(second_call_value, 0.15, places=1)

if __name__ == '__main__':
    unittest.main()