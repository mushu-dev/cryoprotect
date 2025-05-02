import unittest
import requests
from unittest.mock import MagicMock, patch
from chembl.error_handler import ErrorCategory, classify_error, get_recovery_strategy


class TestErrorHandler(unittest.TestCase):
    """Test cases for ChEMBL error handling components."""

    def test_error_categories(self):
        """Test that all expected error categories exist."""
        self.assertTrue(hasattr(ErrorCategory, 'API_RATE_LIMIT'))
        self.assertTrue(hasattr(ErrorCategory, 'CONNECTION_ERROR'))
        self.assertTrue(hasattr(ErrorCategory, 'API_SERVER_ERROR'))
        self.assertTrue(hasattr(ErrorCategory, 'API_CLIENT_ERROR'))
        self.assertTrue(hasattr(ErrorCategory, 'DATA_VALIDATION'))
        self.assertTrue(hasattr(ErrorCategory, 'TRANSFORMATION'))
        self.assertTrue(hasattr(ErrorCategory, 'PARSING_ERROR'))
        self.assertTrue(hasattr(ErrorCategory, 'DATABASE_ERROR'))
        self.assertTrue(hasattr(ErrorCategory, 'UNKNOWN'))
        
        # Verify enum values
        self.assertEqual(ErrorCategory.API_RATE_LIMIT.value, "rate_limit")
        self.assertEqual(ErrorCategory.CONNECTION_ERROR.value, "connection")
        self.assertEqual(ErrorCategory.API_SERVER_ERROR.value, "server_error")
        self.assertEqual(ErrorCategory.API_CLIENT_ERROR.value, "client_error")
        self.assertEqual(ErrorCategory.DATA_VALIDATION.value, "validation")
        self.assertEqual(ErrorCategory.TRANSFORMATION.value, "transformation")
        self.assertEqual(ErrorCategory.PARSING_ERROR.value, "parsing")
        self.assertEqual(ErrorCategory.DATABASE_ERROR.value, "database")
        self.assertEqual(ErrorCategory.UNKNOWN.value, "unknown")

    def test_classify_error_timeout(self):
        """Test classification of timeout errors."""
        timeout_error = requests.exceptions.Timeout()
        category, description = classify_error(timeout_error)
        self.assertEqual(category, ErrorCategory.CONNECTION_ERROR)
        self.assertEqual(description, "Request timed out")

    def test_classify_error_connection(self):
        """Test classification of connection errors."""
        connection_error = requests.exceptions.ConnectionError()
        category, description = classify_error(connection_error)
        self.assertEqual(category, ErrorCategory.CONNECTION_ERROR)
        self.assertEqual(description, "Connection error")

    def test_classify_error_rate_limit(self):
        """Test classification of rate limit errors."""
        http_error = requests.exceptions.HTTPError()
        mock_response = MagicMock()
        mock_response.status_code = 429
        http_error.response = mock_response
        
        category, description = classify_error(http_error)
        self.assertEqual(category, ErrorCategory.API_RATE_LIMIT)
        self.assertEqual(description, "Rate limit exceeded: 429")

    def test_classify_error_client_error(self):
        """Test classification of client errors."""
        for status_code in [400, 401, 403, 404]:
            http_error = requests.exceptions.HTTPError()
            mock_response = MagicMock()
            mock_response.status_code = status_code
            http_error.response = mock_response
            
            category, description = classify_error(http_error)
            self.assertEqual(category, ErrorCategory.API_CLIENT_ERROR)
            self.assertEqual(description, f"Client error: {status_code}")

    def test_classify_error_server_error(self):
        """Test classification of server errors."""
        for status_code in [500, 502, 503, 504]:
            http_error = requests.exceptions.HTTPError()
            mock_response = MagicMock()
            mock_response.status_code = status_code
            http_error.response = mock_response
            
            category, description = classify_error(http_error)
            self.assertEqual(category, ErrorCategory.API_SERVER_ERROR)
            self.assertEqual(description, f"Server error: {status_code}")

    def test_classify_error_json_parsing(self):
        """Test classification of JSON parsing errors."""
        json_error = ValueError("JSON decoding error")
        category, description = classify_error(json_error)
        self.assertEqual(category, ErrorCategory.PARSING_ERROR)
        self.assertTrue("JSON parsing error" in description)

    def test_classify_error_database(self):
        """Test classification of database errors."""
        # Create a mock psycopg2 module
        mock_psycopg2 = MagicMock()
        mock_psycopg2.Error = type('Error', (Exception,), {})
        
        # Create a database error that is an instance of psycopg2.Error
        db_error = MagicMock(spec=Exception)
        db_error.__class__ = mock_psycopg2.Error
        db_error.__str__ = lambda x: "Database connection failed"
        
        # Patch the import to return our mock
        with patch.dict('sys.modules', {'psycopg2': mock_psycopg2}):
            # Force the import to succeed
            with patch('builtins.__import__', return_value=mock_psycopg2):
                category, description = classify_error(db_error)
                self.assertEqual(category, ErrorCategory.DATABASE_ERROR)
                self.assertTrue("Database error" in description)

    def test_classify_error_validation(self):
        """Test classification of validation errors."""
        validation_error = ValueError("Data validation failed")
        category, description = classify_error(validation_error)
        self.assertEqual(category, ErrorCategory.DATA_VALIDATION)
        self.assertTrue("Validation error" in description)

    def test_classify_error_transformation(self):
        """Test classification of transformation errors."""
        transform_error = ValueError("Data transformation failed")
        category, description = classify_error(transform_error)
        self.assertEqual(category, ErrorCategory.TRANSFORMATION)
        self.assertTrue("Transformation error" in description)

    def test_classify_error_unknown(self):
        """Test classification of unknown errors."""
        unknown_error = Exception("Some unexpected error")
        category, description = classify_error(unknown_error)
        self.assertEqual(category, ErrorCategory.UNKNOWN)
        self.assertTrue("Unclassified error" in description)

    def test_recovery_strategy_rate_limit(self):
        """Test recovery strategy for rate limit errors."""
        strategy = get_recovery_strategy(ErrorCategory.API_RATE_LIMIT)
        self.assertEqual(strategy, "RETRY")

    def test_recovery_strategy_connection(self):
        """Test recovery strategy for connection errors."""
        strategy = get_recovery_strategy(ErrorCategory.CONNECTION_ERROR)
        self.assertEqual(strategy, "RETRY")

    def test_recovery_strategy_server_error(self):
        """Test recovery strategy for server errors."""
        strategy = get_recovery_strategy(ErrorCategory.API_SERVER_ERROR)
        self.assertEqual(strategy, "RETRY")

    def test_recovery_strategy_client_error(self):
        """Test recovery strategy for client errors."""
        strategy = get_recovery_strategy(ErrorCategory.API_CLIENT_ERROR)
        self.assertEqual(strategy, "SKIP")

    def test_recovery_strategy_parsing_error(self):
        """Test recovery strategy for parsing errors."""
        strategy = get_recovery_strategy(ErrorCategory.PARSING_ERROR)
        self.assertEqual(strategy, "RETRY")

    def test_recovery_strategy_validation(self):
        """Test recovery strategy for validation errors."""
        strategy = get_recovery_strategy(ErrorCategory.DATA_VALIDATION)
        self.assertEqual(strategy, "SKIP")

    def test_recovery_strategy_transformation(self):
        """Test recovery strategy for transformation errors."""
        strategy = get_recovery_strategy(ErrorCategory.TRANSFORMATION)
        self.assertEqual(strategy, "SKIP")

    def test_recovery_strategy_database(self):
        """Test recovery strategy for database errors."""
        strategy = get_recovery_strategy(ErrorCategory.DATABASE_ERROR)
        self.assertEqual(strategy, "RETRY")

    def test_recovery_strategy_unknown(self):
        """Test recovery strategy for unknown errors."""
        strategy = get_recovery_strategy(ErrorCategory.UNKNOWN)
        self.assertEqual(strategy, "LOG_ONLY")

    def test_recovery_strategy_max_retries(self):
        """Test recovery strategy when max retries is exceeded."""
        # Test with a retryable error but max attempts reached
        strategy = get_recovery_strategy(ErrorCategory.API_RATE_LIMIT, attempt=5)
        self.assertEqual(strategy, "ABORT")
        
        # Test with a retryable error and attempts below max
        strategy = get_recovery_strategy(ErrorCategory.API_RATE_LIMIT, attempt=4)
        self.assertEqual(strategy, "RETRY")


if __name__ == '__main__':
    unittest.main()