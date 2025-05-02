import unittest
from flask import Flask, Response
from security_headers import security_headers, apply_security_headers

class TestSecurityHeaders(unittest.TestCase):
    """Test case for security headers middleware."""

    def setUp(self):
        """Set up test client and app."""
        # Create a simple Flask app for testing
        self.app = Flask(__name__)
        
        # Apply security headers middleware
        apply_security_headers(self.app)
        
        # Add a simple route for testing
        @self.app.route('/test')
        def test_route():
            return 'Test'
        
        self.client = self.app.test_client()

    def test_security_headers_applied(self):
        """Test that security headers are applied to responses."""
        response = self.client.get('/test')
        
        # Verify that all required security headers are present
        self.assertEqual(response.headers.get('Content-Security-Policy'), 
                         "default-src 'self'; script-src 'self' 'unsafe-inline' cdn.jsdelivr.net; style-src 'self' 'unsafe-inline' cdn.jsdelivr.net fonts.googleapis.com; img-src 'self' data:; font-src 'self' fonts.gstatic.com cdn.jsdelivr.net; connect-src 'self' *.supabase.co")
        self.assertEqual(response.headers.get('Strict-Transport-Security'),
                         'max-age=31536000; includeSubDomains; preload')
        self.assertEqual(response.headers.get('X-Content-Type-Options'), 
                         'nosniff')
        self.assertEqual(response.headers.get('X-Frame-Options'), 
                         'SAMEORIGIN')
        self.assertEqual(response.headers.get('X-XSS-Protection'), 
                         '1; mode=block')
        self.assertEqual(response.headers.get('Referrer-Policy'), 
                         'strict-origin-when-cross-origin')
        self.assertEqual(response.headers.get('Permissions-Policy'),
                         "geolocation=(), microphone=(), camera=(), payment=(), usb=(), screen-wake-lock=(), interest-cohort=()")
        self.assertEqual(response.headers.get('Permissions-Policy'),
                         "geolocation=(), microphone=(), camera=(), payment=(), usb=(), screen-wake-lock=(), interest-cohort=()")
        self.assertEqual(response.headers.get('Feature-Policy'),
                         "geolocation 'none'; microphone 'none'; camera 'none'")

    def test_security_headers_function(self):
        """Test the security_headers function directly."""
        # Create a mock response
        app = Flask(__name__)
        with app.test_request_context():
            response = app.make_response('Test')
            
            # Apply security headers
            response = security_headers(response)
            
            # Verify headers
            self.assertEqual(response.headers.get('Content-Security-Policy'), 
                             "default-src 'self'; script-src 'self' 'unsafe-inline' cdn.jsdelivr.net; style-src 'self' 'unsafe-inline' cdn.jsdelivr.net fonts.googleapis.com; img-src 'self' data:; font-src 'self' fonts.gstatic.com cdn.jsdelivr.net; connect-src 'self' *.supabase.co")
            self.assertEqual(response.headers.get('Strict-Transport-Security'),
                             'max-age=31536000; includeSubDomains; preload')
            self.assertEqual(response.headers.get('X-Content-Type-Options'), 
                             'nosniff')
            self.assertEqual(response.headers.get('X-Frame-Options'), 
                             'SAMEORIGIN')
            self.assertEqual(response.headers.get('X-XSS-Protection'), 
                             '1; mode=block')
            self.assertEqual(response.headers.get('Referrer-Policy'), 
                             'strict-origin-when-cross-origin')
            self.assertEqual(response.headers.get('Feature-Policy'), 
                             "geolocation 'none'; microphone 'none'; camera 'none'")


if __name__ == '__main__':
    unittest.main()