from flask import Response, current_app
from functools import wraps
import os

def security_headers(response):
    """Middleware function to add security headers to HTTP responses."""
    # If the response is a Response object
    if isinstance(response, Response):
        # Content Security Policy
        response.headers['Content-Security-Policy'] = "default-src 'self'; script-src 'self' 'unsafe-inline' cdn.jsdelivr.net; style-src 'self' 'unsafe-inline' cdn.jsdelivr.net fonts.googleapis.com; img-src 'self' data:; font-src 'self' fonts.gstatic.com cdn.jsdelivr.net; connect-src 'self' *.supabase.co"
        # HTTP Strict Transport Security
        response.headers['Strict-Transport-Security'] = 'max-age=31536000; includeSubDomains; preload'
        # X-Content-Type-Options
        response.headers['X-Content-Type-Options'] = 'nosniff'
        # X-Frame-Options
        response.headers['X-Frame-Options'] = 'SAMEORIGIN'
        # X-XSS-Protection (deprecated but included for older browsers)
        response.headers['X-XSS-Protection'] = '1; mode=block'
        # Referrer-Policy
        response.headers['Referrer-Policy'] = 'strict-origin-when-cross-origin'
        # Permissions-Policy (formerly Feature-Policy)
        response.headers['Permissions-Policy'] = "geolocation=(), microphone=(), camera=(), payment=(), usb=(), screen-wake-lock=(), interest-cohort=()"
        # Also include the older Feature-Policy header for backward compatibility
        response.headers['Feature-Policy'] = "geolocation 'none'; microphone 'none'; camera 'none'"
    return response

def security_headers_decorator():
    """Decorator version of security headers function."""
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            resp = f(*args, **kwargs)
            return security_headers(resp)
        return decorated_function
    return decorator

def apply_security_headers(app):
    """Apply security headers to all routes of a Flask app."""
    app.after_request(security_headers)