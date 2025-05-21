"""
CORS configuration for the RDKit service.
This module provides proper CORS settings for integrating with the 
main API, frontend on Netlify, and Convex.
"""

import os
from flask_cors import CORS

def configure_cors(app):
    """
    Configure CORS settings for the Flask-based RDKit service to work with
    all CryoProtect services.
    
    Args:
        app: The Flask application instance
    """
    # Default origins to allow in development
    default_origins = [
        "http://localhost:3000",                  # Local Next.js development
        "http://127.0.0.1:3000",                  # Alternative local development
        "https://cryoprotect.netlify.app",        # Netlify deployment
        "https://www.cryoprotect.app",            # Custom domain for frontend
        "https://cryoprotect-8030e4025428.herokuapp.com"  # Main API
    ]
    
    # Add Convex URL if available
    convex_url = os.environ.get("CONVEX_URL", "https://upbeat-parrot-866.convex.cloud")
    if convex_url:
        default_origins.append(convex_url)
        
    # Get additional allowed origins from environment
    additional_origins = os.environ.get("ALLOWED_ORIGINS", "")
    if additional_origins:
        default_origins.extend(additional_origins.split(","))
    
    # Frontend URL from environment
    frontend_url = os.environ.get("FRONTEND_URL", "*")
    if frontend_url != "*" and frontend_url not in default_origins:
        default_origins.append(frontend_url)
    
    # Main API URL from environment
    api_url = os.environ.get("API_URL", "https://api.cryoprotect.app")
    if api_url != "*" and api_url not in default_origins:
        default_origins.append(api_url)
    
    # Configure CORS with all allowed origins
    CORS(app, resources={
        r"/*": {
            "origins": default_origins,
            "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
            "allow_headers": [
                "Content-Type", 
                "Authorization", 
                "X-API-Key", 
                "X-Requested-With"
            ],
            "supports_credentials": True,
            "max_age": 600  # Cache preflight response for 10 minutes
        }
    })
    
    # Add test endpoint for CORS verification
    @app.route('/test-cors', methods=['GET', 'OPTIONS'])
    def test_cors():
        """Endpoint for testing CORS configuration"""
        return {
            "success": True,
            "message": "CORS is configured properly for RDKit service",
            "allowed_origins": default_origins
        }
        
    # Add debug endpoint for CORS configuration
    @app.route('/debug/cors', methods=['GET'])
    def debug_cors():
        """Return the current CORS configuration for debugging"""
        return {
            "configured_origins": default_origins,
            "frontend_url": frontend_url,
            "api_url": api_url,
            "convex_url": convex_url
        }
    
    return app