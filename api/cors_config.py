"""
CORS configuration for the main backend API.
This module provides proper CORS settings for integrating with all frontend
and backend services including Netlify, RDKit service, and Convex.
"""

import os
from flask_cors import CORS

def configure_cors(app):
    """
    Configure CORS settings for the Flask application to work with 
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
    ]
    
    # Add RDKit service URL if available
    rdkit_service_url = os.environ.get("RDKIT_SERVICE_URL", "https://rdkit.cryoprotect.app")
    if rdkit_service_url:
        default_origins.append(rdkit_service_url)
    
    # Add Convex URL if available
    convex_url = os.environ.get("CONVEX_URL", "https://dynamic-mink-63.convex.cloud")
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
            "message": "CORS is configured properly",
            "allowed_origins": default_origins
        }
        
    # Add debug endpoint for CORS configuration
    @app.route('/debug/cors', methods=['GET'])
    def debug_cors():
        """Return the current CORS configuration for debugging"""
        return {
            "configured_origins": default_origins,
            "frontend_url": frontend_url,
            "rdkit_service_url": rdkit_service_url,
            "convex_url": convex_url
        }
    
    return app