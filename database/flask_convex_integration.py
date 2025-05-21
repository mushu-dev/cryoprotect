"""
Flask integration for Convex.

This module provides utilities to integrate Convex with a Flask application,
including authentication, middleware, and database access.
"""

import os
import logging
from typing import Dict, Any, Optional
from flask import Flask, g, request, jsonify

from .enhanced_convex_adapter import ConvexAdapter, create_convex_adapter
from .auth_bridge import AuthBridge, init_auth_bridge

logger = logging.getLogger(__name__)

class FlaskConvexIntegration:
    """
    Flask integration for Convex.
    
    This class provides a convenient wrapper for integrating Convex with a
    Flask application, including authentication, middleware, and database access.
    """
    
    def __init__(self, app: Optional[Flask] = None, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the Flask-Convex integration.
        
        Args:
            app: Optional Flask app to initialize with
            config: Optional configuration dict
        """
        self.app = app
        self.config = config or {}
        self.convex_adapter = None
        self.auth_bridge = None
        
        if app is not None:
            self.init_app(app, config)
    
    def init_app(self, app: Flask, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the integration with a Flask app.
        
        Args:
            app: Flask app instance
            config: Optional configuration dict
        """
        if config:
            self.config.update(config)
        
        # Initialize AuthBridge
        jwt_secret = self.config.get('jwt_secret') or os.environ.get('JWT_SECRET')
        jwt_algorithm = self.config.get('jwt_algorithm') or os.environ.get('JWT_ALGORITHM', 'HS256')
        jwt_expiry = self.config.get('jwt_expiry') or int(os.environ.get('JWT_EXPIRY_SECONDS', 24 * 60 * 60))
        
        self.auth_bridge = init_auth_bridge(app, jwt_secret, jwt_algorithm, jwt_expiry)
        
        # Initialize Convex adapter
        convex_config = {
            'url': self.config.get('convex_url') or os.environ.get('CONVEX_URL', ''),
            'key': self.config.get('convex_key') or os.environ.get('CONVEX_DEPLOYMENT_KEY', ''),
            'timeout': self.config.get('convex_timeout') or int(os.environ.get('CONVEX_TIMEOUT', '30')),
            'retry_count': self.config.get('convex_retry_count') or int(os.environ.get('CONVEX_RETRY_COUNT', '3')),
            'circuit_breaker_threshold': self.config.get('convex_circuit_breaker_threshold') or 
                int(os.environ.get('CONVEX_CIRCUIT_BREAKER_THRESHOLD', '5')),
            'circuit_breaker_timeout': self.config.get('convex_circuit_breaker_timeout') or 
                int(os.environ.get('CONVEX_CIRCUIT_BREAKER_TIMEOUT', '60'))
        }
        
        self.convex_adapter = create_convex_adapter(convex_config, self.auth_bridge)
        
        # Connect to Convex
        if self.convex_adapter.connect():
            logger.info("Connected to Convex")
        else:
            logger.warning("Failed to connect to Convex")
        
        # Add Convex adapter to app context
        app.convex = self.convex_adapter
        
        # Register cleanup function
        @app.teardown_appcontext
        def close_convex_connection(error=None):
            if hasattr(g, 'convex_transaction'):
                if error:
                    # Rollback transaction on error
                    try:
                        self.convex_adapter.rollback_transaction(g.convex_transaction)
                    except Exception as e:
                        logger.error(f"Error rolling back transaction: {str(e)}")
                else:
                    # Commit transaction on success
                    try:
                        self.convex_adapter.commit_transaction(g.convex_transaction)
                    except Exception as e:
                        logger.error(f"Error committing transaction: {str(e)}")
                
                g.convex_transaction = None
        
        # Add middleware for processing JWTs and setting Convex tokens
        @app.before_request
        def set_convex_token():
            # Check if user is authenticated via AuthBridge
            if hasattr(g, 'user_id') and g.user_id and hasattr(g, 'convex_token'):
                # Update Convex adapter with user token
                self.convex_adapter.set_user_token(g.convex_token)
        
        logger.info("Initialized Flask-Convex integration")
    
    def get_convex_adapter(self) -> ConvexAdapter:
        """
        Get the Convex adapter.
        
        Returns:
            ConvexAdapter: Configured ConvexAdapter instance
        """
        return self.convex_adapter
    
    def get_auth_bridge(self) -> AuthBridge:
        """
        Get the AuthBridge.
        
        Returns:
            AuthBridge: Configured AuthBridge instance
        """
        return self.auth_bridge
    
    def transaction(self):
        """
        Context manager for Convex transactions.
        
        Usage:
            with flask_convex.transaction():
                # Do database operations
                # Transaction is automatically committed on success
                # or rolled back on error
                
        Returns:
            A context manager for transactions
        """
        class TransactionContext:
            def __init__(self, integration):
                self.integration = integration
            
            def __enter__(self):
                g.convex_transaction = self.integration.convex_adapter.begin_transaction()
                return g.convex_transaction
            
            def __exit__(self, exc_type, exc_val, exc_tb):
                if exc_type is not None:
                    # Rollback on exception
                    self.integration.convex_adapter.rollback_transaction(g.convex_transaction)
                else:
                    # Commit on success
                    self.integration.convex_adapter.commit_transaction(g.convex_transaction)
                
                g.convex_transaction = None
                return False  # Don't suppress exceptions
        
        return TransactionContext(self)

def init_flask_convex(app: Flask, config: Optional[Dict[str, Any]] = None) -> FlaskConvexIntegration:
    """
    Initialize the Flask-Convex integration.
    
    Args:
        app: Flask app instance
        config: Optional configuration dict
        
    Returns:
        FlaskConvexIntegration: Configured FlaskConvexIntegration instance
    """
    return FlaskConvexIntegration(app, config)