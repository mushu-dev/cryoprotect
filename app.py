#!/usr/bin/env python3
"""
CryoProtect Analyzer API

This is the main Flask application for the CryoProtect Analyzer API.
It provides endpoints for accessing and manipulating data in the Supabase database.
"""

import os
import hashlib
from datetime import datetime, timedelta
from flask import Flask, jsonify, g, render_template, request, redirect, url_for, session
from flask_cors import CORS
from apispec import APISpec
from apispec.ext.marshmallow import MarshmallowPlugin
from flask_apispec.extension import FlaskApiSpec
from flask_mail import Mail

from config import active_config
from api import init_app
from api.docs import init_docs
from api.utils import get_supabase_client, authenticate_user, token_required, release_supabase_connection
from logging_enhanced import setup_enhanced_logging, get_logger, log_with_context
import logging
from api.rate_limiter import configure_rate_limiter, add_rate_limit_headers
import rate_limit_config
from monitoring import init_metrics, PrometheusMiddleware
from monitoring.middleware import update_entity_counts, update_connection_pool_status
from api.observability import init_observability
from backup.backup_manager import BackupManager
from api.lab_verification_resources import LabVerificationResource, VerificationStatsResource
from api import api


# Import authentication configuration
from auth_config import (
    JWT_EXPIRY, JWT_REFRESH_EXPIRY, DEFAULT_ROLE,
    AVAILABLE_ROLES, SESSION_TIMEOUT, REFRESH_TOKEN_ROTATION,
    SECURE_COOKIES, HTTP_ONLY_COOKIES, SAME_SITE_COOKIES,
    USE_JWT_AUTH, USE_SERVICE_ROLE, USER_ID
)

# Import JWT authentication utilities
from api.jwt_auth import (
    jwt_required, role_required, get_current_user,
    extract_user_from_token, get_token_from_request
)


logger = get_logger(__name__)

from security_headers import apply_security_headers
from api.csrf import init_csrf
from api.session_security import init_session_security

def create_app(config_object=None, testing=False):
    """
    Create and configure the Flask application.

    Args:
        config_object: Configuration object to use.
        testing (bool): If True, enables testing mode and overrides config for testing. Default is False.

    Returns:
        Flask: Configured Flask application
    """
    app = Flask(__name__)

    # Setup enhanced logging with ELK integration
    setup_enhanced_logging(app)
    logger.info("Starting CryoProtect Analyzer API application.",
                extra={"event_type": "application_start", "app_version": app.config.get('API_VERSION', '1.0.0')})
    
    # Initialize Prometheus metrics
    init_metrics(app)
    
    # Add Prometheus middleware
    prometheus_middleware = PrometheusMiddleware(app)
    
    # Initialize API observability
    init_observability(app)
    
    # Initialize backup system
    backup_enabled = os.environ.get('BACKUP_ENABLED', '0') == '1'
    if backup_enabled:
        try:
            backup_config_path = os.environ.get('BACKUP_CONFIG_PATH', 'backup/config.json')
            app.backup_manager = BackupManager(config_path=backup_config_path)
            logger.info("Backup system initialized",
                       extra={"event_type": "backup_system_init", "config_path": backup_config_path})
        except Exception as e:
            logger.error(f"Failed to initialize backup system: {str(e)}",
                        extra={"event_type": "backup_system_init_error", "error": str(e)},
                        exc_info=True)

    # Load configuration
    if config_object:
        app.config.from_object(config_object)
    else:
        app.config.from_object(active_config)

    if testing:
        app.config['TESTING'] = True
        # Add other test-specific config overrides here if needed

    # Enable CORS with additional options for documentation access
    CORS(app, resources={
        r"/swagger-ui/*": {"origins": "*"},
        r"/api/docs/*": {"origins": "*"}
    })

    # Initialize API
    # Configure connection pooling
    app.config.update({
        'SUPABASE_MIN_CONNECTIONS': 2,
        'SUPABASE_MAX_CONNECTIONS': 10,
        'SUPABASE_CONNECTION_TIMEOUT': 30
    })
    
    # Configure rate limiting
    app.config.update({
        'RATE_LIMIT_ENABLED': rate_limit_config.RATE_LIMIT_ENABLED,
        'RATE_LIMIT_DEFAULT': rate_limit_config.RATE_LIMIT_DEFAULT,
        'RATE_LIMIT_STORAGE_URL': rate_limit_config.RATE_LIMIT_STORAGE_URL,
        'RATE_LIMIT_STRATEGY': rate_limit_config.RATE_LIMIT_STRATEGY,
        'RATE_LIMIT_HEADERS_ENABLED': rate_limit_config.RATE_LIMIT_HEADERS_ENABLED,
        'RATE_LIMIT_ENDPOINTS': rate_limit_config.RATE_LIMIT_ENDPOINTS,
        'RATE_LIMIT_EXEMPT': rate_limit_config.RATE_LIMIT_EXEMPT,
        'RATE_LIMIT_BY': rate_limit_config.RATE_LIMIT_BY,
        'RATE_LIMIT_ROLES': rate_limit_config.RATE_LIMIT_ROLES,
        'RATE_LIMIT_RETRY_AFTER': rate_limit_config.RATE_LIMIT_RETRY_AFTER
    })
    
    # Configure API documentation settings
    app.config.update({
        'API_TITLE': app.config.get('API_TITLE', 'CryoProtect API'),
        'API_VERSION': app.config.get('API_VERSION', '1.0.0'),
        'OPENAPI_VERSION': app.config.get('OPENAPI_VERSION', '3.0.2'),
        'API_CONTACT_NAME': 'CryoProtect API Support',
        'API_CONTACT_EMAIL': 'support@cryoprotect.com',
        'API_LICENSE_NAME': 'MIT',
    })
    
    # Initialize API
    init_app(app)

    # Register Lab Verification endpoints for CryoProtect v2 workflow
    # api.add_resource(LabVerificationResource, '/api/experiments/<string:experiment_id>/verification', endpoint='experiment_verification')
    # api.add_resource(LabVerificationResource, '/api/verifications/<string:verification_id>', endpoint='verification_update')
    # api.add_resource(VerificationStatsResource, '/api/verification/stats', endpoint='verification_stats')
    
    # Register OpenAPI documentation blueprint
    from api.openapi import register_openapi_blueprint
    register_openapi_blueprint(app)
    
    # Add Swagger UI route
    @app.route('/swagger-ui/')
    def swagger_ui():
        """Render the Swagger UI interface."""
        return render_template('swagger.html')
    
    # Initialize rate limiter
    configure_rate_limiter(app)
    
    # Initialize CSRF protection
    init_csrf(app)
    
    # Initialize session security
    init_session_security(app)
    
    # Apply security headers to all responses
    apply_security_headers(app)
    
    # Configure Flask-Mail
    app.config.update({
        'MAIL_SERVER': os.environ.get('MAIL_SERVER', 'smtp.gmail.com'),
        'MAIL_PORT': int(os.environ.get('MAIL_PORT', 587)),
        'MAIL_USE_TLS': os.environ.get('MAIL_USE_TLS', 'True').lower() in ('true', 'yes', '1'),
        'MAIL_USERNAME': os.environ.get('MAIL_USERNAME', ''),
        'MAIL_PASSWORD': os.environ.get('MAIL_PASSWORD', ''),
        'MAIL_DEFAULT_SENDER': os.environ.get('MAIL_DEFAULT_SENDER', 'noreply@cryoprotect.com')
    })
    mail = Mail(app)
    
    # Configure session
    app.secret_key = os.environ.get('SECRET_KEY', 'dev-secret-key')
    
    # Register error handlers with standardized response format
    from api.api_standards import create_error_response, jsonify_standard_response
    
    @app.errorhandler(400)
    def bad_request(error):
        log_with_context(
            logger, 'warning',
            f"Bad request: {str(error)}",
            context={
                "event_type": "error",
                "error_type": "bad_request",
                "error_details": str(error)
            }
        )
        return jsonify_standard_response(
            *create_error_response(
                error=str(error),
                status_code=400,
                context="Bad request"
            )
        )
    
    @app.errorhandler(401)
    def unauthorized(error):
        log_with_context(
            logger, 'warning',
            f"Unauthorized access: {str(error)}",
            context={
                "event_type": "error",
                "error_type": "unauthorized",
                "error_details": str(error)
            }
        )
        return jsonify_standard_response(
            *create_error_response(
                error=str(error),
                status_code=401,
                context="Unauthorized access"
            )
        )
    
    @app.errorhandler(403)
    def forbidden(error):
        log_with_context(
            logger, 'warning',
            f"Forbidden access: {str(error)}",
            context={
                "event_type": "error",
                "error_type": "forbidden",
                "error_details": str(error)
            }
        )
        return jsonify_standard_response(
            *create_error_response(
                error=str(error),
                status_code=403,
                context="Forbidden access"
            )
        )
    
    @app.errorhandler(404)
    def not_found(error):
        log_with_context(
            logger, 'warning',
            f"Resource not found: {str(error)}",
            context={
                "event_type": "error",
                "error_type": "not_found",
                "error_details": str(error),
                "path": request.path if request else "unknown"
            }
        )
        return jsonify_standard_response(
            *create_error_response(
                error=str(error),
                status_code=404,
                context="Resource not found"
            )
        )
    
    @app.errorhandler(409)
    def conflict(error):
        log_with_context(
            logger, 'warning',
            f"Conflict: {str(error)}",
            context={
                "event_type": "error",
                "error_type": "conflict",
                "error_details": str(error)
            }
        )
        return jsonify_standard_response(
            *create_error_response(
                error=str(error),
                status_code=409,
                context="Resource conflict"
            )
        )
    
    @app.errorhandler(429)
    def too_many_requests(error):
        log_with_context(
            logger, 'warning',
            f"Too many requests: {str(error)}",
            context={
                "event_type": "error",
                "error_type": "rate_limit",
                "error_details": str(error),
                "ip_address": request.remote_addr if request else "unknown"
            }
        )
        return jsonify_standard_response(
            *create_error_response(
                error=str(error),
                status_code=429,
                context="Too many requests"
            )
        )
    
    @app.errorhandler(500)
    def server_error(error):
        log_with_context(
            logger, 'error',
            f"Internal server error: {str(error)}",
            context={
                "event_type": "error",
                "error_type": "server_error",
                "error_details": str(error)
            },
            exc_info=True
        )
        return jsonify_standard_response(
            *create_error_response(
                error=str(error),
                status_code=500,
                context="Internal server error"
            )
        )
    
    @app.errorhandler(503)
    def service_unavailable(error):
        log_with_context(
            logger, 'error',
            f"Service unavailable: {str(error)}",
            context={
                "event_type": "error",
                "error_type": "service_unavailable",
                "error_details": str(error)
            }
        )
        return jsonify_standard_response(
            *create_error_response(
                error=str(error),
                status_code=503,
                context="Service unavailable"
            )
        )
        
    # Add a catch-all exception handler
    @app.errorhandler(Exception)
    def handle_exception(error):
        log_with_context(
            logger, 'error',
            f"Uncaught exception: {str(error)}",
            context={
                "event_type": "error",
                "error_type": "uncaught_exception",
                "error_details": str(error),
                "error_class": error.__class__.__name__
            },
            exc_info=True
        )
        return jsonify_standard_response(
            *create_error_response(
                error=error,
                context="Uncaught exception"
            )
        )
    
    # Register before request handler
    @app.before_request
    def before_request():
        # Initialize Supabase client
        try:
            get_supabase_client()
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Error initializing Supabase client: {str(e)}",
                context={
                    "event_type": "database_error",
                    "error_type": "supabase_init_error",
                    "error_details": str(e)
                },
                exc_info=True
            )
            # Don't abort here, let the request continue and fail gracefully if needed
    
    # Register teardown request handler
    @app.teardown_request
    def teardown_request(exception=None):
        # Clean up resources
        try:
            if hasattr(g, 'supabase'):
                # Release connection back to the pool if using connection pooling
                release_supabase_connection()
                
                # Update connection pool metrics
                if app.config.get('SUPABASE_CONNECTION_POOL_ENABLED', False):
                    active = app.config.get('SUPABASE_ACTIVE_CONNECTIONS', 0)
                    idle = app.config.get('SUPABASE_IDLE_CONNECTIONS', 0)
                    max_conn = app.config.get('SUPABASE_MAX_CONNECTIONS', 10)
                    update_connection_pool_status(active, idle, max_conn)
                
                # For non-pooled connections, just delete the reference
                if hasattr(g, 'supabase'):
                    del g.supabase
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Error in teardown_request: {str(e)}",
                context={
                    "event_type": "resource_cleanup_error",
                    "error_type": "teardown_error",
                    "error_details": str(e)
                },
                exc_info=True
            )
    
    # Register after request handler to add rate limit headers and log responses
    @app.after_request
    def after_request(response):
        # Add rate limit headers to response
        response = add_rate_limit_headers(response)
        
        # Log response details if not a static file or health check
        if not request.path.startswith('/static/') and request.path != '/health':
            status_code = response.status_code
            log_level = 'warning' if status_code >= 400 else 'info'
            
            # Calculate request duration if start_time was set
            duration = None
            if hasattr(g, 'start_time'):
                # Use utcnow() to match the utcnow() used in logging_enhanced.py
                duration = (datetime.utcnow() - g.start_time).total_seconds()
            
            # Get user ID if available
            user_id = getattr(g, 'user_id', 'anonymous')
            
            # Log response with context
            log_with_context(
                logger, log_level,
                f"Response: {request.method} {request.path} {status_code}",
                context={
                    "event_type": "response",
                    "request": {
                        "method": request.method,
                        "path": request.path,
                        "remote_addr": request.remote_addr
                    },
                    "response": {
                        "status_code": status_code,
                        "content_length": response.content_length,
                        "content_type": response.content_type,
                        "duration": duration
                    },
                    "user_id": user_id,
                    "correlation_id": getattr(g, 'correlation_id', 'N/A')
                }
            )
            
            # Update business metrics periodically (every 10 requests)
            if hasattr(app, '_request_count'):
                app._request_count += 1
            else:
                app._request_count = 1
                
            if app._request_count % 10 == 0 and hasattr(g, 'supabase'):
                try:
                    update_entity_counts(g.supabase)
                except Exception as e:
                    logger.error(f"Error updating entity counts: {str(e)}")
        
        return response
    
    # Add health check endpoint
    # Main health check endpoint for overall system health
    @app.route('/health')
    def health_check():
        try:
            # Get deployment color from environment
            deployment_color = os.environ.get('DEPLOYMENT_COLOR', 'production')
            
            # Check database connection
            supabase = get_supabase_client()
            response = supabase.from_("property_types").select("*").limit(1).execute()
            db_status = "connected" if hasattr(response, 'data') else "unknown"

            # Check Redis connection if configured
            redis_status = "not_configured"
            redis_url = (
                app.config.get("CACHE_REDIS_URL")
                or app.config.get("RATE_LIMIT_STORAGE_URL")
                or os.environ.get("REDIS_URL")
            )
            if redis_url and "redis" in redis_url:
                try:
                    import redis as redis_lib
                    from urllib.parse import urlparse

                    # Parse Redis URL
                    parsed = urlparse(redis_url)
                    redis_host = parsed.hostname
                    redis_port = parsed.port or 6379
                    redis_db = int(parsed.path.lstrip("/")) if parsed.path else 0
                    redis_password = parsed.password

                    r = redis_lib.StrictRedis(
                        host=redis_host,
                        port=redis_port,
                        db=redis_db,
                        password=redis_password,
                        socket_connect_timeout=2,
                        socket_timeout=2,
                    )
                    if r.ping():
                        redis_status = "connected"
                    else:
                        redis_status = "unreachable"
                except Exception as re:
                    redis_status = f"error: {str(re)}"
            
            # Check disk space
            disk_status = "unknown"
            try:
                import shutil
                total, used, free = shutil.disk_usage("/")
                disk_free_percent = (free / total) * 100
                disk_status = "ok" if disk_free_percent > 10 else "low"
            except Exception as de:
                disk_status = f"error: {str(de)}"
            
            # Check memory usage
            memory_status = "unknown"
            try:
                import psutil
                memory = psutil.virtual_memory()
                memory_status = "ok" if memory.percent < 90 else "high"
            except Exception as me:
                # psutil might not be available in all environments
                memory_status = "not_available"
            
            # Check for required environment variables
            env_vars = ["FLASK_APP", "FLASK_ENV", "SECRET_KEY"]
            missing_vars = [var for var in env_vars if not os.environ.get(var)]
            env_status = "ok" if not missing_vars else f"missing: {', '.join(missing_vars)}"
            
            # Determine overall status
            overall_status = "ok"
            if db_status != "connected" or (redis_url and redis_status != "connected") or disk_status == "low" or memory_status == "high" or env_status != "ok":
                overall_status = "degraded"

            # Log health check
            log_with_context(
                logger, 'info',
                "Health check completed",
                context={
                    "event_type": "health_check",
                    "status": overall_status,
                    "database_status": db_status,
                    "redis_status": redis_status,
                    "disk_status": disk_status,
                    "memory_status": memory_status,
                    "env_status": env_status,
                    "api_version": app.config['API_VERSION'],
                    "deployment_color": deployment_color
                }
            )

            # Return health status
            return jsonify({
                'status': overall_status,
                'version': app.config['API_VERSION'],
                'deployment': deployment_color,
                'timestamp': datetime.now().isoformat(),
                'services': {
                    'database': db_status,
                    'redis': redis_status,
                    'disk': disk_status,
                    'memory': memory_status,
                    'environment': env_status
                }
            }), 200 if overall_status == "ok" else 207  # 207 Multi-Status for degraded
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Health check failed: {str(e)}",
                context={
                    "event_type": "health_check",
                    "status": "failed",
                    "error_details": str(e),
                    "api_version": app.config['API_VERSION']
                },
                exc_info=True
            )
            return jsonify({
                'status': 'error',
                'version': app.config['API_VERSION'],
                'error': str(e),
                'timestamp': datetime.now().isoformat()
            }), 500
    
    # Liveness probe - lightweight check to verify the application is running
    @app.route('/health/liveness')
    def liveness_check():
        try:
            # Simple check that the app is running
            return jsonify({
                'status': 'alive',
                'timestamp': datetime.now().isoformat()
            }), 200
        except Exception as e:
            return jsonify({
                'status': 'error',
                'error': str(e)
            }), 500
    
    # Readiness probe - check if the application is ready to serve traffic
    @app.route('/health/readiness')
    def readiness_check():
        try:
            # Check database connection
            supabase = get_supabase_client()
            response = supabase.from_("property_types").select("*").limit(1).execute()
            db_ready = hasattr(response, 'data')
            
            # Check if required services are available
            services_ready = db_ready
            
            # Return readiness status
            if services_ready:
                return jsonify({
                    'status': 'ready',
                    'timestamp': datetime.now().isoformat()
                }), 200
            else:
                return jsonify({
                    'status': 'not_ready',
                    'reason': 'Required services not available',
                    'timestamp': datetime.now().isoformat()
                }), 503
        except Exception as e:
            return jsonify({
                'status': 'not_ready',
                'reason': str(e),
                'timestamp': datetime.now().isoformat()
            }), 503
    
    # Startup probe - check if the application has completed startup
    @app.route('/health/startup')
    def startup_check():
        try:
            # Check if all required components are initialized
            all_initialized = True
            
            # Check database connection
            try:
                supabase = get_supabase_client()
                response = supabase.from_("property_types").select("*").limit(1).execute()
                db_initialized = hasattr(response, 'data')
                all_initialized = all_initialized and db_initialized
            except:
                all_initialized = False
            
            # Return startup status
            if all_initialized:
                return jsonify({
                    'status': 'started',
                    'timestamp': datetime.now().isoformat()
                }), 200
            else:
                return jsonify({
                    'status': 'starting',
                    'timestamp': datetime.now().isoformat()
                }), 503
        except Exception as e:
            return jsonify({
                'status': 'error',
                'error': str(e),
                'timestamp': datetime.now().isoformat()
            }), 500
    
    # Add backup health check endpoint
    @app.route('/health/backup')
    def backup_health_check():
        try:
            if not hasattr(app, 'backup_manager'):
                return jsonify({
                    'status': 'disabled',
                    'message': 'Backup system is not enabled',
                    'timestamp': datetime.now().isoformat()
                }), 200
                
            # Check backup system status
            last_backup = app.backup_manager.get_last_backup_info()
            
            if last_backup:
                # Check if the last backup is too old (more than 25 hours for daily backups)
                last_backup_time = datetime.fromisoformat(last_backup.get('timestamp', ''))
                backup_age = datetime.now() - last_backup_time
                
                if backup_age > timedelta(hours=25):
                    status = 'warning'
                    message = f'Last backup is {backup_age.total_seconds() / 3600:.1f} hours old'
                else:
                    status = 'ok'
                    message = 'Backup system is healthy'
                    
                return jsonify({
                    'status': status,
                    'message': message,
                    'last_backup': {
                        'timestamp': last_backup.get('timestamp'),
                        'type': last_backup.get('backup_type'),
                        'path': last_backup.get('backup_path'),
                        'age_hours': backup_age.total_seconds() / 3600
                    },
                    'timestamp': datetime.now().isoformat()
                }), 200
            else:
                return jsonify({
                    'status': 'warning',
                    'message': 'No backups found',
                    'timestamp': datetime.now().isoformat()
                }), 200
                
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Backup health check failed: {str(e)}",
                context={
                    "event_type": "health_check",
                    "component": "backup",
                    "status": "failed",
                    "error_details": str(e)
                },
                exc_info=True
            )
            return jsonify({
                'status': 'error',
                'message': f'Backup health check failed: {str(e)}',
                'timestamp': datetime.now().isoformat()
            }), 500
    
    # Add authentication endpoints
    @app.route('/auth/login', methods=['POST'])
    def login():
        try:
            data = request.get_json()
            email = data.get('email')
            password = data.get('password')
            
            if not email or not password:
                return jsonify({'message': 'Email and password are required'}), 400
            
            # Authenticate user
            user = authenticate_user(email, password)
            if not user:
                return jsonify({'message': 'Invalid credentials'}), 401
            
            # Get Supabase client to access the session
            supabase = get_supabase_client()
            session = supabase.auth.get_session()
            
            if not session:
                return jsonify({'message': 'Failed to create session'}), 500
            
            # Extract tokens
            access_token = session.access_token
            refresh_token = session.refresh_token
            
            # Get client information for session tracking
            ip_address = request.remote_addr
            user_agent = request.headers.get('User-Agent')
            
            # Create session in our session management system
            from api.session_management import get_session_manager
            session_manager = get_session_manager()
            
            # Extract device info if available
            device_info = {}
            if user_agent:
                device_info['user_agent'] = user_agent
            
            # Create session record
            session_manager.create_session(
                user_id=user.id,
                refresh_token=refresh_token,
                ip_address=ip_address,
                user_agent=user_agent,
                device_info=device_info
            )
            
            # Get user roles and permissions
            try:
                from api.rbac import UserRoleManager
                user_roles = UserRoleManager.get_user_roles(user.id)
                user_permissions = UserRoleManager.get_user_permissions(user.id)
            except Exception as e:
                logger.warning(f"Failed to fetch user roles and permissions: {str(e)}")
                user_roles = []
                user_permissions = []
            
            # Create response
            response = jsonify({
                'message': 'Authentication successful',
                'user': {
                    'id': user.id,
                    'email': user.email,
                    'role': getattr(user, 'role', DEFAULT_ROLE),
                    'roles': user_roles,
                    'permissions': user_permissions
                }
            })
            
            # Set cookies if configured to do so
            if HTTP_ONLY_COOKIES:
                from api.session_security import set_secure_cookie
                
                # Set access token cookie with secure attributes
                response = set_secure_cookie(
                    response,
                    'access_token',
                    access_token,
                    max_age=JWT_EXPIRY,
                    httponly=True,
                    secure=SECURE_COOKIES,
                    samesite=SAME_SITE_COOKIES
                )
                
                # Set refresh token cookie with secure attributes
                response = set_secure_cookie(
                    response,
                    'refresh_token',
                    refresh_token,
                    max_age=JWT_REFRESH_EXPIRY,
                    httponly=True,
                    secure=SECURE_COOKIES,
                    samesite=SAME_SITE_COOKIES
                )
                
                # Rotate session on login
                from api.session_security import rotate_session
                rotate_session()
                
                # Store user role in session for privilege change detection
                session['user_role'] = getattr(user, 'role', DEFAULT_ROLE)
                session['user_id'] = user.id
            
            return response, 200
        except Exception as e:
            logger.error(f"Login error: {str(e)}", exc_info=True)
            return jsonify({'message': f'Login error: {str(e)}'}), 500

    @app.route('/auth/logout', methods=['POST'])
    @jwt_required
    def logout():
        try:
            # Get current user
            user = get_current_user()
            if not user:
                return jsonify({'message': 'User not authenticated'}), 401
            
            # Get refresh token from cookies or request body
            refresh_token = None
            if HTTP_ONLY_COOKIES:
                refresh_token = request.cookies.get('refresh_token')
            
            if not refresh_token:
                data = request.get_json() or {}
                refresh_token = data.get('refresh_token')
            
            # Revoke the session in our session management system
            from api.session_management import get_session_manager
            session_manager = get_session_manager()
            
            if refresh_token:
                # Validate and get session
                is_valid, session = session_manager.validate_session(refresh_token)
                if is_valid and session:
                    # Revoke the specific session
                    session_manager.revoke_session(session['id'], reason='logout')
            
            # Get Supabase client
            supabase = get_supabase_client()
            
            # Sign out from Supabase
            supabase.auth.sign_out()
            
            # Create response
            response = jsonify({'message': 'Logout successful'})
            
            # Clear cookies if they were set
            if HTTP_ONLY_COOKIES:
                response.delete_cookie('access_token', path='/', domain=None, secure=SECURE_COOKIES, samesite=SAME_SITE_COOKIES)
                response.delete_cookie('refresh_token', path='/', domain=None, secure=SECURE_COOKIES, samesite=SAME_SITE_COOKIES)
                
                # Clear and rotate session on logout
                session.clear()
                from api.session_security import rotate_session
                rotate_session()
            
            # Log the logout activity
            if user:
                session_manager.log_session_activity(
                    user_id=user.get('id'),
                    action='logout',
                    status='success',
                    ip_address=request.remote_addr,
                    user_agent=request.headers.get('User-Agent')
                )
            
            return response, 200
        except Exception as e:
            logger.error(f"Logout error: {str(e)}", exc_info=True)
            return jsonify({'message': f'Logout error: {str(e)}'}), 500

    @app.route('/auth/refresh', methods=['POST'])
    def refresh():
        try:
            # Get refresh token from request
            refresh_token = None
            
            # Try to get from cookies first
            if HTTP_ONLY_COOKIES:
                refresh_token = request.cookies.get('refresh_token')
            
            # If not in cookies, try to get from request body
            if not refresh_token:
                data = request.get_json() or {}
                refresh_token = data.get('refresh_token')
            
            if not refresh_token:
                return jsonify({'message': 'Refresh token is required'}), 400
            
            # Validate the session in our session management system
            from api.session_management import get_session_manager
            session_manager = get_session_manager()
            
            is_valid, session = session_manager.validate_session(refresh_token)
            if not is_valid:
                return jsonify({'message': 'Invalid or expired refresh token'}), 401
            
            # Get Supabase client
            supabase = get_supabase_client()
            
            # Refresh session
            response = supabase.auth.refresh_session(refresh_token)
            
            if not response or not response.session:
                return jsonify({'message': 'Failed to refresh session'}), 401
            
            # Extract new tokens
            new_access_token = response.session.access_token
            new_refresh_token = response.session.refresh_token
            
            # Rotate the refresh token in our session management system
            session_manager.rotate_refresh_token(
                old_refresh_token=refresh_token,
                new_refresh_token=new_refresh_token,
                ip_address=request.remote_addr,
                user_agent=request.headers.get('User-Agent')
            )
            
            # Create response
            api_response = jsonify({
                'message': 'Session refreshed',
                'token': {
                    'access_token': new_access_token,
                    'refresh_token': new_refresh_token if not HTTP_ONLY_COOKIES else None
                }
            })
            
            # Set cookies if configured to do so
            if HTTP_ONLY_COOKIES:
                from api.session_security import set_secure_cookie
                
                # Set access token cookie with secure attributes
                api_response = set_secure_cookie(
                    api_response,
                    'access_token',
                    new_access_token,
                    max_age=JWT_EXPIRY,
                    httponly=True,
                    secure=SECURE_COOKIES,
                    samesite=SAME_SITE_COOKIES
                )
                
                # Set refresh token cookie with secure attributes
                api_response = set_secure_cookie(
                    api_response,
                    'refresh_token',
                    new_refresh_token,
                    max_age=JWT_REFRESH_EXPIRY,
                    httponly=True,
                    secure=SECURE_COOKIES,
                    samesite=SAME_SITE_COOKIES
                )
                
                # Rotate session on token refresh
                from api.session_security import rotate_session
                rotate_session()
            
            return api_response, 200
        except Exception as e:
            logger.error(f"Token refresh error: {str(e)}", exc_info=True)
            return jsonify({'message': f'Token refresh error: {str(e)}'}), 401

    @app.route('/auth/validate', methods=['GET'])
    @jwt_required
    def validate_token():
        """Validate the current token and return user information"""
        try:
            user = get_current_user()
            if not user:
                return jsonify({'message': 'Invalid token'}), 401
            
            return jsonify({
                'message': 'Token is valid',
                'user': {
                    'id': user.get('id'),
                    'email': user.get('email'),
                    'role': user.get('role', DEFAULT_ROLE)
                }
            }), 200
        except Exception as e:
            logger.error(f"Token validation error: {str(e)}", exc_info=True)
            return jsonify({'message': f'Token validation error: {str(e)}'}), 401

    @app.route('/auth/mfa/initiate', methods=['POST'])
    def mfa_initiate():
        try:
            data = request.get_json()
            email = data.get('email')
            
            if not email:
                return jsonify({'message': 'Email is required'}), 400
            
            # Get Supabase client
            supabase = get_supabase_client()
            
            # Check if user exists
            try:
                # Find user by email
                response = supabase.table("users").select("id").eq("email", email).execute()
                
                if not response.data:
                    return jsonify({'message': 'User not found'}), 404
                
                user_id = response.data[0]['id']
            except Exception as e:
                logger.error(f"Error finding user: {str(e)}")
                return jsonify({'message': 'Error finding user'}), 500
            
            # Generate a random 6-digit code
            import random
            code = str(random.randint(100000, 999999))
            
            # Store the code in the database with expiration time (10 minutes)
            from datetime import datetime, timedelta
            expires_at = datetime.now() + timedelta(minutes=10)
            
            # Generate a unique challenge ID
            import uuid
            challenge_id = str(uuid.uuid4())
            
            # Store MFA challenge in database
            try:
                mfa_data = {
                    "id": challenge_id,
                    "user_id": user_id,
                    "code": code,
                    "expires_at": expires_at.isoformat(),
                    "verified": False,
                    "attempts": 0
                }
                
                response = supabase.table("mfa_challenges").insert(mfa_data).execute()
                
                if response.error:
                    logger.error(f"Error storing MFA challenge: {response.error}")
                    return jsonify({'message': 'Error initiating MFA challenge'}), 500
            except Exception as e:
                logger.error(f"Error storing MFA challenge: {str(e)}")
                return jsonify({'message': 'Error initiating MFA challenge'}), 500
            
            # In a real implementation, send the code to the user's email or phone
            # For this example, we'll just log it
            logger.info(f"MFA code for {email}: {code}")
            
            # Send email with code (mock implementation)
            try:
                from flask_mail import Message
                from flask import current_app
                
                msg = Message(
                    subject="Your CryoProtect MFA Code",
                    recipients=[email],
                    body=f"Your verification code is: {code}\n\nThis code will expire in 10 minutes.",
                    html=f"""
                    <h2>CryoProtect Authentication</h2>
                    <p>Your verification code is: <strong>{code}</strong></p>
                    <p>This code will expire in 10 minutes.</p>
                    """
                )
                
                if hasattr(current_app, 'mail'):
                    current_app.mail.send(msg)
                    logger.info(f"MFA code email sent to {email}")
                else:
                    logger.warning("Mail not configured, skipping email send")
            except Exception as e:
                logger.error(f"Error sending MFA code email: {str(e)}")
                # Continue even if email fails - we'll show the code in the response for testing
            
            return jsonify({
                'message': 'MFA challenge initiated',
                'mfa': {
                    'challenge_id': challenge_id,
                    'expires_at': expires_at.isoformat(),
                    'code': code if current_app.config.get('TESTING', False) else None  # Only include code in testing mode
                }
            }), 200
        except Exception as e:
            logger.error(f"MFA initiation error: {str(e)}", exc_info=True)
            return jsonify({'message': f'MFA initiation error: {str(e)}'}), 500
    
    @app.route('/auth/mfa/verify', methods=['POST'])
    def mfa_verify():
        try:
            data = request.get_json()
            challenge_id = data.get('challenge_id')
            code = data.get('code')
            
            if not challenge_id or not code:
                return jsonify({'message': 'Challenge ID and code are required'}), 400
            
            # Get Supabase client
            supabase = get_supabase_client()
            
            # Get the MFA challenge from the database
            try:
                response = supabase.table("mfa_challenges").select("*").eq("id", challenge_id).execute()
                
                if not response.data:
                    return jsonify({'message': 'Invalid challenge ID'}), 400
                
                challenge = response.data[0]
                
                # Check if challenge is expired
                from datetime import datetime
                expires_at = datetime.fromisoformat(challenge['expires_at'])
                if datetime.now() > expires_at:
                    return jsonify({'message': 'MFA challenge has expired'}), 400
                
                # Check if challenge is already verified
                if challenge['verified']:
                    return jsonify({'message': 'MFA challenge already verified'}), 400
                
                # Check if max attempts reached (5 attempts)
                if challenge['attempts'] >= 5:
                    return jsonify({'message': 'Maximum verification attempts reached'}), 400
                
                # Increment attempts
                supabase.table("mfa_challenges").update({"attempts": challenge['attempts'] + 1}).eq("id", challenge_id).execute()
                
                # Verify code
                if challenge['code'] != code:
                    return jsonify({'message': 'Invalid verification code'}), 400
                
                # Mark challenge as verified
                supabase.table("mfa_challenges").update({"verified": True}).eq("id", challenge_id).execute()
                
                # Get user
                user_id = challenge['user_id']
                user_response = supabase.table("users").select("*").eq("id", user_id).execute()
                
                if not user_response.data:
                    return jsonify({'message': 'User not found'}), 404
                
                user = user_response.data[0]
            except Exception as e:
                logger.error(f"Error verifying MFA challenge: {str(e)}")
                return jsonify({'message': 'Error verifying MFA challenge'}), 500
            
            # Create session for the user
            auth_response = supabase.auth.sign_in_with_password({
                "email": user['email'],
                "password": "dummy_password"  # This is a workaround - in a real implementation, we would use a different approach
            })
            
            if auth_response.error:
                logger.error(f"Error creating session after MFA: {auth_response.error}")
                return jsonify({'message': 'Error creating session after MFA verification'}), 500
            
            # Extract tokens
            access_token = auth_response.session.access_token
            refresh_token = auth_response.session.refresh_token
            
            # Get client information for session tracking
            ip_address = request.remote_addr
            user_agent = request.headers.get('User-Agent')
            
            # Create session in our session management system
            from api.session_management import get_session_manager
            session_manager = get_session_manager()
            
            # Extract device info if available
            device_info = {}
            if user_agent:
                device_info['user_agent'] = user_agent
            
            # Create session record
            session_manager.create_session(
                user_id=user['id'],
                refresh_token=refresh_token,
                ip_address=ip_address,
                user_agent=user_agent,
                device_info=device_info
            )
            
            # Get user roles and permissions
            try:
                from api.rbac import UserRoleManager
                user_roles = UserRoleManager.get_user_roles(user['id'])
                user_permissions = UserRoleManager.get_user_permissions(user['id'])
            except Exception as e:
                logger.warning(f"Failed to fetch user roles and permissions: {str(e)}")
                user_roles = []
                user_permissions = []
            
            # Create response
            api_response = jsonify({
                'message': 'MFA verification successful',
                'user': {
                    'id': user['id'],
                    'email': user['email'],
                    'role': user.get('role', DEFAULT_ROLE),
                    'roles': user_roles,
                    'permissions': user_permissions
                }
            })
            
            # Set cookies if configured to do so
            if HTTP_ONLY_COOKIES:
                api_response.set_cookie(
                    'access_token',
                    access_token,
                    max_age=JWT_EXPIRY,
                    httponly=True,
                    secure=SECURE_COOKIES,
                    samesite=SAME_SITE_COOKIES
                )
                api_response.set_cookie(
                    'refresh_token',
                    refresh_token,
                    max_age=JWT_REFRESH_EXPIRY,
                    httponly=True,
                    secure=SECURE_COOKIES,
                    samesite=SAME_SITE_COOKIES
                )
            
            return api_response, 200
        except Exception as e:
            logger.error(f"MFA verification error: {str(e)}", exc_info=True)
            return jsonify({'message': f'MFA verification error: {str(e)}'}), 500
            
    @app.route('/auth/register', methods=['POST'])
    def register():
        try:
            data = request.get_json()
            email = data.get('email')
            password = data.get('password')
            logger.info("Registration attempt for email: %s", email)
            
            if not email or not password:
                logger.warning("Registration failed: Email and password required.")
                return jsonify({'message': 'Email and password are required'}), 400
                
            supabase = get_supabase_client()
            response = supabase.auth.sign_up({
                "email": email,
                "password": password
            })
            
            if response.error:
                logger.error("Registration error for email %s: %s", email, response.error.message)
                return jsonify({'message': response.error.message}), 400
                
            logger.info("Registration successful for email: %s", email)
            return jsonify({
                'message': 'Registration successful',
                'user': {
                    'id': response.user.id,
                    'email': response.user.email
                }
            }), 201
        except Exception as e:
            logger.exception("Registration error for email %s: %s", email, str(e))
            return jsonify({'message': f'Registration error: {str(e)}'}), 500
            
    @app.route('/auth/reset-password', methods=['POST'])
    def reset_password():
        try:
            data = request.get_json()
            email = data.get('email')
            logger.info("Password reset requested for email: %s", email)
            
            if not email:
                logger.warning("Password reset failed: Email required.")
                return jsonify({'message': 'Email is required'}), 400
                
            supabase = get_supabase_client()
            response = supabase.auth.reset_password_for_email(
                email,
                options={
                    "redirect_to": f"{request.host_url}reset-password"
                }
            )
            
            if response.error:
                logger.error("Password reset error for email %s: %s", email, response.error.message)
                return jsonify({'message': response.error.message}), 400
                
            logger.info("Password reset email sent for: %s", email)
            return jsonify({'message': 'Password reset email sent'}), 200
        except Exception as e:
            logger.exception("Password reset error for email %s: %s", email, str(e))
            return jsonify({'message': f'Password reset error: {str(e)}'}), 500
            
    @app.route('/auth/update-password', methods=['POST'])
    @token_required
    def update_password():
        try:
            data = request.get_json()
            new_password = data.get('password')
            logger.info("Password update requested.")
            
            if not new_password:
                logger.warning("Password update failed: New password required.")
                return jsonify({'message': 'New password is required'}), 400
            
            # Get current user
            user = get_current_user()
            if not user:
                return jsonify({'message': 'User not authenticated'}), 401
            
            supabase = get_supabase_client()
            response = supabase.auth.update_user({
                "password": new_password
            })
            
            if response.error:
                logger.error("Password update error: %s", response.error.message)
                return jsonify({'message': response.error.message}), 400
            
            # Revoke all sessions for the user after password change
            from api.session_management import get_session_manager
            session_manager = get_session_manager()
            session_manager.revoke_all_user_sessions(user.get('id'), reason='password_change')
            
            # Rotate the current session for additional security
            from api.session_security import rotate_session
            rotate_session()
            
            # Log the password change
            session_manager.log_session_activity(
                user_id=user.get('id'),
                action='password_change',
                status='success',
                ip_address=request.remote_addr,
                user_agent=request.headers.get('User-Agent')
            )
            
            logger.info("Password updated successfully.")
            return jsonify({'message': 'Password updated successfully'}), 200
        except Exception as e:
            logger.exception("Password update error: %s", str(e))
            return jsonify({'message': f'Password update error: {str(e)}'}), 500
            
    @app.route('/auth/update-profile', methods=['POST'])
    @token_required
    def update_profile():
        try:
            data = request.get_json()
            user_data = data.get('user_data', {})
            logger.info("Profile update requested.")

            supabase = get_supabase_client()
            response = supabase.auth.update_user({
                "data": user_data
            })
            
            if response.error:
                logger.error("Profile update error: %s", response.error.message)
                return jsonify({'message': response.error.message}), 400
                
            logger.info("Profile updated successfully.")
            return jsonify({
                'message': 'Profile updated successfully',
                'user': {
                    'id': response.user.id,
                    'email': response.user.email,
                    'user_metadata': response.user.user_metadata
                }
            }), 200
        except Exception as e:
            logger.exception("Profile update error: %s", str(e))
            return jsonify({'message': f'Profile update error: {str(e)}'}), 500
    
    # Register session management routes
    from api.session_routes import register_session_routes
    register_session_routes(app)
    
    # Start session cleanup thread
    from api.session_utils import start_session_cleanup_thread
    start_session_cleanup_thread()
    
    apply_security_headers(app)
    return app

# Create the app instance
app = create_app()

# Middleware to check authentication for protected routes
def auth_middleware():
    @app.before_request
    def check_auth():
        # List of routes that require authentication
        protected_routes = [
            '/profile',
            '/molecules',
            '/mixtures',
            '/predictions',
            '/experiments',
            '/comparisons',
            '/protocol-designer',
            '/predictive-models',
            '/teams'
        ]
        
        # API routes that should be excluded from web authentication
        api_routes = ['/api/', '/auth/']

        # Allow unauthenticated GET access to molecules and mixtures pages
        public_get_routes = [
            '/molecules',
            '/molecules/rdkit',
            '/molecules/integrated',
            '/mixtures'
        ]
        if request.method == 'GET' and any(request.path == route for route in public_get_routes):
            return  # Allow public access to these GET routes

        # Check if the current route is protected and not an API route
        if any(request.path.startswith(route) for route in protected_routes) and \
           not any(request.path.startswith(route) for route in api_routes):
            # Check if user is authenticated
            user = get_current_user()
            if not user:
                # If the request accepts JSON, return a JSON response
                if request.headers.get('Accept', '').find('application/json') != -1:
                    return jsonify({'message': 'Authentication required'}), 401
                
                # Otherwise redirect to login page
                return redirect(url_for('login_page',
                                       next=request.path,
                                       error="Please login to access this page"))
    
    # Add context processor to make current_user available in templates
    @app.context_processor
    def inject_user():
        user = get_current_user()
        # Create a user object with is_authenticated property for templates
        class CurrentUser:
            def __init__(self, user_data):
                self.is_authenticated = user_data is not None
                self.data = user_data
                
            def __getattr__(self, name):
                if self.data and name in self.data:
                    return self.data[name]
                return None
                
        return {'current_user': CurrentUser(user)}

# Register the authentication middleware
auth_middleware()

# Add routes for the web interface
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/molecules')
def molecules():
    return render_template('molecules.html')

@app.route('/molecules/rdkit')
def molecules_rdkit():
    return render_template('molecules_rdkit.html')

@app.route('/molecules/integrated')
def molecules_integrated():
    return render_template('molecules_integrated.html')

@app.route('/mixtures')
def mixtures():
    return render_template('mixtures.html')

@app.route('/predictions')
def predictions():
    return render_template('predictions.html')

@app.route('/predictive-models')
def predictive_models():
    return render_template('predictive_models.html')

@app.route('/experiments')
def experiments():
    return render_template('experiments.html')

@app.route('/comparisons')
def comparisons():
    return render_template('comparisons.html')

@app.route('/login')
def login_page():
    # Pass Supabase configuration to the template
    return render_template('login.html',
                          supabase_url=app.config['SUPABASE_URL'],
                          supabase_key=app.config['SUPABASE_KEY'])

@app.route('/register')
def register_page():
    # Pass Supabase configuration to the template
    return render_template('register.html',
                          supabase_url=app.config['SUPABASE_URL'],
                          supabase_key=app.config['SUPABASE_KEY'])

@app.route('/profile')
def profile_page():
    # Pass Supabase configuration to the template
    return render_template('profile.html',
                          supabase_url=app.config['SUPABASE_URL'],
                          supabase_key=app.config['SUPABASE_KEY'])

@app.route('/reset-password')
def reset_password_page():
    # Pass Supabase configuration to the template
    return render_template('reset_password.html',
                          supabase_url=app.config['SUPABASE_URL'],
                          supabase_key=app.config['SUPABASE_KEY'])

@app.route('/protocol-designer')
def protocol_designer_page():
    # Pass Supabase configuration to the template
    return render_template('protocol_designer.html',
                          supabase_url=app.config['SUPABASE_URL'],
                          supabase_key=app.config['SUPABASE_KEY'])

@app.route('/dashboard')
def dashboard_page():
   # Pass Supabase configuration to the template
   return render_template('dashboard.html',
                          supabase_url=app.config['SUPABASE_URL'],
                          supabase_key=app.config['SUPABASE_KEY'])

@app.route('/api/v1/test/rate-limit')
def test_rate_limit():
   """
   Test endpoint to demonstrate rate limiting.
   This endpoint is limited to 5 requests per minute.
   """
   from api.rate_limiter import endpoint_rate_limit
   
   # Apply rate limit to this endpoint
   endpoint_rate_limit("5 per minute")(test_rate_limit)
   
   return jsonify({
       'status': 'success',
       'message': 'Rate limit test endpoint',
       'timestamp': datetime.now().isoformat()
   })

@app.route('/teams')
def teams_page():
   # Pass Supabase configuration to the template
   return render_template('teams.html',
                          supabase_url=app.config['SUPABASE_URL'],
                          supabase_key=app.config['SUPABASE_KEY'])

# Routes for shared items
@app.route('/share/<string:share_id>')
def shared_item(share_id):
    """Display a shared item."""
    try:
        logger.info("Accessing shared item with ID: %s", share_id)
        supabase = get_supabase_client()
        
        # Get share record
        response = supabase.table("shares").select("*").eq("id", share_id).execute()
        
        error_message, status_code = handle_supabase_error(response)
        if error_message:
            logger.warning("Error from Supabase when accessing share: %s", error_message)
            return render_template('error.html', error=error_message), status_code
        
        if not response.data:
            logger.warning("Shared item not found for ID: %s", share_id)
            return render_template('error.html', error="Shared item not found"), 404
        
        share = response.data[0]
        
        # Check if share has expired
        if share.get('expiration'):
            expiration = datetime.fromisoformat(share['expiration'])
            if datetime.now() > expiration:
                logger.info("Shared item expired for ID: %s", share_id)
                return render_template('error.html', error="This shared item has expired"), 410
        
        # For password-protected shares, show password form
        if share.get('password_protected'):
            logger.info("Shared item is password protected for ID: %s", share_id)
            return render_template('shared_item.html',
                                  share_id=share_id,
                                  password_protected=True,
                                  authenticated=False)
        
        # Get the shared item data
        try:
            response = supabase.rpc(
                "get_shared_item_data",
                {"p_share_id": share_id}
            ).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                logger.warning("Error from Supabase RPC for shared item data: %s", error_message)
                return render_template('error.html', error=error_message), status_code
            
            item_data = response.data
            
            # Get user info for "shared by"
            user_response = supabase.table("users").select("email").eq("id", share['user_id']).execute()
            shared_by = user_response.data[0]['email'] if user_response.data else "Unknown User"
            
            # Get item title based on data type
            item_title = "Shared Item"
            if share['data_type'] == 'molecules':
                item_title = item_data.get('name', 'Molecule')
            elif share['data_type'] == 'mixtures':
                item_title = item_data.get('name', 'Mixture')
            elif share['data_type'] == 'predictions':
                item_title = "Predictions"
            elif share['data_type'] == 'experiments':
                item_title = "Experiments"
            elif share['data_type'] == 'comparisons':
                item_title = "Comparisons"
            elif share['data_type'] == 'report':
                item_title = "Report"
            elif share['data_type'] == 'visualization':
                item_title = "Visualization"
            
            # Format expiration date if present
            expiration_formatted = None
            if share.get('expiration'):
                expiration = datetime.fromisoformat(share['expiration'])
                expiration_formatted = expiration.strftime('%Y-%m-%d %H:%M:%S')
            
            logger.info("Successfully accessed shared item with ID: %s", share_id)
            return render_template('shared_item.html',
                                  share_id=share_id,
                                  password_protected=False,
                                  authenticated=True,
                                  data_type=share['data_type'],
                                  item_id=share['item_id'],
                                  item_data=item_data,
                                  item_title=item_title,
                                  shared_by=shared_by,
                                  expiration=expiration_formatted)
            
        except Exception as e:
            logger.exception("Error getting shared item data for ID %s: %s", share_id, str(e))
            return render_template('error.html', error=f"Error getting shared item data: {str(e)}"), 500
        
    except Exception as e:
        logger.exception("Error accessing shared item with ID %s: %s", share_id, str(e))
        return render_template('error.html', error=f"Error accessing shared item: {str(e)}"), 500

@app.route('/share/<string:share_id>', methods=['POST'])
def access_shared_item(share_id):
    """Access a password-protected shared item."""
    try:
        # Get password from form
        password = request.form.get('password')
        if not password:
            return render_template('shared_item.html',
                                  share_id=share_id,
                                  password_protected=True,
                                  authenticated=False,
                                  error="Password is required"), 400
        
        supabase = get_supabase_client()
        
        # Get share record
        response = supabase.table("shares").select("*").eq("id", share_id).execute()
        
        error_message, status_code = handle_supabase_error(response)
        if error_message:
            return render_template('error.html', error=error_message), status_code
        
        if not response.data:
            return render_template('error.html', error="Shared item not found"), 404
        
        share = response.data[0]
        
        # Check if share has expired
        if share.get('expiration'):
            expiration = datetime.fromisoformat(share['expiration'])
            if datetime.now() > expiration:
                return render_template('error.html', error="This shared item has expired"), 410
        
        # Verify password
        if not share.get('password_protected'):
            return redirect(url_for('shared_item', share_id=share_id))
        
        # Simple hash for demo purposes - in production use a proper password verification
        password_hash = hashlib.sha256(password.encode()).hexdigest()
        
        if password_hash != share.get('password_hash'):
            return render_template('shared_item.html',
                                  share_id=share_id,
                                  password_protected=True,
                                  authenticated=False,
                                  error="Invalid password"), 401
        
        # Password is correct, redirect to GET route with a session flag
        session['authenticated_shares'] = session.get('authenticated_shares', [])
        if share_id not in session['authenticated_shares']:
            session['authenticated_shares'].append(share_id)
        
        return redirect(url_for('shared_item', share_id=share_id))
        
    except Exception as e:
        app.logger.error(f"Error accessing shared item: {str(e)}")
        return render_template('error.html', error=f"Error accessing shared item: {str(e)}"), 500

@app.route('/share/<string:share_id>/embed')
def embed_shared_item(share_id):
    """Display a shared item in an embeddable format."""
    # Similar to shared_item but with a simplified template for embedding
    try:
        supabase = get_supabase_client()
        
        # Get share record
        response = supabase.table("shares").select("*").eq("id", share_id).execute()
        
        error_message, status_code = handle_supabase_error(response)
        if error_message:
            return f"<div class='error'>Error: {error_message}</div>", status_code
        
        if not response.data:
            return "<div class='error'>Shared item not found</div>", 404
        
        share = response.data[0]
        
        # Check if share has expired
        if share.get('expiration'):
            expiration = datetime.fromisoformat(share['expiration'])
            if datetime.now() > expiration:
                return "<div class='error'>This shared item has expired</div>", 410
        
        # For password-protected shares, show a message
        if share.get('password_protected'):
            return "<div class='error'>This shared item is password protected. Please view it directly.</div>", 403
        
        # Get the shared item data
        try:
            response = supabase.rpc(
                "get_shared_item_data",
                {"p_share_id": share_id}
            ).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                return f"<div class='error'>Error: {error_message}</div>", status_code
            
            item_data = response.data
            
            # Return embedded content based on data type
            if share['data_type'] == 'visualization':
                # For visualizations, return the image directly
                return f"<img src='data:image/png;base64,{item_data}' alt='Shared Visualization' style='max-width:100%;'>", 200
            else:
                # For other types, return a simplified version of the shared_item template
                return render_template('shared_item_embed.html',
                                      data_type=share['data_type'],
                                      item_id=share['item_id'],
                                      item_data=item_data)
            
        except Exception as e:
            app.logger.error(f"Error getting shared item data: {str(e)}")
            return f"<div class='error'>Error: {str(e)}</div>", 500
        
    except Exception as e:
        app.logger.error(f"Error accessing shared item: {str(e)}")
        return f"<div class='error'>Error: {str(e)}</div>", 500

def get_service_role_auth():
    """
    Authenticate using the service role.
    
    This function creates a user object with the service role credentials,
    which allows bypassing Row Level Security (RLS) policies in Supabase.
    
    Returns:
        User object with service role credentials
    """
    # Create a user object with the service role
    class User:
        def __init__(self):
            self.id = USER_ID
            self.email = "service@cryoprotect.com"
            self.role = "service_role"
    
    return User()

def register_resources_with_docs(app, docs):
    """
    Register API resources with the API documentation system.
    
    Args:
        app: Flask application instance
        docs: FlaskApiSpec instance
    """
    from api.docs import register_resource
    
    # Import all resources
    from api.resources import (
        MoleculeResource, MoleculeListResource,
        MixtureResource, MixtureListResource,
        PredictionResource, ExperimentResource,
        ComparisonResource, PropertyComparisonResource
    )
    
    from api.batch_resources import (
        BatchOperationResource,
        BatchOperationStatusResource,
        BatchOperationCancellationResource
    )
    
    from api.user_profile_resources import UserProfileResource
    
    from api.rdkit_resources import (
        MoleculePropertyResource, MoleculeVisualizationResource,
        SubstructureSearchResource, SimilaritySearchResource,
        MoleculePropertyCalculationResource
    )
    
    from api.scoring_resources import (
        MoleculeScoreResource, MoleculeIdScoreResource,
        MixtureScoreResource, BatchScoringResource
    )
    
    # Register core resources with correct endpoint names
    register_resource(docs, MoleculeListResource, 'api.molecules')
    register_resource(docs, MoleculeResource, 'api.molecules')
    register_resource(docs, MixtureListResource, 'api.mixtures')
    register_resource(docs, MixtureResource, 'api.mixtures')
    register_resource(docs, PredictionResource, 'api.predictions')
    register_resource(docs, ExperimentResource, 'api.experiments')
    register_resource(docs, ComparisonResource, 'api.comparison_resource')
    register_resource(docs, PropertyComparisonResource, 'api.compare-properties')
    
    # Register batch operations
    register_resource(docs, BatchOperationResource, 'api.batch')
    
    # Register user profile
    register_resource(docs, UserProfileResource, 'api.user_profile')
    
    # Register RDKit resources
    register_resource(docs, MoleculePropertyResource, 'api.rdkit/properties')
    register_resource(docs, MoleculeVisualizationResource, 'api.rdkit/visualization')
    register_resource(docs, SubstructureSearchResource, 'api.rdkit/substructure')
    register_resource(docs, SimilaritySearchResource, 'api.rdkit/similarity')
    register_resource(docs, MoleculePropertyCalculationResource, 'api.molecules/calculate-properties')
    
    # Register Scoring resources
    register_resource(docs, MoleculeScoreResource, 'api.scoring/molecules')
    register_resource(docs, MoleculeIdScoreResource, 'api.molecules/score')
    register_resource(docs, MixtureScoreResource, 'api.mixtures/score')
    
    # Register batch scoring resource
    register_resource(docs, BatchScoringResource, 'api.scoring/batch')
    
    # Register batch operation status and cancellation resources
    register_resource(docs, BatchOperationStatusResource, 'api.batch')
    register_resource(docs, BatchOperationCancellationResource, 'api.batch')
    
    app.logger.info("API resources registered with documentation system")

if __name__ == '__main__':
    # Initialize the app without trying to authenticate
    # Authentication will happen during requests
    app.logger.info("Starting CryoProtect Analyzer API")
    
    # Log service role status
    if USE_SERVICE_ROLE:
        app.logger.info(f"Service role authentication is enabled with user ID: {USER_ID}")
    else:
        app.logger.info("Service role authentication is disabled")
    
# Run the application
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)