"""
CryoProtect Analyzer - Rate Limiting Configuration

This module provides configuration settings for the rate limiter.
"""

import os
from distutils.util import strtobool

# Helper function to parse boolean environment variables
def parse_bool(value, default=False):
    if value is None:
        return default
    try:
        return bool(strtobool(value))
    except (ValueError, AttributeError):
        return default

# Enable/disable rate limiting
RATE_LIMIT_ENABLED = parse_bool(os.environ.get('RATE_LIMIT_ENABLED', 'True'), True)

# Default rate limits (applied globally)
# Format: ["X per Y", ...] where X is the number of requests and Y is the time unit (second, minute, hour, day, month, year)
RATE_LIMIT_DEFAULT = [
    "1000 per day",    # 1000 requests per day
    "100 per hour",    # 100 requests per hour
    "20 per minute"    # 20 requests per minute
]

# Endpoint-specific rate limits
# Format: {endpoint: ["X per Y", ...], ...}
RATE_LIMIT_ENDPOINTS = {
    # Authentication endpoints - more permissive to allow login attempts
    "/api/v1/auth/login": ["20 per minute", "100 per hour"],
    "/api/v1/auth/register": ["10 per minute", "50 per hour"],
    
    # API endpoints with higher limits for common operations
    "/api/v1/molecules": ["50 per minute", "200 per hour"],
    "/api/v1/mixtures": ["50 per minute", "200 per hour"],
    
    # Computationally intensive endpoints with lower limits
    "/api/v1/rdkit/properties": ["10 per minute", "100 per hour"],
    "/api/v1/rdkit/visualization": ["10 per minute", "100 per hour"],
    "/api/v1/rdkit/substructure": ["5 per minute", "50 per hour"],
    "/api/v1/rdkit/similarity": ["5 per minute", "50 per hour"],
    
    # Batch operations with very low limits
    "/api/v1/batch": ["2 per minute", "20 per hour"],
}

# Exempt endpoints from rate limiting
# Format: [endpoint, ...]
RATE_LIMIT_EXEMPT = [
    "/api/v1/health",
    "/api/v1/docs",
    "/api/v1/swagger",
    "/api/v1/swagger-ui"
]

# Rate limit by IP address, user ID, or both
# Options: "ip", "user", "hybrid"
RATE_LIMIT_BY = os.environ.get("RATE_LIMIT_BY", "hybrid")

# Storage URL for rate limiting (Redis recommended for production)
# Format: "redis://host:port/db" or "memory://" for in-memory storage
RATE_LIMIT_STORAGE_URL = os.environ.get("RATE_LIMIT_STORAGE_URL", "memory://")

# Rate limiting strategy
# Options: "fixed-window", "moving-window", "fixed-window-elastic-expiry"
RATE_LIMIT_STRATEGY = os.environ.get("RATE_LIMIT_STRATEGY", "fixed-window")

# Enable/disable rate limit headers in responses
RATE_LIMIT_HEADERS_ENABLED = parse_bool(os.environ.get('RATE_LIMIT_HEADERS_ENABLED', 'True'), True)

# Custom rate limit for specific user roles
# Format: {role: ["X per Y", ...], ...}
# Note: This can be overridden by environment variables in a JSON format
# Example: RATE_LIMIT_ROLES='{"admin":["5000 per day"],"premium":["2000 per day"]}'
import json
DEFAULT_ROLE_LIMITS = {
    "admin": ["5000 per day", "500 per hour", "100 per minute"],
    "premium": ["2000 per day", "200 per hour", "40 per minute"],
    "basic": ["1000 per day", "100 per hour", "20 per minute"],
}

# Try to parse role limits from environment variable
try:
    env_role_limits = os.environ.get('RATE_LIMIT_ROLES')
    if env_role_limits:
        RATE_LIMIT_ROLES = json.loads(env_role_limits)
    else:
        RATE_LIMIT_ROLES = DEFAULT_ROLE_LIMITS
except (json.JSONDecodeError, TypeError):
    RATE_LIMIT_ROLES = DEFAULT_ROLE_LIMITS

# Retry-after value for rate limit exceeded responses (in seconds)
try:
    RATE_LIMIT_RETRY_AFTER = int(os.environ.get('RATE_LIMIT_RETRY_AFTER', '60'))
except (ValueError, TypeError):
    RATE_LIMIT_RETRY_AFTER = 60