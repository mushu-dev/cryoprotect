# Authentication Configuration
# Updated to use JWT-based authentication

# JWT Configuration
JWT_EXPIRY = 3600  # 1 hour in seconds
JWT_REFRESH_EXPIRY = 2592000  # 30 days in seconds

# Role Configuration
DEFAULT_ROLE = "user"
AVAILABLE_ROLES = ["user", "admin", "curator", "viewer"]

# RBAC Configuration
RBAC_ENABLED = True
ROLE_HIERARCHY = {
    "admin": ["curator", "user", "viewer"],  # Admin inherits all permissions
    "curator": ["user", "viewer"],           # Curator inherits user and viewer permissions
    "user": ["viewer"],                      # User inherits viewer permissions
    "viewer": []                             # Viewer has no inheritance
}

# Resource types for permissions
RESOURCE_TYPES = [
    "molecules", "mixtures", "experiments", "predictions",
    "projects", "teams", "users", "roles", "permissions",
    "system", "admin"
]

# Actions for permissions
ACTIONS = [
    "create", "read", "update", "delete", "manage", "all"
]

# Session Configuration
SESSION_TIMEOUT = 3600  # 1 hour in seconds
REFRESH_TOKEN_ROTATION = True  # Whether to rotate refresh tokens on use

# Security Configuration
SECURE_COOKIES = True  # Use secure cookies (HTTPS only)
HTTP_ONLY_COOKIES = True  # Use HTTP-only cookies for tokens
SAME_SITE_COOKIES = "Lax"  # SameSite cookie policy (Strict, Lax, None)
HMAC_SECRET_KEY = None  # Set this to a strong random value in production or use SECRET_KEY

# Whether to use JWT-based authentication (should always be True)
USE_JWT_AUTH = True

# Legacy configuration (kept for backward compatibility during transition)
# Will be removed after full migration
USE_SERVICE_ROLE = False
USER_ID = None
