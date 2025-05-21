# CORS and Convex Integration for CryoProtect

This document provides instructions for configuring CORS to support the integration between the CryoProtect frontend, backend services (API and RDKit), and the Convex database.

## Overview

The CryoProtect application now uses several services which need to communicate with each other:

1. **Frontend**: Next.js application deployed on Netlify (https://cryoprotect.netlify.app)
2. **API Backend**: Flask application deployed on Heroku (https://cryoprotect.herokuapp.com)
3. **RDKit Service**: Specialized microservice (https://rdkit.cryoprotect.app)
4. **Convex Database**: New database service (https://dynamic-mink-63.convex.cloud)

For these services to communicate properly, CORS (Cross-Origin Resource Sharing) must be configured correctly on each service.

## Prerequisites

- Heroku CLI installed and authenticated
- Access to the Heroku app `cryoprotect`
- Node.js 18+ installed (for verification scripts)

## Deployment 

### 1. Deploy CORS Configuration to Heroku

We've provided a script to deploy the updated CORS configuration to the Heroku app:

```bash
./deploy-cors-heroku.sh
```

This script:
- Sets the required environment variables (ALLOWED_ORIGINS, CONVEX_URL, etc.)
- Deploys an updated version of the application with CORS configuration
- Updates the Heroku app to use the new configuration

### 2. Verify CORS Configuration

After deployment, verify that CORS is properly configured by running:

```bash
node verify-cors-config.js
```

This script tests CORS between all services and reports any issues.

## Manual Setup (If Needed)

If you prefer to configure CORS manually, follow these steps:

### 1. Set Environment Variables on Heroku

```bash
heroku config:set ALLOWED_ORIGINS="https://cryoprotect.netlify.app,https://rdkit.cryoprotect.app,https://dynamic-mink-63.convex.cloud" --app cryoprotect
heroku config:set CONVEX_URL="https://dynamic-mink-63.convex.cloud" --app cryoprotect
heroku config:set FRONTEND_URL="https://cryoprotect.netlify.app" --app cryoprotect
heroku config:set RDKIT_SERVICE_URL="https://rdkit.cryoprotect.app" --app cryoprotect
```

### 2. Update the API Code

Make sure the API code (`simple_app.py` or `app.py`) uses the CORS configuration from `cors_config.py`:

```python
from api.cors_config import configure_cors

app = Flask(__name__)
app = configure_cors(app)  # Apply CORS configuration
```

## CORS Configuration Details

The CORS configuration in `cors_config.py` is designed to:

1. Allow requests from the frontend, RDKit service, and Convex
2. Support all necessary HTTP methods (GET, POST, PUT, DELETE, OPTIONS)
3. Include required headers like Content-Type and Authorization
4. Support credentials for authenticated requests
5. Provide test endpoints for CORS verification

```python
def configure_cors(app):
    # Default origins to allow
    default_origins = [
        "http://localhost:3000",  # Local development
        "https://cryoprotect.netlify.app",  # Netlify deployment
        "https://www.cryoprotect.app",  # Custom domain
    ]
    
    # Add RDKit service URL
    rdkit_service_url = os.environ.get("RDKIT_SERVICE_URL", "https://rdkit.cryoprotect.app")
    if rdkit_service_url:
        default_origins.append(rdkit_service_url)
    
    # Add Convex URL
    convex_url = os.environ.get("CONVEX_URL", "https://dynamic-mink-63.convex.cloud")
    if convex_url:
        default_origins.append(convex_url)
    
    # Get additional origins from environment
    additional_origins = os.environ.get("ALLOWED_ORIGINS", "")
    if additional_origins:
        default_origins.extend(additional_origins.split(","))
    
    # Configure CORS
    CORS(app, resources={
        r"/*": {
            "origins": default_origins,
            "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
            "allow_headers": ["Content-Type", "Authorization", "X-API-Key", "X-Requested-With"],
            "supports_credentials": True,
            "max_age": 600
        }
    })
    
    # Add test endpoints
    # ...
```

## Netlify Configuration

The Netlify configuration needs to include proper CSP (Content Security Policy) headers to allow connections to all services:

```toml
[[headers]]
  for = "/*"
  [headers.values]
    Content-Security-Policy = "default-src 'self'; connect-src 'self' https://cryoprotect.herokuapp.com https://rdkit.cryoprotect.app https://*.convex.cloud https://dynamic-mink-63.convex.cloud; script-src 'self' 'unsafe-inline' 'unsafe-eval'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; font-src 'self' data:; object-src 'none'; base-uri 'self'; form-action 'self'; frame-ancestors 'self';"
```

## Testing Connection

Use the provided JavaScript test script to verify that all services can communicate properly:

```bash
node test-connection.js
```

This script tests:
1. Basic connectivity between all services
2. CORS configuration for all services
3. Netlify redirects
4. CSP headers

## Troubleshooting

### CORS Issues

If CORS tests fail, check:

1. Environment variables are set correctly on Heroku
2. The CORS configuration is deployed to the app
3. The server is properly configured to use the CORS configuration

### Convex Connection Issues

For Convex-specific issues:

1. Verify that the Convex URL is included in CORS allowed origins
2. Check that the CSP headers on Netlify include the Convex domain
3. Ensure the frontend is properly configured to authenticate with Convex

### Netlify Configuration Issues

If Netlify configuration is problematic:

1. Verify the netlify.toml file is correctly configured
2. Check that CSP headers include all required domains
3. Make sure redirects are properly set up for API endpoints