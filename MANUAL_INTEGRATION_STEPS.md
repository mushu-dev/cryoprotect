# Manual Integration Steps for CryoProtect

This guide provides manual steps to complete the backend-frontend integration for CryoProtect after running the `deploy-minimal-integration.sh` script.

## 1. Update the Heroku App Code

You need to update your Heroku app code to include CORS support:

### Step 1: Add CORS configuration to your Flask app

Create a file named `cors_config.py` in your app:

```python
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
    """
    # Default origins to allow in development
    default_origins = [
        "http://localhost:3000",                 
        "http://127.0.0.1:3000",                 
        "https://cryoprotect.netlify.app",       
        "https://www.cryoprotect.app",           
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
    
    return app
```

### Step 2: Update your main Flask app

Modify your main `app.py` file to use the CORS configuration:

```python
from flask import Flask
from cors_config import configure_cors

app = Flask(__name__)
app = configure_cors(app)

# Rest of your Flask app code...

# Add health check endpoint if not already present
@app.route('/health', methods=['GET'])
def health_check():
    return {"status": "ok", "database": "connected", "timestamp": "now"}

# Ensure you have the API v1 health connectivity endpoint
@app.route('/api/v1/health/connectivity', methods=['GET'])
def connectivity():
    return {"status": "ok", "connectivity": "established", "service": "main-api"}
```

### Step 3: Deploy to Heroku

```bash
# From your app directory
git add .
git commit -m "Add CORS configuration for integration"
git push heroku master
```

## 2. Manually Update Netlify Configuration 

### Step 1: Update Environment Variables

In your Netlify dashboard:

1. Go to Site settings > Build & deploy > Environment
2. Add the following environment variables:
   - `NEXT_PUBLIC_API_URL`: https://cryoprotect.herokuapp.com/v1
   - `NEXT_PUBLIC_RDKIT_API_URL`: https://rdkit.cryoprotect.app
   - `NEXT_PUBLIC_CONVEX_URL`: https://dynamic-mink-63.convex.cloud
   - `NEXT_PUBLIC_USE_CONVEX`: true
   - `NEXT_PUBLIC_ENVIRONMENT`: production
   
### Step 2: Deploy your updated netlify.toml

```bash
# From your project root
cd frontend
netlify deploy --prod
```

## 3. Set Up RDKit Service (When fly.io is available)

### Step 1: Install fly.io CLI

```bash
curl -L https://fly.io/install.sh | sh
```

### Step 2: Log in to fly.io

```bash
fly auth login
```

### Step 3: Deploy the RDKit service

```bash
./deploy-rdkit-service.sh
```

## 4. Test the Integration

### Test API CORS

```bash
curl -I -H "Origin: https://cryoprotect.netlify.app" https://cryoprotect.herokuapp.com/test-cors
```

### Test Netlify Redirects

```bash
curl https://cryoprotect.netlify.app/api/health
```

### Full Test Suite

```bash
./test-connection.sh
```

## 5. Troubleshooting

If you encounter issues, refer to the detailed troubleshooting guide:

```bash
cat INTEGRATION_TROUBLESHOOTING_GUIDE.md
```

## 6. Verifying Convex Integration

Visit the Convex test page to verify the integration:

```
https://cryoprotect.netlify.app/convex-test
```