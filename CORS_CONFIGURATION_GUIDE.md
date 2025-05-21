# CORS Configuration Guide for CryoProtect Services

This guide explains how to configure Cross-Origin Resource Sharing (CORS) for all CryoProtect services to enable proper communication between the frontend, backend API, RDKit service, and Convex database.

## Overview

The CryoProtect application consists of several distinct services that need to communicate with each other:

1. **Frontend**: Next.js application hosted on Netlify (https://cryoprotect.netlify.app)
2. **Main API**: Flask backend hosted on Heroku (https://cryoprotect-8030e4025428.herokuapp.com)
3. **RDKit Service**: Specialized microservice for molecular calculations on fly.io (https://rdkit.cryoprotect.app)
4. **Convex Database**: Real-time database and backend service (https://dynamic-mink-63.convex.cloud)

For these services to work correctly, CORS must be properly configured to allow cross-origin requests between them.

## Main API Configuration (Heroku)

### Step 1: Import the CORS configuration

In your main Flask application file (`app.py`), import and use the CORS configuration module:

```python
from api.cors_config import configure_cors

# Create Flask app
app = Flask(__name__)

# Configure CORS
configure_cors(app)

# Rest of your Flask configuration...
```

### Step 2: Set environment variables

Set the following environment variables on Heroku:

```bash
heroku config:set FRONTEND_URL=https://cryoprotect.netlify.app -a cryoprotect-8030e4025428
heroku config:set RDKIT_SERVICE_URL=https://rdkit.cryoprotect.app -a cryoprotect-8030e4025428
heroku config:set CONVEX_URL=https://dynamic-mink-63.convex.cloud -a cryoprotect-8030e4025428
```

If you need to allow additional origins, you can set them as a comma-separated list:

```bash
heroku config:set ALLOWED_ORIGINS=https://example.com,https://dev.example.com -a cryoprotect-8030e4025428
```

## RDKit Service Configuration (fly.io)

### Step 1: Import the CORS configuration

In the RDKit service's main application file, import and use the CORS configuration module:

```python
from cors_config import configure_cors

# Create Flask app
app = Flask(__name__)

# Configure CORS
configure_cors(app)

# Rest of your Flask configuration...
```

### Step 2: Set environment variables

Set the following environment variables on fly.io:

```bash
fly secrets set FRONTEND_URL=https://cryoprotect.netlify.app
fly secrets set API_URL=https://cryoprotect-8030e4025428.herokuapp.com
fly secrets set CONVEX_URL=https://dynamic-mink-63.convex.cloud
```

## Netlify Configuration

Netlify is configured through the `netlify.toml` file which has already been updated to include:

1. The appropriate redirects for API and RDKit service
2. Content Security Policy headers to allow connections to all services
3. Environment variables to connect to each backend service

## Convex Configuration

Convex doesn't require special CORS configuration as it's designed to work with web applications by default. However, you should ensure that your Convex functions properly handle requests from your application.

## Testing CORS Configuration

Both the Main API and RDKit service have testing endpoints to verify CORS configuration:

```
https://cryoprotect-8030e4025428.herokuapp.com/test-cors
https://rdkit.cryoprotect.app/test-cors
```

Additionally, you can use the debug endpoints to see the current CORS configuration:

```
https://cryoprotect-8030e4025428.herokuapp.com/debug/cors
https://rdkit.cryoprotect.app/debug/cors
```

## Comprehensive Testing Script

You can also use the integration script to test all service connections:

```bash
./integrate-backend-services.sh
```

This script will:
1. Check the availability of all services
2. Verify CORS configuration is working correctly
3. Test redirects on the Netlify site
4. Provide a detailed report of any connection issues

## Troubleshooting

If you encounter CORS issues:

1. **Check browser console**: The browser will show CORS errors in the console
2. **Verify origins**: Ensure all service URLs are correctly listed in the CORS configuration
3. **Check headers**: Make sure the proper headers are being sent and received
4. **Test with curl**: Use curl to manually test CORS preflight requests:

```bash
curl -H "Origin: https://cryoprotect.netlify.app" \
     -H "Access-Control-Request-Method: POST" \
     -H "Access-Control-Request-Headers: X-Requested-With" \
     -X OPTIONS --verbose \
     https://cryoprotect-8030e4025428.herokuapp.com/api/v1/molecules
```

5. **Check environment variables**: Verify that all environment variables are correctly set