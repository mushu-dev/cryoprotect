# Convex Integration Guide for CryoProtect

This document provides a comprehensive guide for integrating Convex with the CryoProtect application, including CORS configuration, frontend setup, and service communication.

## Overview

CryoProtect now uses a microservices architecture with these components:

1. **Frontend**: Next.js application deployed on Netlify (https://cryoprotect.netlify.app)
2. **API Backend**: Flask application deployed on Heroku (https://cryoprotect.herokuapp.com)
3. **RDKit Service**: Specialized microservice (https://rdkit.cryoprotect.app)
4. **Convex Database**: New database service (https://dynamic-mink-63.convex.cloud)

## What's Been Done

### 1. CORS Configuration for API Backend

We've implemented standardized CORS configuration for the API backend to allow cross-origin requests from the frontend, RDKit service, and Convex. The CORS configuration:

- Uses environment variables to define allowed origins
- Allows all necessary HTTP methods and headers
- Supports credentials for authenticated requests
- Includes debug and test endpoints for verification

The implementation is in `api/cors_config.py` and has been deployed to the Heroku app.

### 2. Convex Database Integration

We've updated the Heroku backend to use Convex as the database instead of Supabase:

- Created a Convex adapter that provides a Supabase-like interface for backward compatibility
- Implemented a database factory that can switch between Supabase and Convex based on configuration
- Set environment variables for Convex configuration (URL, deployment key, etc.)
- Added diagnostic endpoints to verify the database connection

### 3. Netlify Configuration

We've updated the Netlify configuration in `frontend/netlify.toml` to:

- Point to the correct Heroku app URL (https://cryoprotect-8030e4025428.herokuapp.com)
- Include Convex URL in environment variables
- Set up redirects for API endpoints
- Configure Content Security Policy headers to allow Convex connections

## Implementation Status

### 1. API Backend CORS Configuration ✅

The CORS configuration has been implemented and can be deployed with:

```bash
./deploy-cors-heroku.sh
```

This script:
- Sets the required environment variables on Heroku
- Deploys the updated app with CORS configuration
- Verifies the deployment

### 2. Convex API Functions ⏳

The Convex API functions have been created but deployment requires additional setup:

```bash
./deploy-convex-api.sh
```

Note: Before running this script, a proper Convex development environment needs to be set up and TypeScript issues resolved.

### 3. Convex Integration on Heroku ✅

The Convex adapter and database factory have been implemented and can be deployed with:

```bash
./deploy-convex-heroku.sh
```

This script:
- Sets the Convex environment variables on Heroku (URL, deployment key, etc.)
- Deploys the Convex adapter and database factory
- Updates the app to use Convex instead of Supabase
- Adds diagnostic endpoints to verify the database connection

### 4. Netlify Configuration ✅

The Netlify configuration has been updated and can be deployed with:

```bash
./deploy-netlify-config.sh
```

This script:
- Updates the netlify.toml file with Convex configuration
- Deploys the updated Netlify configuration
- Verifies the deployment

### 5. Complete Integration ⏳

The complete integration script will run all of the above steps:

```bash
./complete-integration.sh
```

Note: Currently this script is blocked by the Convex API function deployment.

### 6. RDKit Service (To Be Done) ⏳

The RDKit service needs to be deployed with CORS configuration:

```bash
# Install fly.io CLI if needed
curl -L https://fly.io/install.sh | sh

# Deploy RDKit service
cd rdkit-service
flyctl deploy
```

## Testing and Verification

### 1. Quick Verification

```bash
./verify-integration.sh
```

This script performs a quick check of all services and reports any issues.

### 2. Test Convex Adapter

```bash
python test-convex-adapter.py
```

This script tests the Convex adapter functionality:
- Basic CRUD operations
- Authentication methods
- Error handling

### 3. Full Integration Test

```bash
node test-full-integration.js
```

This script tests all aspects of the integration:
- Basic connectivity between all services
- CORS configuration
- Netlify redirects
- CSP headers
- End-to-end data flow

## Remaining Tasks

1. **Deploy RDKit Service**: 
   - Install fly.io CLI
   - Deploy the RDKit service with CORS configuration

2. **Complete Convex Authentication**:
   - Implement authentication in the Convex client provider
   - Test authentication flow
   
3. **Optimize Convex API Functions**:
   - Improve performance of query operations
   - Add indexing for frequently accessed fields
   - Implement caching for common queries

4. **Build Frontend Convex Integration**:
   - Update frontend components to use Convex directly where needed
   - Implement real-time updates using Convex subscriptions

## Troubleshooting

### CORS Issues

If CORS tests fail:
- Verify environment variables on Heroku match what's expected
- Check that CORS configuration is properly deployed
- Use the browser's developer tools to inspect CORS errors

### Connection Issues

If services cannot connect:
- Check that all services are running
- Verify DNS resolution for all domains
- Check network firewalls and security settings

### Authentication Issues

If authentication fails:
- Verify session handling in the Convex client provider
- Check that the frontend has the correct Convex URL
- Ensure proper credentials are being passed

### Database Errors

If database operations fail:
- Check the Heroku logs: `heroku logs --tail --app cryoprotect`
- Verify the Convex API functions work directly: `curl -X POST https://dynamic-mink-63.convex.cloud/api/query -H "Content-Type: application/json" --data '{"table":"molecules","limit":1}'`
- Ensure the Convex deployment key is set correctly

## References

- [CORS Documentation](https://developer.mozilla.org/en-US/docs/Web/HTTP/CORS)
- [Convex Auth Documentation](https://docs.convex.dev/auth)
- [Netlify Redirects](https://docs.netlify.com/routing/redirects/)
- [Content Security Policy](https://developer.mozilla.org/en-US/docs/Web/HTTP/CSP)