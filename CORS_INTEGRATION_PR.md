# CORS Configuration for Convex Integration

## Summary
This PR adds CORS configuration to support integration between the CryoProtect frontend, backend services, and the new Convex database. It provides standardized CORS handling across services to enable seamless communication.

## Changes
- Added standardized CORS configuration module (`cors_config.py`) for all services
- Updated Heroku app to use dynamic CORS origins from environment variables
- Added test endpoints for verifying CORS configuration
- Created deployment and verification scripts for CORS setup
- Added comprehensive documentation for the integration

## Test Plan
1. Run the deployment script: `./deploy-cors-heroku.sh`
2. Verify CORS configuration: `node verify-cors-config.js`
3. Test full connectivity: `node test-connection.js`
4. Verify frontend-Convex authentication with real-world operations

## Implementation Details

### CORS Configuration
The `cors_config.py` module provides a standardized approach to CORS that:
- Uses environment variables for configuration
- Includes all necessary origins (frontend, RDKit service, Convex)
- Supports all required HTTP methods and headers
- Provides test endpoints for verification

### Deployment Script
The `deploy-cors-heroku.sh` script automates:
- Setting required environment variables on Heroku
- Deploying updated app code with CORS configuration
- Verifying the deployment

### Environment Variables
The following environment variables control CORS behavior:
- `ALLOWED_ORIGINS`: Comma-separated list of allowed origins
- `FRONTEND_URL`: URL of the frontend application
- `RDKIT_SERVICE_URL`: URL of the RDKit service
- `CONVEX_URL`: URL of the Convex database

## Documentation
See `CORS_CONVEX_INTEGRATION.md` for detailed instructions on:
- CORS configuration for all services
- Deployment procedures
- Testing and verification
- Troubleshooting

## Security Considerations
- CORS is configured to only allow specific domains, not wildcard origins
- Content Security Policy headers are updated on Netlify to allow Convex connections
- Authentication flows are preserved with credential support
- No sensitive information is exposed in CORS headers