# CryoProtect Integration Troubleshooting Guide

This guide provides solutions for common issues encountered when integrating the frontend (Netlify), backend services (Heroku and fly.io), and Convex database.

## Connection Test Results

The connection test identified several integration issues:

1. **RDKit Service Connectivity** - The RDKit service at rdkit.cryoprotect.app is not responding
2. **Convex Connection** - The Convex endpoint is returning 404 errors
3. **CORS Configuration** - The CORS test endpoints are not available
4. **Netlify Redirects** - Redirects for API and health endpoints are not working
5. **Content Security Policy** - CSP headers are not being properly applied

## Step-by-Step Solutions

### 1. Fix RDKit Service Issues

The RDKit service needs to be deployed with proper CORS configuration:

```bash
# First, install the fly.io CLI
curl -L https://fly.io/install.sh | sh

# Make the script executable
chmod +x ./deploy-rdkit-service.sh

# Deploy the RDKit service
./deploy-rdkit-service.sh
```

After deployment, verify the service is working:

```bash
curl https://rdkit.cryoprotect.app/health
curl -H 'Origin: https://cryoprotect.netlify.app' https://rdkit.cryoprotect.app/test-cors
```

### 2. Fix Main API CORS Configuration

Deploy the updated CORS configuration to your Heroku app:

```bash
# Make the script executable
chmod +x ./deploy-integration-changes.sh

# Deploy only to Heroku
DEPLOY_NETLIFY=false ./deploy-integration-changes.sh
```

Verify the API is working with CORS:

```bash
curl https://cryoprotect.herokuapp.com/health
curl -H 'Origin: https://cryoprotect.netlify.app' https://cryoprotect.herokuapp.com/test-cors
```

### 3. Fix Netlify Configuration

Deploy the updated netlify.toml with correct redirects and CSP headers:

```bash
# Deploy only to Netlify
DEPLOY_HEROKU=false ./deploy-integration-changes.sh
```

Verify the Netlify redirects are working:

```bash
curl https://cryoprotect.netlify.app/api/health
curl https://cryoprotect.netlify.app/rdkit-api/health
```

### 4. Fix Convex Connection Issues

The Convex endpoint should not be directly accessible via HTTP GET requests. Instead, you should:

1. Ensure your Convex project is properly set up:

```bash
cd frontend
npm run convex:codegen
```

2. Verify the Convex URL in your environment variables:

```bash
netlify env:set NEXT_PUBLIC_CONVEX_URL "https://dynamic-mink-63.convex.cloud" -s cryoprotect
```

3. Test the integration via the Convex test page after deployment:

```bash
npm run build:with-convex
netlify deploy --prod
```

Then visit: https://cryoprotect.netlify.app/convex-test

### 5. Deploy Everything at Once

For a complete setup, deploy all services:

```bash
# Deploy RDKit service
./deploy-rdkit-service.sh

# Deploy Main API and Netlify frontend
./deploy-integration-changes.sh
```

## Common Issues and Solutions

### The RDKit service is still not responding

1. Verify the service is running on fly.io:
   ```bash
   fly status -a cryoprotect-rdkit
   ```

2. Check if the DNS record is correctly configured:
   ```bash
   dig rdkit.cryoprotect.app
   ```

3. Ensure the fly.io certificate is created:
   ```bash
   fly certs create rdkit.cryoprotect.app
   ```

### CORS tests still fail

1. Verify that both services have the CORS config applied:
   ```bash
   curl -I -H "Origin: https://cryoprotect.netlify.app" https://cryoprotect-8030e4025428.herokuapp.com/test-cors
   curl -I -H "Origin: https://cryoprotect.netlify.app" https://rdkit.cryoprotect.app/test-cors
   ```

2. Check for CORS headers in the response:
   ```
   Access-Control-Allow-Origin: https://cryoprotect.netlify.app
   Access-Control-Allow-Methods: GET, POST, PUT, DELETE, OPTIONS
   Access-Control-Allow-Headers: Content-Type, Authorization
   ```

3. If headers are missing, ensure the CORS config is correctly integrated into the main application code.

### Netlify redirects are not working

1. Verify the netlify.toml file has been correctly deployed:
   ```bash
   netlify sites:list
   ```

2. Check redirect rules in the Netlify dashboard.

3. Try clearing the Netlify cache:
   ```bash
   netlify plugins:install netlify-purge-cache-plugin
   netlify purge-cache
   ```

### Content Security Policy headers are missing

1. Verify the headers in the deployed site:
   ```bash
   curl -I https://cryoprotect.netlify.app
   ```

2. Check the headers configuration in netlify.toml.

3. Ensure that the netlify.toml file has been correctly deployed.

## Verifying Complete Integration

After fixing all issues, run the full connection test:

```bash
./test-connection.sh
```

All tests should pass. If any test fails, refer to the corresponding section in this guide.

## Additional Resources

- [CORS_CONFIGURATION_GUIDE.md](./CORS_CONFIGURATION_GUIDE.md)
- [CONVEX_README.md](./CONVEX_README.md)
- [BACKEND_FRONTEND_INTEGRATION_SUMMARY.md](./BACKEND_FRONTEND_INTEGRATION_SUMMARY.md)
- [fly.io Documentation](https://fly.io/docs/)
- [Heroku Documentation](https://devcenter.heroku.com/)
- [Netlify Documentation](https://docs.netlify.com/)
- [Convex Documentation](https://docs.convex.dev/)