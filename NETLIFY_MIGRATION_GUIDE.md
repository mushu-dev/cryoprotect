# Netlify Migration Guide

This document outlines the steps taken to migrate the CryoProtect frontend from Vercel to Netlify while maintaining connectivity to the Heroku backend.

## Overview

The migration process involved:

1. Creating a Netlify configuration file (`netlify.toml`)
2. Modifying Next.js configuration for compatibility
3. Setting up build scripts compatible with Netlify
4. Configuring redirects to maintain API connectivity with Heroku
5. Adding API connectivity testing functionality
6. Resolving React dependency conflicts

## Configuration Files

### 1. netlify.toml

The main Netlify configuration file provides:

- Build settings (command, publish directory)
- Environment variables
- Redirect rules for API endpoints
- Security headers
- Client-side routing support

```toml
[build]
  base = "frontend"
  command = "chmod +x ./minimal-build.sh && ./minimal-build.sh"
  publish = ".next"

[build.environment]
  NEXT_PUBLIC_API_URL = "https://cryoprotect-8030e4025428.herokuapp.com/v1"
  NEXT_PUBLIC_USE_MOCK_DATA = "false"
  NEXT_PUBLIC_ENABLE_API_LOGGING = "true"
  NEXT_PUBLIC_ENVIRONMENT = "production"
  NEXT_TELEMETRY_DISABLED = "1"
  NODE_VERSION = "18"

[functions]
  node_bundler = "esbuild"
  external_node_modules = ["@netlify/next-runtime"]

# API redirects for connectivity to Heroku backend
[[redirects]]
  from = "/api/v1/health/connectivity"
  to = "https://cryoprotect-8030e4025428.herokuapp.com/api/v1/health/connectivity"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

[[redirects]]
  from = "/api/*"
  to = "https://cryoprotect-8030e4025428.herokuapp.com/api/:splat"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}
```

### 2. Next.js Configuration

The Next.js configuration (`next.config.js`) was updated to:

- Remove deprecated options
- Configure API rewrites
- Support Netlify deployment

```javascript
/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    domains: ['localhost', 'api.cryoprotect.app'],
    unoptimized: process.env.NODE_ENV !== 'production',
  },
  // API routes are handled by netlify.toml redirects
  async rewrites() {
    return [
      {
        source: '/api/:path*',
        destination: process.env.NEXT_PUBLIC_API_URL + '/:path*',
      },
    ]
  },
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  },
  // Keep commented out for now
  // output: 'export', 
}

module.exports = nextConfig
```

### 3. Build Scripts

Custom build scripts were created to handle the build process on Netlify:

**minimal-build.sh**:
```bash
#!/bin/bash
# Minimal build script for Netlify that uses yarn to handle dependency conflicts
set -ex

# Install yarn
npm install -g yarn

# Show versions
echo "Node version: $(node -v)"
echo "Yarn version: $(yarn -v)"

# Install dependencies with yarn
echo "Installing dependencies with yarn..."
yarn install

# Build Next.js
echo "Building Next.js app..."
yarn build

echo "Build completed successfully!"
```

## API Connectivity

Two important changes were made to ensure API connectivity:

1. Added an API connectivity endpoint to the backend:
```python
@app.route('/api/v1/health/connectivity')
def api_connectivity_check():
    """Simple API connectivity test endpoint for frontend-backend connection."""
    try:
        heroku_app_name = os.environ.get('HEROKU_APP_NAME', 'cryoprotect-8030e4025428')
        frontend_url = os.environ.get('VERCEL_FRONTEND_URL', 
                       os.environ.get('NETLIFY_URL', 'https://cryoprotect.netlify.app'))
        
        # Return connectivity information
        return jsonify({
            'status': 'connected',
            'api_version': os.environ.get('API_VERSION', '1.0.0'),
            'environment': os.environ.get('FLASK_ENV', 'production'),
            'deployment': {
                'backend': f"https://{heroku_app_name}.herokuapp.com",
                'frontend': frontend_url
            },
            'timestamp': str(datetime.now().isoformat()),
            'cors_enabled': True
        }), 200
    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': str(e),
            'timestamp': str(datetime.now().isoformat())
        }), 500
```

2. Enhanced the frontend connectivity test component to handle both Netlify and direct API connectivity:
```typescript
// In ApiConnectivityTest.tsx
useEffect(() => {
  const checkConnectivity = async () => {
    try {
      // Try the main endpoint
      const response = await apiClient.get('/api/v1/health/connectivity');
      
      setStatus({
        connected: response.data.status === 'connected',
        loading: false,
        error: null,
        details: response.data,
      });
    } catch (error: any) {
      // Try the fallback endpoint
      try {
        const netlifyApiUrl = process.env.NEXT_PUBLIC_API_URL || 
          'https://cryoprotect-8030e4025428.herokuapp.com/v1';
        
        const directResponse = await fetch(`${netlifyApiUrl}/health/connectivity`);
        
        if (directResponse.ok) {
          const data = await directResponse.json();
          setStatus({
            connected: data.status === 'connected',
            loading: false,
            error: null,
            details: data,
          });
        } else {
          throw new Error(`HTTP error! status: ${directResponse.status}`);
        }
      } catch (fallbackError: any) {
        setStatus({
          connected: false,
          loading: false,
          error: error.message || 'Failed to connect to API',
          details: null,
        });
      }
    }
  };

  checkConnectivity();
}, []);
```

## Dependency Management

React dependency conflicts were resolved by:

1. Using yarn instead of npm for building
2. Adding a `resolutions` field to package.json to pin React versions:

```json
"resolutions": {
  "react": "18.2.0",
  "react-dom": "18.2.0",
  "@types/react": "18.2.0",
  "@types/react-dom": "18.2.0"
}
```

## Testing

A simple bash script was created to test API connectivity:

```bash
#!/bin/bash
# Test API connectivity to the Heroku backend

echo "Testing API connectivity to Heroku backend..."
curl -s https://cryoprotect-8030e4025428.herokuapp.com/api/v1/health/connectivity | json_pp

echo "Testing API connectivity using simplified endpoint..."
curl -s https://cryoprotect-8030e4025428.herokuapp.com/v1/health/connectivity | json_pp
```

## Setup Instructions

### 1. Install Netlify CLI

```bash
npm install -g netlify-cli
netlify login
```

### 2. Initialize Your Netlify Site

```bash
cd /home/mushu/Projects/cryoprotect
netlify init
```

Follow the prompts to:
- Create & configure a new site
- Specify build command: `cd frontend && chmod +x ./minimal-build.sh && ./minimal-build.sh`
- Specify publish directory: `frontend/.next`

### 3. Set Up Environment Variables

```bash
cd /home/mushu/Projects/cryoprotect
netlify env:set NEXT_PUBLIC_API_URL "https://cryoprotect-8030e4025428.herokuapp.com/v1"
netlify env:set NEXT_PUBLIC_USE_MOCK_DATA "false"
netlify env:set NEXT_PUBLIC_ENABLE_API_LOGGING "true"
netlify env:set NEXT_PUBLIC_ENVIRONMENT "production"
```

### 4. Deploy Manually

```bash
cd /home/mushu/Projects/cryoprotect
netlify deploy --prod
```

### 5. Verify API Connectivity

```bash
# Run the API connectivity test script
chmod +x ./test-api-connectivity.sh
./test-api-connectivity.sh

# Visit the deployed Netlify site and check the API connectivity page
```

## Troubleshooting

### Common Issues

1. **CORS errors**:
   - Check headers in the netlify.toml redirects
   - Verify API endpoints have proper CORS headers

2. **API connectivity issues**:
   - Test direct API connectivity with the test script
   - Check redirection rules in netlify.toml
   - Verify environment variables are correctly set

3. **Build failures**:
   - Try using yarn instead of npm
   - Check for React dependency conflicts
   - Verify Node.js version compatibility (use Node 18)

4. **Next.js routing issues**:
   - Ensure redirects are correctly configured in netlify.toml
   - Check for client-side routing conflicts

### Differences from Vercel

Netlify handles Next.js apps differently from Vercel:

1. **Redirects**: API redirects need to be explicitly defined in `netlify.toml`
2. **Build Settings**: We use a custom build script for React dependency compatibility
3. **Environment Variables**: Set through Netlify UI or CLI
4. **Deployment Process**: Direct GitHub integration or via Netlify CLI

## Rollback Plan

If needed, you can roll back to Vercel by:

1. Re-enabling the Vercel integration in GitHub
2. Pushing to the Vercel-connected branch
3. Updating any backend services to point to the Vercel URL

## References

- [Netlify Docs for Next.js](https://docs.netlify.com/frameworks/next-js/)
- [Next.js on Netlify](https://www.netlify.com/with/nextjs/)
- [Netlify CLI Documentation](https://cli.netlify.com/)
- [Next.js Documentation](https://nextjs.org/docs)