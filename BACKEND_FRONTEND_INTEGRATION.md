# Backend to Frontend Integration Plan

This document outlines the integration plan between our backend services (Heroku, fly.io), Netlify frontend, and Convex database.

## 1. Current Architecture

- **Frontend**: Netlify-hosted Next.js app (https://cryoprotect.netlify.app/)
- **Main Backend**: Heroku API (https://api.cryoprotect.app) - Points to cryoprotect-8030e4025428.herokuapp.com
- **RDKit Service**: fly.io (https://rdkit.cryoprotect.app) - Points to cryoprotect-rdkit.fly.dev
- **Database**: Transitioning to Convex (https://dynamic-mink-63.convex.cloud)

## 2. DNS Configuration

| Subdomain | Record Type | Points To |
|-----------|-------------|-----------|
| www.cryoprotect.app | CNAME | cryoprotect.netlify.app |
| api.cryoprotect.app | CNAME | cryoprotect-8030e4025428.herokuapp.com |
| rdkit.cryoprotect.app | CNAME | cryoprotect-rdkit.fly.dev |

## 3. Environment Variables

| Variable | Value | Service |
|----------|-------|---------|
| NEXT_PUBLIC_API_URL | https://api.cryoprotect.app/v1 | Netlify |
| NEXT_PUBLIC_RDKIT_API_URL | https://rdkit.cryoprotect.app | Netlify |
| NEXT_PUBLIC_CONVEX_URL | https://dynamic-mink-63.convex.cloud | Netlify |
| NEXT_PUBLIC_USE_CONVEX | true | Netlify |
| NEXT_PUBLIC_NETLIFY | true | Netlify |

## 4. Cross-Origin Resource Sharing (CORS)

### 4.1 Heroku Backend (simple_app.py)

Update CORS configuration to allow requests from all relevant domains:

```python
CORS(app, resources={
    r"/*": {"origins": [
        "http://localhost:3000",
        "https://cryoprotect.netlify.app",
        "https://www.cryoprotect.app",
        "https://rdkit.cryoprotect.app",
        os.environ.get("FRONTEND_URL", "*")
    ], "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"], 
       "allow_headers": ["Content-Type", "Authorization"]}
})
```

### 4.2 RDKit Service (fly.io)

Ensure the fly.io service has similar CORS settings:

```python
CORS(rdkit_app, resources={
    r"/*": {"origins": [
        "http://localhost:3000",
        "https://cryoprotect.netlify.app",
        "https://www.cryoprotect.app",
        "https://api.cryoprotect.app",
        os.environ.get("FRONTEND_URL", "*")
    ], "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"], 
       "allow_headers": ["Content-Type", "Authorization"]}
})
```

## 5. Netlify Configuration (netlify.toml)

Update or confirm the netlify.toml configuration:

```toml
[build]
  command = "npm run build:with-convex"
  publish = "out"

[build.environment]
  NEXT_PUBLIC_ENABLE_API_LOGGING = "true"
  NEXT_PUBLIC_USE_CONVEX = "true"
  NEXT_PUBLIC_NETLIFY = "true"
  NEXT_PUBLIC_API_URL = "https://api.cryoprotect.app/v1"
  NEXT_PUBLIC_RDKIT_API_URL = "https://rdkit.cryoprotect.app"
  NEXT_PUBLIC_CONVEX_URL = "https://dynamic-mink-63.convex.cloud"

# API redirects to Heroku backend
[[redirects]]
  from = "/api/*"
  to = "https://api.cryoprotect.app/:splat"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

# RDKit service redirects
[[redirects]]
  from = "/rdkit-api/*"
  to = "https://rdkit.cryoprotect.app/:splat"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

# Handle dynamic routes for molecules
[[redirects]]
  from = "/molecules/*"
  to = "/molecules/[id].html"
  status = 200
  force = false

# Handle dynamic routes for mixtures
[[redirects]]
  from = "/mixtures/*"
  to = "/mixtures/[id].html"
  status = 200
  force = false

# This ensures that all other client-side routes work correctly
[[redirects]]
  from = "/*"
  to = "/index.html"
  status = 200
```

## 6. Unified API Client

Create a unified API client that can work with all three services:

```typescript
// src/services/api-client.ts
import axios from 'axios';
import { convex } from '../convex/client';

// Main backend API client (Heroku)
export const apiClient = axios.create({
  baseURL: process.env.NEXT_PUBLIC_API_URL || 'https://api.cryoprotect.app/v1',
  timeout: 15000,
  headers: { 'Content-Type': 'application/json' },
  withCredentials: true,
});

// RDKit service client (fly.io)
export const rdkitClient = axios.create({
  baseURL: process.env.NEXT_PUBLIC_RDKIT_API_URL || 'https://rdkit.cryoprotect.app',
  timeout: 30000,
  headers: { 'Content-Type': 'application/json' },
});

// Export convex client for direct use in components
export { convex };

// Health check function for all services
export async function checkAllServices() {
  const results = {
    mainApi: false,
    rdkitService: false,
    convexDb: false
  };
  
  try {
    const apiResponse = await apiClient.get('/health');
    results.mainApi = apiResponse.status === 200;
  } catch (error) {
    console.error('Main API check failed:', error);
  }
  
  try {
    const rdkitResponse = await rdkitClient.get('/health');
    results.rdkitService = rdkitResponse.status === 200;
  } catch (error) {
    console.error('RDKit service check failed:', error);
  }
  
  try {
    // Use a simple Convex query to check connection
    const convexResponse = await convex.query({ name: 'healthCheck' });
    results.convexDb = convexResponse === 'ok';
  } catch (error) {
    console.error('Convex check failed:', error);
  }
  
  return results;
}
```

## 7. Convex Integration

### 7.1 Create Health Check Function

```typescript
// convex/health.ts
import { query } from './_generated/server';

export const healthCheck = query({
  args: {},
  handler: async () => {
    return 'ok';
  },
});
```

### 7.2 Data Service with Conditional Fetching

```typescript
// src/services/data-service.ts
import { apiClient, convex } from './api-client';

// Determine if we should use Convex based on environment
const useConvex = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';

export const dataService = {
  getMolecules: async (params) => {
    if (useConvex) {
      return convex.query({ name: 'molecules:list', args: params });
    } else {
      const response = await apiClient.get('/molecules', { params });
      return response.data;
    }
  },
  
  getMolecule: async (id) => {
    if (useConvex) {
      return convex.query({ name: 'molecules:get', args: { id } });
    } else {
      const response = await apiClient.get(`/molecules/${id}`);
      return response.data;
    }
  },
  
  // Add other data methods for mixtures, etc.
};
```

## 8. Implementation Checklist

1. [ ] Update CORS configuration on Heroku backend
2. [ ] Update CORS configuration on RDKit service
3. [ ] Confirm Netlify environment variables match DNS configuration
4. [ ] Verify netlify.toml configuration
5. [ ] Implement unified API client
6. [ ] Create Convex health check function
7. [ ] Implement data service with conditional fetching
8. [ ] Test connectivity between all services
9. [ ] Verify custom domain works with all services
10. [ ] Update frontend components to use the data service

## 9. Testing

### 9.1 Connection Testing Script

Create a simple testing script to verify connections between services:

```javascript
// scripts/test-connections.js
const axios = require('axios');

async function testConnections() {
  const endpoints = [
    { name: 'Frontend', url: 'https://cryoprotect.netlify.app/' },
    { name: 'API', url: 'https://api.cryoprotect.app/health' },
    { name: 'RDKit', url: 'https://rdkit.cryoprotect.app/health' },
    { name: 'Convex', url: 'https://dynamic-mink-63.convex.cloud' }
  ];

  for (const endpoint of endpoints) {
    try {
      const response = await axios.get(endpoint.url, { timeout: 5000 });
      console.log(`✅ ${endpoint.name}: Connected successfully (${response.status})`);
    } catch (error) {
      console.error(`❌ ${endpoint.name}: Connection failed - ${error.message}`);
    }
  }
}

testConnections();
```

### 9.2 CORS Testing

Add a CORS testing endpoint to both backend services:

```python
@app.route('/test-cors')
def test_cors():
    """Test CORS configuration."""
    origin = request.headers.get('Origin', 'Unknown')
    return jsonify({
        'status': 'success',
        'message': 'CORS test successful',
        'origin': origin,
        'cors_enabled': True
    })
```

## 10. Deployment Steps

1. Update backend services first:
   ```bash
   # Deploy Heroku changes
   git push heroku master
   
   # Deploy fly.io changes
   fly deploy
   ```

2. Update Convex schema if needed:
   ```bash
   npm run convex:deploy
   ```

3. Update and deploy frontend:
   ```bash
   cd frontend
   npm run convex:setup
   npm run deploy:netlify
   ```

4. Test connectivity between all services.

---

This integration plan provides a structured approach to ensure seamless connectivity between our Netlify frontend, Heroku and fly.io backend services, and Convex database.