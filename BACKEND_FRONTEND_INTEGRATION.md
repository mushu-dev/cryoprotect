# Backend and Frontend Integration Plan

> **Status: Completed!** The Convex integration has been successfully completed. See [CONVEX_INTEGRATION_COMPLETE.md](CONVEX_INTEGRATION_COMPLETE.md) for the final report.

This document outlines the integration plan between the CryoProtect frontend (Netlify), backend services (Heroku and fly.io), and the Convex database.

## Current Architecture

The CryoProtect application consists of several components:

1. **Frontend**: Next.js React application deployed on Netlify
2. **Main Backend**: Flask application deployed on Heroku
3. **RDKit Service**: Specialized Flask microservice for molecular calculations deployed on fly.io
4. **Database**: Currently using Supabase/PostgreSQL, transitioning to Convex

### Service URLs and Domains

| Service | Current URL | Custom Domain |
|---------|------------|---------------|
| Frontend | https://cryoprotect.netlify.app | https://www.cryoprotect.app |
| Main API | https://cryoprotect-8030e4025428.herokuapp.com | https://api.cryoprotect.app |
| RDKit Service | https://cryoprotect-rdkit.fly.dev | https://rdkit.cryoprotect.app |
| Convex DB | https://dynamic-mink-63.convex.cloud | N/A |

## Integration Requirements

1. **DNS Configuration**:
   - Frontend: www.cryoprotect.app → Netlify
   - API: api.cryoprotect.app → Heroku
   - RDKit: rdkit.cryoprotect.app → fly.io

2. **CORS Configuration**:
   - All services must accept requests from all other services
   - Proper headers for handling credentials and methods

3. **API Client Design**:
   - Unified client for all backend services
   - Support for both REST APIs and Convex
   - Fallback mechanisms for offline operation

4. **Authentication Flow**:
   - Consistent token handling across services
   - Proper session management

## Implementation Plan

### 1. DNS and Environment Configuration

#### Update Netlify Environment Variables:

```
NEXT_PUBLIC_API_URL=https://api.cryoprotect.app/v1
NEXT_PUBLIC_RDKIT_API_URL=https://rdkit.cryoprotect.app
NEXT_PUBLIC_CONVEX_URL=https://dynamic-mink-63.convex.cloud
NEXT_PUBLIC_ENVIRONMENT=production
NEXT_PUBLIC_USE_CONVEX=true
```

#### Update netlify.toml:

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
```

### 2. Update CORS Configuration

#### Main Backend (Heroku):

```python
# Enable CORS with specific configuration for all services
CORS(app, resources={
    r"/*": {"origins": [
        "http://localhost:3000",                  # Local development
        "https://cryoprotect.netlify.app",        # Netlify site
        "https://www.cryoprotect.app",            # Custom domain
        "https://rdkit.cryoprotect.app",          # RDKit service
        os.environ.get("FRONTEND_URL", "*")       # Dynamic frontend URL from env
    ], "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"], 
       "allow_headers": ["Content-Type", "Authorization"]}
})
```

#### RDKit Service (fly.io):

```python
# Enable CORS with specific configuration for all services
CORS(app, resources={
    r"/*": {"origins": [
        "http://localhost:3000",                  # Local development
        "https://cryoprotect.netlify.app",        # Netlify site
        "https://www.cryoprotect.app",            # Custom domain
        "https://api.cryoprotect.app",            # API service
        os.environ.get("FRONTEND_URL", "*")       # Dynamic frontend URL from env
    ], "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"], 
       "allow_headers": ["Content-Type", "Authorization"]}
})
```

### 3. Unified API Client Implementation

Create a unified API client to handle both REST and Convex:

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
    const apiResponse = await apiClient.get('/health/connectivity');
    results.mainApi = apiResponse.status === 200;
  } catch (error) {
    console.error('Main API health check failed:', error);
  }
  
  try {
    const rdkitResponse = await rdkitClient.get('/health');
    results.rdkitService = rdkitResponse.status === 200;
  } catch (error) {
    console.error('RDKit service health check failed:', error);
  }
  
  try {
    // For Convex, we can check if the client is initialized
    results.convexDb = !!convex;
  } catch (error) {
    console.error('Convex connection check failed:', error);
  }
  
  return results;
}
```

### 4. Data Service with Conditional Fetching

Create a data service that can conditionally use either REST API or Convex:

```typescript
// src/services/data-service.ts
import { apiClient, rdkitClient, convex } from './api-client';
import { useQuery } from 'convex/react';
import { api } from '../../convex/_generated/api';

// Determine if we should use Convex
const useConvex = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';

export function useMolecules({ limit = 10, offset = 0 }) {
  // Use Convex if enabled
  if (useConvex) {
    return useQuery(api.molecules.query.list, { limit, skip: offset });
  }
  
  // Otherwise use REST API
  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  
  useEffect(() => {
    const fetchData = async () => {
      try {
        setLoading(true);
        const response = await apiClient.get('/molecules', {
          params: { limit, offset }
        });
        setData(response.data);
        setError(null);
      } catch (err) {
        setError(err);
        setData(null);
      } finally {
        setLoading(false);
      }
    };
    
    fetchData();
  }, [limit, offset]);
  
  return { data, loading, error };
}

export function useRDKitPropertyCalculation() {
  const calculateProperties = async (moleculeData, inputFormat = 'smiles') => {
    try {
      const response = await rdkitClient.post('/api/calculate-properties', {
        molecule_data: moleculeData,
        input_format: inputFormat
      });
      return response.data;
    } catch (error) {
      console.error('RDKit property calculation failed:', error);
      throw error;
    }
  };
  
  return { calculateProperties };
}
```

### 5. Testing and Troubleshooting

Create a script to test connections between all services:

```javascript
// scripts/test-connections.js
const axios = require('axios');

const config = {
  timeout: 10000,
  frontend: 'https://www.cryoprotect.app',
  api: 'https://api.cryoprotect.app',
  rdkit: 'https://rdkit.cryoprotect.app',
  convex: 'https://dynamic-mink-63.convex.cloud'
};

async function testConnection(name, url, endpoint = '', options = {}) {
  const fullUrl = `${url}${endpoint}`;
  console.log(`Testing connection to ${name}: ${fullUrl}`);
  
  try {
    const response = await axios.get(fullUrl, { 
      timeout: config.timeout,
      ...options
    });
    
    const status = response.status;
    const contentType = response.headers['content-type'] || '';
    
    console.log(`✅ ${name}: Connected successfully (${status})`);
    console.log(`   Content-Type: ${contentType}`);
    
    if (contentType.includes('application/json')) {
      console.log(`   Response data: ${JSON.stringify(response.data).substring(0, 100)}...`);
    }
    
    return true;
  } catch (error) {
    console.error(`❌ ${name}: Connection failed`);
    
    if (error.response) {
      console.error(`   Status: ${error.response.status}`);
      console.error(`   Message: ${error.message}`);
    } else if (error.request) {
      console.error(`   Network error: No response received`);
    } else {
      console.error(`   Error: ${error.message}`);
    }
    
    return false;
  }
}

async function runTests() {
  console.log('🧪 Starting connection tests...\n');
  
  // Test direct connections
  await testConnection('Frontend', config.frontend);
  await testConnection('API', config.api, '/health');
  await testConnection('RDKit Service', config.rdkit, '/health');
  
  // Test CORS
  await testConnection('API CORS Test', config.api, '/test-cors', {
    headers: { 'Origin': config.frontend }
  });
  
  await testConnection('RDKit CORS Test', config.rdkit, '/test-cors', {
    headers: { 'Origin': config.frontend }
  });
  
  // Test Netlify redirects
  await testConnection('Netlify API Redirect', config.frontend, '/api/health');
  await testConnection('Netlify RDKit Redirect', config.frontend, '/rdkit-api/health');
  
  console.log('\n🏁 Connection tests completed');
}

runTests();
```

### 6. Integration Checklist

- [x] Update DNS settings for all domains
- [x] Configure environment variables in Netlify
- [x] Update netlify.toml with proper redirects
- [x] Implement CORS fixes in both backend services
- [x] Create unified API client
- [x] Implement conditional data fetching service
- [x] Test connections between all services
- [x] Deploy updates to all platforms
- [x] Verify integrations in production

## Common Issues and Solutions

### CORS Issues

**Symptoms**: Browser console shows CORS errors when making requests.

**Solutions**:
1. Verify that the CORS configuration on all servers includes all domains.
2. Check that the correct headers are being sent (especially for authenticated requests).
3. Ensure Netlify redirects include appropriate headers.

### Connection Timeouts

**Symptoms**: Requests hang or time out.

**Solutions**:
1. Check that all services are running.
2. Verify DNS records are pointing to the correct services.
3. Confirm that Heroku and fly.io services have not scaled to zero.

### Authentication Problems

**Symptoms**: 401 Unauthorized errors.

**Solutions**:
1. Verify token handling in API client.
2. Check that cookies or tokens are being properly sent.
3. Ensure CORS is configured to allow credentials.

## Convex Transition Plan ✅

The Convex transition has been successfully completed. The following phases were executed:

1. **Phase 1: Dual Operation** ✅
   - Configured frontend to use both REST API and Convex
   - Implemented feature flag to control which one is used
   - Built data migration utilities

2. **Phase 2: Data Migration** ✅
   - Implemented data copying from PostgreSQL to Convex
   - Verified data integrity between systems
   - Ran validation tests

3. **Phase 3: Full Transition** ✅
   - Switched feature flag to use Convex exclusively
   - Maintained REST API for backward compatibility
   - Updated documentation and team training

## Deployment Sequence ✅

1. ✅ Updated DNS records
2. ✅ Deployed CORS fixes to backend services
3. ✅ Updated Netlify configuration
4. ✅ Tested connections
5. ✅ Deployed frontend with updated API client
6. ✅ Completed Convex transition

For full details on the implementation, see [CONVEX_INTEGRATION_COMPLETE.md](CONVEX_INTEGRATION_COMPLETE.md).