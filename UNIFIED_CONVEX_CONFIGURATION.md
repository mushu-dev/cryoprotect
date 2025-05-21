# Unified Convex Configuration Strategy

This document outlines the strategy for unifying Convex configuration across both the main and minimal frontends of the CryoProtect application. The goal is to create a consistent approach to Convex integration that works in both development and production environments.

## Current State

### Main Frontend (`/frontend`)
- Uses a more complex Convex configuration in `src/convex/client.ts`
- Conditionally wraps application with `ConvexClientProvider` in `src/pages/_app.js`
- Uses environment variable `NEXT_PUBLIC_USE_CONVEX` to determine if Convex is enabled
- References Convex URL from `NEXT_PUBLIC_CONVEX_URL` with fallback to "https://upbeat-parrot-866.convex.cloud"
- Includes advanced resilience features like timeouts and presence state

### Minimal Frontend (`/minimal-frontend`)
- Uses a simpler Convex client in `src/convex/client.js`
- Does not wrap application with ConvexProvider in `pages/_app.js`
- Explicitly disables Convex in production with `NEXT_PUBLIC_USE_CONVEX="false"`
- References Convex URL from `NEXT_PUBLIC_CONVEX_URL` with fallback to "https://primary-meerkat-478.convex.cloud"
- Includes health check but lacks advanced resilience features

### Backend
- Uses a Convex adapter (`database/convex_adapter.py`) that provides a Supabase-like interface
- Connection pooling with resilience features (`database/convex_pool.py`)
- Database factory that can switch between Supabase and Convex
- Uses environment variable `USE_CONVEX` to determine whether to use Convex

## Unified Configuration Approach

### 1. Standardize Environment Variables

Use consistent environment variables across all components:

- `NEXT_PUBLIC_USE_CONVEX`: Controls whether Convex is enabled in frontend components
- `NEXT_PUBLIC_CONVEX_URL`: URL of the Convex deployment to use
- `USE_CONVEX`: Controls whether Convex is enabled in backend components
- `CONVEX_URL`: URL of the Convex deployment for backend components
- `CONVEX_DEPLOYMENT_KEY`: Deployment key for authenticated backend operations

### 2. Create Unified Convex Client

Create a unified Convex client that can be used by both frontends:

```typescript
// /shared/convex/client.ts
import { ConvexReactClient } from "convex/react";

// Determine if Convex is enabled based on environment variable
const isConvexEnabled = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';

// Get Convex URL from environment variable with fallback to production
const convexUrl = process.env.NEXT_PUBLIC_CONVEX_URL || "https://upbeat-parrot-866.convex.cloud";

// Create and export the Convex client instance with production configuration
export const convex = new ConvexReactClient(convexUrl, {
  // Production configuration options
  unsavedChangesWarning: false, // Disable unsaved changes warning
  networkTimeout: 10000, // 10 seconds timeout for network requests
  initialPresence: {}, // Start with empty presence state
});

// Add health check method for diagnostics
convex.health = async function() {
  try {
    await this.getAuth();
    return true;
  } catch (error) {
    console.error("Convex health check failed:", error);
    return false;
  }
};

// Export convenience function to check if Convex is enabled
export const isEnabled = () => isConvexEnabled;
```

### 3. Create Unified ConvexProvider Component

Create a shared ConvexProvider component:

```tsx
// /shared/components/ConvexProvider.tsx
import React from 'react';
import { ConvexProvider } from "convex/react";
import { convex, isEnabled } from '../convex/client';

type ConvexWrapperProps = {
  children: React.ReactNode;
};

export function ConvexWrapper({ children }: ConvexWrapperProps) {
  // Only wrap with ConvexProvider if Convex is enabled
  if (!isEnabled()) {
    return <>{children}</>;
  }

  return (
    <ConvexProvider client={convex}>
      {children}
    </ConvexProvider>
  );
}
```

### 4. Update Frontend Application Entry Points

Update both frontends to use the shared ConvexWrapper:

```jsx
// Both frontends' _app.js
import { ConvexWrapper } from '../shared/components/ConvexProvider';

function MyApp({ Component, pageProps }) {
  return (
    <>
      <Head>
        <meta name="viewport" content="width=device-width, initial-scale=1" />
      </Head>
      <ConvexWrapper>
        <Component {...pageProps} />
      </ConvexWrapper>
    </>
  );
}
```

### 5. Implement Dynamic Data Fetching

Create a unified approach to data fetching that can work with either Convex or API backend:

```jsx
// /shared/hooks/useData.js
import { useQuery } from "convex/react";
import { api } from "../convex/api";
import { isEnabled } from "../convex/client";
import { useState, useEffect } from "react";

export function useData(collection, options = {}) {
  const { id, filters, limit } = options;
  
  // State for API data
  const [apiData, setApiData] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  
  // Try to use Convex if enabled
  const convexData = isEnabled() 
    ? useQuery(api[collection].get, { id, filters, limit }) 
    : undefined;
  
  // Fall back to API if Convex is not enabled or not available
  useEffect(() => {
    if (!isEnabled()) {
      setIsLoading(true);
      
      // Create API endpoint based on collection
      const endpoint = id 
        ? `/api/${collection}/${id}` 
        : `/api/${collection}`;
      
      // Fetch from API
      fetch(endpoint)
        .then(res => res.json())
        .then(data => {
          setApiData(data);
          setIsLoading(false);
        })
        .catch(err => {
          setError(err);
          setIsLoading(false);
        });
    }
  }, [collection, id, isEnabled()]);
  
  return {
    data: isEnabled() ? convexData : apiData,
    isLoading: isEnabled() ? convexData === undefined && !error : isLoading,
    error: isEnabled() ? undefined : error
  };
}
```

### 6. Shared Environment Configuration for Netlify Deployment

Create a shared Netlify configuration that works for both frontends:

```toml
# netlify.toml for both frontends
[build]
  publish = "out"
  command = "npm run build"

[build.environment]
  NEXT_PUBLIC_USE_CONVEX = "false"  # Disable Convex by default in production
  NODE_VERSION = "18"

# Enable Convex in specific contexts if needed
[context.develop.environment]
  NEXT_PUBLIC_USE_CONVEX = "true"
  NEXT_PUBLIC_CONVEX_URL = "https://hallowed-malamute-424.convex.cloud"

# Add necessary security headers for Convex
[[headers]]
  for = "/*"
  [headers.values]
    Content-Security-Policy = "default-src 'self'; connect-src 'self' https://cryoprotect-8030e4025428.herokuapp.com https://*.convex.cloud; script-src 'self' 'unsafe-inline' 'unsafe-eval'; style-src 'self' 'unsafe-inline';"
```

### 7. Configuration Script for Development

Create a script to set up development environment with the correct Convex configuration:

```bash
#!/bin/bash
# setup-convex-dev.sh

# Check for environment
if [ "$1" == "main" ]; then
  export USE_CONVEX=true
  export CONVEX_URL=https://hallowed-malamute-424.convex.cloud
  export NEXT_PUBLIC_USE_CONVEX=true
  export NEXT_PUBLIC_CONVEX_URL=https://hallowed-malamute-424.convex.cloud
  cd frontend
elif [ "$1" == "minimal" ]; then
  export USE_CONVEX=true
  export CONVEX_URL=https://hallowed-malamute-424.convex.cloud
  export NEXT_PUBLIC_USE_CONVEX=true
  export NEXT_PUBLIC_CONVEX_URL=https://hallowed-malamute-424.convex.cloud
  cd minimal-frontend
else
  echo "Please specify 'main' or 'minimal' as an argument"
  exit 1
fi

# Start development server
npm run dev
```

## Implementation Plan

1. Create the shared Convex configuration files
2. Implement the ConvexWrapper component
3. Update both frontends to use the shared configuration
4. Update Netlify configuration for production deployment
5. Test the unified configuration in development and production
6. Document the final approach for future maintenance

## Database Population Approach

The existing database population scripts should continue to work with minimal changes:

1. `direct_convex_population.py` - No changes needed; already uses environment variables
2. `database/migrations/convex_bridge.py` - No changes needed; uses the factory pattern
3. `unified_chembl_import.py` - May need updates to export to Convex

The key is ensuring that the environment variables are consistently set across all environments.