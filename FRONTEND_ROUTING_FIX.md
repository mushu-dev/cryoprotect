# Frontend Routing Fix - Implementation Plan

## Issue Summary

The CryoProtect application is experiencing 404 errors when navigating to core routes like `/molecules` and `/mixtures`. Our investigation has identified several key issues:

1. **Missing Page Routes**: Navigation components reference routes that don't exist in the codebase
2. **Incomplete Static Export Configuration**: The `exportPathMap` in `next.config.js` is missing critical routes
3. **Mixed Routing Architecture**: The project contains both `/pages` and `/app` directories, suggesting a transition between routing systems
4. **Deployment Configuration**: Netlify deployment may not be correctly handling the static export

## Implementation Plan

### 1. Add Missing Page Routes

Create the following missing page files:

```jsx
// /frontend/src/pages/molecules/index.js
import React from 'react';
import Head from 'next/head';
import { useEffect, useState } from 'react';
import MoleculesList from '../../features/molecules/components/molecules-list';

export default function MoleculesPage() {
  const [isLoading, setIsLoading] = useState(true);
  const [molecules, setMolecules] = useState([]);
  const [error, setError] = useState(null);

  useEffect(() => {
    async function fetchMolecules() {
      try {
        setIsLoading(true);
        const response = await fetch(`${process.env.NEXT_PUBLIC_API_URL}/molecules`);
        if (!response.ok) {
          throw new Error(`API error: ${response.status}`);
        }
        const data = await response.json();
        setMolecules(data);
      } catch (err) {
        console.error('Error fetching molecules:', err);
        setError('Unable to load molecules. Please try again later.');
      } finally {
        setIsLoading(false);
      }
    }

    fetchMolecules();
  }, []);

  return (
    <>
      <Head>
        <title>Molecules - CryoProtect</title>
        <meta name="description" content="Browse and search for cryoprotectant molecules" />
      </Head>

      <div className="container mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold mb-6">Molecules</h1>
        
        {error && (
          <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded mb-4">
            {error}
          </div>
        )}
        
        {isLoading ? (
          <div className="text-center py-8">
            <p>Loading molecules...</p>
          </div>
        ) : (
          <MoleculesList molecules={molecules} />
        )}
      </div>
    </>
  );
}
```

```jsx
// /frontend/src/pages/mixtures/index.js
import React from 'react';
import Head from 'next/head';
import { useEffect, useState } from 'react';

export default function MixturesPage() {
  const [isLoading, setIsLoading] = useState(true);
  const [mixtures, setMixtures] = useState([]);
  const [error, setError] = useState(null);

  useEffect(() => {
    async function fetchMixtures() {
      try {
        setIsLoading(true);
        const response = await fetch(`${process.env.NEXT_PUBLIC_API_URL}/mixtures`);
        if (!response.ok) {
          throw new Error(`API error: ${response.status}`);
        }
        const data = await response.json();
        setMixtures(data);
      } catch (err) {
        console.error('Error fetching mixtures:', err);
        setError('Unable to load mixtures. Please try again later.');
      } finally {
        setIsLoading(false);
      }
    }

    fetchMixtures();
  }, []);

  return (
    <>
      <Head>
        <title>Mixtures - CryoProtect</title>
        <meta name="description" content="Browse and search for cryoprotectant mixtures" />
      </Head>

      <div className="container mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold mb-6">Mixtures</h1>
        
        {error && (
          <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded mb-4">
            {error}
          </div>
        )}
        
        {isLoading ? (
          <div className="text-center py-8">
            <p>Loading mixtures...</p>
          </div>
        ) : (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {mixtures.map((mixture) => (
              <div key={mixture.id} className="border rounded-lg p-4 shadow-sm hover:shadow-md transition-shadow">
                <h2 className="text-xl font-semibold mb-2">{mixture.name}</h2>
                <p className="text-gray-600 mb-3">{mixture.description}</p>
                <div className="text-sm text-gray-500">
                  Components: {mixture.components?.length || 0}
                </div>
              </div>
            ))}
          </div>
        )}
      </div>
    </>
  );
}
```

```jsx
// /frontend/src/pages/properties/index.js
import React from 'react';
import Head from 'next/head';

export default function PropertiesPage() {
  return (
    <>
      <Head>
        <title>Properties - CryoProtect</title>
        <meta name="description" content="Explore cryoprotectant properties and measurements" />
      </Head>

      <div className="container mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold mb-6">Properties</h1>
        
        <div className="prose max-w-none">
          <p>
            This page provides access to the physical and chemical properties of cryoprotectant molecules 
            in our database. You can explore different property categories and their relationships.
          </p>
          
          <h2 className="text-2xl font-semibold mt-8 mb-4">Property Categories</h2>
          
          <ul className="list-disc pl-5 space-y-2">
            <li>Physical Properties (melting point, boiling point, density)</li>
            <li>Chemical Properties (molecular weight, LogP, pKa)</li>
            <li>Cryoprotective Properties (glass transition temperature, vitrification)</li>
            <li>Toxicity Indicators (cytotoxicity, genotoxicity)</li>
            <li>Experimental Measurements (cell viability, membrane permeability)</li>
          </ul>
          
          <div className="bg-blue-50 border-l-4 border-blue-500 p-4 my-6">
            <p className="text-blue-700">
              Detailed property exploration features are coming soon. 
              Check back for interactive property visualization and comparison tools.
            </p>
          </div>
        </div>
      </div>
    </>
  );
}
```

### 2. Update Next.js Export Configuration

Modify `next.config.js` to include all routes in the `exportPathMap`:

```js
exportPathMap: async function() {
  const paths = {
    '/': { page: '/' },
    '/molecules': { page: '/molecules' },
    '/mixtures': { page: '/mixtures' },
    '/experiments': { page: '/experiments' },
    '/protocols': { page: '/protocols' },
    '/properties': { page: '/properties' },
  };
  
  // Generate experiment detail pages (1-6)
  for (let i = 1; i <= 6; i++) {
    paths[`/experiments/${i}`] = { 
      page: '/experiments/[id]',
      query: { id: i.toString() } 
    };
  }
  
  // Generate protocol detail pages (1-5)
  for (let i = 1; i <= 5; i++) {
    paths[`/protocols/${i}`] = { 
      page: '/protocols/[id]',
      query: { id: i.toString() } 
    };
  }
  
  // Add protocol create page
  paths['/protocols/create'] = { page: '/protocols/create' };
  
  return paths;
}
```

### 3. Consolidate Routing Architecture

To avoid confusion, we should commit to one routing approach. Since the project is already using the Pages Router (`/pages`) more extensively, we should:

1. Move any App Router (`/app`) functionality to Pages Router
2. Update the `next.config.js` to disable the App Router if needed:

```js
const nextConfig = {
  // Existing configuration...
  
  // Disable App Router
  experimental: {
    appDir: false,
  },
  
  // Export path map configuration...
};
```

### 4. Optimize Netlify Deployment Configuration

Update the `netlify.toml` file to improve static site deployment:

```toml
[build]
  publish = "frontend/out"
  command = "cd frontend && npm run build && npm run export"
  
[build.environment]
  NEXT_PUBLIC_ENVIRONMENT = "production"
  NEXT_PUBLIC_ENABLE_API_LOGGING = "true"
  NEXT_PUBLIC_NETLIFY = "true"
  NEXT_PUBLIC_API_URL = "https://cryoprotect-8030e4025428.herokuapp.com/v1"
  NEXT_PUBLIC_RDKIT_API_URL = "https://cryoprotect-rdkit.fly.dev"
  NEXT_PUBLIC_CONVEX_URL = "https://upbeat-parrot-866.convex.cloud"
  NEXT_PUBLIC_USE_CONVEX = "false"

# Add specific redirects for spa routing
[[redirects]]
  from = "/molecules/*"
  to = "/molecules/index.html"
  status = 200
  
[[redirects]]
  from = "/mixtures/*"
  to = "/mixtures/index.html"
  status = 200
  
[[redirects]]
  from = "/experiments/*"
  to = "/experiments/index.html"
  status = 200
  
[[redirects]]
  from = "/protocols/*"
  to = "/protocols/index.html"
  status = 200
  
[[redirects]]
  from = "/properties/*"
  to = "/properties/index.html"
  status = 200

# Existing API redirects...

# This handles client-side routing for all other routes
[[redirects]]
  from = "/*"
  to = "/index.html"
  status = 200
```

### 5. Add Export Script to package.json

Add a static export script to `frontend/package.json`:

```json
"scripts": {
  // Existing scripts...
  "export": "next export",
  "build-export": "next build && next export"
}
```

### 6. Create a Deployment Verification Test

Create a script to verify deployment and check all routes:

```js
// /frontend/scripts/verify-deployment.js
const fetch = require('node-fetch');

const BASE_URL = 'https://cryoprotect.app';
const ROUTES = [
  '/',
  '/molecules',
  '/mixtures',
  '/experiments',
  '/protocols',
  '/properties',
  '/api/health'
];

async function verifyRoutes() {
  console.log('Starting deployment verification...');
  
  for (const route of ROUTES) {
    const url = `${BASE_URL}${route}`;
    try {
      console.log(`Checking ${url}`);
      const response = await fetch(url);
      const status = response.status;
      
      if (status === 200) {
        console.log(`✅ ${route} - OK (${status})`);
      } else {
        console.error(`❌ ${route} - Failed (${status})`);
      }
    } catch (error) {
      console.error(`❌ ${route} - Error: ${error.message}`);
    }
  }
  
  console.log('Deployment verification complete');
}

verifyRoutes();
```

## Implementation Steps

1. Create missing page files for `/molecules`, `/mixtures`, and `/properties`
2. Update Next.js configuration to include all routes in `exportPathMap`
3. Consolidate routing approach to use Pages Router exclusively
4. Update Netlify configuration for improved static site handling
5. Add export script to package.json
6. Create and run deployment verification test
7. Deploy the changes
8. Verify all routes are accessible

## Expected Outcomes

After implementation, the application should:

1. Have functioning navigation to all core sections
2. Display appropriate content on each page
3. Properly handle static generation of all routes
4. Maintain proper API integration
5. Show no 404 errors for main navigation links

## Future Improvements

1. **API Fallback**: Implement proper loading and error states for when API is unavailable
2. **Route Monitoring**: Add monitoring to detect routing issues in production
3. **Type Safety**: Add TypeScript interfaces for API responses
4. **Pre-rendering Optimization**: Explore Incremental Static Regeneration (ISR) for pages with dynamic data
5. **Unified Routing Architecture**: Consider complete migration to App Router when ready