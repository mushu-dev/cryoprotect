# Handling Dynamic Routes with Netlify Static Export

This document explains how we've configured Next.js dynamic routes to work with Netlify's static export deployment strategy.

## The Challenge

When using Next.js with `output: 'export'` (static site generation) and deploying to Netlify, dynamic routes like `/molecules/[id]` and `/mixtures/[id]` present a challenge because:

1. Static exports require all dynamic routes to be pre-rendered at build time
2. We can't pre-render all possible routes since we have thousands of molecules
3. Next.js requires the `generateStaticParams()` function to specify which routes to pre-render

## Our Solution

We've implemented a hybrid approach that combines:

1. Static pre-rendering for a small set of common routes
2. Client-side rendering for the actual dynamic content
3. Netlify redirects to handle non-pre-rendered routes

### 1. Modified `next.config.js`

We've configured Next.js to use static export with settings optimized for Netlify:

```javascript
// next.config.js
const nextConfig = {
  output: 'export',
  images: {
    unoptimized: true, // Required for static export
  },
  trailingSlash: true, // For cleaner URLs with Netlify
  experimental: {
    // Settings to improve compatibility with static exports
    staticPageGenerationTimeout: 120,
    outputFileTracingRoot: './',
    outputFileTracingExcludes: {
      '*': [
        'node_modules/@swc/core-linux-x64-gnu',
        'node_modules/@swc/core-linux-x64-musl',
        'node_modules/@esbuild/linux-x64',
      ],
    },
  }
}
```

### 2. Pre-generating Key Routes

We've updated the `generateStaticParams()` functions to pre-generate a small set of common routes:

```typescript
// molecules/[id]/generateStaticParams.ts
export function generateStaticParams() {
  return [
    { id: 'placeholder' },
    { id: '962' },   // Glycerol
    { id: '176' },   // DMSO
    { id: '6276' },  // Ethylene Glycol 
    { id: '8857' }   // Propylene Glycol
  ]
}
```

### 3. Client-Side Data Fetching

Our dynamic pages use client-side data fetching with custom hooks:

```typescript
// Inside the page component
const params = useParams()
const id = params.id as string
const { data: molecule, isLoading, isError } = useMolecule(id)
```

This means the actual data is always fetched client-side, regardless of whether the page was pre-rendered or not.

### 4. Netlify Redirects

We've configured Netlify redirects to handle routes that weren't pre-rendered:

```toml
# netlify.toml
[[redirects]]
  from = "/molecules/*"
  to = "/molecules/[id]/index.html"
  status = 200
  force = false

[[redirects]]
  from = "/mixtures/*"
  to = "/mixtures/[id]/index.html"
  status = 200
  force = false
```

These redirects ensure that requests for any molecule or mixture ID are served the pre-rendered template, which then loads the correct data client-side.

## Validation

We've created a validation script (`validate-netlify-build.js`) that checks:

1. The correct configuration in `next.config.js`
2. The presence of `generateStaticParams()` functions in dynamic routes
3. The netlify.toml configuration
4. Required environment variables
5. Analytics implementation

Run this script before deploying to ensure everything is set up correctly:

```bash
node frontend/validate-netlify-build.js
```

## Deployment Process

1. Run the validation script to check configuration
2. Build the Next.js app with `minimal-build.sh`
3. Deploy the `out` directory to Netlify
4. The Netlify build plugin will handle the redirects and hosting

## Troubleshooting

If you encounter issues with dynamic routes:

1. Check the network tab in browser dev tools to see if the data is being fetched
2. Verify that the redirects in netlify.toml are correct
3. Ensure all dynamic route pages have `'use client'` at the top
4. Make sure the `generateStaticParams()` function is exporting the expected format
5. Verify that the client-side hooks are working correctly

## Testing Routes

After deployment, test these routes to verify everything works:

1. Pre-rendered routes (e.g., `/molecules/962`)
2. Non-pre-rendered routes (e.g., `/molecules/12345`)
3. Invalid routes should redirect to your 404 page