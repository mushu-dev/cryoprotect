# Netlify Deployment Fix for CryoProtect with Convex

## Current Issue

The CryoProtect application deployed on Netlify is showing a static fallback page instead of the full Next.js application with Convex integration. The page displays:

```
This is a static fallback page. Please enable JavaScript for the full experience.
```

This is occurring because the current configuration is using a static export approach, which doesn't support the dynamic functionality needed for the Convex integration.

## Previous Issue (Resolved)

Previously, the Netlify deployment was failing with the error:

```
Error message
Command failed with exit code 1: npm run build
```

This was fixed by:
1. Creating a root `package.json` with redirect scripts
2. Updating the `netlify.toml` configuration 
3. Creating a deployment script

## New Solution for Dynamic App

We've created a new deployment script that fixes the static fallback page issue:

1. **Updates `next.config.js`**:
   - Removes `output: 'export'` to enable server-side rendering
   - Keeps other configuration options unchanged

2. **Updates `netlify.toml` configuration**:
   - Changes `publish = "out"` to `publish = ".next"`
   - Adds the necessary redirects for Next.js server-side rendering
   - Updates CORS and Content-Security-Policy headers for Convex and RDKit service

3. **Updates package.json build script**:
   - Changes from `next build && next export -o out` to just `next build`
   - Installs the Netlify Next.js plugin (`@netlify/plugin-nextjs`)

4. **Deploys with the Netlify Next.js integration**:
   - Uses the Netlify CLI to deploy with the updated configuration
   - Enables server-side rendering for dynamic functionality

## How to Deploy with the New Fix

Run the following commands to deploy:

```bash
cd /home/mushu/Projects/cryoprotect
chmod +x ./deploy-netlify-fix.sh
./deploy-netlify-fix.sh
```

## Key Configuration Changes

### Updated netlify.toml

```toml
[build]
  command = "npm run build"
  publish = ".next"  # Changed from "out" to ".next"

# This redirects all requests to the Next.js handler
[[redirects]]
  from = "/*"
  to = "/.netlify/functions/next"
  status = 200
```

### Updated next.config.js

```js
// Removed this line to enable server-side rendering
// output: 'export',

const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    unoptimized: true,
  },
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  }
};
```

### Updated RDKit Service URL

The configuration now correctly references the fly.io deployed RDKit service:

```
NEXT_PUBLIC_RDKIT_API_URL = "https://cryoprotect-rdkit.fly.dev"
```

## Why This Works

Next.js applications with dynamic functionality (like Convex integration) need server-side rendering, which the Netlify Next.js plugin provides. By removing the static export configuration and using the proper Next.js directory structure, we enable full functionality while still benefiting from Netlify's hosting and CI/CD capabilities.

## Testing the Deployment

After deploying, you should see the full dynamic application at https://cryoprotect.netlify.app instead of the static fallback page. All Convex functionality should work correctly.