# Vercel to Netlify Migration Guide

This guide documents the process of migrating the CryoProtect frontend from Vercel to Netlify deployment.

## Why Migrate?

While Vercel provides excellent support for Next.js applications, we've moved to Netlify for:

1. Better integration with our broader CI/CD pipeline
2. More flexible domain and team management
3. Improved analytics capabilities
4. Support for our backend services
5. Cost optimization for our usage patterns

## Migration Process

### 1. Project Configuration

We've made the following changes to support Netlify deployment:

#### 1.1 Static Export Configuration

Next.js is now configured for static export in `next.config.js`:

```javascript
const nextConfig = {
  output: 'export',
  images: {
    unoptimized: true, // Required for static export
  },
  trailingSlash: true, // For cleaner URLs with Netlify
  // Other settings...
}
```

#### 1.2 Netlify Configuration

Created a `netlify.toml` file for deployment settings:

```toml
[build]
  command = "chmod +x ./minimal-build.sh && ./minimal-build.sh"
  publish = "out"
  
[build.environment]
  NEXT_PUBLIC_ENVIRONMENT = "production"
  NEXT_PUBLIC_ENABLE_API_LOGGING = "true"
  NEXT_PUBLIC_NETLIFY = "true"
  NEXT_PUBLIC_API_URL = "https://cryoprotect-8030e4025428.herokuapp.com/v1"
```

#### 1.3 Redirects & Rewrites

API proxying and client-side routing are handled by Netlify redirects:

```toml
# Handle API redirects to Heroku backend
[[redirects]]
  from = "/api/*"
  to = "https://cryoprotect-8030e4025428.herokuapp.com/api/:splat"
  status = 200
  force = true
  
# Handle client-side routing
[[redirects]]
  from = "/*"
  to = "/index.html"
  status = 200
```

#### 1.4 Dynamic Routes

Special handling for dynamic routes (see [NETLIFY_DYNAMIC_ROUTES.md](./NETLIFY_DYNAMIC_ROUTES.md) for details):

```toml
# Handle dynamic routes for molecules
[[redirects]]
  from = "/molecules/*"
  to = "/molecules/[id]/index.html"
  status = 200
  force = false
```

### 2. Code Adaptations

#### 2.1 Environment Detection

Updated environment detection to support both platforms:

```javascript
const isVercel = process.env.NEXT_PUBLIC_VERCEL === 'true';
const isNetlify = process.env.NEXT_PUBLIC_NETLIFY === 'true';
```

#### 2.2 Analytics Implementation

Implemented a dual-track analytics approach (see [NETLIFY_ANALYTICS_SETUP.md](./NETLIFY_ANALYTICS_SETUP.md)):

```javascript
// Use Vercel Analytics when on Vercel, Netlify Analytics when on Netlify
const analyticsId = isVercel ? true : undefined;
```

#### 2.3 Build Script Adaptation

Created `minimal-build.sh` with environment variable setup:

```bash
#!/bin/bash
# Build script for Netlify
yarn install
echo "Creating .env.local file..."
cat > .env.local << EOL
NEXT_PUBLIC_API_URL=${NEXT_PUBLIC_API_URL:-https://cryoprotect-8030e4025428.herokuapp.com/v1}
NEXT_PUBLIC_NETLIFY=true
NEXT_PUBLIC_USE_MOCK_DATA=false
EOL
yarn build
```

### 3. Deployment Steps

1. **Install Netlify CLI**:
   ```bash
   npm install -g netlify-cli
   ```

2. **Link to Netlify**:
   ```bash
   netlify login
   netlify link
   ```

3. **Test Deployment**:
   ```bash
   cd frontend
   netlify deploy --build
   ```

4. **Production Deployment**:
   ```bash
   netlify deploy --prod
   ```

### 4. Migration Verification

After deployment, verify:

1. All pages load correctly
2. Dynamic routes work properly
3. API requests succeed
4. Analytics is tracking properly
5. Authentication flows work
6. Images and assets display correctly

Use the validation script to check configuration:

```bash
node frontend/validate-netlify-build.js
```

### 5. DNS Configuration

Update DNS configuration to point to Netlify:

1. Log in to your domain registrar
2. Update the DNS records:
   - Add a CNAME record for `www` pointing to your Netlify site
   - Set ALIAS/ANAME record for the apex domain

### 6. Post-Migration Tasks

1. **Enable Netlify Analytics**:
   - Navigate to your site dashboard in Netlify
   - Go to Analytics and click "Enable Netlify Analytics"

2. **Set up Branch Deployments**:
   - Configure branch deployments for development/staging

3. **Set up Build Hooks**:
   - Add build hooks for automated deployments

4. **Monitor Performance**:
   - Use Netlify Analytics to monitor site performance
   - Compare performance metrics with Vercel deployment

## Troubleshooting Common Issues

### Dynamic Routes Not Working

If dynamic routes return 404 errors:
- Verify the redirects in `netlify.toml`
- Check that `generateStaticParams()` is implemented
- See [NETLIFY_DYNAMIC_ROUTES.md](./NETLIFY_DYNAMIC_ROUTES.md) for details

### Build Failures

Common build failures:
- TypeScript errors (resolved with `ignoreBuildErrors: true`)
- Missing environment variables
- Dynamic route export errors

### API Connectivity Issues

If API requests fail:
- Check CORS configuration in `netlify.toml`
- Verify API URL environment variables
- Inspect network requests in browser devtools

## Resources

- [Next.js Static Export Documentation](https://nextjs.org/docs/pages/building-your-application/deploying/static-exports)
- [Netlify Next.js Plugin](https://github.com/netlify/next-runtime)
- [Netlify Redirects Documentation](https://docs.netlify.com/routing/redirects/)
- [Netlify Analytics Documentation](https://docs.netlify.com/monitor-sites/analytics/)