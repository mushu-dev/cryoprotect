# CryoProtect Frontend Deployment Guide

This document provides detailed instructions for deploying the CryoProtect frontend application to Vercel. It includes configuration steps, environment variables, and troubleshooting tips.

## Deployment Checklist

Before deploying, ensure you have the following ready:

- [ ] A Vercel account with access to the CryoProtect team
- [ ] The Vercel CLI installed globally (`npm i -g vercel`)
- [ ] Access to the production API endpoint
- [ ] A secure `NEXTAUTH_SECRET` value (either generated or provided)

## Environment Variables

The application requires the following environment variables to function properly:

| Variable | Description | Example Value |
|----------|-------------|---------------|
| `NEXT_PUBLIC_API_URL` | URL of the backend API | `https://api.cryoprotect.app/v1` |
| `NEXTAUTH_URL` | URL of the deployed frontend | `https://www.cryoprotect.app` |
| `NEXTAUTH_SECRET` | Secret key for NextAuth.js | (32+ character random string) |

> **Important:** The `NEXTAUTH_SECRET` must remain consistent across all deployments to maintain session compatibility. If you change this value, all user sessions will be invalidated.

## Deployment Methods

### Method 1: Using the Deployment Script

We've created a deployment script that handles all the necessary configuration:

1. Ensure you're logged in to Vercel CLI:
   ```bash
   vercel login
   ```

2. Run the deployment script:
   ```bash
   ./deploy-to-vercel.sh
   ```

3. You'll be prompted to choose a build method:
   - **Production build (recommended)**: Uses custom configs to avoid SSR errors
   - **Static export**: Generates a fully static site with client-side rendering only
   - **Default Next.js build**: Standard Next.js build (may have SSR errors)

4. The script will:
   - Generate a secure `NEXTAUTH_SECRET` if none is provided
   - Display this secret in the terminal (save it for future reference)
   - Build the application using your selected method
   - Configure all required environment variables
   - Deploy the application to Vercel

### Method 2: Manual Deployment

If you prefer to deploy manually:

1. Build the application using one of the provided scripts:
   ```bash
   # Recommended production build
   npm run build:prod
   
   # Or static export (fully client-side)
   npm run build:static
   
   # Or standard Next.js build (may have SSR errors)
   npm run build
   ```

2. Deploy to Vercel:
   ```bash
   vercel deploy --prod
   ```

3. Configure the following environment variables:
   - `NEXT_PUBLIC_API_URL=https://api.cryoprotect.app/v1`
   - `NEXTAUTH_URL=https://www.cryoprotect.app`
   - `NEXTAUTH_SECRET=<your-secure-random-secret>`

## Build Methods

The CryoProtect frontend offers several build methods to address different deployment scenarios and potential issues:

### 1. Production Build (`build-prod.sh`)

This is the recommended build method for production deployment. It:
- Forces dynamic rendering for authenticated pages to avoid SSR errors
- Skips prerendering for protected routes like profile and settings
- Uses client-side rendering for pages that access session data
- Includes optimizations for production performance
- Works reliably with Vercel deployment

### 2. Static Export (`build-static.sh`)

This method generates a fully static site that can be deployed to any static hosting:
- Generates all pages as static HTML, CSS, and JavaScript
- All data fetching happens client-side
- No server-side rendering occurs at all
- Can be hosted on any static file server
- May have slower initial page loads but no build errors

### 3. Standard Next.js Build

The default Next.js build process:
- Uses server-side rendering where possible
- May encounter errors with pages that use authentication
- Not recommended unless you've fixed all SSR issues
- Faster initial page loads when working correctly

## Troubleshooting Common Issues

### Build Failures

If your build fails, check for these common issues:

1. **SSR Errors with Authentication**: The most common error is trying to access `session` or other browser-specific objects during server-side rendering. Solutions:
   - Use the `build-prod.sh` script instead of the standard build
   - Add the `'use client'` directive to all pages that use authentication
   - Use optional chaining for all session access: `session?.user?.name`
   - Check for window/document access: `typeof window !== 'undefined'`

2. **"Cannot read properties of undefined"**: This usually appears during prerendering of protected pages:
   - Use the client version of pages (`page.client.tsx`) that defer all data loading to client-side
   - Configure pages for dynamic rendering only
   - Use the Static Export build method

3. **TypeScript Errors**: We've configured the build to ignore TypeScript errors. If you're still getting them:
   - Check that settings in `next.config.js` has `ignoreBuildErrors: true`
   - Use the force flag when building: `next build --force`
   - Create proper type declarations for external libraries

### Runtime Errors

1. **Authentication Issues**: If users can't log in, check:
   - The `NEXTAUTH_URL` matches your actual deployment URL
   - The `NEXTAUTH_SECRET` is properly set and consistent with previous deployments
   - The API endpoints are correctly configured and accessible

2. **API Connection Problems**: If the application can't connect to the API:
   - Verify the `NEXT_PUBLIC_API_URL` is correct
   - Check API CORS settings allow your deployment domain
   - Check network requests in browser console for specific errors

## Deployment Domains

When deploying to Vercel, you'll get a default domain like `cryoprotect-abc123.vercel.app`. 

For production, we use custom domains:
- Production: `www.cryoprotect.app`
- API: `api.cryoprotect.app`

## Post-Deployment Verification

After deploying, verify these key features:

1. **Authentication**: Test sign-in and profile access
2. **API Integration**: Ensure molecule and mixture data loads correctly
3. **Visualization**: Check that 3D molecule viewer and charts render properly
4. **Theme Switching**: Verify dark/light mode toggle works

## Continuous Deployment

With our GitHub integration, Vercel will automatically deploy:
- Production branch (`master`) changes to the production environment
- Other branches to preview environments

To promote a preview to production, merge your branch to `master` or use Vercel's "Promote to Production" feature.

## Security Considerations

- Always use HTTPS (enabled by default on Vercel)
- Keep the `NEXTAUTH_SECRET` secure and consistent
- Use the Content Security Policy headers defined in `vercel.json`
- Regularly rotate API keys (if applicable)

## Support

If you encounter issues with your deployment, please contact the development team or file an issue in the GitHub repository.