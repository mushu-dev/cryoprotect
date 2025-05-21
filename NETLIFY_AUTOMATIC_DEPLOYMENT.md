# Netlify Automatic Deployment Guide

This document explains how the automatic deployment to Netlify is configured and how to troubleshoot common issues.

## Deployment Strategy

The project uses GitHub Actions for automatic deployment to Netlify whenever changes are pushed to the main branch or when pull requests are created.

### Key Components

1. **GitHub Actions Workflow**: `.github/workflows/netlify-deploy.yml`
   - Automatically triggered on pushes to main/master or pull requests
   - Builds the frontend and deploys to Netlify
   - Updates Heroku with the Netlify URL

2. **Pages Router Configuration**: `frontend/use-pages-router.sh`
   - Temporarily hides the App Router during build
   - Ensures compatibility with Netlify static export

3. **Next.js Configuration**: `frontend/next.config.js`
   - Configured for static export using `output: 'export'`
   - Set up to properly handle dynamic routes

## How It Works

1. When you push changes to the main/master branch:
   - GitHub Actions workflow is triggered
   - The frontend is built with the Pages Router approach
   - Static files are generated in the `out` directory
   - These files are deployed to Netlify production

2. When you create a pull request:
   - GitHub Actions creates a deploy preview
   - A comment is added to the PR with the preview URL
   - You can test changes before merging

## Vulnerability Fixes

The workflow includes steps to automatically fix vulnerabilities:

1. Runs `npm audit fix --force` during the build process
2. Updates key dependencies to secure versions

## Troubleshooting

### Issue: Build fails with App Router error

**Solution**: The workflow temporarily disables the App Router during build. Make sure the `use-pages-router.sh` script exists and is working correctly.

### Issue: UI components not rendering

**Solution**: Make sure your UI components are properly imported in the Pages Router files, especially in `/pages/_app.js`.

### Issue: API endpoints not working

**Solution**: Check the `_redirects` file in the deployed site. It should contain proper redirects to the backend services.

## Local Testing

To test the Netlify deployment locally:

1. Run the `deploy-to-netlify-fixed.sh` script:
   ```bash
   cd frontend
   ./deploy-to-netlify-fixed.sh
   ```

2. This will:
   - Configure the environment for Pages Router
   - Build the static export
   - Prepare all necessary Netlify configuration files

## Manual Deployment

If you need to deploy manually:

1. Install the Netlify CLI:
   ```bash
   npm install -g netlify-cli
   ```

2. Authenticate:
   ```bash
   netlify login
   ```

3. Deploy from the frontend directory:
   ```bash
   cd frontend
   netlify deploy --prod --dir=out
   ```

## Environment Variables

Make sure these GitHub secrets are configured:

- `NETLIFY_AUTH_TOKEN`: Your Netlify authentication token
- `NETLIFY_SITE_ID`: Your Netlify site ID
- `NEXTAUTH_SECRET`: Secret for Next.js authentication
- `PROTECTION_BYPASS`: Token to bypass frontend protection
- `HEROKU_API_KEY`: Only needed if updating Heroku configs