# Netlify Deployment Guide

This guide explains how to deploy this minimal frontend to Netlify.

## Automatic Deployment via GitHub

This app is set up for automatic deployment when changes are pushed to the `netlify-autodeploy` branch.

### How Auto-deployment Works

1. GitHub Actions workflow is configured in `.github/workflows/minimal-frontend-deploy.yml`
2. When code is pushed to the `netlify-autodeploy` branch, the workflow:
   - Builds the Next.js app
   - Deploys it to Netlify
   - Sets up proper redirects for dynamic routes

### Requirements for Auto-deployment

- GitHub repository must have these secrets:
  - `NETLIFY_AUTH_TOKEN`: Your Netlify personal access token
  - `NETLIFY_SITE_ID_MINIMAL`: The site ID for your Netlify site

## Manual Deployment

You can also deploy manually to Netlify:

1. Build the site:
   ```
   npm run build
   ```

2. Deploy using Netlify CLI:
   ```
   npm run deploy:netlify
   ```
   
   Or directly:
   ```
   npx netlify deploy --prod --dir=out
   ```

## Configuration Files

### netlify.toml

The `netlify.toml` file contains the Netlify deployment configuration:
- Build settings
- Environment variables
- Redirects for client-side routing
- Security headers

### next.config.js

The Next.js configuration file is set up for static export with:
- `output: 'export'` - enables static HTML export
- `trailingSlash: false` - prevents trailing slashes in URLs
- `images.unoptimized: true` - required for static exports with images

## Handling Dynamic Routes

Since this is a static export, dynamic routes like `/molecules/[id]` need special handling.

For this we:
1. Generate pages during the build process
2. Use redirects in Netlify to handle client-side routing

The `_redirects` file is created during deployment with rules to direct dynamic routes properly.

## Troubleshooting

### Build Failures

If your build fails, check:
- Build logs in GitHub Actions
- Netlify deploy logs
- Run `npm run build` locally to see if there are any errors

### Routing Issues

If routes aren't working properly:
- Check the `_redirects` file in your published site
- Verify the Next.js configuration (`next.config.js`)
- Test navigation in development mode (`npm run dev`) first