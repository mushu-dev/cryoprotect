# Netlify Deployment Guide for CryoProtect Frontend

This guide walks through the process of deploying the CryoProtect frontend to Netlify while maintaining connectivity with the Heroku backend.

## Prerequisites

- [Node.js](https://nodejs.org/) 18 or later
- [Netlify CLI](https://docs.netlify.com/cli/get-started/) (`npm install -g netlify-cli`)
- [Heroku CLI](https://devcenter.heroku.com/articles/heroku-cli) (optional, for backend updates)

## Deployment Steps

### Step 1: Initial Setup

1. Navigate to the frontend directory:
   ```bash
   cd frontend
   ```

2. Install dependencies:
   ```bash
   npm run install-deps
   ```

3. Login to Netlify (if not already):
   ```bash
   netlify login
   ```

### Step 2: Configure Netlify Site

1. Initialize Netlify site (if not already linked):
   ```bash
   netlify init
   ```
   
   This will guide you through creating a new site or linking to an existing one.

2. Migrate environment variables from Vercel to Netlify:
   ```bash
   npm run migrate-to-netlify
   ```

   This script automatically sets up all necessary environment variables in Netlify.

### Step 3: Update API Endpoints

1. Update any hardcoded API endpoints in the code:
   ```bash
   npm run update-api-endpoints
   ```

   This ensures that all API calls use environment variables instead of hardcoded URLs.

### Step 4: Deploy to Netlify

1. Deploy using our custom script:
   ```bash
   npm run deploy:netlify
   ```

   This script will:
   - Update API endpoints
   - Set environment variables
   - Build the project
   - Deploy to Netlify
   - Verify the deployment
   - Update Heroku with the new frontend URL (if HEROKU_API_KEY is set)

### Step 5: Verify Deployment

1. After deployment, verify the connectivity:
   ```bash
   npm run verify-netlify
   ```

   This script checks both the Netlify deployment and API connectivity.

2. Alternatively, use Playwright to test the site visually:
   ```bash
   # From project root
   ./mcp-playwright-final.sh browser_navigate https://your-site-name.netlify.app
   ./mcp-playwright-final.sh browser_take_screenshot https://your-site-name.netlify.app screenshot.png
   ```

## Continuous Deployment

To set up continuous deployment with GitHub:

1. Connect your GitHub repository in the Netlify UI
2. Configure build settings:
   - Build command: `npm run build`
   - Publish directory: `.next`

3. Set up environment variables in the Netlify UI or through the CLI

## Troubleshooting

### API Connectivity Issues

If you encounter issues with API connectivity:

1. Verify CORS settings on the Heroku backend are updated for the Netlify domain
2. Check the Netlify redirects in `netlify.toml` are properly configured
3. Ensure environment variables are correctly set:
   ```bash
   netlify env:list
   ```

### Build Errors

For build errors:

1. Check the build logs in the Netlify UI
2. Try a local build to debug:
   ```bash
   npm run build
   ```

3. Verify Next.js configuration in `next.config.js`

## Switching Back to Vercel

If needed, you can easily switch back to Vercel:

```bash
npm run deploy
```

## Additional Resources

- [Netlify documentation for Next.js](https://docs.netlify.com/configure-builds/common-configurations/next-js/)
- [Netlify CLI documentation](https://docs.netlify.com/cli/get-started/)
- [Next.js on Netlify](https://www.netlify.com/with/nextjs/)