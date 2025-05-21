# Netlify Deployment Quickstart

This guide provides direct steps to deploy the CryoProtect application to Netlify without any complications.

## Setup Instructions

### 1. Clone the Repository

```bash
git clone https://github.com/mushu-dev/cryoprotect.git
cd cryoprotect
```

### 2. Install the Netlify CLI

```bash
npm install -g netlify-cli
netlify login
```

### 3. Initialize Your Netlify Site

```bash
netlify init
```

Follow the prompts to create a new site.

### 4. Deploy Directly Using the Config File

The repository includes a `.netlify.toml` file that overrides any Netlify UI settings. This file:

- Configures the build command to install dependencies with `--legacy-peer-deps`
- Sets up the correct Node.js version and environment variables
- Configures Next.js runtime settings
- Sets up API redirects to your Heroku backend

Deploy using:

```bash
netlify deploy --prod
```

This will:
1. Use the build command: `cd frontend && npm install --legacy-peer-deps && npm run build`
2. Publish from the `frontend/.next` directory
3. Apply all redirects and headers from the configuration file

### 5. Check for Deployment Success

After deployment, Netlify will provide a URL where your application is hosted. Verify that:

1. The frontend loads correctly
2. API requests to the Heroku backend work properly
3. Authentication functions correctly

## Troubleshooting

If you encounter issues:

1. **Build failures**:
   - Check the build logs for specific errors
   - Make sure you're using Node.js v20
   - Ensure that `.netlify.toml` file is in the root directory

2. **API connectivity issues**:
   - Verify that the redirects are properly set up
   - Check that CORS is configured correctly on your Heroku backend
   - Test API endpoints directly

3. **Force a clean rebuild**:
   ```bash
   netlify deploy --prod --clear-cache
   ```

## Next Steps

- Set up continuous deployment with GitHub
- Configure custom domain settings
- Set up environment-specific variables for staging/production