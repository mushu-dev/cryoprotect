# Experimental Data Enhancement Deployment Guide

This document outlines the process for deploying the Experimental Data Enhancement feature to Netlify.

## Prerequisites

Before deploying, ensure you have:

1. Netlify CLI installed and authenticated
2. Node.js and npm installed
3. Playwright installed for verification tests
4. Git access to the repository
5. Access to the Netlify site settings

## Deployment Process

The deployment process is automated through the `deploy-experimental-data-enhancement.sh` script, which handles:

1. Running validation tests
2. Setting environment variables
3. Building the application
4. Deploying to Netlify
5. Running verification tests
6. Generating a deployment report

### Step 1: Validate the Feature Locally

Before deployment, run the validation tests to ensure everything's working as expected:

```bash
./validate-experimental-data-enhancement.sh
```

### Step 2: Deploy using the Script

Run the deployment script:

```bash
./deploy-experimental-data-enhancement.sh
```

You can customize the deployment by setting environment variables:

```bash
NETLIFY_SITE_NAME=your-site-name \
API_URL=https://your-api-url.com/v1 \
PROTECTION_BYPASS=your-bypass-token \
./deploy-experimental-data-enhancement.sh
```

### Step 3: Verify the Deployment

After deployment, you can run a quick verification:

```bash
./verify-experimental-data-deployment.sh
```

Or for more comprehensive verification:

```bash
NETLIFY_URL=https://your-site.netlify.app \
./validate-experimental-data-enhancement.sh
```

## Deployment Configuration

### Environment Variables

The following environment variables are set during deployment:

| Variable | Description | Default |
|----------|-------------|---------|
| `NEXT_PUBLIC_API_URL` | Backend API URL | https://cryoprotect-8030e4025428.herokuapp.com/v1 |
| `NEXT_PUBLIC_ENABLE_EXPERIMENTAL_FEATURES` | Enable experimental features | true |
| `NEXT_PUBLIC_USE_MOCK_DATA` | Use mock data when API is unavailable | false |
| `NEXT_PUBLIC_ENABLE_API_LOGGING` | Enable API request logging | true |
| `NEXT_PUBLIC_ENVIRONMENT` | Environment name | production |
| `NEXT_PUBLIC_NETLIFY` | Indicates Netlify deployment | true |

### Deployment Report

After each deployment, a report is generated in the `test-results` directory with information about:

- Deployment timestamp
- Deployed URL
- Feature status
- Verification results

The report is available at `./test-results/deployment-report-{timestamp}.md`.

## Troubleshooting

### Common Issues

1. **Deployment Fails with 401 Error**
   - Ensure you're logged in to Netlify CLI
   - Run `netlify login` to authenticate

2. **Validation Tests Fail**
   - Check browser dependencies: `npx playwright install --with-deps`
   - Try running specific tests to identify the issue

3. **Environment Variables Not Applied**
   - Verify variables in Netlify dashboard
   - Run `netlify env:list` to confirm settings

4. **Feature Not Visible After Deployment**
   - Verify `NEXT_PUBLIC_ENABLE_EXPERIMENTAL_FEATURES` is set to `true`
   - Check browser console for errors
   - Clear browser cache

### Getting Help

If you encounter issues with the deployment:

1. Check the Netlify logs: `netlify sites:list`
2. Review the full test output: `npx playwright show-report`
3. Run individual tests to pinpoint issues: 
   ```bash
   npx playwright test tests/playwright/experimental-data-ui.spec.js --debug
   ```

## Current Implementation Status

| Feature | Status |
|---------|--------|
| Basic experiment listing | ✅ Complete |
| Experiment detail view | ✅ Complete |
| Enhanced visualization | ⚠️ In progress |
| Filtering | ⚠️ In progress |
| Data sorting | ⚠️ In progress |
| Data comparison | ⚠️ Planned |
| Data export | ⚠️ Planned |
| Mobile responsiveness | ✅ Complete |