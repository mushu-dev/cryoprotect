# Netlify AutoDeploy Setup

This branch (`netlify-autodeploy`) is configured to automatically deploy to Netlify when changes are pushed.

## How It Works

1. The GitHub Action workflow (`.github/workflows/minimal-frontend-deploy.yml`) is triggered on pushes to this branch that modify files in the `minimal-frontend` directory.

2. The workflow:
   - Checks out the code
   - Sets up Node.js
   - Installs dependencies
   - Builds the Next.js app
   - Deploys to Netlify using the Netlify CLI

## Required Secrets

For this workflow to function properly, the following secrets must be set in your GitHub repository:

- `NETLIFY_AUTH_TOKEN`: Your Netlify personal access token
- `NETLIFY_SITE_ID_MINIMAL`: The site ID for your minimal-frontend Netlify site

## How to Use

1. Work on your changes in the main development branch
2. When you want to deploy to Netlify, merge your changes into this branch:

```bash
# From your development branch
git push origin mybranch

# Create a PR from your branch to netlify-autodeploy
# Or directly merge/cherry-pick changes to this branch

# Push this branch to trigger the deploy
git push origin netlify-autodeploy
```

3. The GitHub Action will run automatically and deploy your changes to Netlify

## Checking Deployment Status

You can check the status of your deployments in:
- GitHub Actions tab in your repository
- Netlify dashboard for your site

## Troubleshooting

If your deployment fails:
1. Check the GitHub Actions logs for errors
2. Verify your Netlify API token hasn't expired
3. Ensure your NETLIFY_SITE_ID_MINIMAL is correct
4. Check if the build process is failing (run it locally to debug)

## Notes on Static Exports

- The minimal-frontend Next.js app is configured for static export with `output: 'export'`
- The build generates static HTML/CSS/JS files in the `out` directory
- For route handling on Netlify, redirects are set up in `out/_redirects`