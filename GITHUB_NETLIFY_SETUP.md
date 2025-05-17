# GitHub and Netlify Integration Setup

This guide explains how to set up the necessary secrets and configuration for automatic deployment from GitHub to Netlify.

## Required GitHub Secrets

To enable the GitHub Actions workflow for Netlify deployment, you need to add the following secrets to your GitHub repository:

1. **NETLIFY_AUTH_TOKEN**: Your Netlify personal access token
2. **NETLIFY_SITE_ID**: The ID of your Netlify site
3. **NEXTAUTH_SECRET**: A random string used for NextAuth.js session encryption
4. **PROTECTION_BYPASS**: The protection bypass token for your application
5. **HEROKU_API_KEY** (optional): Your Heroku API key to update the backend with the Netlify URL

## Steps to Set Up GitHub Secrets

1. **Get your Netlify Auth Token**:
   - Log in to Netlify
   - Go to User Settings > Applications > Personal Access Tokens
   - Create a new access token and copy it

2. **Get your Netlify Site ID**:
   - From the Netlify dashboard, select your site
   - Go to Site Settings > General > Site details
   - Copy the Site ID

3. **Add secrets to GitHub**:
   - Go to your GitHub repository
   - Navigate to Settings > Secrets and variables > Actions
   - Click "New repository secret"
   - Add each of the required secrets with their values

## Manual Site Setup in Netlify

Before using the GitHub Actions workflow, you should manually set up your site in Netlify at least once:

1. **Create a new site in Netlify**:
   ```bash
   cd frontend
   netlify init
   ```

2. **Set up the environment variables**:
   ```bash
   npm run migrate-to-netlify
   ```

3. **Deploy once manually**:
   ```bash
   npm run deploy:netlify
   ```

## GitHub Actions Workflow

Our GitHub Actions workflow (`.github/workflows/netlify-deploy.yml`) will:

1. Automatically deploy all pushes to the main/master branch to Netlify production
2. Create preview deployments for all pull requests
3. Add a comment to pull requests with the preview URL
4. Update the Heroku backend with the Netlify URL (production deployments only)

## Testing the Setup

After configuring all secrets and pushing the workflow file to GitHub:

1. Make a small change to the frontend code
2. Commit and push to GitHub
3. Check the Actions tab in your GitHub repository to monitor the workflow
4. Verify the deployment in Netlify

## Troubleshooting

If you encounter issues with the GitHub Actions workflow:

1. **Check workflow logs** in the GitHub Actions tab for specific error messages
2. **Verify secret values** are correctly set up
3. **Try a manual deployment** to see if the issue is with Netlify or with the workflow
4. **Check Netlify build logs** in the Netlify dashboard

## Local Workflow Testing

You can test the GitHub Actions workflow locally using [act](https://github.com/nektos/act):

```bash
# Install act
brew install act  # macOS
# OR
curl https://raw.githubusercontent.com/nektos/act/master/install.sh | sudo bash  # Linux

# Run the workflow locally
act push -s NETLIFY_AUTH_TOKEN=your_token -s NETLIFY_SITE_ID=your_site_id
```