# GitHub Organization Setup and Deployment Integration

This guide outlines how to properly set up and maintain the CryoProtect repository within a GitHub organization and ensure all deployment platforms (Vercel, Heroku, and Fly.io) are correctly integrated.

## Table of Contents

1. [GitHub Organization Repository](#github-organization-repository)
2. [Vercel Integration](#vercel-integration)
3. [Heroku Integration](#heroku-integration)
4. [Fly.io Integration](#flyio-integration)
5. [GitHub Actions Workflows](#github-actions-workflows)
6. [Troubleshooting](#troubleshooting)

## GitHub Organization Repository

### Repository Configuration

The CryoProtect repository is hosted under the `blueprint-house` GitHub organization. This provides several advantages:
- Centralized ownership and access control
- Team-based permissions
- Organization-level secrets and settings

### Branch Protection Rules

To ensure code quality, the `master` branch has the following protection rules:
- Require pull request reviews before merging
- Require status checks to pass (CI/CD workflows)
- Restrict who can push to the branch

### Clone and Push Commands

When working with the organization repository, use the following commands:

```bash
# Clone the repository
git clone https://github.com/blueprint-house/CryoProtect.git

# Add your changes, commit them
git add .
git commit -m "Your commit message"

# Push to the remote
git push origin <branch-name>
```

## Vercel Integration

### Connecting Vercel to GitHub Organization

1. Log in to your Vercel account and navigate to "Settings" > "Git"
2. Under "GitHub", click "Connect" or "Configure" if already connected
3. Grant access to the `blueprint-house` organization
4. Import the CryoProtect project

### Deployment Configuration

The Vercel deployment settings are:
- **Framework Preset**: Next.js
- **Root Directory**: `frontend`
- **Build Command**: `npm run build`
- **Output Directory**: `.next`
- **Install Command**: `npm ci`

### Environment Variables

Required environment variables for Vercel:
- `NEXT_PUBLIC_API_URL`: The URL of the backend API
- `NEXT_PUBLIC_RDKIT_SERVICE_URL`: The URL of the RDKit service
- `NEXTAUTH_URL`: Authentication URL for NextAuth
- `NEXTAUTH_SECRET`: Secret for NextAuth

### GitHub Actions for Vercel

The `.github/workflows/vercel-deploy-enhanced.yml` workflow is set up to:
1. Run code quality checks
2. Generate environment-specific configurations
3. Deploy to Vercel
4. Verify the deployment

## Heroku Integration

### Connecting Heroku to GitHub Organization

1. Log in to your Heroku account and go to "Dashboard"
2. Select the "cryoprotect" application
3. Go to "Deploy" tab > "Deployment method" > "GitHub"
4. Connect to GitHub and select the `blueprint-house` organization
5. Search for and connect to the "CryoProtect" repository

### Deployment Configuration

The Heroku deployment uses the following:
- **Buildpack**: Python
- **Procfile**: Defines the web process command
- **requirements.txt**: Lists all Python dependencies

### Environment Variables

Required environment variables for Heroku:
- `DATABASE_URL`: Supabase database connection URL
- `SECRET_KEY`: Secret key for Flask sessions and security
- `RDKIT_SERVICE_URL`: URL of the RDKit service on Fly.io
- `API_KEY`: API key for external service communication

### GitHub Actions for Heroku

The `.github/workflows/heroku-deploy.yml` workflow is set up to:
1. Run tests
2. Deploy to production Heroku app from master branch
3. Create review apps for pull requests

## Fly.io Integration

### Connecting Fly.io to GitHub Organization

Fly.io doesn't offer direct GitHub integration like Vercel and Heroku, so we use GitHub Actions to deploy:

1. Install the Fly.io CLI: `curl -L https://fly.io/install.sh | sh`
2. Log in to Fly.io: `fly auth login`
3. Add a Fly.io API token to your GitHub organization's secrets

### Deployment Configuration

The RDKit service deployment uses:
- **Configuration**: `fly.toml` in the `rdkit-service-files` directory
- **Docker**: Dockerfile in the `rdkit-service-files` directory

### Environment Variables

Required environment variables for Fly.io:
- `PORT`: Port for the application to listen on
- `ALLOWED_ORIGINS`: Comma-separated list of allowed CORS origins

### GitHub Actions for Fly.io

The `.github/workflows/fly-deploy.yml` workflow is set up to:
1. Run tests
2. Deploy to Fly.io from the master branch
3. Verify the deployment

## GitHub Actions Workflows

The repository includes several GitHub Actions workflows for CI/CD:

1. **Vercel Deployment** (`vercel-deploy-enhanced.yml`): Deploys the frontend to Vercel
2. **Heroku Deployment** (`heroku-deploy.yml`): Deploys the backend API to Heroku
3. **Fly.io Deployment** (`fly-deploy.yml`): Deploys the RDKit service to Fly.io
4. **CI/CD** (`ci-cd.yml`): Runs tests and other checks

## Troubleshooting

### Vercel Deployment Issues

- **Problem**: Vercel deployment fails to connect to GitHub organization
- **Solution**: Ensure Vercel has access to the GitHub organization in Vercel's Git Integration settings

### Heroku Deployment Issues

- **Problem**: Heroku deployment fails with "No app found"
- **Solution**: 
  1. Verify Heroku app name in settings
  2. Check that GitHub integration is properly set up
  3. Ensure Heroku API key is valid

### Fly.io Deployment Issues

- **Problem**: GitHub Actions deployment to Fly.io fails
- **Solution**:
  1. Check that `FLY_API_TOKEN` secret is set in GitHub
  2. Verify the fly.toml configuration
  3. Ensure Docker build is successful

### GitHub Organization Issues

- **Problem**: Cannot push to organization repository
- **Solution**:
  1. Verify you have the correct permissions in the organization
  2. Use a personal access token with the appropriate scopes
  3. Ensure your local git config is correct