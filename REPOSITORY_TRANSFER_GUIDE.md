# Repository Transfer Guide: blueprint-house → mushu-dev

This guide outlines the complete process for transferring the CryoProtect repository from the blueprint-house organization to the mushu-dev personal account, along with all necessary updates to maintain functionality with deployment platforms.

## Table of Contents

1. [Repository Transfer Process](#repository-transfer-process)
2. [Local Repository Updates](#local-repository-updates)
3. [Deployment Platform Updates](#deployment-platform-updates)
4. [GitHub Actions Workflows](#github-actions-workflows)
5. [Testing and Verification](#testing-and-verification)
6. [Troubleshooting](#troubleshooting)

## Repository Transfer Process

### 1. Transfer Repository on GitHub

#### Prerequisites
- You must have admin permissions on both the source (blueprint-house organization) and destination (mushu-dev account)
- The mushu-dev account cannot have an existing repository with the same name
- You should back up any repository data before proceeding

#### Steps

1. Log in to GitHub as the account with admin access to blueprint-house
2. Navigate to the CryoProtect repository: https://github.com/blueprint-house/CryoProtect
3. Go to "Settings" > "Options" (scroll down to the bottom)
4. In the "Danger Zone" section, click "Transfer"
5. Enter "blueprint-house/CryoProtect" as the repository name for confirmation
6. Enter "mushu-dev" as the new owner
7. Click "I understand, transfer this repository"
8. If prompted, enter your GitHub password to confirm

### 2. Repository Visibility Settings

Ensure the transferred repository has the correct visibility settings:

1. Go to the newly transferred repository: https://github.com/mushu-dev/CryoProtect
2. Go to "Settings" > "Options" 
3. Scroll down to "Danger Zone"
4. If you want the repository to be public, click "Change visibility" and select "Public"
5. If you want to keep it private, ensure it's set to "Private"

## Local Repository Updates

After the repository transfer is complete, update your local repository:

```bash
# Navigate to your local repository
cd /home/mushu/Projects/CryoProtect

# View current remote
git remote -v  # Should show blueprint-house URLs

# Update remote URL
git remote set-url origin https://github.com/mushu-dev/CryoProtect.git

# Verify the update
git remote -v  # Should now show mushu-dev URLs

# Fetch from new remote
git fetch origin
```

## Deployment Platform Updates

### Vercel

1. Log in to Vercel dashboard: https://vercel.com/dashboard
2. Go to the CryoProtect project
3. Go to "Settings" > "Git"
4. Disconnect the current GitHub repository (blueprint-house)
5. Connect to the new repository (mushu-dev/CryoProtect)
6. Update environment variables if necessary
7. Deploy the project to verify the connection

**Using Vercel CLI:**

```bash
# Go to frontend directory
cd /home/mushu/Projects/CryoProtect/frontend

# Remove current project link
rm -rf .vercel

# Link to your new project
vercel link

# Deploy to production
vercel --prod
```

### Heroku

1. Log in to Heroku dashboard: https://dashboard.heroku.com/apps
2. Go to the CryoProtect app
3. Go to "Deploy" tab
4. Under "Deployment method", select "GitHub"
5. If needed, disconnect the current GitHub repository
6. Connect to the new repository (mushu-dev/CryoProtect)
7. Enable automatic deploys for the desired branch
8. Deploy the app to verify the connection

**Using Heroku CLI:**

```bash
# View current remotes
git remote -v

# Remove the current Heroku remote
git remote remove heroku

# Add the new Heroku remote
heroku git:remote -a cryoprotect

# Verify the change
git remote -v
```

### Fly.io

For Fly.io, update the GitHub Actions workflow to reference the new repository:

1. Edit `.github/workflows/fly-deploy.yml`
2. Update any references to blueprint-house/CryoProtect to mushu-dev/CryoProtect
3. Ensure the FLY_API_TOKEN is set in the GitHub secrets for the new repository

## GitHub Actions Workflows

1. Update all GitHub Actions workflow files in `.github/workflows/`:
   - `vercel-deploy-enhanced.yml`
   - `heroku-deploy.yml`
   - `fly-deploy.yml`
   - Any other workflow files

2. For each file, check:
   - Repository references
   - Commented URLs
   - GitHub organization/owner references

3. Transfer GitHub Secrets:
   - Go to Settings > Secrets and variables > Actions in both repositories
   - Copy all secrets from the blueprint-house repository to the mushu-dev repository

## Update URLs in Project Files

Search for all references to blueprint-house in the codebase and update them:

```bash
# Find all files containing blueprint-house
grep -r "blueprint-house" /home/mushu/Projects/CryoProtect --include="*.md" --include="*.js" --include="*.ts" --include="*.tsx" --include="*.py" --include="*.yml" --include="*.json"

# For each file found, edit it to replace blueprint-house with mushu-dev
```

Common places to check:
- README.md
- Documentation
- Package.json
- Configuration files
- Deployment scripts
- GitHub Actions workflows

## Testing and Verification

After completing all updates, verify that everything works:

1. **Local Development:**
   - Clone the repository from the new location
   - Install dependencies
   - Run the development server
   - Verify functionality

2. **CI/CD:**
   - Push a small change to trigger CI/CD workflows
   - Verify that GitHub Actions run successfully
   - Check that deployments complete successfully

3. **Deployments:**
   - Verify Vercel deployment
   - Verify Heroku deployment
   - Verify Fly.io deployment
   - Check that cross-service communication works

## Troubleshooting

### GitHub Transfer Issues

- **Problem:** Transfer option is not available
- **Solution:** Ensure you have admin access to both the source and target accounts

### Deployment Connection Issues

- **Problem:** Vercel cannot connect to the new repository
- **Solution:** Check repository visibility settings and that Vercel has access to it

- **Problem:** Heroku deployment fails
- **Solution:** Verify GitHub integration credentials and repository link

- **Problem:** GitHub Actions fail after transfer
- **Solution:** Ensure all secrets are transferred to the new repository

### Local Repository Issues

- **Problem:** Cannot push to new repository
- **Solution:** 
  ```bash
  git remote set-url origin https://github.com/mushu-dev/CryoProtect.git
  git push -u origin master
  ```

- **Problem:** Pull requests to old repository
- **Solution:** Update default remote and branch
  ```bash
  git config --local branch.master.remote origin
  git config --local branch.master.merge refs/heads/master
  ```