#!/bin/bash
# Script to update Heroku deployment settings for new repository

# Check if the Heroku CLI is installed
if ! command -v heroku &> /dev/null; then
    echo "Heroku CLI is not installed. Please install it."
    exit 1
fi

# Ensure user is logged in to Heroku
heroku whoami &> /dev/null
if [ $? -ne 0 ]; then
    echo "Please log in to Heroku first using: heroku login"
    exit 1
fi

# Set variables
NEW_OWNER="mushu-dev"
REPO_NAME="CryoProtect"
APP_NAME="cryoprotect"  # Use your main Heroku app name

echo "Updating Heroku deployment for new repository location..."

# Check if we're in the git repository
if [ ! -d ".git" ]; then
    echo "Not in a git repository. Please run this script from the root of your project."
    exit 1
fi

# View current git remotes
echo "Current git remotes:"
git remote -v

# Remove the existing Heroku remote if it exists
if git remote | grep -q "heroku"; then
    echo "Removing existing Heroku remote..."
    git remote remove heroku
fi

# Add the new Heroku remote
echo "Adding Heroku remote for app: $APP_NAME"
heroku git:remote -a "$APP_NAME"

# Verify the remote was added
echo "Updated git remotes:"
git remote -v

# Update Heroku GitHub integration via API
echo "Updating Heroku GitHub integration..."
echo "This step needs to be done manually in the Heroku dashboard:"
echo "1. Go to https://dashboard.heroku.com/apps/$APP_NAME/deploy/github"
echo "2. If necessary, disconnect the current GitHub repository"
echo "3. Connect to the new repository: $NEW_OWNER/$REPO_NAME"
echo "4. Set up automatic deploys for your preferred branch"

# Instructions for the GitHub Actions workflow
echo ""
echo "Don't forget to update your GitHub Actions workflow:"
echo "1. Ensure the Heroku API key is set in the GitHub secrets for the new repository"
echo "2. Update any references to the old GitHub organization in .github/workflows/heroku-deploy.yml"

# Manually force a deploy
read -p "Do you want to attempt a manual deployment now? (y/N) " DEPLOY
if [[ "$DEPLOY" =~ ^[Yy]$ ]]; then
    echo "Deploying to Heroku..."
    git push heroku master
    heroku open
else
    echo "Skipping manual deployment. You can deploy later with: git push heroku master"
fi

echo "Heroku repository update complete. Please verify that deployments work correctly."