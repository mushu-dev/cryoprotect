#!/bin/bash
# Script to update Vercel deployment settings for new repository

# Check if the Vercel CLI is installed
if ! command -v vercel &> /dev/null; then
    echo "Vercel CLI is not installed. Please install it using: npm install -g vercel"
    exit 1
fi

# Ensure user is logged in to Vercel
vercel whoami &> /dev/null
if [ $? -ne 0 ]; then
    echo "Please log in to Vercel first using: vercel login"
    exit 1
fi

# Set variables
NEW_OWNER="mushu-dev"
REPO_NAME="CryoProtect"
FRONTEND_DIR="frontend"
PROJECT_NAME="cryoprotect1"  # Use the project name from vercel project ls

echo "Updating Vercel project for new repository location..."

# Go to frontend directory
cd "$FRONTEND_DIR" || (echo "Frontend directory not found"; exit 1)

# Check if project is already linked
if [ -f ".vercel/project.json" ]; then
    echo "Project is already linked. Removing link to reconnect with new repository..."
    vercel unlink
    echo "Project unlinked."
fi

# Link the project
echo "Linking Vercel project to new GitHub repository..."
vercel link -p "$PROJECT_NAME"

# Connect to the new GitHub repository
echo "Configuring project to use new GitHub repository..."
vercel git connect "github" "$NEW_OWNER/$REPO_NAME"

# Pull environment variables
echo "Pulling environment variables..."
vercel env pull .env.local

# Deploy to ensure everything is working
echo "Deploying project to verify configuration..."
vercel --prod

echo "Vercel project has been updated to use the new GitHub repository: $NEW_OWNER/$REPO_NAME"
echo "Please check the Vercel dashboard to verify the configuration:"
echo "https://vercel.com/dashboard/projects/$PROJECT_NAME/settings"