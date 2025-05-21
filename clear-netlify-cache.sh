#!/bin/bash
# Script to clear Netlify cache and redeploy

# Check if Netlify CLI is installed
if ! command -v netlify &> /dev/null; then
    echo "Netlify CLI not found. Installing..."
    npm install -g netlify-cli
fi

# Clear Netlify cache
echo "Clearing Netlify cache..."
netlify build --clear-cache

# Deploy to Netlify with production flag
echo "Deploying to Netlify with clean cache..."
netlify deploy --prod --clear-cache --message "Deployment with cleared cache"

echo "Done! Check the Netlify dashboard for the deployment status."