#!/bin/bash
# Script to deploy to Netlify directly using their Next.js runtime plugin

# Change to the frontend directory
cd frontend || { echo "Frontend directory not found"; exit 1; }

# Check if netlify CLI is installed
if ! command -v netlify &> /dev/null; then
    echo "Netlify CLI not found. Installing..."
    npm install -g netlify-cli
fi

# Make sure Next.js plugin is installed
npm list @netlify/plugin-nextjs || npm install --save-dev @netlify/plugin-nextjs

# Initialize site if needed
netlify status || netlify init

# Deploy to Netlify without static export (using Next.js runtime)
echo "Deploying to Netlify with Next.js runtime..."
netlify deploy --prod

echo "Deployment completed!"