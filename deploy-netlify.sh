#!/bin/bash
# Script to deploy to Netlify from the root directory

# Make script executable
chmod +x frontend/netlify-build.sh

# Check if netlify CLI is installed
if ! command -v netlify &> /dev/null; then
    echo "Netlify CLI not found. Installing..."
    npm install -g netlify-cli
fi

# Initialize site if needed
netlify status || netlify init

# Deploy to Netlify
netlify deploy --prod