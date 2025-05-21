#!/bin/bash
# Script to deploy to Netlify from the root directory - Fixed version

# Change to the frontend directory to run build
cd frontend || { echo "Frontend directory not found"; exit 1; }

# Run the build
echo "Running build..."
npm run build-export

# Check if build was successful
if [ ! -d "out" ]; then
    echo "Build failed - 'out' directory not found"
    exit 1
fi

# Check if netlify CLI is installed
if ! command -v netlify &> /dev/null; then
    echo "Netlify CLI not found. Installing..."
    npm install -g netlify-cli
fi

# Initialize site if needed
netlify status || netlify init

# Deploy to Netlify with the correct publish directory
netlify deploy --prod --dir=out

echo "Deployment completed!"