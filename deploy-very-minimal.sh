#!/bin/bash
# Script to deploy a very minimal version of the site to Netlify

echo "Building and deploying very minimal version of the site..."

# Change to the minimal-frontend directory
cd minimal-frontend || { echo "minimal-frontend directory not found"; exit 1; }

# Install dependencies
echo "Installing dependencies..."
npm install

# Run the build 
echo "Running build..."
npm run build

# Check if build was successful
if [ ! -d "out" ]; then
    echo "Build failed - 'out' directory not found"
    echo "Running export manually..."
    npm run export
    if [ ! -d "out" ]; then
        echo "Export failed - 'out' directory still not found"
        exit 1
    fi
fi

# Deploy to Netlify
echo "Deploying to Netlify..."
netlify deploy --prod

echo "Deployment completed!"