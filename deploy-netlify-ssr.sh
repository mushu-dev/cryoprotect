#!/bin/bash
# Script to deploy to Netlify using server-side rendering with Next.js runtime

echo "Deploying to Netlify with server-side rendering..."

# Change to the frontend directory
cd frontend || { echo "Frontend directory not found"; exit 1; }

# Clean any previous build artifacts
echo "Cleaning previous builds..."
rm -rf .next out

# Run the build 
echo "Running build..."
npm run build

# Check if build was successful
if [ ! -d ".next" ]; then
    echo "Build failed - '.next' directory not found"
    exit 1
fi

# Deploy to Netlify with the NextJS runtime
echo "Deploying to Netlify..."
netlify deploy --prod

echo "Deployment completed!"