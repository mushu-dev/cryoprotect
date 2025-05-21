#!/bin/bash
# Script to deploy Convex functions and setup the integration

# Make sure we have Convex CLI
if ! command -v convex &> /dev/null; then
    echo "Convex CLI is not installed. Installing..."
    npm install -g convex
fi

echo "Setting up Convex deployment..."

# Go to the frontend directory
cd frontend

# Generate Convex types
echo "Generating Convex types..."
npx convex codegen

# Deploy Convex functions
echo "Deploying Convex functions..."
npx convex deploy

# Ensure the Convex URL is set in .env.local
if [ ! -f .env.local ]; then
    echo "Creating .env.local with Convex URL..."
    echo "NEXT_PUBLIC_USE_CONVEX=true" > .env.local
    echo "NEXT_PUBLIC_CONVEX_URL=https://upbeat-parrot-866.convex.cloud" >> .env.local
fi

echo "Convex deployment complete! Your application is now configured to use Convex."
echo "Run 'npm run dev:with-convex' to test locally."
echo "Netlify has been updated to use Convex in production."