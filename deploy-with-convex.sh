#!/bin/bash

# Exit on error
set -e

echo "=== CryoProtect Deployment with Convex Integration ==="
echo ""

# Check if Convex CLI is installed
if ! command -v convex &> /dev/null; then
  echo "Convex CLI not found. Installing..."
  npm install -g convex
fi

# Ensure we have the latest dependencies
echo "Installing dependencies..."
npm install

# Enter frontend directory to install its dependencies
echo "Installing frontend dependencies..."
cd frontend
npm run install-deps
cd ..

# Set up environment variables
echo "Setting up environment variables..."
npm run convex:setup

# Generate Convex TypeScript types
echo "Generating Convex TypeScript types..."
npm run convex:codegen

# Deploy Convex backend
echo "Deploying Convex backend..."
npm run convex:deploy

# Build and deploy the frontend
echo "Building and deploying frontend with Convex integration..."
cd frontend
NEXT_PUBLIC_USE_CONVEX=true npm run build:static
npm run deploy:netlify
cd ..

echo ""
echo "✅ Deployment complete!"
echo "Your application is now deployed with Convex integration."
echo ""
echo "To view your application, visit: https://cryoprotect.netlify.app"
echo "To test the Convex integration, visit: https://cryoprotect.netlify.app/convex-test"