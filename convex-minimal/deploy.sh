#!/bin/bash
# Script to deploy minimal Convex API functions

echo "Deploying minimal Convex API functions..."

# Set Convex deployment
export CONVEX_DEPLOYMENT=upbeat-parrot-866

# Initialize Convex project if needed
if [ ! -d "convex/_generated" ]; then
  echo "Initializing Convex project..."
  npx convex init
fi

# Deploy
npm run deploy

echo "Deployment complete!"