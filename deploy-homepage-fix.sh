#!/bin/bash
# Script to deploy the homepage fix

set -e

echo "Deploying homepage fix..."

# Navigate to the project root
cd "$(dirname "$0")"

# Navigate to the frontend directory
cd frontend

# Build the project
echo "Building the project..."
npm run build

# Deploy to Netlify
echo "Deploying to Netlify..."
npx netlify deploy --prod

echo "Deployment complete!"
echo "The homepage fix should now be visible at https://cryoprotect.app"