#!/bin/bash
# Script to deploy the homepage fix

set -e

echo "Deploying homepage fix..."

# Navigate to the project root
cd "$(dirname "$0")"

# Navigate to the frontend directory
cd frontend

# Install dependencies if needed
if [ ! -d "node_modules" ]; then
  echo "Installing dependencies..."
  npm install
fi

# Build the project with static export
echo "Building the project with static export..."
npm run build
npm run export

# Deploy to Netlify using the out directory
echo "Deploying to Netlify..."
npx netlify deploy --prod --dir=out

echo "Deployment complete!"
echo "The homepage fix should now be visible at https://cryoprotect.app"