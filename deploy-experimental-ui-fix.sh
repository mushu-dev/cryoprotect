#!/bin/bash
# Script to deploy the experimental data enhancement UI fix

set -e

echo "Deploying experimental data enhancement UI fix..."

# Navigate to the project root
cd "$(dirname "$0")"

# Navigate to the frontend directory
cd frontend

# Install dependencies if needed
if [ ! -d "node_modules" ]; then
  echo "Installing dependencies..."
  npm install
fi

# Build the project
echo "Building the project..."
npm run build

# Deploy to Netlify (assuming Netlify CLI is installed)
echo "Deploying to Netlify..."
npx netlify deploy --prod

echo "Deployment complete!"
echo "The experimental data enhancement UI should now be fully functional."
echo "Visit https://cryoprotect.app to see the changes."