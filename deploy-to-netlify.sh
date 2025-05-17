#!/bin/bash
# Master script to deploy the frontend to Netlify from project root

set -e  # Exit on any error

echo "Deploying frontend to Netlify..."

# Navigate to frontend directory
cd frontend

# Ensure dependencies are installed
echo "Installing dependencies..."
npm run install-deps

# Update API endpoints
echo "Updating API endpoints..."
npm run update-api-endpoints

# Run the Netlify deployment script
echo "Running deployment script..."
npm run deploy:netlify

# Return to original directory
cd ..

echo "Deployment process completed!"
echo "Check NETLIFY_DEPLOYMENT_GUIDE.md for more information and troubleshooting tips."