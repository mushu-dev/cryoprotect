#!/bin/bash
# Script to deploy the updated site with fixed build settings

set -e

echo "Deploying site with fixed Netlify configuration..."

# Navigate to the project root
cd "$(dirname "$0")"

# Deploy to Netlify using the root directory to pick up netlify.toml
echo "Deploying to Netlify..."
npx netlify deploy --prod

echo "Deployment complete!"
echo "The site should now be correctly deployed at https://cryoprotect.app"
echo "Please verify the homepage and experimental data enhancement pages are working."