#!/bin/bash
# Deploy updated Netlify configuration
# This script updates the Netlify site with the new configuration

set -e

# Configuration variables with defaults
NETLIFY_SITE_ID=${NETLIFY_SITE_ID:-b3069bed-40fd-4f61-82bc-e5f2fc010ccd}

echo "==== CryoProtect Netlify Configuration Deployment ===="
echo "Target Netlify site ID: $NETLIFY_SITE_ID"
echo

# Verify Netlify CLI is available
if ! command -v netlify &> /dev/null; then
    echo "❌ Netlify CLI is not installed or not available in PATH"
    echo "Please install the Netlify CLI: npm install -g netlify-cli"
    exit 1
fi

# Verify logged in to Netlify
echo "Verifying Netlify authentication..."
netlify status || {
    echo "❌ Not logged in to Netlify"
    echo "Please run 'netlify login' to authenticate"
    exit 1
}

# Navigate to the frontend directory
cd frontend

# Make sure we have the right configuration files
if [ ! -f "netlify.toml" ]; then
    echo "❌ netlify.toml not found in frontend directory"
    exit 1
fi

# Update netlify.toml with Convex configuration
echo "Updating netlify.toml with Convex configuration..."

# Update the netlify.toml file to use the correct Heroku app name
# and ensure Convex is properly configured
sed -i 's/cryoprotect-8030e4025428\.herokuapp\.com/cryoprotect.herokuapp.com/g' netlify.toml

# Ensure the correct Convex URL is set
sed -i 's|NEXT_PUBLIC_CONVEX_URL = ".*"|NEXT_PUBLIC_CONVEX_URL = "https://upbeat-parrot-866.convex.cloud"|g' netlify.toml

# Ensure Convex environment variables are set
grep -q "NEXT_PUBLIC_CONVEX_URL" netlify.toml || 
  sed -i '/NEXT_PUBLIC_ENVIRONMENT/a \  NEXT_PUBLIC_CONVEX_URL = "https://dynamic-mink-63.convex.cloud"' netlify.toml

grep -q "NEXT_PUBLIC_USE_CONVEX" netlify.toml || 
  sed -i '/NEXT_PUBLIC_ENVIRONMENT/a \  NEXT_PUBLIC_USE_CONVEX = "true"' netlify.toml

# Ensure build command includes Convex
sed -i 's/npm run build/npm run build:with-convex/g' netlify.toml

# Add a direct Convex redirect if needed
grep -q "from = \"/convex/" netlify.toml || 
  sed -i '/from = "\/api\/\*"/a \
[[redirects]]\
  from = "/convex/*"\
  to = "https://dynamic-mink-63.convex.cloud/:splat"\
  status = 200\
  force = true\
  headers = {Access-Control-Allow-Origin = "*"}' netlify.toml

echo "Deploying Netlify configuration..."
netlify deploy --site="$NETLIFY_SITE_ID" --prod

echo "✅ Netlify configuration deployed successfully"
echo "Run the test-connection.js script to verify the configuration."