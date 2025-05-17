#!/bin/bash
# Script to migrate environment variables from Vercel to Netlify

# Check if netlify CLI is installed
if ! command -v netlify &> /dev/null; then
    echo "Netlify CLI not found. Please install it with: npm install -g netlify-cli"
    exit 1
fi

# Initialize Netlify site if not already linked
netlify status || netlify init

# Get the Netlify site ID
NETLIFY_SITE_ID=$(netlify status | grep "Site Id" | awk '{print $3}')

if [ -z "$NETLIFY_SITE_ID" ]; then
    echo "Failed to get Netlify site ID. Please run netlify init first."
    exit 1
fi

echo "Netlify site ID: $NETLIFY_SITE_ID"

# Heroku backend URL
HEROKU_BACKEND_URL="https://cryoprotect-8030e4025428.herokuapp.com"

# Set up environment variables
echo "Setting up environment variables in Netlify..."

# Generate a random NEXTAUTH_SECRET if not provided
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

# Set API URL
netlify env:set NEXT_PUBLIC_API_URL "${HEROKU_BACKEND_URL}/v1"

# Common environment variables
netlify env:set NEXT_PUBLIC_USE_MOCK_DATA "false"
netlify env:set NEXT_PUBLIC_ENABLE_API_LOGGING "true"
netlify env:set NEXT_PUBLIC_ENVIRONMENT "production"

# Authentication related variables
netlify env:set NEXTAUTH_URL "https://${NETLIFY_SITE_ID}.netlify.app"
netlify env:set NEXTAUTH_SECRET "$NEXTAUTH_SECRET"

# Protection bypass
PROTECTION_BYPASS=${PROTECTION_BYPASS:-TAt23KbtFE8dkZobJU3hpgTP4L5ja07V}
netlify env:set PROTECTION_BYPASS "$PROTECTION_BYPASS"
netlify env:set NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS "$PROTECTION_BYPASS"

echo "Environment variables have been set up in Netlify."
echo "Next steps:"
echo "1. Run 'npm run build' to build the project"
echo "2. Run 'netlify deploy --prod' to deploy to Netlify"
echo "3. Run 'npm run verify-netlify' to verify the deployment"