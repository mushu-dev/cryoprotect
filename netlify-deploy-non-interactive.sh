#!/bin/bash
# Deploy to Netlify non-interactively

set -e  # Exit on error

# Change to the frontend directory
cd frontend

# Install dependencies
echo "Installing dependencies..."
npm run install-deps

# Build the project
echo "Building the project..."
npm run build

# Create a new Netlify site (if not already created)
SITE_NAME="cryoprotect-$(date +%s)"
echo "Creating new Netlify site: $SITE_NAME"

# Create a new site and extract site ID
SITE_ID=$(netlify sites:create --name "$SITE_NAME" --json | jq -r '.id')

if [ -z "$SITE_ID" ]; then
  echo "Failed to create Netlify site"
  exit 1
fi

echo "Site created with ID: $SITE_ID"

# Link the site
echo "Linking to site..."
netlify link --id "$SITE_ID" --force

# Set environment variables
echo "Setting environment variables..."
netlify env:set NEXT_PUBLIC_API_URL "https://cryoprotect-8030e4025428.herokuapp.com/v1"
netlify env:set NEXT_PUBLIC_USE_MOCK_DATA "false"
netlify env:set NEXT_PUBLIC_ENABLE_API_LOGGING "true"
netlify env:set NEXT_PUBLIC_ENVIRONMENT "production"

# Generate a random NEXTAUTH_SECRET
NEXTAUTH_SECRET=$(openssl rand -base64 32)
netlify env:set NEXTAUTH_SECRET "$NEXTAUTH_SECRET"

# Set up protection bypass
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V"
netlify env:set PROTECTION_BYPASS "$PROTECTION_BYPASS"
netlify env:set NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS "$PROTECTION_BYPASS"

# Deploy to Netlify
echo "Deploying to Netlify..."
netlify deploy --prod --dir=.next

# Print site URL
echo "Deployment complete!"
echo "Your site is available at: https://$SITE_NAME.netlify.app"