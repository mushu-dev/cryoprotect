#!/bin/bash
# Script to deploy the frontend to Netlify with environment variables

# Generate a random NEXTAUTH_SECRET if not provided
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

# Get Heroku app name from environment or default
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect}

# Construct the backend API URL
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"

# Get Netlify site name from environment or default
NETLIFY_SITE_NAME=${NETLIFY_SITE_NAME:-cryoprotect}
FRONTEND_URL="https://${NETLIFY_SITE_NAME}.netlify.app"

# Get protection bypass token from environment or use default
PROTECTION_BYPASS=${PROTECTION_BYPASS:-TAt23KbtFE8dkZobJU3hpgTP4L5ja07V}

echo "Deploying frontend to Netlify with the following configuration:"
echo "Backend API URL: $API_URL"
echo "Frontend URL: $FRONTEND_URL"
echo "Protection Bypass: Configured"

# First, make sure any API endpoints in the code are properly using environment variables
echo "Updating API endpoints in code..."
npm run update-api-endpoints

# Set environment variables in Netlify
echo "Setting environment variables in Netlify..."
netlify env:set NEXT_PUBLIC_API_URL "$API_URL"
netlify env:set NEXT_PUBLIC_USE_MOCK_DATA false
netlify env:set NEXT_PUBLIC_ENABLE_API_LOGGING true
netlify env:set NEXT_PUBLIC_ENVIRONMENT production
netlify env:set NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS "$PROTECTION_BYPASS"
netlify env:set NEXTAUTH_URL "$FRONTEND_URL"
netlify env:set NEXTAUTH_SECRET "$NEXTAUTH_SECRET"
netlify env:set PROTECTION_BYPASS "$PROTECTION_BYPASS"

# Build the project
echo "Building the project..."
npm run build

# Deploy to Netlify
echo "Deploying to Netlify..."
netlify deploy --prod

# Once deployment is complete, verify and notify the Heroku backend
if [ $? -eq 0 ]; then
  echo "Netlify deployment successful!"
  
  # Verify the deployment
  echo "Verifying the deployment..."
  NETLIFY_URL="$FRONTEND_URL" NEXT_PUBLIC_API_URL="$API_URL" npm run verify-netlify
  
  echo "Notifying Heroku backend about frontend URL..."
  
  # Only run this if HEROKU_API_KEY is set
  if [ -n "$HEROKU_API_KEY" ]; then
    curl -X PATCH \
      -H "Content-Type: application/json" \
      -H "Accept: application/vnd.heroku+json; version=3" \
      -H "Authorization: Bearer $HEROKU_API_KEY" \
      -d "{\"config\": {\"NETLIFY_FRONTEND_URL\": \"$FRONTEND_URL\"}}" \
      "https://api.heroku.com/apps/${HEROKU_APP_NAME}/config-vars"
      
    echo "Heroku configuration updated with frontend URL."
  else
    echo "HEROKU_API_KEY not set. Skipping Heroku configuration update."
    echo "Please manually set NETLIFY_FRONTEND_URL config var on Heroku to: $FRONTEND_URL"
  fi
  
  echo "Deployment complete! Your app is now live at: $FRONTEND_URL"
else
  echo "Netlify deployment failed."
fi